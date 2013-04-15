#!/usr/bin/perl
#TODO: to get the real fusion candidates

use strict;
use Data::Dumper;
use Getopt::Long;

my %opts = (
            'score'=>20,
            'identity'=>95,
            'overlap'=>20,
            'unique'=>10,
            'extension'=>20,
            'transcripts'=>'',
           );

GetOptions (
               "score|s=f"    => \$opts{'score'},
               "identity|i=f" => \$opts{'identity'},
               "overlap|o=i"  => \$opts{'overlap'},
               "unique|u=i"   => \$opts{'unique'},
               "extension|e=i"=> \$opts{'extension'},
               "transcripts|t=s" => \$opts{'transcripts'},
               "help|h" => sub{
                                print "usage: $0 [options]\n\nOptions:\n\t--transcripts\tthe assembled transcript file in fasta format.\n";
                                print "\t--score\tminimum score cutoff for BLAT hits, default: $opts{'score'}.\n";
                                print "\t--identity\tminimum percentage identity for a valid segment hits, default: $opts{'identity'}.\n";
                                print "\t--overlap\tmaximum overlapping bases between two connecting segments, default: $opts{'overlap'}.\n";
                                print "\t--unique\tminimum unique bases a segment has to have after masking hits, default: $opts{'unique'}.\n";
                                print "\t--help\t\tprint this help message\n\n";
                                exit 0;
                              }
           );


open TRANS, "$opts{'transcripts'}";
my %transcripts_fasta;
my $transcript_fasta;
while ( <TRANS> ){
   chomp;
   if ($_ =~ /^>(.+)$/){
      $transcript_fasta = $1;
   }
   else {
      $transcripts_fasta{$transcript_fasta} .= "$_\n";
   }
}
close TRANS;


my %blat;
while ( <> ) {
   chomp;
   next unless /^\d/;

   my @cols = split /\t/;

   my $blat;
   (
    $blat->{'matches'},
    $blat->{'misMatches'},
    $blat->{'repMatches'},
    $blat->{'nCount'},
    $blat->{'qCountGap'},
    $blat->{'qBasesGap'},
    $blat->{'tCountGap'},
    $blat->{'tBasesGap'},
    $blat->{'strand'},
    $blat->{'qName'},
    $blat->{'qSize'},
    $blat->{'qStart'},
    $blat->{'qEnd'},
    $blat->{'tName'},
    $blat->{'tSize'},
    $blat->{'tStart'},
    $blat->{'tEnd'},
    $blat->{'blockCount'},
    $blat->{'blockSizes'},
    $blat->{'qStarts'},
    $blat->{'tStarts'}
   ) =  @cols;

   $blat->{'qName'} =~ /^Locus_\d+_Transcript_\d+\/\d+_Confidence_([^_]+)/;
   $blat->{'confidence'} = $1;

   $blat->{'perc_id'} = &calc_percent_identity(@cols);
   $blat->{'score'}   = ($blat->{'matches'} + ($blat->{'repMatches'} >> 1))-$blat->{'misMatches'}-$blat->{'qCountGap'}-$blat->{'tCountGap'};
   $blat->{'ratio'}   = ($blat->{'qEnd'} - $blat->{'qStart'} + 1)/$blat->{'qSize'};

   unless ( $blat->{'misMatches'}/$blat->{'matches'} >= 0.1 or $blat->{'misMatches'} >= 20 ) {  #if the mismatches are too much fraction
      push (@{$blat{$blat->{'qName'}}}, $blat);
   }

}

#print Dumper (\@{$blat{'Locus_1_Transcript_3/3_Confidence_0.714_Length_1008_id_51705'}});

my $count;

#for each transcript, there is an array of blat segments
foreach my $transcript (keys %blat) {

  my @blat = sort { $a->{'qStart'} <=> $b->{'qStart'} or $b->{'matches'} <=> $a->{'matches'} } @{$blat{$transcript}};
  my $transLength;

  my %baseHitCount;
  my $max_ratio = 0;
  foreach my $blat (@{$blat{$transcript}}) {
    if ( $blat->{'perc_id'} > $opts{'identity'} ) {
      for (my $i = $blat->{'qStart'}+1; $i<= $blat->{'qEnd'}; $i++) {
        $baseHitCount{$i} += $blat->{'perc_id'};
      }
    }
    $transLength = $blat->{'qSize'} if( $transLength < $blat->{'qSize'} );
    $max_ratio = $blat->{'ratio'} if ( $max_ratio < $blat->{'ratio'} );
  } # each blat

  next if $max_ratio > 0.90;

#=pod
  #mask segments with too many non-unique hits
  my @blat_uniq;
  foreach my $blat (@blat) {
    my $numUniqBase=0;
    for (my $i = $blat->{'qStart'}+1; $i <= $blat->{'qEnd'}; $i++) {
      if( !defined $baseHitCount{$i} || $baseHitCount{$i} <= 100 ){
        $numUniqBase++;
      }
    }
    push @blat_uniq, $blat if($numUniqBase >= $opts{'unique'});
  }
  @blat = @blat_uniq;
#=cut
  #print Dumper (\@blat) if $transcript eq 'Locus_1_Transcript_3/3_Confidence_0.714_Length_1008_id_51705';

  next if ( scalar(@blat) < 2 );

  #get connections
  for (my $i=0; $i<=$#blat; $i++) {
    $blat[$i]->{'id'}=$i;
    my $a = $blat[$i];
    next unless ($a->{'perc_id'} >= $opts{'identity'});

    #pair-wise connection
    for(my $j=$i+1; $j<=$#blat; $j++){
      my $b=$blat[$j];
      if(    $b->{'qStart'} <= $a->{'qEnd'} + 1                   #contiguous: no unmapped base in between
         and $b->{'qEnd'}   >  $a->{'qEnd'} + $opts{'extension'}  #Must extend at least
         and $b->{'qStart'} >= $a->{'qEnd'} - $opts{'overlap'}    #Allow at most $opts{'overlap'} bp overlap
         and $b->{'perc_id'}>= $opts{'identity'}                  #minimum percent_identify
        ){
        push @{$a->{'next'}},$j;  #connect
        push @{$b->{'prev'}},$i;  #connect
      }
    } #foreach j
  } #foreach i

  #print "$transcript\n";
  #print Dumper (\@blat);

  my @paths=&GetAllPaths(@blat);
  my ($PathBest, $Path2Best) =(0,0);
  my ($ScoreBest, $Score2Best)=(0,0);
  my ($leftMarginBest,$rightMarginBest);
  foreach my $path (@paths) {
    my $score=0;
    my $a;
    my ($leftMargin,$rightMargin);
    my $tNameSwitch=0;
    foreach my $id(split /\./,$path){
      my $b = $blat[$id];
      if(defined $a){
        my $overlap = $a->{'qEnd'} - $b->{'qStart'};
        $score += $b->{'score'}-1-(($overlap>0)?$overlap:0);
        $tNameSwitch++ if ($b->{'tName'} ne $a->{'tName'});
      }
      else{
        $score += $b->{'score'};
        $leftMargin = $b->{'qStart'};
      }
      $rightMargin= $transLength - $b->{'qEnd'};
      $a=$b;
    }
    next if($tNameSwitch > 2);

    if( $ScoreBest < $score ){
      $ScoreBest = $score;
      $PathBest = $path;
      ($leftMarginBest,$rightMarginBest) = ($leftMargin,$rightMargin);
    }
    elsif($Score2Best < $score) {
      $Score2Best = $score;
      $Path2Best = $path;
    }
  }

  #print "$PathBest\t$ScoreBest\n";
  #print "$Path2Best\t$Score2Best\n";

  my $qRealSize = $transLength;
  $qRealSize -= $leftMarginBest if(defined $leftMarginBest && $leftMarginBest <= $opts{'s'}+1);
  $qRealSize -= $rightMarginBest if(defined $rightMarginBest && $rightMarginBest <= $opts{'s'}+1);
  my $chimeric_score = exp(($ScoreBest-$qRealSize)/10)-exp(($Score2Best-$qRealSize)/10);

  my $pb;
  my @breakpoints;
  my $PathBestSize;

  foreach my $idx (split /\./,$PathBest) {   # Assessing breakpoints between nodes and blats

    my $a = $blat[$idx];
    my $newHit = 1;

    #print Dumper (\$a) if $transcript eq 'Locus_1_Transcript_6/9_Confidence_0.381_Length_566_id_10615';;

    $PathBestSize += $a->{'qEnd'} - $a->{'qStart'};

    my @qs = split( /\,/, $a->{'qStarts'} );  #0-based
    my @ts = split( /\,/, $a->{'tStarts'} );  #0-based
    my @bsize = split( /\,/, $a->{'blockSizes'} );

    my $qbstart = $a->{'qStart'}; #query breakpoint start
    my $blatStart = $qbstart + 1;
    my $blatGenomeStart;
    my $bpInner = 0;

    foreach my $i ( 0..($a->{'blockCount'}-1) ) {
      my $b;
      #in PSL, negative strand alignments starts from 3' (reverse complemented)
      if ( $a->{'strand'} eq '-' ) {
        $b->{'qStart'} = $qbstart + 1;  #1-based
        $b->{'qEnd'} = $b->{'qStart'} + $bsize[$a->{'blockCount'} - 1 - $i] - 1;  #1-based
        $qbstart = $b->{'qEnd'};
        $b->{'tEnd'} = $ts[$a->{'blockCount'} - 1 - $i] + 1;  #1-based
        $b->{'tStart'} = $b->{'tEnd'} + $bsize[$a->{'blockCount'} - 1 - $i] - 1;   #1-based
      }
      else {
        $b->{'qStart'} = $qs[$i] + 1;  #1-based
        $b->{'qEnd'} = $qs[$i] + $bsize[$i];  #1-based

        $b->{'tStart'} = $ts[$i] + 1;  #1-based
        $b->{'tEnd'} = $ts[$i] + $bsize[$i];   #1-based
      }

      $b->{'tName'} = $a->{'tName'};
      $b->{'strand'} = $a->{'strand'};
      $b->{'perc_id'} = $a->{'perc_id'};
      $b->{'blatEnd'} = $b->{'qEnd'};
      $blatGenomeStart = $b->{'tStart'} if (!defined $blatGenomeStart);
      $b->{'blatGenomeEnd'} = $b->{'tEnd'};

      my $bktype;
      if( defined $pb ) {
        if( $b->{'tName'} ne $pb->{'tName'} ){
          $bktype='inter_C';
          $blatStart = $a->{'qStart'}+1;             #set new blatStart
          $blatGenomeStart = $b->{'tStart'};
        }
        else {
          $bktype='intra_C';
          if ( abs($b->{'tStart'} - $pb->{'tEnd'}) > 230000 ) {
            $blatStart = ($newHit == 1)?($a->{'qStart'}+1):($b->{'qStart'});         #set new blatStart
            $blatGenomeStart = $b->{'tStart'};
          }
          elsif (scalar(@breakpoints) > 0){
            $breakpoints[$#breakpoints]->{'bk2'}->{'blatEnd'} = $b->{'qEnd'};
            $breakpoints[$#breakpoints]->{'bk2'}->{'blatGenomeEnd'} = $b->{'tEnd'};
          }
        }
        if( defined $bktype ) {
          my $bk;
          #if( &ChromGT($pb->{'tName'}, $b->{'tName'}) ) {
            #$pb->{'strand'} = ( $pb->{'strand'} eq '-' )?'+':'-';
            #$b->{'strand'} = ($b->{'strand'} eq '-' )?'+':'-';

            #my $tmp = $b->{'tStart'};
            #$b->{'tStart'} = $b->{'tEnd'};
            #$b->{'tEnd'} = $tmp;

            #$tmp = $pb->{'tStart'};
            #$pb->{'tStart'} = $pb->{'tEnd'};
            #$pb->{'tEnd'} = $tmp;

            #$bk->{'bk1'} = $b;
            #$bk->{'bk2'} = $pb;

            #$bk->{'reverse_complement'}=1;
          #}
          #else {
            $bk->{'bk1'} = $pb;
            $bk->{'bk2'} = $b;
          #}

          $bk->{'bktype'}=$bktype;
          $bk->{'score'} = $chimeric_score;
          $bk->{'perc_id'} = ($b->{'perc_id'} + $pb->{'perc_id'})/2;
          unless ( $bktype eq 'intra_C' and abs($b->{'tStart'}-$pb->{'tEnd'}) <= 230000 ){ #now you see a "breakpoint"
             $bpInner = 1 if ($i > 0);
             push (@breakpoints, $bk) if $i == 0;
          }
        }
      }
      else {
        $blatGenomeStart = $b->{'tStart'};
      }

      $b->{'blatStart'} = $blatStart;
      $b->{'blatGenomeStart'} = $blatGenomeStart;
      if ($i == ($a->{'blockCount'}-1) and $bpInner == 0) {  #last block and there is no inner breakpoint
        $b->{'blatEnd'} = $a->{'qEnd'} if $b->{'blatEnd'} != $a->{'qEnd'};
        $b->{'blatGenomeEnd'} = $a->{'tEnd'} if ($b->{'blatGenomeEnd'} != $a->{'tEnd'} and $b->{'strand'} eq '+');
        $b->{'blatGenomeEnd'} = $a->{'tStart'}+1 if ($b->{'blatGenomeEnd'} != ($a->{'tStart'}+1) and $b->{'strand'} eq '-');
      }
      $pb = $b;
      $newHit = 0 if ($newHit == 1);

    } #foreach block of a hit
  } #get breakpoint foreach path node (blat)

  next if ($PathBestSize/$transLength < 0.5);
  if (scalar(@breakpoints) > 0) {
    $count++;
  }

  next if ($chimeric_score == 0);
  my $print_fasta = 0;
  my $breakpoint_index = 0;
  foreach my $breakpoint (@breakpoints) {
    $breakpoint_index++;
    print "$transcript\t$breakpoint_index\t$transLength\t";
    $print_fasta = 1 if ($print_fasta == 0);

    my $bp1 = $breakpoint->{'bk1'};
    my $bp2 = $breakpoint->{'bk2'};
    my $bpPos1 = $bp1->{'blatEnd'};
    my $bpPos2 = $bp2->{'blatStart'};
    my ($blat1, $blat2);
    if ($bp1->{'strand'} eq '+'){
       $blat1 = $bp1->{'blatStart'}.'-'.$bp1->{'blatEnd'}.'('.$bp1->{'strand'}.')'.$bp1->{'tName'}.':'.$bp1->{'blatGenomeStart'}.'-'.$bp1->{'blatGenomeEnd'};
    } else {
       $blat1 = $bp1->{'blatStart'}.'-'.$bp1->{'blatEnd'}.'('.$bp1->{'strand'}.')'.$bp1->{'tName'}.':'.$bp1->{'blatGenomeEnd'}.'-'.$bp1->{'blatGenomeStart'};
    }
    if ($bp2->{'strand'} eq '+'){
       $blat2 = $bp2->{'blatStart'}.'-'.$bp2->{'blatEnd'}.'('.$bp2->{'strand'}.')'.$bp2->{'tName'}.':'.$bp2->{'blatGenomeStart'}.'-'.$bp2->{'blatGenomeEnd'};
    } else {
       $blat2 = $bp2->{'blatStart'}.'-'.$bp2->{'blatEnd'}.'('.$bp2->{'strand'}.')'.$bp2->{'tName'}.':'.$bp2->{'blatGenomeEnd'}.'-'.$bp2->{'blatGenomeStart'};
    }

    print "$bpPos1\.\.$bpPos2\t";
    print "$blat1\t$blat2\t";
    printf "%.3f\n", $chimeric_score;

    #if ($transcript eq 'Locus_3_Transcript_2/2_Confidence_0.667_Length_459_id_20430'){
    #  print Dumper ( \($breakpoint->{'bk1'}) );
    #  print Dumper ( \($breakpoint->{'bk2'}) );
    #}
  }

  if ($print_fasta == 1){
     my $seq = $transcripts_fasta{$transcript};
     print STDERR "\>$transcript\n$seq";
  }

} #each transcript

#print "$count\n";

sub calc_percent_identity {
  my @cols = @_;
  my $perc_id = (100.0 - (&pslCalcMilliBad(@cols) * 0.1));
  return $perc_id;
}

sub pslCalcMilliBad { #this function is borrowed and modified from blat wesite

  my @cols = @_;

  # sizeMul set to 1 for nucl sequence (3 for protein)
  my $sizeMul = 1;

  my $milliBad = 0;

  my $qAliSize = $sizeMul * ($cols[12] - $cols[11]); #qEnd-qStart
  my $tAliSize = $cols[16] - $cols[15];  #tEnd-tStart

  # aliSize is the minimum of qAliSize and tAliSize
  my $aliSize;
  $qAliSize < $tAliSize ? $aliSize = $qAliSize : $aliSize =$tAliSize;

  # return 0 if AliSize == 0
  return 0 if ($aliSize <= 0);

  # size difference between q and t
  my $sizeDif = $qAliSize - $tAliSize;
  $sizeDif = 0 if ($sizeDif < 0); #mRNA

  # insert Factor
  my $insertFactor = $cols[4]; #qNumInsert

  my $total =  ($sizeMul * ($cols[0] + $cols[2]+ $cols[1])); #mathces + repMatches + misMatches

  if ($total != 0) {
    $milliBad = (1000 * ($cols[1]*$sizeMul + $insertFactor + &round(3*log(1 + $sizeDif)))) / $total;
  }

  return $milliBad;
}

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub GetAllPaths{
  my @blat = @_;
  my @allpaths;
  foreach my $blat (@blat) {
    next if( defined $blat->{'prev'} );  #not the first segment
    my ($paths, $scores) = &GetPath($blat, @blat);
    foreach my $p (@{$paths}) {
      push(@allpaths, $p);
    }
  }
  return @allpaths;
}

sub GetPath {
  my ($node,@blat)=@_;
  my @paths;
  my @scores;
  if( !defined $node->{'next'} ){   #the last segment
    push @paths,  $node->{'id'};
    push @scores, $node->{'score'};
  }
  else {
    foreach my $id (@{$node->{'next'}}) {

      my ($ps,$ss) = &GetPath($blat[$id], @blat);

      my $diff = ($#{$ss}==0)?0:1;
      my $s1;
      foreach my $s2 ( @{$ss} ){
	next if ( !defined $s1 || $s2 ne $s1 );
	$diff = 1;
      }

      if( $diff ) {  #Extension will help distinguish pathes
	for( my $i=0; $i<=$#{$ps}; $i++ ) {
	  my $p=$$ps[$i];
	  my $s=$$ss[$i];
	  my $path=join('.',$node->{'id'},$p);
	  push @paths,  $path;
	  push @scores, $node->{'score'}+$s;
	}
      }
      else {  #Futher extension makes no difference, terminating
	push @paths, join('.',$node->{'id'},$id);
	push @scores,$node->{'score'}+$blat[$id]->{'score'};
      }
    } #foreach next id
  } #not the last element
  return (\@paths,\@scores);
}

sub ChromGT{
  my ($a,$b)=@_;
  my $GT=0;
  $a =~ s/chr//;
  $b =~ s/chr//;
  if($a=~/\D/ || $b=~/\D/){
    $GT=1 if ($a gt $b);
  }
  else{
    $GT=1 if ($a > $b);
  }
  return $GT;
}

