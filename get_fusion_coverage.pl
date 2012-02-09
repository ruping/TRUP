#!/usr/bin/perl

######
#TODO: 1)transforming the title of fusion candidates (addingthe 5', 3')
###### 2)get fusion coverage from the pair-end mapping and the assembly

use strict;
use Getopt::Long;

my $type;
my $mapping;
my $read_length;
my $accepthits;
my $encomcov;
my $gene_annotation;
my $ensemble_ref_name;
my $locnamefile;
my $ucsc_refgene;
my $repeatmasker;

GetOptions (
              "type|t=s"              => \$type,
              "mappingfile|m=s"       => \$mapping,
	      "readlength|r=i"        => \$read_length,
	      "accepthits|ah=s"       => \$accepthits,
	      "encomcov=s"            => \$encomcov,
	      "geneanno|ga=s"         => \$gene_annotation,
              "ensrefname|ern=s"      => \$ensemble_ref_name,
              "locname|loc=s"         => \$locnamefile,
              "refgene|ref=s"         => \$ucsc_refgene,
              "repeatmasker=s"        => \$repeatmasker,
	      "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--type\t\tthe type of the reads, either pair or single\n";
                             print "\t--mappingfile\tthe mappingfile of the reads onto the assembled fusion candidates\n";
			     print "\t--readlength\tthe length of the read\n";
			     print "\t--geneanno\tthe ensemble gene anotation file\n";
                             print "\t--ensrefname\tthe ensemble, refseq, gene name list table\n";
                             print "\t--locname\tloc_gene_name2location file\n";
                             print "\t--refgene\tdownloaded from ucsc the refgene annotation file\n";
			     print "\t--help\t\tprint this help message\n";
                             exit 0;
			     }
           );

#repeat mask
my %repeatmask;
my @rs;       #start array for each chr
my $old_chr;  #checking the chr
my $ptr;      #pointer for repeatmask sub
open REPEAT, "$repeatmasker";
while ( <REPEAT> ){
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    my $chr = $tmp[0];
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
}
close REPEAT;


open ENS, "$ensemble_ref_name";
my %name2ens;
while ( <ENS> ){
   chomp;
   next if /^Ensembl/;
   my ($ensemble_id, $refseq_id, $gene_name) = split /\t/;
   $name2ens{$gene_name} = $ensemble_id;
}
close ENS;

open LOC, "$locnamefile";
my %loc;
while (<LOC>){
   next if /^Gene/;
   chomp;
   my ($gene, $ens, $chr, $start, $end, $strand) = split /\t/;
   unless ($chr eq ''){
      if ($chr !~ /^chr/){ 
         $chr = 'chr'.$chr;
      }
      if ($strand == 1){
         $strand = '+';
      }
      if ($strand == -1){
         $strand = '-';
      }
      $loc{$gene}{'chr'}        = $chr;
      $loc{$gene}{'start'}      = $start;
      $loc{$gene}{'end'}        = $end;
      $loc{$gene}{'strand'}     = $strand;
   }
}
close LOC;

open REFGENE, "$ucsc_refgene";
my %refgenes;
while ( <REFGENE> ){
   chomp;
   next if /^#/;
   my @cols = split /\t/;
   my $refid   = $cols[1];
   my $chr     = $cols[2];
   my $strand  = $cols[3];
   my $start   = $cols[4];
   my $end     = $cols[5];
   $refgenes{$refid}{'chr'}    = $chr;
   $refgenes{$refid}{'start'}  = $start;
   $refgenes{$refid}{'end'}    = $end;
   $refgenes{$refid}{'strand'} = $strand;
}
close REFGENE;


open GA, "$gene_annotation";
my %gene;
while ( <GA> ){
   next if /^#/;
   chomp;
   my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$tag) = split /\t/;
   if ($type eq 'gene'){
       $tag =~ /^ID=(.+?);Name=(.+?);/;
       my $ensemble = $1;
       my $gene = $2;
       $gene{$ensemble}{'chr'}    = $chr;
       $gene{$ensemble}{'start'}  = $start;
       $gene{$ensemble}{'end'}    = $end;
       $gene{$ensemble}{'strand'} = $strand;
       $gene{$gene}{'chr'}        = $chr;
       $gene{$gene}{'start'}      = $start;
       $gene{$gene}{'end'}        = $end;
       $gene{$gene}{'strand'}     = $strand;
   }
}
close GA;


open IN, "$mapping";
my %for_rm;
my %coverage;
my %for_encompass;
while ( <IN> ){

   chomp;
   my @cols = split /\t/;

   if ($type =~ /pair/){ #bowtie
      my ($read, $strand, $candidate, $start, $read_seq, $read_qual, $multi, $mismatch) = @cols;

      $candidate =~ /^(.+?Confidence_([10]\.\d+).*?)\|(\d+)\|(.+?)\+(.+?)\|(.+?)\|(\d+)\.\.(\d+)\|(.+?)\|(.+?)\|(.+)$/;
      my $transcript = $1;
      my $confidence = $2;
      my $length = $3;
      my $part1 = $4;
      my $part2 = $5;
      my $orientation = $6;
      my $bp_s = $7;
      my $bp_e = $8;
      my $ron  = $9;
      my $blat1 = $10;
      my $blat2 = $11;

      $part1 =~ /^(.+?)\((.+?)\)$/;
      my $message1      = $1;
      my $message_ref1  = $2;

      $part2 =~ /^(.+?)\((.+?)\)$/;
      my $message2      = $1;
      my $message_ref2  = $2;

      if ($coverage{$transcript}{'info'} eq '') {

      my ($mg1_chr, $mg1_start, $mg1_end, $mg1_strand, $mg2_chr, $mg2_start, $mg2_end, $mg2_strand);

      #message1
      if ($gene{$message1} ne '') {
         $mg1_chr   = $gene{$message1}{'chr'};
         $mg1_start = $gene{$message1}{'start'};
         $mg1_end   = $gene{$message1}{'end'};
         $mg1_strand= $gene{$message1}{'strand'};
      }
      elsif ($refgenes{$message_ref1} ne '') {
         $mg1_chr   = $refgenes{$message_ref1}{'chr'};
         $mg1_start = $refgenes{$message_ref1}{'start'};
         $mg1_end   = $refgenes{$message_ref1}{'end'};
         $mg1_strand= $refgenes{$message_ref1}{'strand'};
      }
      elsif ($loc{$message1} ne '') {
         $mg1_chr    = $loc{$message1}{'chr'};
         $mg1_start  = $loc{$message1}{'start'};
         $mg1_end    = $loc{$message1}{'end'};
         $mg1_strand = $loc{$message1}{'strand'};
      }
      elsif ($name2ens{$message1} ne '') {
         my $ensemble_id1 = $name2ens{$message1};
         if ($gene{$ensemble_id1} ne ''){
              $mg1_chr   = $gene{$ensemble_id1}{'chr'};
              $mg1_start = $gene{$ensemble_id1}{'start'};
              $mg1_end   = $gene{$ensemble_id1}{'end'};
              $mg1_strand= $gene{$ensemble_id1}{'strand'};
         }
      }

      #message2
      if ($gene{$message2} ne '') {
         $mg2_chr   = $gene{$message2}{'chr'};
         $mg2_start = $gene{$message2}{'start'};
         $mg2_end   = $gene{$message2}{'end'};
         $mg2_strand= $gene{$message2}{'strand'};
      }
      elsif ($refgenes{$message_ref2} ne '') {
         $mg2_chr   = $refgenes{$message_ref2}{'chr'};
         $mg2_start = $refgenes{$message_ref2}{'start'};
         $mg2_end   = $refgenes{$message_ref2}{'end'};
         $mg2_strand= $refgenes{$message_ref2}{'strand'};
      }
      elsif ($loc{$message2} ne '') {
         $mg2_chr    = $loc{$message2}{'chr'};
         $mg2_start  = $loc{$message2}{'start'};
         $mg2_end    = $loc{$message2}{'end'};
         $mg2_strand = $loc{$message2}{'strand'};
      }
      elsif (exists $name2ens{$message2}) {
         my $ensemble_id2 = $name2ens{$message2};
         if ($gene{$ensemble_id2} ne ''){
              $mg2_chr   = $gene{$ensemble_id2}{'chr'};
              $mg2_start = $gene{$ensemble_id2}{'start'};
              $mg2_end   = $gene{$ensemble_id2}{'end'};
              $mg2_strand= $gene{$ensemble_id2}{'strand'};
         }
      }

      if ($mg1_chr ne '' and $mg2_chr ne ''){
         $for_encompass{$transcript}{'chr1'}   = $mg1_chr;
         $for_encompass{$transcript}{'start1'} = $mg1_start;
         $for_encompass{$transcript}{'end1'}   = $mg1_end;
         $for_encompass{$transcript}{'chr2'}   = $mg2_chr;
         $for_encompass{$transcript}{'start2'} = $mg2_start;
         $for_encompass{$transcript}{'end2'}   = $mg2_end;
      }

      $blat1 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
      my $blat1_qs  = $1;
      my $blat1_qe  = $2;
      my $blat1_st  = $3;
      my $blat1_chr = $4;
      my $blat1_ts  = $5;
      my $blat1_te  = $6;
      $blat2 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
      my $blat2_qs  = $1;
      my $blat2_qe  = $2;
      my $blat2_st  = $3;
      my $blat2_chr = $4;
      my $blat2_ts  = $5;
      my $blat2_te  = $6;

      #for repeat masker to get the breakpoint genomic positions
      if ($blat1_qs < $blat2_qs){
	  if ($blat1_st eq '+' and $blat2_st eq '+'){
	      $for_rm{$blat1_chr}{$blat1_te}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_ts}{$transcript} = '';
          }
          if ($blat1_st eq '+' and $blat2_st eq '-'){
              $for_rm{$blat1_chr}{$blat1_te}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_te}{$transcript} = '';
          }
          if ($blat1_st eq '-' and $blat2_st eq '+'){
              $for_rm{$blat1_chr}{$blat1_ts}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_ts}{$transcript} = '';
          }
          if ($blat1_st eq '-' and $blat2_st eq '-'){
              $for_rm{$blat1_chr}{$blat1_ts}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_te}{$transcript} = '';
          }
      }
      if ($blat1_qs > $blat2_qs){
          if ($blat1_st eq '+' and $blat2_st eq '+'){
              $for_rm{$blat1_chr}{$blat1_ts}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_te}{$transcript} = '';
          }
          if ($blat1_st eq '+' and $blat2_st eq '-'){
              $for_rm{$blat1_chr}{$blat1_ts}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_ts}{$transcript} = '';
          }
          if ($blat1_st eq '-' and $blat2_st eq '+'){
              $for_rm{$blat1_chr}{$blat1_te}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_te}{$transcript} = '';
          }
          if ($blat1_st eq '-' and $blat2_st eq '-'){
              $for_rm{$blat1_chr}{$blat1_te}{$transcript} = '';
              $for_rm{$blat2_chr}{$blat2_ts}{$transcript} = '';
          }
      }


      my ($range_s, $range_e);
      if ($bp_s > $bp_e){$range_s = $bp_e; $range_e = $bp_s;}
      if ($bp_s <= $bp_e){$range_s = $bp_s; $range_e = $bp_e;}

      my $flag;
      if ($blat1_chr eq $blat2_chr){
         if (abs($blat2_ts - $blat1_ts) < 300000){
             $flag = 'Read_Th';
         }
         else {
             $flag = 'Intra_C';
         }
      }
      else{
         $flag = 'Inter_C';
      }

      my %newinfo = ('con'=>$confidence, 'length'=>$length, 'gene1'=>$part1, 'gene2'=>$part2, 'ron'=>$ron, 'ori'=>$orientation, 'bps'=>$range_s, 'bpe'=>$range_e, 'blat1'=>$blat1, 'blat2'=>$blat2, 'flag'=>$flag);

      $coverage{$transcript}{'info'} = \%newinfo;
      }

      #$read =~ /(.+?)\/([12])/;
      $read =~ /(.+?)[\/\s]([12])/;
      my $read_root = $1;
      my $read_end = $2;

      $start += 1;
      $coverage{$transcript}{$read_root}{$read_end}{$start} = $_;

   }
}

close IN;


#do the repeatmasker for the breakpoint
my %bprepeat;
foreach my $chr (sort keys %for_rm){
    foreach my $pos (sort {$a<=>$b} keys %{$for_rm{$chr}}){
	my $rmflag = repeatmask($chr, $pos);
        if ($rmflag == 1){
	    foreach my $transcript (keys %{$for_rm{$chr}{$pos}}){
                $bprepeat{$transcript} .= 'R';
            }
        }
        if ($rmflag == 0){
            foreach my $transcript (keys %{$for_rm{$chr}{$pos}}){
                $bprepeat{$transcript} .= 'N';
            }
        }
    }
}


#generate the coverage encompassing first
open AH, "samtools view $accepthits |";
my (%buffer1, %buffer2);

#sort out the encompassing chr and positions
my %encoposA;
my %encoposB;
foreach my $transcript (keys %for_encompass) {

       my $chr1    = $for_encompass{$transcript}{'chr1'};
       my $start1  = $for_encompass{$transcript}{'start1'};
       my $end1    = $for_encompass{$transcript}{'end1'};
       my $chr2    = $for_encompass{$transcript}{'chr2'};
       my $start2  = $for_encompass{$transcript}{'start2'};
       my $end2    = $for_encompass{$transcript}{'end2'};

       push (@{$encoposA{$chr1}{$start1}{$end1}}, [$transcript, $start2, $end2]);
       push (@{$encoposB{$chr2}{$start2}{$end2}}, [$transcript, $start1, $end1]);

}

my ($old_chrA, $old_chrB);
my (@ssA, @ssB);
my ($ptrA, $ptrB);
while ( <AH> ) {

    next if /^@/;                #ignore comments
    next if ($_ !~ /NH\:i\:1$/); #ignore multiple mappable reads
    chomp;

    my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;
    next unless ($mateRname eq '=');
    next if ($FLAG & 2);         #ignore correct pair

    if ($Rname ne $old_chrA) {
       @ssA  = sort {$a <=> $b} keys %{$encoposA{$Rname}};
       $ptrA = 0;
    }
    if ($Rname ne $old_chrB) {
       @ssB  = sort {$a <=> $b} keys %{$encoposB{$Rname}};
       $ptrB = 0;
    }

    my @endsA = keys %{$encoposA{$Rname}{$ssA[$ptrA]}};
    while (($ptrA<=$#ssA) and ($endsA[0] < $Pos)){
      $ptrA++;
      @endsA = keys %{$encoposA{$Rname}{$ssA[$ptrA]}};
    }
    if ($ssA[$ptrA] <= $Pos){
       foreach my $array_ref (@{$encoposA{$Rname}{$ssA[$ptrA]}{$endsA[0]}}){
          if ($matePos >= $array_ref->[1] and $matePos <= $array_ref->[2]){
             $buffer1{$array_ref->[0]}{$Qname} = $_;
          }
       }
    }

    my @endsB = keys %{$encoposB{$Rname}{$ssB[$ptrB]}};
    while (($ptrB<=$#ssB) and ($endsB[0] < $Pos)){
       $ptrB++;
       @endsB = keys %{$encoposB{$Rname}{$ssB[$ptrB]}};
    }
    if ($ssB[$ptrB] <= $Pos){
       foreach my $array_ref (@{$encoposB{$Rname}{$ssB[$ptrB]}{$endsB[0]}}){
          if ($matePos >= $array_ref->[1] and $matePos <= $array_ref->[2]){
             $buffer2{$array_ref->[0]}{$Qname} = $_;
          }
       }
    }

    $old_chrA = $Rname;
    $old_chrB = $Rname;

}
close AH;

foreach my $transcript (keys %buffer1){
  foreach my $read (keys %{$buffer1{$transcript}}){
     if (exists $buffer2{$transcript}{$read}){
         push (@{$for_encompass{$transcript}{'hits'}}, $buffer1{$transcript}{$read}, $buffer2{$transcript}{$read});         
         my @cols = split (/\t/, $buffer1{$transcript}{$read});
         my $mapo = join ("_", sort {$a <=> $b} ($cols[3],$cols[7]));
         $for_encompass{$transcript}{'cov'}{$mapo}++;
     }
  }
}

#end of encompassing cov

open ENCOMCOV, ">$encomcov";

foreach my $transcript (sort { $coverage{$b}{info}->{con} <=> $coverage{$a}{info}->{con} } keys %coverage){

   my $length = $coverage{$transcript}{info}->{length};
   my $gene1 = $coverage{$transcript}{info}->{gene1};
   my $gene2 = $coverage{$transcript}{info}->{gene2};
   my $bps = $coverage{$transcript}{info}->{bps}; 
   my $bpe = $coverage{$transcript}{info}->{bpe};
   my $blat1 = $coverage{$transcript}{info}->{blat1};
   my $blat2 = $coverage{$transcript}{info}->{blat2};
   my $flag  = $coverage{$transcript}{info}->{flag};
   my $ori   = $coverage{$transcript}{info}->{ori};
   my $ron   = $coverage{$transcript}{info}->{ron};
   my $bpr;
   $bpr = 'R' if ($bprepeat{$transcript} =~ /R/);
   $bpr = 'N' if ($bprepeat{$transcript} !~ /R/);

   foreach my $read_root (keys %{$coverage{$transcript}}) {

      my @pair = keys %{$coverage{$transcript}{$read_root}};

      if (scalar(@pair) == 1){ #only one end map
         my @starts = keys %{$coverage{$transcript}{$read_root}{$pair[0]}};
         next if (scalar (@starts) >= 2);
         my $start = $starts[0];
         my $end = $start+$read_length-1;
         if ($start < $bps and $end > $bpe) { #spanning
            next if ( ($bps-$start) < 5 or ($end-$bpe) < 5 );    #requiring minimal anchoring length
            $coverage{$transcript}{'spanning'}{$start}++;
            push (@{$coverage{$transcript}{reads}},  $coverage{$transcript}{$read_root}{$pair[0]}{$start});
         }
         elsif (($start < $bps and $end > $bps) or ($start >= $bps and $start <= $bpe)) { #other overlapping with bp, not interesting
            #push (@{$coverage{$transcript}{reads}},  $coverage{$transcript}{$read_root}{$pair[0]}{$start});
         }
      }

      if (scalar(@pair) == 2) { #pair mapped
         my @starts1 = keys %{$coverage{$transcript}{$read_root}{1}};
         my @starts2 = keys %{$coverage{$transcript}{$read_root}{2}};
         next if (scalar (@starts1) >= 2 or scalar (@starts2) >= 2);
         my $start1 = $starts1[0];
         my $start2 = $starts2[0];
         my $end1 = $start1+$read_length-1;
         my $end2 = $start2+$read_length-1;

         my $frag_s;
         my $frag_e;
         if ($start1 <  $start2){$frag_s=$start1;$frag_e=$end2;}
         if ($start1 >= $start2){$frag_s=$start2;$frag_e=$end1;}

         my $flag = 0;

         if ($start1 < $bps and $end1 > $bpe) { # 5'end spanning
           if ( ($bps-$start1) >= 5 or ($end1-$bpe) >= 5 ) {   #minimal anchoring length
             $coverage{$transcript}{'spanning'}{$start1}++;
             $flag = 1;
           }
         }
         elsif (($start1 < $bps and $end1 > $bps) or ($start1 >= $bps and $start1 <= $bpe)){ #5'end overlapping
            $flag = 2;
         }

         if ($flag == 0 and ($start2 < $bps and $end2 > $bpe)) { # 3'end spanning
           if ( ($bps-$start2) >= 5 or ($end2-$bpe) >= 5 ) { #minimal anchoring length
             $coverage{$transcript}{'spanning'}{$start2}++;
             $flag = 3;
           }
         }
         elsif ($flag == 0 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){ #3'end overlapping
            $flag = 4;
         }

         elsif ($flag == 1 and ($start2 < $bps and $end2 > $bpe)) { # 5' end and 3' end both spanning
            if ( ($bps-$start2) >= 5 or ($end2-$bpe) >= 5 ) { #minimal anchoring length
              $flag = 5;
            }
         }

         elsif ($flag == 1 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){  #3'end overlapping and 5'end spanning
            $flag = 6;
         }

         elsif ($flag == 2 and ($start2 < $bps and $end2 > $bpe)) { # 5' end overlapping and 3' end spanning
            if ( ($bps-$start2) >= 5 or ($end2-$bpe) >= 5 ) { #minimal anchoring length
              $coverage{$transcript}{'spanning'}{$start2}++;
              $flag = 7;
            }
         }

         elsif ($flag == 2 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){  #both 5' and 3' end overlapping
            $flag = 8;
         }


         if ($flag =~ /[0248]/ and ($frag_s < $bps and $frag_e > $bpe)){
            if ( ($bps-$frag_s) >= 5 or ($frag_e-$bpe) >= 5 ) { #minimal anchoring length
              $coverage{$transcript}{'encompass'}{$frag_s}++;
              $flag = 9;
            }
         }

         if ($flag =~ /[135679]/) {
            push (@{$coverage{$transcript}{'reads'}},  $coverage{$transcript}{$read_root}{1}{$start1});
            push (@{$coverage{$transcript}{'reads'}},  $coverage{$transcript}{$read_root}{2}{$start2});
         }
      }
   }


   my $span = scalar (keys %{$coverage{$transcript}{'spanning'}}); #pile-up reads only count once
   my $enco = scalar (keys %{$coverage{$transcript}{'encompass'}}); #pile-up frags only count once
   my $cov  = $span.'('.$span.'+'.$enco.')'; 
   my $real_enco = scalar(keys (%{$for_encompass{$transcript}{'cov'}}));
   my $all  = $span + $real_enco;

   next if ($all == 0);

   my $newtitle = join("\t", $transcript, $length, $gene1.'-'.$gene2, $ori, $bps.'..'.$bpe, $bpr, $flag, $ron, $blat1, $blat2, $cov, $real_enco, $all);
   print "#$newtitle\n";
   print ENCOMCOV "#$newtitle\n";
   foreach my $razers (@{$coverage{$transcript}{reads}}){
      print "$razers\n";
   }
   foreach my $hits (@{$for_encompass{$transcript}{'hits'}}){
      print ENCOMCOV "$hits\n";
   }
}

close ENCOMCOV;


exit;

sub repeatmask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr){
	@rs = sort {$a <=> $b} keys %{$repeatmask{$chr}};
	$ptr = 0;
    }
    while (($ptr<=$#rs) and ($repeatmask{$chr}{$rs[$ptr]} < $coor)){
	$ptr++;
    }
    if ($rs[$ptr] <= $coor){
	$flag = 1;
    }
    $old_chr = $chr;
    return $flag;
}

