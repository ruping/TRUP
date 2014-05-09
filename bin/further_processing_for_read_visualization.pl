#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

my $fusionseqfile;
my $coveragefile;
my $read_length;
my $genomeBlatPred;
my $seqType;

GetOptions (
              "fusionseqfile|f=s"  => \$fusionseqfile,
              "coveragefile|c=s"   => \$coveragefile,
              "readlength|r=i"     => \$read_length,
              "genomeBlatPred"     => \$genomeBlatPred,
              "seqType=s"          => \$seqType,
              "help|h" => sub{
                             print "usage: $0 [options]\n\nOptions:\n\t--fusionseqfile\t\tthe file fusion after filteration seq\n";
                             print "\t--coveragefile\tthe coverage file of fusion candidates\n";
                             print "\t--readlength\tthe length of the read\n";
                             print "\t--genomeBlatPred\tusing genome Blat.\n";
                             print "\t--seqType\tthe sequencing type (p)aired-end or (s)ingle-end.\n";
                             print "\t--help\t\tprint this help message\n\n";
                             exit 0;
                            }
           );



open SEQ, "$fusionseqfile";
my $candidate;
my %seq;
while ( <SEQ> ) {
  chomp;
  if ($genomeBlatPred) {
    if ($_ =~ /^>(.+?)$/) {
      $candidate = $1;
    } else {
      $seq{$candidate} .= $_;
    }
  } else {
    if ($_ =~ /^>(.+?)\|/) {
      $candidate = $1;
    } else {
      $seq{$candidate} .= $_;
    }
  }
}
close SEQ;

open IN, "$coveragefile";
my $transcript;
my %mapping;
my %fusion_point_cov;
my %spanScore;
my %alignInDel;  #align insertion deletion
my ($confidence, $length, $part1, $part2, $ori, $bp_s, $bp_e, $ron, $blat1, $blat2);

while ( <IN> ){
  chomp;
  if ($_ =~ /^#/) {
     $transcript = $_;
     my ($candidate, $idx, $cScore, $fusion_genes, $breakpoint, $rep, $sc, $type, $cov1, $cov2, $cov3, $spanAll, $spanScore);
     if ($genomeBlatPred) {
       ($candidate, $idx, $cScore, $length, $fusion_genes, $ori, $breakpoint, $rep, $sc, $type, $ron, $blat1, $blat2, $cov1, $cov2, $cov3, $spanAll, $spanScore) = split (/\t/, $transcript);
     } else {
       ($candidate, $cScore, $length, $fusion_genes, $ori, $breakpoint, $rep, $sc, $type, $ron, $blat1, $blat2, $cov1, $cov2, $cov3, $spanAll, $spanScore) = split (/\t/, $transcript);
     }
     $candidate =~ s/^#//;
     $breakpoint =~ /^(\d+)\.\.(\d+)$/;
     my $bps = $1; $bp_s = $bps;
     my $bpe = $2; $bp_e = $bpe;
     $fusion_genes =~ /^(.+?\(.+?\))\-(.+?\(.+?\))/; $part1 = $1; $part2 = $2;
     my $seq = $seq{$candidate};
     my $s1 = substr($seq, 0, $bps-1);
     my $s2 = substr($seq, $bps-1, $bpe-$bps+1);
     $s2 =~ tr/ACGTN/acgtn/;
     my $s3 = substr($seq, $bpe);
     $seq{$candidate} = $s1.$s2.$s3;
     $fusion_point_cov{$transcript} = $cov3;
     $spanScore{$transcript} = $spanScore;
  }
  else {
     my ($read, $flag, $candidate, $start, $mapQ, $cigar, $mateR, $matePos, $insert, $read_seq, $read_qual, @tags) = split /\t/;
     my $strand = '+';
     my $read_end = 1;
     if ($flag & 16) {$strand = '-';} #the read strand -
     if ($seqType =~ /^p/){
       if ($flag & 64) {$read_end = '1';} #the read is the first in the pair
       elsif ($flag & 128) {$read_end = '2';} #the read is the second in the pair
     }
     $read .= '/'.$read_end;

     my $end;
     $candidate =~ /Confidence_([10]\.\d+)/;
     $confidence = $1;
     $end = $start+$read_length-1;

     my ($range_s, $range_e);
     if ($bp_s > $bp_e){$range_s = $bp_e; $range_e = $bp_s;}
     if ($bp_s <= $bp_e){$range_s = $bp_s; $range_e = $bp_e;}

     my $new_read_seq = $read_seq;

     if ($end >= $range_e and $start <= $range_s){
        my $s1 = substr($read_seq, 0, $range_s-$start);
        my $s2 = substr($read_seq, $range_s-$start, $range_e-$range_s+1);
        $s2 =~ tr/ACGTN/acgtn/;
        my $s3 = substr($read_seq, $range_e-$start+1, $end-$range_e);
        $new_read_seq = $s1.$s2.$s3;
     }

     elsif ($end < $range_e and $end >= $range_s){
        my $s1 = substr($read_seq, 0, $range_s-$start);
        my $s2 = substr($read_seq, $range_s-$start);
        $s2 =~ tr/ACGTN/acgtn/;
        $new_read_seq = $s1.$s2;
     }

     elsif ($start > $range_s and $start <= $range_e){
        my $s1 = substr($read_seq, 0, $range_e-$start+1);
        $s1 =~ tr/ACGTN/acgtn/;
        my $s2 = substr($read_seq, $range_e-$start+1);
        $new_read_seq = $s1.$s2;
     }

     my $cigarPos = 0;
     my $refPos = $start - 1;

     #processing cigar
     while ($cigar =~ /(\d+)([MID])/g) {
        my $blockSize = $1;
        my $blockType = $2;
        if ($blockType eq 'M'){ #match
          $cigarPos += $blockSize;
          $refPos += $blockSize;
        } elsif ($blockType eq 'I'){
          $cigarPos += $blockSize;
          $alignInDel{$transcript}{$refPos}{'isize'} = $blockSize;  #insertion pos in the ref
          $alignInDel{$transcript}{$refPos}{$read} = 1;
        } elsif ($blockType eq 'D'){ #deletion
          $refPos += $blockSize;
          substr($new_read_seq, $cigarPos, 0) = '-'x$blockSize;
          $cigarPos += $blockSize;
        }
     }

     my %hash_tmp = ('READ'=>$read, 'START'=>$start, 'END'=>$refPos, 'STRAND'=>$strand, 'READSEQ'=>$new_read_seq);
     $mapping{$transcript}{$read} = \%hash_tmp;

  }
}
close IN;


foreach my $transcript (sort {($spanScore{$b}<=>$spanScore{$a}) or ($fusion_point_cov{$b}<=>$fusion_point_cov{$a})} keys %mapping){

  my $first_start;
  my $indi = 0;
  print "$transcript\n";

  $transcript =~ /^#(.+?)\t/;
  my $seq = $seq{$1};

  my %currentREADS;
  foreach my $read ( keys %{$mapping{$transcript}} ) {
     my $start = $mapping{$transcript}{$read}->{'START'};
     my $end = $mapping{$transcript}{$read}->{'END'};
     $currentREADS{$start}{$read} = $end;
  }
  my @readStarts = sort {$a <=> $b} keys (%currentREADS);
  my $ptr = 0;

  my %increment;
  my $incre = 0;
  foreach my $iPos (sort {$a <=> $b} keys(%{$alignInDel{$transcript}})) {
    my $reiPos = $iPos + $incre;
    my $iSize = $alignInDel{$transcript}{$iPos}{'isize'};
    substr($seq, $reiPos, 0) = '-'x$iSize;   #add insert part for the reference
    $incre += $iSize;

    while ( $ptr <= $#readStarts ) {
      my @readsOrderEnd = sort { $currentREADS{$readStarts[$ptr]}{$b} <=> $currentREADS{$readStarts[$ptr]}{$a} } keys %{$currentREADS{$readStarts[$ptr]}};
      if ( $currentREADS{$readStarts[$ptr]}{$readsOrderEnd[0]} < $iPos ){
        $ptr++;
      } else {
        last;
      }
    } #ptr

    my $off = 0;
    while ( ($ptr+$off) <= $#readStarts ) {
      my $readStart = $readStarts[$ptr+$off];
      if ( $readStart <= $iPos ) {
        foreach my $read (keys %{$currentREADS{$readStart}}) {
          my $readEnd   = $currentREADS{$readStart}{$read};
          if ( $readEnd > $iPos ) {
            #whether this read contain this insert?
            if ($alignInDel{$transcript}{$iPos}{$read} ne '') {
              $mapping{$transcript}{$read}->{'END'} += $iSize;
            } else {
              #find the iPos position in the read and change it to '-'
              my $readIpos = $iPos - $readStart + 1;
              substr($mapping{$transcript}{$read}->{'READSEQ'}, $readIpos, 0) = '-'x$iSize;   #add the insert part for the normal reads
              $mapping{$transcript}{$read}->{'END'} += $iSize;
            }
          } #end > ipos
        } #foreach read
      } #if start < ipos
      else { #start > ipos
        $increment{$readStart} = $incre;
        #last;
      }
      $off++;
    } #ptr+off
  } #each insert pos

  my $maxIncre = 0;
  my $maxStartwithIncre = -1;
  foreach my $readStart (@readStarts) {
    my $increment;
    if (exists $increment{$readStart}){
      $increment = $increment{$readStart};
      foreach my $read (keys %{$currentREADS{$readStart}}) {
       $mapping{$transcript}{$read}->{'START'} += $increment;
       $mapping{$transcript}{$read}->{'END'} += $increment;
      }
      $maxIncre = $increment;
      $maxStartwithIncre = $readStart;
    } else {
      if ($maxStartwithIncre != -1 and $readStart > $maxStartwithIncre and $maxIncre > 0){
        foreach my $read (keys %{$currentREADS{$readStart}}) {
          $mapping{$transcript}{$read}->{'START'} += $maxIncre;
          $mapping{$transcript}{$read}->{'END'} += $maxIncre;
        }
      }
    } #not recoreded, may be left out
  }


  my @readEnds = sort { $mapping{$transcript}{$b}->{'END'} <=> $mapping{$transcript}{$a}->{'END'} } keys %{$mapping{$transcript}};
  my $lastEnd = $mapping{$transcript}{$readEnds[0]}->{'END'};
  my $cutlength = $lastEnd - $readStarts[0] + 1;
  my $cutseq = substr($seq, $readStarts[0] - 1, $cutlength);

  if ($1 eq "Locus_1_Transcript_8/8_Confidence_0.500_Length_329_id_43505"){
    print "the sequence is: $lastEnd\t$readStarts[0]\t$cutlength\n";
  }

  print "$cutseq\n";

  foreach my $read ( sort { $mapping{$transcript}{$a}->{'START'} <=> $mapping{$transcript}{$b}->{'START'} } keys %{$mapping{$transcript}} ) {
      my $start = $mapping{$transcript}{$read}->{'START'};
      if ($indi == 0){$first_start = $start; $indi = 1;}
      my $read_seq = $mapping{$transcript}{$read}->{'READSEQ'};

      my $space   = $start - $first_start;
      my $prefix  = ' 'x$space;
      my $printer = $prefix.$read_seq;
      print "$printer\t$read\n";
  }
}

exit;
