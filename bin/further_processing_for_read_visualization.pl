#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

my $fusionseqfile;
my $coveragefile;
my $read_length;
my $genomeBlatPred;

GetOptions (
              "fusionseqfile|f=s"  => \$fusionseqfile,
              "coveragefile|c=s"   => \$coveragefile,
              "readlength|r=i"     => \$read_length,
              "genomeBlatPred"     => \$genomeBlatPred,
              "help|h" => sub{
                             print "usage: $0 [options]\n\nOptions:\n\t--fusionseqfile\t\tthe file fusion after filteration seq\n";
                             print "\t--coveragefile\tthe coverage file of fusion candidates\n";
                             print "\t--readlength\tthe length of the read\n";
                             print "\t--genomeBlatPred\tusing genome Blat.\n";
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
     my ($read, $strand, $candidate, $start, $read_seq, $read_qual, $multi, $mismatch) = split /\t/;
     my $end;
     $candidate =~ /Confidence_([10]\.\d+)/;
     $confidence = $1;
     $start += 1;
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

     my %hash_tmp = ('READ'=>$read, 'START'=>$start, 'END'=>$end, 'STRAND'=>$strand, 'READSEQ'=>$new_read_seq);
     push (@{$mapping{$transcript}}, \%hash_tmp);
  }
}
close IN;

foreach my $transcript (sort {($spanScore{$b}<=>$spanScore{$a}) or ($fusion_point_cov{$b}<=>$fusion_point_cov{$a})} keys %mapping){
  my $first_start;
  my $indi = 0;
  print "$transcript\n";

  $transcript =~ /^#(.+?)\t/;
  my $seq = $seq{$1};
  my @starts = map { $_->{'START'} } sort { $a->{'START'} <=> $b->{'START'} } @{$mapping{$transcript}};
  my $cutlength = $starts[$#starts]+$read_length-$starts[0];
  my $cutseq = substr($seq, $starts[0]-1, $cutlength);
  print "$cutseq\n";

  foreach my $mapping (sort { $a->{'START'} <=> $b->{'START'} } @{$mapping{$transcript}}){
      my $start = $mapping->{'START'};
      if ($indi == 0){$first_start = $start; $indi = 1;}
      my $end = $mapping->{'END'};
      my $strand = $mapping->{'STRAND'};
      my $read_seq = $mapping->{'READSEQ'};
      my $readname = $mapping->{'READ'};

      my $space = $start-$first_start;
      my $prefix = ' 'x$space;
      my $printer = $prefix.$read_seq;
      print "$printer\t$readname\n";
  }
}

exit;
