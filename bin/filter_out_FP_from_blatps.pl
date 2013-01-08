#!/usr/bin/perl
#TODO: read the blast ps result to filter out false positives in fusion transcripts

use strict;
use Data::Dumper;
use Getopt::Long;

my ($fusion_bfseq, $fusion_bf, $fusion_blatps);

GetOptions (
               "fusion_bf_seq|fbs=s"  => \$fusion_bfseq,
               "fusion_bf|fb=s"       => \$fusion_bf,
               "fusion_bf_blat|fbb=s" => \$fusion_blatps,
               "help|h" => sub{
                                print "usage: $0 [options]\n\nOptions:\n\t--fusion_bf_seq\t\tthe fasta sequence file of  fusion candidates before filtration\n";
                                print "\t--fusion_bf\tthe fusion candidates file before filtration\n";
                                print "\t--fusion_bf_blat\tthe genome blat file of the fusion candidates\n";
                                print "\t--help\t\tprint this help message\n";
                                exit 0;
                              }


           );

open IN1, "$fusion_bfseq";
my %fusion_bfseq;
my $transcript;
while ( <IN1> ){
   chomp;
   if ($_ =~ /^>(.+)$/){
      $transcript = $1;
   }
   else {
      $fusion_bfseq{$transcript} .= "$_\n";
   }
}
close IN1;


open IN2, "$fusion_bf";
my %title;
my %breakpoint;
my %length;
my %symbol;
while ( <IN2> ){
    chomp;
    my ($transcript, $length, $message1, $message2, $orientation, $breakpoint) = split /\t/;

    $length{$transcript} = $length;

    $breakpoint =~ /^(\d+)\.\.(\d+)$/;

    my $distance = abs($1-$2);
    @{$breakpoint{$transcript}} = sort {$a<=>$b} ($1,$2);

    #-----------------------------------+/-
    $message1 =~ /^(.+?)\(.+?\)$/;
    my $symbol1 = $1;
    $message2 =~ /^(.+?)\(.+?\)$/;
    my $symbol2 = $1;
    my $symbol_combined = $symbol1.'_'.$symbol2;
    $symbol{$symbol_combined} = 1;
    #-----------------------------------+/-

    if ($distance < 10) {
       $title{$transcript} = $transcript.'|'.$length.'|'.$message1.'+'.$message2.'|'.$orientation.'|'.$breakpoint;
    }
}
close IN2;

foreach my $transcript (keys %title){
    $title{$transcript} =~ /^.+?\|\d+\|(.+?)\+(.+?)\|.+?\|\d+\.\.\d+$/;
    my $message1 = $1;
    my $message2 = $2;
    $message1 =~ /^(.+?)\(.+?\)$/;
    my $symbol1 = $1;
    $message2 =~ /^(.+?)\(.+?\)$/;
    my $symbol2 = $1;
    my $symbol_combined = $symbol2.'_'.$symbol1;

    my $ron;
    if (exists $symbol{$symbol_combined}){
      $ron = 'st_both';
    }
    else {
      $ron = 'st_sing';
    }

    $title{$transcript} .= '|'.$ron;
}


open IN3, "$fusion_blatps";
my %blat;
my %trans_confidence;
while ( <IN3> ){
   chomp;
   next unless /^\d/;
   my @cols = split /\t/;
   my $match = $cols[0];
   my $mis = $cols[1];
   my $strand = $cols[8];
   my $transcript = $cols[9];

   $transcript =~ /^Locus_\d+_Transcript_\d+\/\d+_Confidence_([^_]+)/;
   my $confidence = $1;
   my $length = $cols[10];
   my $qs = $cols[11];
   my $qe = $cols[12];
   my $ta = $cols[13];
   my $ts = $cols[15];
   my $te = $cols[16];

   my $blat_len = $qe-$qs+1;
   my $ratio = $blat_len/$length;
   unless ($mis/$match >= 0.03 or $mis >= 10) {
     my %tmp_hash = ('ratio'=>$ratio, 'qs'=>$qs, 'qe'=>$qe, 'ta'=>$ta, 'ts'=>$ts, 'te'=>$te, 'strand'=>$strand);
     push (@{$blat{$transcript}}, \%tmp_hash);
     $trans_confidence{$transcript} = $confidence;
   }
}
close IN3;

#print Dumper(\%trans_confidence);

foreach my $transcript (sort { $trans_confidence{$b} <=> $trans_confidence{$a} } keys %trans_confidence){
  my @ratio = sort { $b->{ratio} <=> $a->{ratio} } @{$blat{$transcript}};

  if ($ratio[0]->{ratio} >=0.99) {  #skip if all mapped
     next;
  }

  if ( scalar(@ratio) != 2 ) {

    unless ( scalar(@ratio) == 3 ){
      next;
    }

    pop(@ratio);

  }



  #filter out some strange blat
  my $blat_start = $ratio[0]->{qs} < $ratio[1]->{qs}? $ratio[0]->{qs} : $ratio[1]->{qs};
  my $blat_end   = $ratio[0]->{qs} < $ratio[1]->{qs}? $ratio[1]->{qe} : $ratio[0]->{qe};
  my $blat_dis   = $blat_end - $blat_start;
  next if ($blat_dis/$length{$transcript} < 0.85);

  #filter out the blat covering the breakpoint
  if ($ratio[0]->{qs} < $breakpoint{$transcript}[0] and $ratio[0]->{qe} > $breakpoint{$transcript}[1]) {

     my $dis1 = $breakpoint{$transcript}[0] - $ratio[0]->{qs};
     my $dis2 = $ratio[0]->{qe} - $breakpoint{$transcript}[1];

     next if ($dis1 > 10 and $dis2 > 10 and $ratio[0]->{ratio} >= 0.9);
     next if ($dis1 > 20 and $dis2 > 20);
  }

  if ($ratio[1]->{qs} < $breakpoint{$transcript}[0] and $ratio[1]->{qe} > $breakpoint{$transcript}[1]){

     my $dis1 = $breakpoint{$transcript}[0] - $ratio[1]->{qs};
     my $dis2 = $ratio[1]->{qe} - $breakpoint{$transcript}[1];

     next if ($dis1 > 10 and $dis2 > 10 and $ratio[1]->{ratio} >= 0.9);
     next if ($dis1 > 20 and $dis2 > 20);
  }

  if ($title{$transcript} ne '') {

      my $title = $title{$transcript};

      my $qs1  = $ratio[0]->{'qs'};
      my $qe1  = $ratio[0]->{'qe'};
      my $chr1 = $ratio[0]->{'ta'};
      my $ts1  = $ratio[0]->{'ts'};
      my $te1  = $ratio[0]->{'te'};
      my $strand1  = $ratio[0]->{'strand'};

      my $qs2  = $ratio[1]->{'qs'};
      my $qe2  = $ratio[1]->{'qe'};
      my $chr2 = $ratio[1]->{'ta'};
      my $ts2  = $ratio[1]->{'ts'};
      my $te2  = $ratio[1]->{'te'};
      my $strand2  = $ratio[1]->{'strand'};

      my $dis_overlap;
      $dis_overlap = abs($qe1-$qs2) if ($qs1 < $qs2);
      $dis_overlap = abs($qe2-$qs1) if ($qs1 > $qs2);
      next if ($dis_overlap > 17);   #further filter out those overlapping of far apart stuff

      my $addinfo = '|'.$qs1.'-'.$qe1.'('.$strand1.')'.$chr1.':'.$ts1.'-'.$te1.'|'.$qs2.'-'.$qe2.'('.$strand2.')'.$chr2.':'.$ts2.'-'.$te2;
      $title .= $addinfo;

      my $chr2 = $ratio[0]->{'ta'};
      my $seq = $fusion_bfseq{$transcript};
      print ">$title\n$seq";
  }

}

exit;
