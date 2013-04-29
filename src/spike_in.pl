#!/usr/bin/perl

my $bowtie2 = shift;

open IN, "$bowtie2";
my $last_frag;
my $last_start;
my %fragl;
while ( <IN> ) {

   next if /^\@/;
   chomp;

   my ($frag, $flag, $candidate, $start, $mapQ, $cigar, $mateR, $matePos, $insert, $read_seq, $read_qual, @tags) = split /\t/;
   next unless $insert > 0;
   next if $insert > 1000;
   #my $strand = '+';
   #my $read_end;
   #if ($flag & 16) {$strand = '-';} #the read strand -
   #if ($flag & 64) {$read_end = '1';} #the read is the first in the pair
   #elsif ($flag & 128) {$read_end = '2';} #the read is the second in the pair

   $fragl{$frag} = $insert;
}
close IN;

foreach my $frag (keys %fragl){
    print "$frag\t$fragl{$frag}\n";
}

exit;
