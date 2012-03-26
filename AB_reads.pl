#!/usr/bin/perl

#accept compressed read files in .gz format

use strict;
use File::Basename;

my $mate1 = shift;
my $mate2 = shift;
my $half = shift;

open M1, "gzip -d -c $mate1 |";
open M2, "gzip -d -c $mate2 |";

my $q_line = 1;

my %reads = ();

while ( <M1> ){

   chomp;
   if ($q_line == 1) {
     if ($_ =~ /^(.+)\/1/) {
       $reads{'title1A'} = $1."A/1";
       $reads{'title1B'} = $1."B/1";
     } elsif ( $_ =~ /^(.+?)\s+(.+)$/ ) {
       $reads{'title1A'} = $1."A ".$2;
       $reads{'title1B'} = $1."B ".$2;
     }
     else {
       #print STDERR "the reads identifier is wird, check it and re-run the pipeline 1\n";
       exit 1;
     }
     $_ = <M2>;
     chomp;
     if ($_ =~ /^(.+)\/2/) {
       $reads{'title2A'} = $1."A/2";
       $reads{'title2B'} = $1."B/2";
     }
     elsif ($_ =~ /^(.+?)\s+(.+)$/) {
       $reads{'title2A'} = $1."A ".$2;
       $reads{'title2B'} = $1."B ".$2;
     }
     else {
       #print STDERR "the reads identifier is wird, check it and re-run the pipeline 2\n";
       exit 1;
     }
     $q_line++;
     next;
   }
   if ($q_line == 2) {
      $half = length($_)/2 if ($half == 0);
      $reads{'read1A'} = substr($_, 0, $half);
      $reads{'read1B'} = substr($_, $half, $half);
      $_ = <M2>;
      chomp;
      $reads{'read2A'} = substr($_, $half, $half);
      $reads{'read2B'} = substr($_, 0, $half);
      $q_line++;
      next;
   }
   if ($q_line == 3) {
     if ($_ =~ /^(.+)\/1/) {
       $reads{'third1A'} = $1."A/1";
       $reads{'third1B'} = $1."B/1";
     }
     elsif ($_ =~ /^(.+?)\s+(.+)$/){
       $reads{'third1A'} = $1."A ".$2;
       $reads{'third1B'} = $1."B ".$2;
     }
     elsif ($_ =~ /^(\+)$/ ) {
       $reads{'third1A'} = $1;
       $reads{'third1B'} = $1;
     }
     else {
       #print STDERR "the reads identifier is wird, check it and re-run the pipeline 3\n";
       exit 1;
     }
     $_ = <M2>;
     chomp;
     if ($_ =~ /^(.+)\/2/) {
       $reads{'third2A'} = $1."A/2";
       $reads{'third2B'} = $1."B/2";
     }
     elsif ($_ =~ /^(.+?)\s+(.+)$/) {
       $reads{'third2A'} = $1."A ".$2;
       $reads{'third2B'} = $1."B ".$2;
     }
     elsif ($_ =~ /^(\+)$/ ) {
       $reads{'third2A'} = $1;
       $reads{'third2B'} = $1;
     }
     else {
       #print STDERR "the reads identifier is wird, check it and re-run the pipeline 4\n";
       exit 1;
     }
     $q_line++;
     next;
   }
   if ($q_line == 4) {
      $reads{'qual1A'} = substr($_, 0, $half);
      $reads{'qual1B'} = substr($_, $half, $half);
      $_ = <M2>;
      chomp;
      $reads{'qual2A'} = substr($_, $half, $half);
      $reads{'qual2B'} = substr($_, 0, $half);
      print "$reads{'title1A'}\n$reads{'read1A'}\n$reads{'third1A'}\n$reads{'qual1A'}\n";
      print "$reads{'title1B'}\n$reads{'read1B'}\n$reads{'third1B'}\n$reads{'qual1B'}\n";
      print STDERR "$reads{'title2A'}\n$reads{'read2A'}\n$reads{'third2A'}\n$reads{'qual2A'}\n";
      print STDERR "$reads{'title2B'}\n$reads{'read2B'}\n$reads{'third2B'}\n$reads{'qual2B'}\n";
      $q_line = 1;
      %reads = ();
      next;
   }

}

close M1;
close M2;
