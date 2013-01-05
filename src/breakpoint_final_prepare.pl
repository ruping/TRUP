#!/usr/bin/perl

use strict;

open IN, shift;

my %remember;
my %forget;
my %redundancy;


while ( <IN> ) {
   chomp;
   my ($id, $type, $chr, $coor, $support, $repornot) = split /\t/;
   $remember{$id} .= "$_\n";

   my $bp = $chr."\t".$coor;
   if ($redundancy{$bp}{'id'} ne '') {
     if ($support <= $redundancy{$bp}{'su'}){
        $forget{$id} = '';
     }
     else {
        $forget{$redundancy{$bp}{'id'}} = '';
        $redundancy{$bp}{'id'} = $id;
        $redundancy{$bp}{'su'} = $support;
     }
   }
   else {
     $redundancy{$bp}{'id'} = $id;
     $redundancy{$bp}{'su'} = $support;
   }


   if ($repornot eq 'R') {
     $forget{$id} = '';
   }

   if ($support < 3) {
     $forget{$id} = '';
   }
}

close IN;


foreach my $id (keys %remember){
  if (! exists($forget{$id})) {
     print "$remember{$id}";
  }
}

exit 0;
