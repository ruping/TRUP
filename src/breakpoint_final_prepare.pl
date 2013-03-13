#!/usr/bin/perl

use strict;

open IN, shift;

my %remember;
my %forget;
my %redundancy;

#for each coordinate, only remember one

while ( <IN> ) {
   chomp;
   my ($id, $type, $chr, $coor, $support, $pw, $repornot) = split /\t/;
   $remember{$id} .= "$_\n";

   my $bp = $chr."\t".$coor;
   if ( $redundancy{$bp}{'id'} ne '' ) {
     if ( ($type eq 's') and ($redundancy{$bp}{'type'} eq 'p') ){
        $forget{$id} = '';
     }
     elsif ( ($type eq $redundancy{$bp}{'type'}) and ($pw < $redundancy{$bp}{'pw'}) ){
        $forget{$id} = '';
     }
     elsif ( ($type eq $redundancy{$bp}{'type'}) and ($pw == $redundancy{$bp}{'pw'}) and ($support <= $redundancy{$bp}{'su'}) ) {
        $forget{$id} = '';
     }
     else {
        $forget{$redundancy{$bp}{'id'}} = '';
        $redundancy{$bp}{'id'} = $id;
        $redundancy{$bp}{'su'} = $support;
        $redundancy{$bp}{'type'} = $type;
        $redundancy{$bp}{'pw'} = $pw;
     }
   }
   else {
     $redundancy{$bp}{'id'} = $id;
     $redundancy{$bp}{'su'} = $support;
     $redundancy{$bp}{'type'} = $type;
     $redundancy{$bp}{'pw'} = $pw;
   }


   if ($repornot eq 'R') {
     $forget{$id} = '';
   }

   if ( $type eq 'p' and $support < 3 ) {
     $forget{$id} = '';
   }

   if ( $type eq 's' and ($support < 5 or $pw < 3) ){
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
