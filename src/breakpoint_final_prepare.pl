#!/usr/bin/perl

use strict;
use Data::Dumper;

open IN, shift;

my %remember;
my @data;
my %pair;

#store data
while ( <IN> ) {

  chomp;
  my ($id, $type, $chr, $coor, $support, $pw, $repornot) = split /\t/;

  $remember{$id} .= "$_\n";

  my $data;
  $data->{'id'} = $id;
  $data->{'type'} = $type;
  $data->{'chr'} = $chr;
  $data->{'coor'} = $coor;
  $data->{'support'} = $support;
  $data->{'pw'} = $pw;
  $data->{'rep'} = $repornot;
  push @data, $data;


  if ($type eq 'p'){
    my $pid;
    $pid->{'chr'} = $chr;
    $pid->{'coor'} = $coor;
    $pid->{'su'} = $support;
    $pid->{'pw'} = $pw;
    $pid->{'rep'} = $repornot;
    push @{$pair{$id}}, $pid;
  }

}
close IN;

my %forget;
my %redundancy;

foreach my $data (@data) {

   my $id = $data->{'id'};
   my $type = $data->{'type'};
   my $chr = $data->{'chr'};
   my $coor = $data->{'coor'};
   my $support = $data->{'support'};
   my $pw = $data->{'pw'};
   my $repornot = $data->{'rep'};

   my $bp = $chr."\t".$coor;
   $redundancy{$bp}{'count'}++;      #remember how many times this $bp occur

   if ( $redundancy{$bp}{'id'} ne '' ) {

     #here is for SINGLE type
     if ( $type eq 's' ) {
       if ( $redundancy{$bp}{'type'} eq 'p' ) {    #the previous one is paired
         $forget{$id} = '';
       } else {                                    #both are single type
         if ($pw < $redundancy{$bp}{'pw'}) {
           $forget{$id} = '';
         } elsif ($pw == $redundancy{$bp}{'pw'} and $support <= $redundancy{$bp}{'su'}) {
           $forget{$id} = '';
         } else {                                  #the current one seems better
            $forget{$redundancy{$bp}{'id'}} = '';  #forget the old one
            $redundancy{$bp}{'id'} = $id;
            $redundancy{$bp}{'su'} = $support;
            $redundancy{$bp}{'type'} = $type;
            $redundancy{$bp}{'pw'} = $pw;
            $redundancy{$bp}{'rep'} = $repornot;
         }
       }
     }    #single type

     #here is for PAIRED type
     elsif ( $type eq 'p' ) {

       if ( $redundancy{$bp}{'type'} eq 'p' ) {      #both are paired ,compare now

          my $old_id = $redundancy{$bp}{'id'};
          my ($old_supportA, $supportA);
          my ($old_pwA, $pwA);
          foreach my $old_bp ( @{$pair{$old_id}} ){
             $old_supportA += $old_bp->{'su'};
             $old_pwA += $old_bp->{'pw'};
          }
          foreach my $new_bp ( @{$pair{$id}} ){
             $supportA += $new_bp->{'su'};
             $pwA += $new_bp->{'pw'};
          }

          if ($supportA < $old_supportA){
             $forget{$id} = '';
          } elsif ($supportA == $old_supportA and $pwA <= $old_pwA) {
             $forget{$id} = '';
          } else {                                     #the current one is better
             $forget{$redundancy{$bp}{'id'}} = '';     #forget the old one
             $redundancy{$bp}{'id'} = $id;
             $redundancy{$bp}{'su'} = $support;
             $redundancy{$bp}{'type'} = $type;
             $redundancy{$bp}{'pw'} = $pw;
             $redundancy{$bp}{'rep'} = $repornot;
          }

       } else {                                        #the old one is single
         $forget{$redundancy{$bp}{'id'}} = '';
         $redundancy{$bp}{'id'} = $id;
         $redundancy{$bp}{'su'} = $support;
         $redundancy{$bp}{'type'} = $type;
         $redundancy{$bp}{'pw'} = $pw;
         $redundancy{$bp}{'rep'} = $repornot;
       }
     }   #paired type
   }     #defined
   else {  #the bp is not defined
     $redundancy{$bp}{'id'} = $id;
     $redundancy{$bp}{'su'} = $support;
     $redundancy{$bp}{'type'} = $type;
     $redundancy{$bp}{'pw'} = $pw;
     $redundancy{$bp}{'rep'} = $repornot;
   }

   if ($repornot eq 'R') {
     $forget{$id} = '';
   }

   if ($chr eq 'chrM') {
     $forget{$id} = '';
   }

   if ( $type eq 'p' and $support < 3 and $pw < 5 ) {
     $forget{$id} = '';
   }

   if ( $type eq 's' and ($support < 5 or $pw < 3) ) {
     $forget{$id} = '';
   }
}


#scan id
my %printed;
foreach my $id (keys %remember) {
  if (! exists($forget{$id})) {
    print "$remember{$id}";
    $printed{$id} = '';
  }
}

#scan coordinates
my @save;
my %meet;
foreach my $bp (keys %redundancy) {
   $bp =~ /^(.+?)\t(\d+)$/;
   my $chr = $1;
   my $coor = $2;
   my $id = $redundancy{$bp}{'id'};        #the best id for this bp
   my $type = $redundancy{$bp}{'type'};
   my $support = $redundancy{$bp}{'su'};
   my $pw = $redundancy{$bp}{'pw'};
   my $count = $redundancy{$bp}{'count'};
   my $repornot = $redundancy{$bp}{'rep'};
   next if ($type eq 's');

   my ($chr2, $coor2, $bp2, $support2, $pw2, $count2, $repornot2);

   if ( ${$pair{$id}}[0]->{'chr'} eq $chr and ${$pair{$id}}[0]->{'coor'} eq $coor ){
      $chr2 = ${$pair{$id}}[1]->{'chr'};
      $coor2 = ${$pair{$id}}[1]->{'coor'};
      $bp2 = $chr2."\t".$coor2;
      $support2 = ${$pair{$id}}[1]->{'su'};
      $pw2 = ${$pair{$id}}[1]->{'pw'};
      $repornot2 = ${$pair{$id}}[1]->{'rep'};
   } else {
      $chr2 = ${$pair{$id}}[0]->{'chr'};
      $coor2 = ${$pair{$id}}[0]->{'coor'};
      $bp2 = $chr2."\t".$coor2;
      $support2 = ${$pair{$id}}[0]->{'su'};
      $pw2 = ${$pair{$id}}[0]->{'pw'};
      $repornot2 = ${$pair{$id}}[0]->{'rep'};
   }
   $count2 = $redundancy{$bp2}{'count'};

   my $supportA = $support + $support2;
   my $pwA = $pw + $pw2;

   next if ($repornot eq 'R' or $repornot2 eq 'R');
   next if ($chr eq 'chrM' or $chr2 eq 'chrM');
   next if ( $supportA < 6 or $pwA < 10 );
   next if ($count < 3);
   if (! exists($printed{$id})) {              #it wasn't printed
     #my $info;
     #$info->{'supportA'} = $supportA;
     #$info->{'pwA'} = $pwA;
     #$meet{$bp}{$id} = $info;
     #$meet{$bp2}{$id} = $info;
     #print "$remember{$id}";
     #$printed{$id} = '';
   } #save it
}



exit 0;
