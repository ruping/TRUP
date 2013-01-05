#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;

my $breakpoint_file = shift;  # a sorted breakpoint file

open BP, "$breakpoint_file";

my %breakpoints;          #storing data
my %breakpoints_coor;     #storing all breakpoints

my $last_chrom = 'SRP';
my $last_coor = 0;

while ( <BP> ) {
  chomp;
  my ($chr, $bp, $tag) = split /\t/;

  if ( $chr ne $last_chrom ) {
    $breakpoints{$chr}{$bp}->{$tag} = '';
    $last_chrom = $chr;
    $last_coor = $bp;
    my $bp_coor = $chr."\t".$bp;
    $breakpoints_coor{$bp_coor} = '';
  }
  else {
    my $distance = abs($last_coor - $bp);
    if ($distance > 222) {
      $breakpoints{$chr}{$bp}->{$tag} = '';
      $last_coor = $bp;
      my $bp_coor = $chr."\t".$bp;
      $breakpoints_coor{$bp_coor} = '';
    }
    else {
      $breakpoints{$chr}{$last_coor}->{$tag} = '';
    }
  }
}
close BP;

#print Dumper (\%breakpoints_coor);

#now we have a coordinate hash containing a hash with keys indicating tags
#now reverse the data structure

my %tag_bp;
foreach my $chr (sort keys %breakpoints){
  foreach my $bp (sort {$a <=> $b} keys %{$breakpoints{$chr}}){
    my $bpc = $chr."\t".$bp;
    foreach my $tag (keys %{$breakpoints{$chr}{$bp}}){
       push (@{$tag_bp{$tag}}, $bpc);
    }
  }
}

#print Dumper (\%tag_bp);
%breakpoints = ();

#check pair
my %bp_pair;
my %bp_single;


foreach my $tag (keys %tag_bp) {
  if (@{$tag_bp{$tag}} > 2){
    #print STDERR "this tag: $tag contains more than 2 breakpoints... ";
    #print STDERR Dumper (\@{$tag_bp{$tag}});
    next;
  }
  if (@{$tag_bp{$tag}} == 1) {
    if ($bp_single{${$tag_bp{$tag}}[0]}{'c'} ne ''){
     $bp_single{${$tag_bp{$tag}}[0]}{'c'}++;
    }
    else{
     $bp_single{${$tag_bp{$tag}}[0]}{'c'} = 1;
    }
    next;
  }

  #should check the pair whether far apart
  ${$tag_bp{$tag}}[0] =~ /^(chr\w+)\t(\d+)$/;
  my $chr1 = $1;
  my $bp1 = $2;
  ${$tag_bp{$tag}}[1] =~ /^(chr\w+)\t(\d+)$/;
  my $chr2 = $1;
  my $bp2 = $2;
  next if (($chr1 eq $chr2) and (abs($bp1-$bp2) <= 5000000 ));

  if ($bp_pair{${$tag_bp{$tag}}[0]}{'p'} ne ''){
     $bp_pair{${$tag_bp{$tag}}[0]}{'c'}++;
  }
  else {
    $bp_pair{${$tag_bp{$tag}}[0]}{'p'} = ${$tag_bp{$tag}}[1];
    $bp_pair{${$tag_bp{$tag}}[0]}{'c'} = 1;

  }
  if ($bp_pair{${$tag_bp{$tag}}[1]}{'p'} ne ''){
     $bp_pair{${$tag_bp{$tag}}[1]}{'c'}++;
  }
  else {
    $bp_pair{${$tag_bp{$tag}}[1]}{'p'} = ${$tag_bp{$tag}}[0];
    $bp_pair{${$tag_bp{$tag}}[1]}{'c'} = 1;
  }
}

#print Dumper (\%bp_pair);
#%tag_bp = ();


#generate final breakpoint list
my %printed_pair;
my $id = 1;
foreach my $breakpoint_coor ( keys %breakpoints_coor ){
  if ( $bp_pair{$breakpoint_coor}{'p'} ne '' ) {
     my $breakpoint_pair = $bp_pair{$breakpoint_coor}{'p'};
     my $support_c = $bp_pair{$breakpoint_coor}{'c'};

     my $print_pair = "$breakpoint_coor\t$breakpoint_pair";
     my $print_pair2 = "$breakpoint_pair\t$breakpoint_coor";
     unless (exists $printed_pair{$print_pair} or exists $printed_pair{$print_pair2}){
       print "$id\tp\t$breakpoint_coor\t$support_c\n";
       print "$id\tp\t$breakpoint_pair\t$support_c\n";
       $printed_pair{$print_pair} = '';
       $printed_pair{$print_pair2} = '';
       $id++;
     }
  }
  elsif ( $bp_single{$breakpoint_coor}{'c'} ne '' ) {
     my $support_c = $bp_single{$breakpoint_coor}{'c'};
     print "$id\ts\t$breakpoint_coor\t$support_c\n";
     $id++;
  }
}
