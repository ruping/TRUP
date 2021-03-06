#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;

my $breakpoint_file = shift;  # a sorted breakpoint file
my $type = shift;
my $maxIntron = shift;

open BP, "gzip -d -c $breakpoint_file |";

my %breakpoints;          #storing data
my %breakpoints_coor;     #storing all breakpoints
my %breakpoints_ct;

my $last_chrom = 'SRP';
my $last_coor = 0;

while ( <BP> ) {
  chomp;
  my ($chr, $bp, $tag, $pw, $ct) = split /\t/;    # pw:proper or wrong; ct:clip type
  my $add = 0;
  $add = 1 if ($pw eq 'w');
  $ct = ($ct eq 'H')?1:0;

  if ( $chr ne $last_chrom ) {
    $breakpoints{$chr}{$bp}->{$tag} = '';
    $last_chrom = $chr;
    $last_coor = $bp;
    my $bp_coor = $chr."\t".$bp;
    $breakpoints_coor{$bp_coor} += $add;
    $breakpoints_ct{$bp_coor} += $ct;
  }
  else {
    my $distance = abs($last_coor - $bp);
    if ($distance > 200) {
      $breakpoints{$chr}{$bp}->{$tag} = '';
      $last_coor = $bp;
      my $bp_coor = $chr."\t".$bp;
      $breakpoints_coor{$bp_coor} += $add;
      $breakpoints_ct{$bp_coor} += $ct;
    }
    else {
      $breakpoints{$chr}{$last_coor}->{$tag} = '';
      my $last_bp = $chr."\t".$last_coor;
      $breakpoints_coor{$last_bp} += $add;
      $breakpoints_ct{$last_bp} += $ct;
    }
  }
}
close BP;


#now we have a coordinate hash containing a hash with keys indicating tags
#now reverse the data structure

my %tag_bp;
foreach my $chr (sort keys %breakpoints){
  foreach my $bp (sort {$a <=> $b} keys %{$breakpoints{$chr}}) {
    my $bpc = $chr."\t".$bp;
    foreach my $tag (keys %{$breakpoints{$chr}{$bp}}) {
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

  if (@{$tag_bp{$tag}} > 2) {
    #print STDERR "this tag: $tag contains more than 2 breakpoints... ";
    next;
  }

  if (@{$tag_bp{$tag}} == 1) {
    my $bp_now = ${$tag_bp{$tag}}[0];
    if ( ($type =~ /^p/ and $breakpoints_coor{$bp_now} > 0) or $type =~ /^s/ ) {
      if ( exists ($bp_single{$bp_now}) ) {
        $bp_single{$bp_now}++;
      }
      else {
        $bp_single{$bp_now} = 1;
      }
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

  next if (($chr1 eq $chr2) and (abs($bp1-$bp2) <= $maxIntron ));

  $bp_pair{${$tag_bp{$tag}}[0]}{${$tag_bp{$tag}}[1]} += 1;
  $bp_pair{${$tag_bp{$tag}}[1]}{${$tag_bp{$tag}}[0]} += 1;

}


#generate final breakpoint list
my %printed_pair;
my $id = 1;
foreach my $breakpoint_coor ( keys %breakpoints_coor ){

  my $wrong_add = $breakpoints_coor{$breakpoint_coor};
  my $ct = $breakpoints_ct{$breakpoint_coor};

  if ( $bp_pair{$breakpoint_coor} ne '' ) {

     my @pairs = sort keys %{$bp_pair{$breakpoint_coor}};

     foreach my $breakpoint_pair (@pairs) {
        my $print_pair = "$breakpoint_coor\t$breakpoint_pair";
        my $print_pair2 = "$breakpoint_pair\t$breakpoint_coor";
        my $support_c = $bp_pair{$breakpoint_coor}{$breakpoint_pair};
        my $wrong_add2 = $breakpoints_coor{$breakpoint_pair};
        my $ct2 = $breakpoints_ct{$breakpoint_pair};

        unless (exists $printed_pair{$print_pair} or exists $printed_pair{$print_pair2}){
          print "$id\tp\t$breakpoint_coor\t$support_c\t$wrong_add\t$ct\n";
          print "$id\tp\t$breakpoint_pair\t$support_c\t$wrong_add2\t$ct2\n";
          $printed_pair{$print_pair} = '';
          $printed_pair{$print_pair2} = '';
          $id++;
        }
     }
  }

  if ( exists ($bp_single{$breakpoint_coor}) ) {
      my $support_c = $bp_single{$breakpoint_coor};
      print "$id\ts\t$breakpoint_coor\t$support_c\t$wrong_add\t$ct\n" unless ($support_c > 1000000);
      $id++;
  }

}
