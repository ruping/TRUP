#!/usr/bin/perl
#this is a script to get the sequence of the reference genome for a list of regions provided by breakpointer finder

use strict;

my $breakpoints = shift;
my $genome = shift;

open IN, "$breakpoints";
my %bp_regions;
while ( <IN> ){
   chomp;
   my ($id, $type, $chr, $coor, $support, $pw, $rep) = split /\t/;

   my $start = $coor - 200;
   my $end = $coor + 200;

   if (! exists $bp_regions{$id} ) {
     $bp_regions{$id}{'chr1'} = $chr;
     $bp_regions{$id}{'start1'} = $start;
     $bp_regions{$id}{'end1'} = $end;
   }
   else {
     $bp_regions{$id}{'chr2'} = $chr;
     $bp_regions{$id}{'start2'} = $start;
     $bp_regions{$id}{'end2'} = $end;
   }

   $bp_regions{$id}{'bps'} .= "$id\t$chr\t$coor\t$support\t$pw\t$type\n";

}
close IN;

my %genome = ();

open IN, "$genome";

undef $/;
$/ = '>';

while ( <IN> )
 {
   next if (/^>$/);
   s/\n>//;
   /^(chr\w+).*?\n(.*)/s;
   my $chr = $1;
   my $seq = $2;
   $seq =~ s/\n//g;
   $genome{$chr} = $seq;
 }
close IN;    #read Genome into memory


foreach my $bp_id (sort {$a<=>$b} keys %bp_regions){

  my @pors = keys %{$bp_regions{$bp_id}};

  my $chr1;
  my $start1;
  my $end1;
  my $chr2;
  my $start2;
  my $end2;
  my $seq1;
  my $seq2;

  if (scalar(@pors) == 3){
    $chr1   = $bp_regions{$bp_id}{'chr1'};
    $start1 = $bp_regions{$bp_id}{'start1'};
    $end1   = $bp_regions{$bp_id}{'end1'};
    $seq1   = substr ($genome{$chr1}, $start1-1, $end1-$start1+1);
    print ">$chr1\:$start1\-$end1\n$seq1\n";
  }

  elsif (scalar(@pors) == 6){
    $chr1   = $bp_regions{$bp_id}{'chr1'};
    $start1 = $bp_regions{$bp_id}{'start1'};
    $end1   = $bp_regions{$bp_id}{'end1'};

    $chr2   = $bp_regions{$bp_id}{'chr2'};
    $start2 = $bp_regions{$bp_id}{'start2'};
    $end2   = $bp_regions{$bp_id}{'end2'};

    $seq1 = substr ($genome{$chr1}, $start1-1, $end1-$start1+1);
    $seq2 = substr ($genome{$chr2}, $start2-1, $end2-$start2+1);

    print ">$chr1\:$start1\-$end1\n$seq1\n";
    print ">$chr2\:$start2\-$end2\n$seq2\n";
  }

}

exit;
