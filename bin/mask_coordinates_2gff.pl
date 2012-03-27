#!/usr/bin/perl

#TODO: use another file to mask the record in this file(both gff file)
#      by the comparison of coordinates in 2 gff files
#      written by Ruping Sun

use strict;

my $maskfile;
my $original;
my $t=0; #tolerant
my %printer;

while ($ARGV[0]){
  my $arg = shift @ARGV;
  if ($arg =~ /-m/){$maskfile = shift @ARGV;}
  elsif ($arg =~ /-t/){$t     = shift @ARGV;}
  elsif ($arg =~ /-o/){$original = shift @ARGV;}
  elsif ($arg =~ /-h/){print "useage: mask_coordinates_2gff.pl -o: origninal_filename -m: maskfile\n";}  
}

#the original file is loaded in a hash
open ORI, "$original";
my %original;
while ( <ORI> ){
  next if /^#/;
  chomp;
  my @col =split /\t/;
  my $chr = $col[0];
  $chr =~ s/chr//;
  $chr = 'MT' if ($chr eq 'M');
  my $start = $col[3];
  my $end = $col[4];
  push (@{$original{$chr}{$start}{$end}}, $_);
}
close ORI;

#the mask file is loaded in a hash as well
open MASK, "$maskfile";
my %mask;
while ( <MASK> ){
  next if /^#/;
  chomp;
  my @col =split /\t/;
  my $chr = $col[0];
  $chr =~ s/chr//;
  $chr = 'MT' if ($chr eq 'M');
  my $start = $col[3];
  my $end = $col[4];
  push (@{$mask{$chr}{$start}{$end}}, $_);
}
close MASK;


my $counter;
foreach my $chr (sort {$a cmp $b} keys %original){
  my @sortedpos = sort {$a <=> $b} keys %{$mask{$chr}};  #each chromosome sorted mask record into a array
  my $ptr = 0;
  foreach my $start (sort {$a <=> $b} keys %{$original{$chr}}){
    foreach my $end (sort {$a <=> $b} keys %{$original{$chr}{$start}}){


       while (($ptr<$#sortedpos) and ((maxmaskend($sortedpos[$ptr],$chr)+$t) < $start)){$ptr++;}  #set pointer close to the mask_start
       #add 1 by 1 to test
       my $off=0;
       while (($ptr+$off<$#sortedpos) and (($sortedpos[$ptr+$off]-$t) < $end)){
         my $mstart = $sortedpos[$ptr+$off];
	 foreach my $mend (sort {$a <=> $b} keys %{$mask{$chr}{$mstart}}){
	  #test if overlapped
	  if ( $end >= ($mstart-$t) and $start <= ($mend+$t) ){
	     foreach my $mask (@{$mask{$chr}{$mstart}{$mend}}){
	       push (@{$printer{$chr}{$start}{$end}}, $mask);
	     }
	   }
	 }
         $off++;
       }

       foreach my $record (@{$original{$chr}{$start}{$end}}) {
         print "$record";
         foreach my $printer (@{$printer{$chr}{$start}{$end}}) {
           #print "$printer\n";
         }
         if (@{$printer{$chr}{$start}{$end}} != 0){
            print "\tyes\n";
         }
         else{
            print "\tno\n";
         }
       }

       %printer = undef;
       $counter++;
       print STDERR "$counter\n" if ($counter%100000 == 0);

    }
  }
}

=pod
#now print the result
foreach my $chr (sort {$a cmp $b} keys %original){
  foreach my $start (sort {$a <=> $b} keys %{$original{$chr}}){
    foreach my $end (sort {$a <=> $b} keys %{$original{$chr}{$start}}){
      foreach my $record (@{$original{$chr}{$start}{$end}}) {
        print "$record\n";
	foreach my $printer (@{$printer{$chr}{$start}{$end}}) {
	  print "$printer\n";
	}
      }
    }
  }
}
=cut

sub maxmaskend{
  my ($mstart, $chr) = @_;
  my @tmp = sort {$b <=> $a} keys %{$mask{$chr}{$mstart}};
  return $tmp[0];
}
