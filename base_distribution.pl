#!/usr/bin/perl
#TODO: calculate the base distribution around the read
#

use Getopt::Long;

my $read_length;
GetOptions(
              "length|l=i" => \$read_length,
     
          );

open HS, "/scratch/local/ElandSquashedGenomes/RazerS/Hs.dna.fa";
my %genome = ();
my $chr = undef;
while ( <HS> )
 {
   if (/^>(\w+).*?\n$/) {$chr = $1; $chr =~ s/^chr//;}
   else { s/\n//g; s/\s//g; $genome{$chr}.=$_;}
 }
close HS;
print STDERR "genome loaded\n";

my %pos;
while ( <> ){
  chomp;
  my @col = split /\t/;
  my $chr = $col[0];
  $chr =~ s/^chr//;
  my $start = $col[3];
  my $end = $col[4];
  
  next unless (($end-$start) == ($read_length-1)); #skip the splitly mapped reads
  
  my $strand = $col[6]; 
  $col[8] =~ /^Pileup=(\d+)$/;
  my $times = $1;
  my $seq;
  if ($strand eq '+'){
    $seq = substr($genome{$chr},$start-21,$read_length+20);
  }
  if ($strand eq '-'){
    my $seq_o = substr($genome{$chr},$start-1,$read_length+20);
    $seq = reverse $seq_o;
  }
  $pos{$chr}{$start}{p} = $times;
  $pos{$chr}{$start}{s} = $seq;
}


my %base;
for(my $i=1; $i<=($read_length+20); $i++){
  foreach my $chr (keys %pos){
    foreach my $start (keys %{$pos{$chr}}){
      my $base = substr($pos{$chr}{$start}{s},$i-1,1);
      $base{$i}{$base} += $pos{$chr}{$start}{p};
    }
  }
  #here I have at $i-something base kind
  my $total = $base{$i}{A}+$base{$i}{C}+$base{$i}{G}+$base{$i}{T};
  my $ratio_A = $base{$i}{A}/$total;
  my $ratio_C = $base{$i}{C}/$total;
  my $ratio_G = $base{$i}{G}/$total;
  my $ratio_T = $base{$i}{T}/$total;
  print "$i\t$ratio_A\t$ratio_C\t$ratio_G\t$ratio_T\n";
}

exit 0;
