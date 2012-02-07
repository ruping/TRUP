#!/usr/bin/perl

my ($bowtie, $readlen) = @ARGV;

open IN, "$bowtie";
my $last_frag;
my $last_start;
my %fragl;
while ( <IN> ){
   chomp;
   my @cols = split /\t/;
   my $read = $cols[0];
   my $strand = $cols[1];
   my $chr = $cols[2];
   my $start = $cols[3];
   my $frag;
   if ( $read !~ /\s+/ ){
     $read =~ /^(.+?)\//;
     $frag = $1;
   }
   else {
     $read =~ /^(.+?)\s+/;
     $frag = $1;
   }
   if ($frag eq $last_frag){
     $fragl{$frag} = $start - $last_start + $readlen;
   }
   $last_frag = $frag;
   $last_start = $start;
}
close IN;

foreach my $frag (keys %fragl){
    print "$frag\t$fragl{$frag}\n";
}

exit;
