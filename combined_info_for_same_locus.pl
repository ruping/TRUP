#!/usr/bin/perl
use strict;

open IN, "$ARGV[0]";
my %combined;
while ( <IN> ){
    chomp;
    my ($transcript, $length, $fusion, $breakpoint, $coverage) = split /\t/;
    $transcript =~ /^#Locus_(\d+)_Transcript/;
    my $locus = $1;
    $_ =~ s/#//;
    my %hash_tmp = ('info'=>$_, 'cov'=>$coverage); 
    push (@{$combined{$locus}}, \%hash_tmp);
}

close IN;


foreach my $locus (sort { $combined{$b}[0]->{cov} <=> $combined{$a}[0]->{cov} } keys %combined){
    print "\@Locus\_$locus\n";
    foreach my $record (@{$combined{$locus}}){
       print "$record->{info}\n";
    }
}

exit;
