#!/usr/bin/perl 
use strict;
use Getopt::Long;

my $type;
GetOptions ("type|y=s" => \$type);



if ($type eq 'qa'){
my $line = 1;
my $id = undef;
my $seq = undef;


while ( <> ){
  chomp;
  if ($line == 1){
    $_ =~ /^@(.*)$/;
    $id = $1;
  }
  if ($line == 2){
    $seq = $_;
  }
  if ($line == 3){
    #doing nothing
  }
  if ($line == 4){
    #my $av = 0;
    #my $quality = $_;
    #foreach my $q (split //,$quality){
    #  my $qual = ord($q)-64;
    #  $av += $qual;
    #}
    #$av = $av/95; 
    #printf ">$id\|%.1f\|%s\n$seq\n", $av, $quality;
    print ">$id\n$seq\n";
    $line = 0;
  }
  $line++;
}
}

if ($type eq 'aq'){
my $line1;
my $seq;
my $line3;
my $qual;
my $id;
while ( <> ){
  chomp;
  if ($_ =~ /^>(.+)\|(.+)\|(.+)$/){
    $id = $1;
    $qual = $3;
    $line1 = '@'.$id;
    $line3 = '+'.$id;
    print "$line1\n";
  }
  else{
    $seq = $_;
    print "$seq\n";
    print "$line3\n";
    print "$qual\n";
  }
}
}

exit;
