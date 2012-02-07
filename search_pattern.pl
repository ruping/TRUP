#!/usr/bin/perl

use strict;
use Getopt::Long;

my $pattern;
GetOptions ("pattern|p=s" => \$pattern,);

my $title;
my %data;
while ( <> ){
   chomp;
   if (/^>/){
     $title = $_;
   }
   else{
     s/\n//g; 
     s/\s//g;
     $data{$title} .= $_;
   }
}

foreach my $title (keys %data){
  if ($data{$title} =~ /$pattern/){
    print "$title\n$data{$title}\n";
  }
}
exit;
