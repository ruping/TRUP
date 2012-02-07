#!/usr/bin/perl
use strict;

use Getopt::Long;

my $format;
my $trim;
my $read_length;

GetOptions (
              "format|f=s" => \$format,
	      "trim|t=i"   => \$trim,
	      "read|r=i"   => \$read_length,
              "help|h"     => sub {print "usage: perl trim_reads.pl -f=[format:fastq|fasta] -t=[length after trimming] -r=[read length] <input>\n"; exit 22;}
           );
	   
my $trim_end = $read_length - $trim;
my $q_line = 1;  #line indicator for fastq file; 
  
while ( <> ){
   chomp;
   #fastq
   if ($format eq 'fastq'){
     if ($q_line == 1){
       print "$_\n";
       $q_line++;
       next;
     }
     if ($q_line == 2){
       my $read = substr($_, 0, $trim);
       print "$read\n";
       $q_line++;
       next;
     }
     if ($q_line == 3){
       print "$_\n";
       $q_line++;
       next;
     }
     if ($q_line == 4){
       my $qual = substr($_, 0, $trim);
       print "$qual\n";
       $q_line = 1;
       next;
     }
   }
   
   #fasta   
   if ($format eq 'fasta'){
      if ($_ =~ /^>/){
        my $keep_len = length($_) - $trim_end;
	my $keep = substr($_, 0, $keep_len);
        print "$keep\n";
      }
      else {
        my $read = substr($_, 0, $trim);
	print "$read\n";
      }
   }
}

exit;
