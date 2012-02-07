use strict;
use Data::Dumper;

my ($task, $file) = @ARGV;


if ($task  =~ /refseq/){
open REFSEQ, "$file";
while ( <REFSEQ> ){
    chomp;
    if (/^>(gi\|.+?\|ref\|(.+?)\|).+\((.+?)\).+?$/){
        my $id = $1;
        my $ref = $2;
        my $name = $3;
        #if ($name =~ /^LOC/){
          print "$name\n";
        #}
    }
}
close REFSEQ;
}

if ($task  =~ /ensemble/){
  my %ensemble;
  open ENS, "$file";
  while ( <ENS> ){
      chomp;
      next if /^Ensembl/;
      my ($ensemble_id, $refseq_id, $gene_name) = split /\t/;
      $ensemble{$gene_name} = $ensemble_id;
  }
  close ENS;
  print Dumper (\%ensemble);
  
}