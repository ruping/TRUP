use strict;
use Data::Dumper;
use File::Glob ':glob';

my $dir = shift;
my $toptags = shift;
my $ifDE = shift;
my $spaired = shift;


my %IFDE;

open IFDE, "$dir/$ifDE";
while ( <IFDE> ){
  chomp;
  next if /^logFC/;
  if ($spaired == 1){
    my ($ens, $logFC, $logCPM, $LR, $pvalue, $ifDE) = split /\t/;
    $IFDE{$ens} = $ifDE;
  } else {
    my ($ens, $logFC, $logCPM, $pvalue, $ifDE) = split /\t/;
    $IFDE{$ens} = $ifDE;
  }
}
close IFDE;

=pod
my %anno;
open ANNO, "$ens_annotation";
while ( <ANNO> ){
  chomp;
  my @cols = split /\t/;
  $cols[$#cols] =~ /ID\=(.+?)\;Name\=(.+?)\;/;
  my $ensid = $1;
  my $ensname = $2;
  $anno{$ensid} = $ensname;
}
close ANNO;

my %entrz;
open ENTRZ, "$ens_entrz";
while ( <ENTRZ> ){
  chomp;
  my @cols = split /\t/;
  my $ensid = $cols[0];
  my $entrz = $cols[1];
  my $desc = $cols[2];
  $entrz{$ensid} = $entrz."\t".$desc;
}
close ENTRZ;
=cut

if ($spaired == 1) {
  printf("%s\n", join("\t", "ens-id", "logFC", "logCPM", "LR", "P-value", "FDR", "ifDE"));
} else {
  printf("%s\n", join("\t", "ens-id", "logFC", "logCPM", "P-value", "FDR", "ifDE"));
}

open TOP, "$dir/$toptags";
while ( <TOP> ) {
  chomp;
  next if /^logFC/;
  if ($spaired == 1){
    my ($ens, $logFC, $logCPM, $LR, $PValue, $FDR) = split /\t/;
    printf("%s\n", join("\t", $ens, $logFC, $logCPM, $LR, $PValue, $FDR, $IFDE{$ens}));
  } else {
    my ($ens, $logFC, $logCPM, $PValue, $FDR) = split /\t/;
    printf("%s\n", join("\t", $ens, $logFC, $logCPM, $PValue, $FDR, $IFDE{$ens}));
  }
}
close TOP;
