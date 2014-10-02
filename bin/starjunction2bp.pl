use strict;

my $breakpoints = shift;
#my $lastbpid = shift;
my $starjunction = shift;

open IN, "$breakpoints";
my %breakpoints;
while ( <IN> ){
  chomp;
  my @cols = split /\t/;
  $breakpoints{$cols[2]} = '';
}
close IN;

open IN, "$starjunction";
while ( <IN> ) {
  chomp;
  my ($chr1, $pos1, $strand1, $chr2, $pos2, $strand2, $junctionType, $repLeft, $repRight, $read, $loc1, $cigar1, $loc2, $cigar2) = split /\t/;
  my $pow = 'w';
  if (($chr1 eq $chr2) and (abs($pos1-$pos2) < 230000)) {
     $pow = 'p';
  }
  next if (exists($breakpoints{$read}));
  while ( $cigar1 =~ /(\d+)S/g ) {
    if ( $1 > 20 ) {
       printf("%s\n", join("\t", $chr1, $pos1, $read, $pow, 'S'));
       last;
    }
  }
  while ( $cigar2 =~ /(\d+)S/g ) {
    if ( $1 > 20 ) {
       printf("%s\n", join("\t", $chr2, $pos2, $read, $pow, 'S'));
       last;
    }
  }
}
close IN;
