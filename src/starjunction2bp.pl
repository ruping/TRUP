use strict;

#my $breakpoints = shift;
my $lastbpid = shift;
my $starjunction = shift;
my $maxIntron = shift;

if ($maxIntron eq ''){
  print STDERR "error: the maxIntron for starjunction2bp.pl is not set!!!!\n";
  exit 22;
}

#open IN, "$breakpoints";
#my %breakpoints;
#while ( <IN> ){
#  chomp;
#  my @cols = split /\t/;
#  $breakpoints{$cols[2]} = '';
#}
#close IN;

open IN, "$starjunction";
my %starjunction;
while ( <IN> ) {
  chomp;
  my ($chr1, $pos1, $strand1, $chr2, $pos2, $strand2, $junctionType, $repLeft, $repRight, $read, $loc1, $cigar1, $loc2, $cigar2) = split /\t/;

  #my $pow = 'w';

  if (($chr1 eq $chr2) and (abs($pos1-$pos2) < $maxIntron)) {
     #$pow = 'p';
     next;
  }

  #16202   p       chr1    138804  3       4       0       N       20
  #id     type     chr      pos   support  pw     ct      rep     discordant

  #next if (exists($breakpoints{$read}));

  #while ( $cigar1 =~ /(\d+)S/g ) {
  #  if ( $1 > 20 ) {
  #     printf("%s\n", join("\t", $chr1, $pos1, $read, $pow, 'S'));
  #     last;
  #  }
  #}
  #while ( $cigar2 =~ /(\d+)S/g ) {
  #  if ( $1 > 20 ) {
  #     printf("%s\n", join("\t", $chr2, $pos2, $read, $pow, 'S'));
  #     last;
  #  }
  #}

  $starjunction{$chr1."\t".$pos1}{$chr2."\t".$pos2} += 1;

}
close IN;

my $id = $lastbpid+1;
foreach my $bp1 (sort {$a =~ /^(chr\w+)\t(\d+)$/; my $ca=$1; my $pa=$2; $b =~ /^(chr\w+)\t(\d+)$/; my $cb=$1; my $pb=$2; ($ca cmp $cb) || ($pa <=> $pb)} keys %starjunction) {
  foreach my $bp2 (sort {$a =~ /^(chr\w+)\t(\d+)$/; my $ca=$1; my $pa=$2; $b =~ /^(chr\w+)\t(\d+)$/; my $cb=$1; my $pb=$2; ($ca cmp $cb) || ($pa <=> $pb)} keys %{$starjunction{$bp1}}) {
    my $support = $starjunction{$bp1}{$bp2};
    printf("%s\n", join("\t", $id, 'p', $bp1, $support, $support, 0, 'A', $support));
    printf("%s\n", join("\t", $id, 'p', $bp2, $support, $support, 0, 'A', $support));
    $id += 1;
  }
}
