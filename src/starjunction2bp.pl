use strict;

#my $breakpoints = shift;
my $lastbpid = shift;
my $starjunction = shift;
my $maxIntron = shift;
my $repeatmasker = shift;

if ($maxIntron eq ''){
  print STDERR "error: the maxIntron for starjunction2bp.pl is not set!!!!\n";
  exit 22;
}

#repeat mask
my %repeatmask;
my @rs;       #start array for each chr
my $old_chr;  #checking the chr
my $ptr;      #pointer for repeatmask sub
open REPEAT, "$repeatmasker";
while ( <REPEAT> ) {
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    my $chr = $tmp[0];
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
  }
close REPEAT;


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

  if (($chr1 eq $chr2) and (abs($pos1-$pos2) <= $maxIntron)) {
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

  $starjunction{$chr1}{$pos1} += 1;
  $starjunction{$chr2}{$pos2} += 1;

}
close IN;

my $id = $lastbpid+1;
foreach my $chr (sort {$a cmp $b} keys %starjunction) {

  my $last_coor = 0;
  my %starbp;

  foreach my $pos (sort {$a <=> $b} keys %{$starjunction{$chr}}) {

    my $support = $starjunction{$chr}{$pos};
    my $distance = abs($last_coor - $pos);
    if ($distance > 200) {
      if ($last_coor != 0) {  #print out
         my $supall = $starbp{$last_coor};
         my $rmflag = &repeatmask($chr, $last_coor);
         my $rep = ($rmflag == 1)? 'R':'N';
         printf("%s\n", join("\t", $id, 's', $chr, $last_coor, $supall, $supall, 0, $rep, $supall));
         delete $starbp{$last_coor};
         $id++;
      }
      $starbp{$pos} = $support;
      $last_coor = $pos;
    } else {
      $starbp{$last_coor} += $support;
    }

  }

}

exit 0;


sub repeatmask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr){
        @rs = sort {$a <=> $b} keys %{$repeatmask{$chr}};
        $ptr = 0;
      }
    while (($ptr<=$#rs) and ($repeatmask{$chr}{$rs[$ptr]} < $coor)){
        $ptr++;
      }
    if ($rs[$ptr] <= $coor){
        $flag = 1;
      }
    $old_chr = $chr;
    return $flag;
}
