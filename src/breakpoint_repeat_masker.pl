use strict;

my $breakpoints = shift;
my $repeatmasker = shift;

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

open IN, "$breakpoints";

while ( <IN> ){
   chomp;
   my ($id, $type, $chr, $coor, $support, $pw, $ct) = split /\t/;
   my $rmflag = repeatmask($chr, $coor);
   if ($rmflag == 1) {
      print "$_\tR\n";
   }
   else {
      print "$_\tN\n";
   }
}
close IN;


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
