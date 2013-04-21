use strict;

my $dismate = shift;
my $repeatmasker = shift;
my $breakpoints = shift;

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

#breakpoints_mask
my %breakmask;
my @rs2;       #start array for each chr
my $old_chr2;  #checking the chr
my $ptr2;      #pointer for repeatmask sub
open BPM, "$breakpoints";
while ( <BPM> ) {
    chomp;
    my @tmp = split /\t/;
    my $chr = $tmp[2];
    my $bpm_s = $tmp[3] - 250;
    my $bpm_e = $tmp[3] + 250;
    $breakmask{$chr}{$bpm_s} = $bpm_e;
}
close BPM;


open IN, "$dismate";
while ( <IN> ){
   chomp;
   my ($id, $type, $chr, $coor, $support, $pw, $ct, $rep, $dc) = split /\t/;
   my $rmflag = &repeatmask($chr, $coor);
   my $bpmflag = &breakmask($chr, $coor);

   if ($rmflag == 1) {
      $rep = 'R';
   } else {
      $rep = 'N';
   }

   if ($bpmflag != 1 and $rep eq 'N'){
     printf "%s\n", join("\t", $id, $type, $chr, $coor, $support, $pw, $ct, $rep, $dc);
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

sub breakmask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr2){
        @rs2 = sort {$a <=> $b} keys %{$breakmask{$chr}};
        if ($#rs2 < 0) {
           return $flag;
        }
        $ptr2 = 0;
    }
    while (($ptr2<=$#rs2) and ($breakmask{$chr}{$rs2[$ptr2]} < $coor)){
        $ptr2++;
    }
    if ($ptr2 <= $#rs2 and $rs2[$ptr2] <= $coor){
        $flag = 1;
    }
    $old_chr2 = $chr;
    return $flag;
}
