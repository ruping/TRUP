use strict;

my $vis = shift;

my $readlength = shift;

open VIS, "$vis";
my $title = 'SRP';
while ( <VIS> ){
  chomp;
  if ($_ =~ /^#/){    #fusion title
    $title = $_;
  } else {
    if ($title ne 'SRP'){  #do things here
      my ($candidate, $idx, $cScore, $length, $fusion_genes, $ori, $breakpoint, $rep, $sc, $type, $ron, $blat1, $blat2, $cov1, $cov2, $cov3, $spanAll, $spanScore) = split(/\t/,$title);
      $candidate =~ s/^#//;
      my $identifier = $candidate.'|'.$idx.'|'.$blat1.'|'.$blat2;
      $_ =~ /[ACGTN\-]+([acgtn]+)[ACGTN\-]+/;
      my $bps = $-[1];
      my $bpe = $+[1];
      my $bpl = $bpe - $bps;
      my $restl = $readlength - $bpl;
      my $needl = ($restl % 2 == 0)? ($restl/2):(($restl/2)+0.5);
      my $needr = $needl;
      my $pos;
      if ($bps < $needl) {
        $needl = $bps;
        $needr += $needl - $bps;
        $pos = 0;
      } elsif ((length($_) - $bpe + 1) < $needr) {
        $needr = length($_) - $bpe + 1;
        $needl += $needr - (length($_) - $bpe + 1);
        $pos = $bps - $needl;
      } else {
        $pos = $bps - $needl;
      }
      my $seq = substr($_, $pos, $readlength);
      $seq =~ s/[\-]+//g;
      print "\>$identifier\n$seq\n";
      $title = 'SRP'; #reset title
    } #the first line after title
  }
}
close VIS;
