

my $mapfile  = shift;
my %mapping;

open MAP, "<$mapfile" || die "can't open $mapfile";

while ( <MAP> ) {
  chomp;
  my @cols = split /\t/;
  my $chr = $cols[2];
  my $start = $cols[3];
  $mapping{$chr}{$start}{$end}{$strand}++;
}
close MAP;

print STDERR "$mapfile loaded\n";


#now we got all the starts#######################################################

foreach my $chr (sort keys %mapping) {

  my $ptr = 0;
  my @starts = sort {$a <=> $b} keys %{$mapping{$chr}};

  my $total_starts = scalar @starts;
  print STDERR "there are total $total_starts starts on $chr\n";

  my $speed_count=0;
  foreach my $start (@starts) {

    my $winstart = $start;
    my $winend   = $start+$winsize;


    my $off = 0;
    while ((($ptr+$off) <= $#starts) and ($starts[$ptr+$off] <= $winend)) {

      my @ends = sort {$a <=> $b} keys %{$mapping{$chr}{$starts[$ptr+$off]}};

      foreach my $end (@ends) {
        foreach my $strand (keys %{$mapping{$chr}{$starts[$ptr+$off]}{$end}}) {
          if ($starts[$ptr+$off] >= $winstart) {
            my $dis1 = $starts[$ptr+$off]-$start;
            my $dis2 = $start-$starts[$ptr+$off];
            $ss{$dis1}++;
            $ss{$dis2}++;
          }
        }
      }                         #foreach

      $off++;

    }                           #while
    $ptr++;
    print STDERR "now processed $ptr\n" if ($ptr%10000==0);
  }

}                               #for each chr to get the edepth

foreach my $dis (sort {$a <=> $b} keys %ss) {
  print OUT1 "$dis\t$ss{$dis}\n";
}
close OUT1;
}


