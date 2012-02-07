my @files = @ARGV;

my %fusion;
my $flag = 1;

foreach my $file (@files){
  open IN, "$file";
  while ( <IN> ){
    chomp;
    my ($candidate,$length,$fusion_genes,$breakpoint,$coverage) = split /\t/;
    $fusion{$fusion_genes}{file} .= $flag;
    $coverage =~ /^(.+?)\(/;
    $fusion{$fusion_genes}{cov} .= "\t$1";
  }
  close IN;
  $flag++;
}


$flag -= 1;

foreach my $fusion (keys %fusion){
  my $cov = $fusion{$fusion}{cov};
  my $indi = 'YES';
  for (1..$flag){
     if ($fusion{$fusion}{file} !~ /$_/){
        $indi = 'NO';
     }
  }
  print "$fusion$cov\n" if $indi eq 'YES';
}

exit;
