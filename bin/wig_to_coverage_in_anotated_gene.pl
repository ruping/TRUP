use strict;
use Data::Dumper;

my ($wig, $annotation, $path, $lanename) = @ARGV;


my %chr;
for (1..22){
  my $c = 'chr'.$_;
  $chr{$c} = '';
}
$chr{'chrX'}  = '';
$chr{'chrY'}  = '';
$chr{'chrMT'} = '';


#rep mask
my %covmask;
my @covs;      #start array for each chr
my $old_chr; #checking the chr
my $ptr;     #pointer for repeatmask sub
my $old_run;

open WIG, "$wig";
  while ( <WIG> ){
    next if ($_ !~ /^chr/);
    chomp;
    my ($chr, $start, $end, $coverage) = split /\t/;
    $chr = 'chr'.$chr if ($chr !~ /^chr/ );
    $chr = 'chrMT' if ($chr =~ /chrM$/);

    $start += 1;
    if ($covmask{$chr} eq ''){
      #print STDERR "wig $chr\n";
    }
    $covmask{$chr}{$start}{'END'} = $end;
    $covmask{$chr}{$start}{'COV'} = $coverage;
  }
close WIG;
print STDERR "wig loaded\n";


open ANNO, "$annotation";
my %transcriptome;
while ( <ANNO> ) {
  next if /^#/;
  chomp;
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $tag) = split /\t/;
  $chr = 'chr'.$chr if ($chr !~ /^chr/ );
  next if (! exists($chr{$chr}));

  if ($type eq 'CDS' or $type eq 'exon' or $type eq 'three_prime_UTR' or $type eq 'five_prime_UTR') {

    #if ($source =~ /RNA/) {

      if ($transcriptome{$chr} eq '') {
        #print STDERR "transcriptome $chr\n";
      }

      if ($transcriptome{$chr}{$start} eq '') {
        $transcriptome{$chr}{$start} = $end;
      } elsif ($transcriptome{$chr}{$start} ne '') {
        my $old_end = $transcriptome{$chr}{$start};
        if ($old_end < $end) {
          $transcriptome{$chr}{$start} = $end;
        }
      }

    #}

  }
}
close ANNO;
print STDERR "Transcriptome loaded\n";


#sort the start of transcriptome and merge overlapping region into big non-overlapping region
my %tr_region;
foreach my $chr (sort keys %transcriptome){
  my $last_start = 0;
  my $last_end = 0;
  foreach my $start (sort {$a <=> $b} keys %{$transcriptome{$chr}}){
     my $end = $transcriptome{$chr}{$start};
     my $start_now = $start;
     my $end_now = $end;

     if ($start > $last_end){
        $tr_region{$chr}{$start} = $end;
     }

     elsif ($start <= $last_end){
        if ($end > $last_end){
          $tr_region{$chr}{$last_start} = $end;
          $start_now = $last_start;
        }
        elsif ($end <= $last_end) {
          $start_now = $last_start;
          $end_now = $last_end;
        }
     }

     $last_start = $start_now;
     $last_end = $end_now;

  } #start
} #$chr
print STDERR "non-overlapping region generated\n";

%transcriptome = ();

my $total_bases;
my $total_covered;
my %intr;
open XCOV, ">$path/$lanename\.coverage_in_transcriptome\.perbase";
print XCOV "#coverage in transcriptome per base\n";
foreach my $chr (sort keys %tr_region){
  #print STDERR "$chr coverage calculating\n";
  foreach my $start (sort {$a <=> $b} keys %{$tr_region{$chr}}){
    my $end = $tr_region{$chr}{$start};
    for ($start..$end){
      my $cov = covmask(1, $chr, $_);
      $total_bases++;
      $total_covered += $cov;
      $intr{$chr}{bases}++;
      $intr{$chr}{covered} += $cov;
      print XCOV "$chr\t$_\t$cov\n"; 
    }
  }
}
close XCOV;

my $ratio = sprintf ("%.1f", $total_covered/$total_bases);
print "ALL\t$total_bases\t$total_covered\t$ratio\n";

foreach my $chr (sort keys %intr){
  my $bases = $intr{$chr}{bases};
  my $covered = $intr{$chr}{covered};
  my $ratio = sprintf ("%.1f", $covered/$bases);
  print "$chr\t$bases\t$covered\t$ratio\n";
}


exit 22;

#subroutine Region
sub covmask {
    my ($run, $chr, $coor) = @_;
    my $cov = 0;
    if (($chr ne $old_chr) or ($run ne $old_run)){
       @covs = sort {$a <=> $b} keys %{$covmask{$chr}};
       $ptr = 0;
    }
    while (($ptr<=$#covs) and ($covmask{$chr}{$covs[$ptr]}{END} < $coor)){
      $ptr++;
    }
    if ($covs[$ptr] <= $coor){
      $cov = $covmask{$chr}{$covs[$ptr]}{COV};
    }
    $old_chr = $chr;
    $old_run = $run;
    return $cov;
}


