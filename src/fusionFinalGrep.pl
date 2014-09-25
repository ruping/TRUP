use strict;

my $filteredFusions = shift;
my $vis = shift;

open FF, "$filteredFusions";
my %fusion;
while ( <FF> ){
  chomp;
  my ($identifier, $idx, $length, $bp, $blat1, $blat2, $score) = split /\t/;
  $identifier =~ /(.+?)\|[12]\|\d+\-\d+\([+-]\)(chr\w+)\:(\d+)\-(\d+)\|\d+\-\d+\([+-]\)(chr\w+)\:(\d+)\-(\d+)$/;
  my $filtered = $1;
  my $blat1Chr = $2;
  my $blat1Start = $3;
  my $blat1End = $4;
  my $blat2Chr = $5;
  my $blat2Start = $6;
  my $blat2End = $7;
  $blat1 =~ /\d+\-\d+\([+-]\)(chr\w+)\:(\d+)\-(\d+)/;
  next unless ($1 eq $blat1Chr and ($2 <= $blat1End and $3 >= $blat1Start));
  $blat2 =~ /\d+\-\d+\([+-]\)(chr\w+)\:(\d+)\-(\d+)/;
  next unless ($1 eq $blat2Chr and ($2 <= $blat2End and $3 >= $blat2Start));
  $fusion{$filtered} = '';
}
close FF;

open VIS, "$vis";
while ( <VIS> ){
   chomp;
   my @cols = split /\t/;
   if ($cols[0] =~ /^#(.+)$/){
     my $candidate = $1;
     if (exists ($fusion{$candidate})){
        print "$_\n";
     }
   } else {
     next;
   }
}
close VIS;
