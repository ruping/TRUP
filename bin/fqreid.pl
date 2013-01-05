use strict;

my $fastq = shift;

open IN, "gzip -dc $fastq|";

while ( <IN> ){

  chomp;
  $_ =~ /^(\@.+?\/[12])\S+$/;
  my $identifier = $1;
  print "$identifier\n";

  $_ = <IN>;
  print "$_";

  $_ = <IN>;
  chomp;
  if ($_ =~ /^(\+.+?\/[12])\S+$/){
    $identifier = $1;
    print "$identifier\n";
  }
  else {
    print "$_\n";
  }

  $_ = <IN>;
  print "$_";

}

close IN;
