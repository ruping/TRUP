use strict;

#automatically decide the problem of the fastq identifier
#1) @fastqname/1(N)?
#2) @fastqname/(4)?
#more to ba added in the future...

my $goodFlag = 0;

my $fastq = shift;

my $decompress = "";
if ($fastq =~ /\.gz$/){
  $decompress = "gzip -dc"
} elsif ($fastq =~ /\.bz2$/){
  $decompress = "bzip2 -dc";
}

open IN, "$decompress $fastq|";

while ( <IN> ){

  chomp;

  if ($_ =~ /^(\@.+?\/[123])$/) {
    $goodFlag = 1;   #it is a good identifier
    last;
  }

  if ($_ =~ /^\@.+?\s+[\d]/) {
    $goodFlag = 1;   #it is a good identifier
    last;
  }


  if ($_ =~ /^(\@.+?\/[1234])[\S]?$/) {
    my $identifier = $1;          #problem 1)
    $identifier =~ s/4$/3/;       #problem 2)
    print "$identifier\n";
  }

  $_ = <IN>;
  print "$_";

  $_ = <IN>;
  chomp;
  if ($_ =~ /^(\+.+?\/[1234])[\S]?$/){
    my $identifier = $1;          #problem 1)
    $identifier =~ s/4$/3/;       #problem 2)
    print "$identifier\n";
  } else {
    print "$_\n";
  }

  $_ = <IN>;
  print "$_";

}

close IN;
