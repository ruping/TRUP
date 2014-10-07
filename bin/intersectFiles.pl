#  (c) 2014 - Sun Ruping
#  rs3412@columbia.edu
#
#  TODO: compare two files of coordinates, with a minimal consumption of memory and running-time

use strict;
use Data::Dumper;
use File::Basename;

my $maskfile;
my $original;
my $bgzip = 0;
my $gzip = 0;
my $snp = 0;
my $t=0;                                                               #tolerant
my $count = 0;
my $overlap = 0;
my $nonoverlap = 0;
my $column='';

my $michr = 0;
my $mistart = 1;
my $miend = 1;
my $oichr = 0;
my $oistart = 1;
my $oiend = 1;


my %printer;

while ($ARGV[0]) {
  my $arg = shift @ARGV;
  if ($arg eq '-m'){$maskfile    = shift @ARGV;}
  elsif ($arg eq '-t'){$t        = shift @ARGV;}
  elsif ($arg eq '-o'){$original = shift @ARGV;}
  elsif ($arg eq '-bgz'){$bgzip = 1;}
  elsif ($arg eq '-gz'){$gzip = 1;}
  elsif ($arg eq '-s'){$snp = 1;}
  elsif ($arg eq '-michr'){$michr = shift @ARGV;}
  elsif ($arg eq '-mistart'){$mistart = shift @ARGV;}
  elsif ($arg eq '-miend'){$miend = shift @ARGV;}
  elsif ($arg eq '-oichr'){$oichr = shift @ARGV;}
  elsif ($arg eq '-oistart'){$oistart = shift @ARGV;}
  elsif ($arg eq '-oiend'){$oiend = shift @ARGV;}
  elsif ($arg eq '-count'){$count = 1;}
  elsif ($arg eq '-nonoverlap'){$nonoverlap = 1;}
  elsif ($arg eq '-overlap'){$overlap = 1;}
  elsif ($arg eq '-column'){$column = shift @ARGV;}
  elsif ($arg eq '-h'){print "useage: intersectFiles.pl -o: origninal_filename -m: maskfile [-s:snpdb]  [-gz/-bgz] [-count/overlap/nonoverlap] [-column:get that column from mask file] [-t:tolerant]\n"; exit 0;}
  else {print "useage: intersectFiles.pl -o: origninal_filename -m: maskfile [-s:snpdb] [-gz/-bgz] [-t:tolerant]\n"; exit 0;}
}


#the original file is loaded in a hash
my $oopen;
if ($bgzip == 1) {
  $oopen = "bgzip -dc $original |";
} elsif ($gzip == 1) {
  $oopen = "gzip -dc $original |";
} else {
  $oopen = "$original";
}

my @columnIndex = ();
my @columnNames = ();
my %maskChrJumper = &getchrpos($maskfile);
print STDERR Dumper(\%maskChrJumper);
print STDERR "column set: $column\tindex:\n";
print STDERR Dumper(\@columnIndex);
print STDERR Dumper(\@columnNames);

my $line;
my $old_chr = "SRP";

#regions for the input of region file
my @variants;  #it is an Array like a deque

my $isComment = 0;
open ORI, "$original";
$line = <ORI>;
$isComment = eatline($line, \@variants);


while ($isComment == 1) {
   #print comment
   if ($count == 1) {
     chomp($line);
     my $maskadd = basename($maskfile);
     print "$line\t$maskadd\n";
   } elsif ($column ne '') {
     chomp($line);
     my $maskadd = basename($maskfile);
     if (scalar(@columnIndex) == 1){
       print "$line\t$maskadd\.$column\n";
     } else {
       chomp($columnNames[$#columnNames]);
       $maskadd = join("\t", @columnNames);
       print "$line\t$maskadd\n";
     }
   } else {
     print "$line";
   }
   $line = <ORI>;
   $isComment = eatline($line, \@variants);                            # eat a new line
}                                                                      # is comment

my $it = 0;


while ( $variants[$it]->{'chr'} ne $old_chr ) {

    $old_chr = $variants[$it]->{'chr'};                                # set the current chr as old chr

    my $chr_id  = $variants[$it]->{'chr'};

    if ( ! exists($maskChrJumper{$chr_id}) ) {                         # reference not found

      while ( $it <= $#variants and $variants[$it]->{'chr'} eq $old_chr ) {
        &var_processing($variants[$it]);                               # print the old region info
        splice(@variants, $it, 1);                                     # erase the current region ########################################################################
      }

      while ( $#variants == -1 ) {
        $line = <ORI>;
        if (!defined $line) {
          print STDERR "finished: end of region file, zone 0\n";
          exit 0;
        }
        eatline($line, \@variants);
        $it = 0;
        if ( $variants[$it]->{'chr'} eq $old_chr ) {
          &var_processing($variants[$it]);
          @variants = ();
          next;
        }
      }
      next;
    } #reference not found

    my $jumper = $maskChrJumper{$chr_id};
    print STDERR "$chr_id\tjumper: $jumper\t$old_chr\t$variants[$it]->{'chr'}\t$variants[$it]->{'start'}\t$variants[$it]->{'end'}\n";

    open MASK, "$maskfile";
    seek(MASK, $jumper, 0);
    while ( <MASK> ) {

      chomp;
      my @cols = split /\t/;
      my $chr_now  =  $cols[$michr];
      $chr_now = 'chr'.$chr_now unless ($chr_now =~ /^chr/);
      $chr_now = 'chrM' if ($chr_now eq 'chrMT');

      my $dbStart  =  $cols[$mistart];
      my $dbEnd    =  $cols[$miend];

      if ($chr_now ne $old_chr) {                                 # jump out if reaching the next chromosome
         last;
      }

      my $iter = 0;

      #skip wierd thing with empty chr and information##############################################################
      if ($variants[$iter]->{'chr'} eq '') {
         splice(@variants, $iter, 1);
         if ( $#variants == -1 ) {
             $line = <ORI>;
             if (!defined $line) {
                print STDERR "finished: end of region file, zone 2.2\n";
                 exit 0;
             }
             eatline($line, \@variants);                              # eat a line and put it into the duque
             $iter = 0;
         }
      }
      ##############################################################################################################


      if ( $variants[$iter]->{'start'} > $dbEnd ) {
         next;                                                         # skip reads not overlapping with the first region
      }


      while ( $variants[$iter]->{'chr'} eq $old_chr and $variants[$iter]->{'start'} <= $dbEnd and $iter <= $#variants ) {


        if ( $variants[$iter]->{'end'} < $dbStart ) {                  # the region end is beyond the alignmentStart

          &var_processing($variants[$iter]);                           # processing
          splice(@variants, $iter, 1);                                 # this region should be removed
          if ( $#variants == -1 ) {
              $line = <ORI>;
              if (!defined $line) {
                 print STDERR "finished: end of region file, zone 3\n";
                 exit 0;
              }
              eatline($line, \@variants);                              # eat a line and put it into the duque
              $iter = 0;
          }
          next;
        }


        if ( $variants[$iter]->{'end'} >= $dbStart and $variants[$iter]->{'start'} <= $dbEnd ) {  #overlapping, should take action

           push(@{$variants[$iter]->{'printer'}}, $_);

        }                                                              # overlapping take action!


        if ( $iter != $#variants ) {
          $iter++;                                                     # if this region is not the last element in the deque
        } else {                                                       # the last element
            $line = <ORI>;                                             # get a line of region file
            if (!defined $line){
               &var_processing($variants[$iter-1]);
               print STDERR "finished: end of region file, zone 4\n";
               exit 0;
            }
            eatline($line, \@variants);                                # eat a line and put it into the duque
            $iter = $#variants;
        }   #the last element

      }     #while

    }       # read a mask file line

    # somehow to loop back/
    $it = 0;                                                                    # reset to begin
    while ( $it <= $#variants and $variants[$it]->{'chr'} eq $old_chr ) {
      &var_processing($variants[$it]);                                          # print the old region info
      splice(@variants, $it, 1);                                                # erase the current region
    }

    while ( $#variants == -1  ) {

      $line = <ORI>;
      if (!defined $line) {
         print STDERR "finished: end of region file, zone 5\n";
         exit 0;
      }
      eatline($line, \@variants);
      $it = 0;
      if ( $variants[$it]->{'chr'} eq $old_chr ) {
        &var_processing($variants[$it]);
        @variants = ();
        next;
      }

    }

} # region chr != old chr



sub eatline {

  my $line = shift;
  my $variants = shift;

  my $isComment = 0;
  if ($line =~ /^#/ || $line =~ /^[cC][hH][rR]\t/) {
    $isComment = 1;
    return $isComment;
  }

  chomp($line);
  my @cols = split (/\t/, $line);

  my %variant;

  $variant{'chr'} = $cols[$oichr];
  $variant{'chr'} = 'chr'.$variant{'chr'} if ($variant{'chr'} !~ /^chr/);
  $variant{'chr'} = 'chrM' if ($variant{'chr'} eq 'chrMT');
  $variant{'start'} = $cols[$oistart] - $t;
  $variant{'end'} = $cols[$oiend] + $t;
  $variant{'info'} = $line;

  push(@{$variants}, \%variant);
  return $isComment;

}


sub getchrpos {

  my $dbfile = shift;

  my $chr_old = "UNDEF";

  my %chr_start;

  my $jumper = 0;

  open DBFILE, "$dbfile" or die "The db file read error!";

  while ( <DBFILE> ) {

      if ($_ =~ /^[\#\@]/ or $_ =~ /^[cC][hH][rR]\t/ or $_ =~ /^FID/) {
        $jumper = tell DBFILE;
        if ($column ne '' and ($_ =~ /^#/ or $_ =~ /^FID/)) {
          $_ =~ s/^#//;
          my @colnames = split(/\t/, $_);
          if ($column =~ /^(\d+)\-(\d+)$/) {
            foreach my $i ($1..$2){
              push(@columnIndex, $i);
            }
          } elsif ($column =~ /\,/) {
            my @colneed = split(",", $column);
            foreach my $i (@colneed) {
              push(@columnIndex, $i);
            }
          } else {
            for (my $i = 0; $i <= $#colnames; $i++) {
              if ($colnames[$i] eq $column) {
                push(@columnIndex, $i); #this is a global variable
              }
            }
          }

          foreach my $index (@columnIndex){
             push(@columnNames, $colnames[$index]);
          }
        } #find column index
        next;
      }
      chomp;

      my @cols = split /\t/;
      my $chr  = $cols[$michr];

      $chr = 'chr'.$chr unless ($chr =~ /^chr/);
      $chr = 'chrM' if ($chr eq 'chrMT');

      if ($chr ne $chr_old) {
        $chr_start{$chr} = $jumper;
      }

      $chr_old = $chr;
      $jumper = tell DBFILE;
  }

  close DBFILE;
  print STDERR "$dbfile chr_start loaded\n";

  return %chr_start;
}


sub var_processing {
   my $variant = shift;
   my $original = $variant->{'info'};
   return if $original eq '';
   if ($count == 0 and $overlap == 0 and $nonoverlap == 0 and $column eq '') {
     print "###\t$original\n";
     foreach my $printer ( @{$variant->{'printer'}} ) {
        print "$printer\n";
     }
   } elsif ($count == 1) {
     my $nump = 0;
     foreach my $printer ( @{$variant->{'printer'}} ) {
        $nump++;
     }
     print "$original\t$nump\n" if ($count == 1 and $nonoverlap == 0 and $overlap == 0);
     print "$original\t$nump\n" if ($count == 1 and $nonoverlap == 1 and $nump == 0);
     print "$original\t$nump\n" if ($count == 1 and $overlap == 1 and $nump > 0);
   } elsif ($column ne ''){
     print "$original\t";
     my $columnAdd;
     foreach my $printer ( @{$variant->{'printer'}} ) {
       my @tmp = split (/\t/, $printer);
       my @columnAdd;
       foreach my $index (@columnIndex){
         push(@columnAdd, $tmp[$index]);
       }
       if (scalar(@columnAdd) == 1){
         $columnAdd = $columnAdd[0];
       } else {
         $columnAdd = join("\t", @columnAdd);
       }
     }
     if ($columnAdd eq ''){
       $columnAdd = "NA";
     }
     print "$columnAdd\n";
   }
}
