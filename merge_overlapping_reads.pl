         #!/usr/bin/perl

######
#TODO: merge the overlapping pair-end reads into one read
###### 
      
use strict;
use Getopt::Long;

my $read_file1;
my $read_file2;
my $read_length;
my $anchor_size = 3;
my $mis_rate = 2;

GetOptions (
              "readfile1=s"      =>   \$read_file1,
              "readfile2=s"      =>   \$read_file2,
	      "readlength|r=i"   =>   \$read_length,
              "anchor=i"         =>   \$anchor_size,
              "mismatch=f"       =>   \$mis_rate,     
	      "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--readfile1\tthe first read file of a pair-end sample\n";
                             print "\t--readfile2\tthe second read file of a pair-end sample\n";
			     print "\t--readlength\tthe length of the read\n";
                             print "\t--anchor\tthe anchor length for first matching, default 3\n";
                             print "\t--mismatch\tthe mismatch rate threshold for the matching string\n";
			     print "\t--help\t\tprint this help message\n";
                             exit 0;
			     }
           );

open READF1, "< $read_file1";
open READF2, "< $read_file2";

while ( <READF1> ){
     my $title1 = $_;
     my $title2 = <READF2>;
     my $read1  = <READF1>;
     my $read2  = <READF2>;
     my $third1 = <READF1>;
     my $third2 = <READF2>;
     my $qual1  = <READF1>;
     my $qual2  = <READF2>;
     chomp ($title1,$title2,$read1,$read2,$third1,$third2,$qual1,$qual2);
     
     my $title;
     $title1 =~ /^@(.+?)\//;
     $title = $1;
     
     my $count_N1;
     my $count_N2;
     while ($read1 =~ /N/g){$count_N1++;}
     while ($read2 =~ /N/g){$count_N2++;}
     
     next if ($count_N1 >= 1 or $count_N2 >= 1);
     
     
     #reverse complementary the sequence
     $read2 = reverse($read2);
     $read2 =~ tr/ACGT/TGCA/;
     $qual2 = reverse($qual2);
     
     my $flag = 0;
     
     
     my $anchor = substr($read2, 0, $anchor_size);
     while ($read1 =~ /$anchor/g){
       my $pos = pos($read1);
       my $toend1 = substr($read1, $pos-$anchor_size);
       my $leng = length $toend1;
       my $toend2 = substr($read2, 0, $leng);     
       my $edis = levenshtein($toend1, $toend2);
       my $ratio = $edis/$leng;
       
       next if ($edis > 2 or $ratio > 0.05);
       
       my $add_string  = substr($read2, -($read_length-$leng));
       my $add_qual    = substr($qual2, -($read_length-$leng));
       my $merged_read = $read1.$add_string;
       my $merged_qual = $qual1.$add_qual;
       print "\@$title\+\n";
       print "$merged_read\n";
       print "\+$title\+\n";
       print "$merged_qual\n";       
       $flag = 1;
       last;
     }
     
     next if $flag == 1;
     $anchor = substr($read1, 0, $anchor_size);
     while ($read2 =~ /$anchor/g){
       my $pos = pos($read2);
       my $toend1 = substr($read2, $pos-$anchor_size);
       my $leng = length $toend1;
       my $toend2 = substr($read1, 0, $leng);     
       my $edis = levenshtein($toend1, $toend2);
       my $ratio = $edis/$leng;
       
       next if ($edis > 2 or $ratio > 0.05);
       
       my $add_string  = substr($read1, -($read_length-$leng));
       my $add_qual    = substr($qual1, -($read_length-$leng));
       my $merged_read = $read2.$add_string;
       my $merged_qual = $qual2.$add_qual;
       print "\@$title\-\n";
       print "$merged_read\n";
       print "\+$title\-\n";
       print "$merged_qual\n";       
       $flag = 2;
       last;
     }
     
     print STDERR "$title\n" if $flag == 0;
}

close READF1;
close READF2;

#edit distance of two equal length string
sub levenshtein($$){
  my @A=split //, lc shift; # lower case
  my @B=split //, lc shift;
  my @W=(0..@B);
  my ($i, $j, $cur, $next);
  for $i (0..$#A){ #
      $cur=$i+1;
      for $j (0..$#B){ #
          $next=min(
              $W[$j+1]+1,
              $cur+1,
              ($A[$i] ne $B[$j])+$W[$j]
          );
		$W[$j]=$cur;
		$cur=$next;
	}
	$W[@B]=$next;
  }
  return $next;
}

sub min($$$){
  if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
  return $_[0] < $_[1]? $_[0]:$_[1];
}

exit;