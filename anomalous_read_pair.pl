use strict;
use Getopt::Long;

my $mapping_file;
my $AB;
my $read_file_1;
my $read_file_2;

GetOptions (
              "mappingfile|m=s"  => \$mapping_file,
	      "readfile1|r1=s"   => \$read_file_1,
              "readfile2|r2=s"   => \$read_file_2,
              "AB"               => \$AB,
	      "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--type\t\tthe type of the mapping, either genome or refseq\n";
                             print "\t--mappingfile\tthe mapping file of genome mapping\n";
			     print "\t--readfile1\tthe 5' end reads\n";
			     print "\t--readfile2\tthe 3' end reads\n";
			     print "\t--help\t\tprint this help message\n";
                             exit 0;
			    }
	   );



my %ARP;  #hash to remember the anomalous read pairs
my %UMR;  #hash to remember the mapping information

open MAPPING, "$mapping_file";
my $old_chr;
my %chr_start;
while ( <MAPPING> ) {

      next if /^@/;
      chomp;

      my @cols = split /\t/;
      my $chr = $cols[2];

      $chr = 'chr'.$chr unless ($chr =~ /^chr/);
      $chr .= 'T' if ($chr eq 'chrM');

      if ($chr ne $old_chr){
         $chr_start{$chr} = tell MAPPING;
      }

      $old_chr = $chr;

}
close MAPPING;

#select anomalous read pairs of both mapped pairs

foreach my $chr (sort keys %chr_start){

  load_mapping($chr);
  print STDERR "$chr UMR loading for searching ARP\n";

  foreach my $read (keys %UMR){

    #there is only one mapping record for this UMR, meaning that the other end mapped to another chr
    if (scalar(@{$UMR{$read}}) == 1){
       $ARP{$read} = '';
       next;
    }

    my $pos1 = $UMR{$read}[0]->{'pos'};
    my $pos2 = $UMR{$read}[1]->{'pos'};
    my $dis = abs($pos1-$pos2);
    if ($dis >= 100000){
      $ARP{$read} = '';
    }
  }
}

open R1, "$read_file_1";
(my $a_R1 = $read_file_1) =~ s/fq$/ARP\.fq/;
open AR1, ">$a_R1";
while ( <R1> ){
  if ($_ =~ /^@(.+?)[\/\s]/){
     if (exists $ARP{$1}){
        print AR1 "$_";
        $_ = <R1>;
        print AR1 "$_";
        $_ = <R1>;
        print AR1 "$_";
        $_ = <R1>;
        print AR1 "$_";
     }
     else {
        $_ = <R1>;
        $_ = <R1>;
        $_ = <R1>;
     }
  }
}
close R1;
close AR1;

open R2, "$read_file_2";
(my $a_R2 = $read_file_2) =~ s/fq$/ARP\.fq/;
open AR2, ">$a_R2";
while ( <R2> ){
  if ($_ =~ /^@(.+?)[\/\s]/){
     if (exists $ARP{$1}){
        print AR2 "$_";
        $_ = <R2>;
        print AR2 "$_";
        $_ = <R2>;
        print AR2 "$_";
        $_ = <R2>;
        print AR2 "$_";
     }
     else {
        $_ = <R2>;
        $_ = <R2>;
        $_ = <R2>;
     }
  }
}
close R2;
close AR2;

#subregion############################################################################################################

sub load_mapping {

    my $chr_w = shift;

    %UMR = ();

    open MAPPING, "<$mapping_file";
    seek(MAPPING, $chr_start{$chr_w}, 0);
    while( <MAPPING> ){

       next if ($_ !~ /\tNH:i:1$/);   #filter out multiple mapping

       my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;
       $Qname =~ s/[AB]$// if ($AB);
       my $chr = $Rname;
       $chr = 'chr'.$chr unless ($chr =~ /^chr/);
       $chr .= 'T' if ($chr eq 'chrM');

       last if ($chr ne $chr_w);

       my $code = decode_FLAG($FLAG);

       if ($code =~ /U/){
         $ARP{$Qname} = '';
         next;
       }

       my %hash_tmp = ('chr'=>$chr, 'pos'=>$Pos);
       push (@{$UMR{$Qname}}, \%hash_tmp);
    }
    close MAPPING;
}

sub decode_FLAG {
   my $compound = $_[0];
   my $code;
   if ($compound & 1)  {$code .= 'p';} #the read is paired in sequencing
   if ($compound & 2)  {$code .= 'P';} #the read is mapped in a proper pair
   if ($compound & 4)  {$code .= 'u';} #the read is unmapped
   if ($compound & 8)  {$code .= 'U';} #the mate is unmapped

   if ($compound & 16) {$code .= '-';} #the read strand -
   else                {$code .= '+';} #the read strand +
   if ($compound & 32) {$code .= '-';} #the mate strand -
   else                {$code .= '+';} #the mate strand +

   if ($compound & 64) {$code .= '1';} #the read is the first in the pair
   if ($compound & 128){$code .= '2';} #the read is the second in the pair

   if ($compound & 256){$code .= 'N';} #not primary record

   return $code;
}
