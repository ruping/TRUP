use strict;
use Getopt::Long;

my $g1;
my $g2;
my $task = 1;
my $read_file_1 = '';
my $read_file_2 = '';

GetOptions(
             "g1=s"   => \$g1,
             "g2=s"   => \$g2,
             "task=i" => \$task,
             "r1=s"   => \$read_file_1,
             "r2=s"   => \$read_file_2,
             "help|h" => sub{
	                     print "\nusage: $0 [options] samfile\n\nOptions:\n\t--g1\tthe position of fusion gene1, like chrX:1-100\n";
                             print "\t--g2\tthe position of fusion gene2, like chrX:1-100\n";
			     print "\t--help\tprint this help message\n";
                             exit 0;
			     },
          );


$g1 =~ /^(chr\w+)\:(.+?)\-(.+)$/;
my $chr1 = $1;
my $gene1start = $2; 
my $gene1end = $3;
$gene1start =~ s/,//g;
$gene1end =~ s/,//g;

$g2 =~ /^(chr\w+)\:(.+?)\-(.+)$/;
my $chr2 = $1;
my $gene2start = $2; 
my $gene2end = $3;
$gene2start =~ s/,//g;
$gene2end =~ s/,//g;



if ($task ==1){


print "#$chr1\:$gene1start\.\.$gene1end\n";
print "#$chr2\:$gene2start\.\.$gene2end\n";
  
  my %buffer1;
  my %buffer2;
 
  while ( <> ){
     
     next if /^@/;                #ignore comments
     next if ($_ !~ /NH\:i\:1$/); #ignore multiple mappable reads
     chomp;
     my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;

     next unless ($mateRname eq '=');
     next if ($FLAG & 2);         #ignore correct pair
     
     if ($Rname eq $chr1){
        if (($Pos >= $gene1start and $Pos <= $gene1end) and ($matePos >= $gene2start and $matePos <= $gene2end)){
            $buffer1{$Qname} = $_;
        }
     }
     
     if ($Rname eq $chr2){
        if (($Pos >= $gene2start and $Pos <= $gene2end) and ($matePos >= $gene1start and $matePos <= $gene1end)){
            $buffer2{$Qname} = $_;
        }
     }
  }

  foreach my $qname (keys %buffer1){
     if (exists $buffer2{$qname}){
        print "$buffer1{$qname}\n";
        print "$buffer2{$qname}\n";
     }
  }
}

if ($task == 2){
   
   my %buffer; 
   
   while ( <> ){
      next if /^@/;                #ignore comments
      next if ($_ !~ /NH\:i\:1$/); #ignore multiple mappable reads
      chomp;
      my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;
      
      if ($Rname eq $chr1 and $Pos >= $gene1start and $Pos <= $gene1end){
         $buffer{$Qname} = '';
      }
      
      if ($Rname eq $chr2 and $Pos >= $gene2start and $Pos <= $gene2end){
         $buffer{$Qname} = '';
      }
   }
   
   open R1, "$read_file_1";
   my $flag1 = 0;
   while ( <R1> ){
     if ($_ =~ /^@(.+?)\//){
        if (exists $buffer{$1}){
           print "$_";
           $flag1 = 1;
        }
        else {
           $flag1 = 0;
        }
     }  
     elsif($flag1 == 1){
       print "$_";
     }
   }
   close R1;
   
   
   open R2, "$read_file_2";
   my $flag2 = 0;
   while ( <R2> ){
     if ($_ =~ /^@(.+?)\//){
        if (exists $buffer{$1}){
           print STDERR "$_";
           $flag2 = 1;
        }
        else {
           $flag2 = 0;
        }
     }  
     elsif($flag2 == 1){
       print STDERR "$_";
     }
   }
   close R2;
   
}

exit;
