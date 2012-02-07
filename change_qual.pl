
open IN, "$ARGV[0]";
while ( <IN> ){
 chomp;
 $_=$_-31;
 print "$_\n";
}
close IN;
