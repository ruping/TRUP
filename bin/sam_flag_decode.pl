my $flag = $ARGV[0];

my $code = decode_FLAG($flag);

print "$code\n";

sub decode_FLAG {
   my $compound = $_[0];
   my $code;
   if ($compound & 1)  {$code .= 'p';} #the read is paired in sequencing
   if ($compound & 2)  {$code .= 'P';} #the read is mapped in a proper pair
   if ($compound & 4)  {$code .= 'u';} #the read is unmapped
   if ($compound & 8)  {$code .= 'U';} #the mate is unmapped
   
   if ($compound & 16) {$code .= '-';} #the read strand -
   else                {$code .= '+';} #the read strand +
   
   if ($compound & 32)   {$code .= '-';} #the mate strand -
   elsif (! $compound & 8) {$code .= '+';} #the mate strand +
   
   if ($compound & 64) {$code .= '1';} #the read is the first in the pair
   if ($compound & 128){$code .= '2';} #the read is the second in the pair
   
   if ($compound & 256){$code .= 'N';} #not primary record (a read having split hits may have multiple primary alignment records)
   
   return $code;
}
