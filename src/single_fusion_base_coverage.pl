#!/usr/bin/perl
use strict;
use Getopt::Long;

my $fusion;
my $mapping;
my $read_length;

GetOptions (
              "fusion|f=s"         => \$fusion,
              "mappingfile|m=s"    => \$mapping,
              "readlength|r=i"     => \$read_length,
              "help|h" => sub{
                             print "usage: $0 [options]\n\nOptions:\n\t--fusion\t\tthe name of the fusion transcript\n";
                             print "\t--mappingfile\tthe mappingfile of the reads onto the assembled fusion candidates\n";
                             print "\t--readlength\tthe length of the read\n";
                             print "\t--help\t\tprint this help message\n";
                             exit 0;
                            }
           );


open MAPPING, "$mapping";
my %coverage;
while ( <MAPPING> ){

      chomp;
      my ($read, $strand, $candidate, $start, $read_seq, $read_qual, $multi, $mismatch) = split /\t/;
      $candidate =~ /^(Locus_\d+_Transcript_\d+\/\d+_Confidence_([10]\.\d+))\|(\d+)\|(.+?)\_(.+?)\|(\d+)\.\.(\d+)\|(.+?)\|(.+)$/;
      my $transcript = $1;
      my $confidence = $2;
      my $length = $3;
      my $message1 = $4;
      my $message2 = $5;
      my $bp_s = $6;
      my $bp_e = $7;
      my $blat1 = $8;
      my $blat2 = $9;
      
      my ($range_s, $range_e);
      if ($bp_s > $bp_e){$range_s = $bp_e; $range_e = $bp_s;}
      if ($bp_s <= $bp_e){$range_s = $bp_s; $range_e = $bp_e;}

      my %newinfo = ('con'=>$confidence, 'length'=>$length, 'gene1'=>$message1, 'gene2'=>$message2, 'bps'=>$range_s, 'bpe'=>$range_e, 'blat1'=>$blat1, 'blat2'=>$blat2);

      if ($transcript =~ /$fusion/){

          $read =~ /(.+?)\/([12])/;
          my $read_root = $1;
          my $read_end = $2;
          $start += 1;
          $coverage{$fusion}{info} = \%newinfo;
          $coverage{$fusion}{$read_root}{$read_end}{$start} = $_;

      }
}
close MAPPING;


my $length = $coverage{$fusion}{info}->{'length'};
my $gene1 = $coverage{$fusion}{info}->{'gene1'};
my $gene2 = $coverage{$fusion}{info}->{'gene2'};
my $bps = $coverage{$fusion}{info}->{'bps'};
my $bpe = $coverage{$fusion}{info}->{'bpe'};
my $blat1 = $coverage{$fusion}{info}->{'blat1'};
my $blat2 = $coverage{$fusion}{info}->{'blat2'};


my $newtitle = join("\t", $fusion, $length, $gene1.'-'.$gene2, $bps.'..'.$bpe, $blat1, $blat2);
print "#$newtitle\n";


my %base;
foreach my $read_root (keys %{$coverage{$fusion}}){

     my @pair = keys %{$coverage{$fusion}{$read_root}};
     
     if (scalar(@pair) == 1){ #only one end map 
         
         my @starts = keys %{$coverage{$fusion}{$read_root}{$pair[0]}};
         next if (scalar (@starts) >= 2);
         my $start = $starts[0];
         my $end = $start+$read_length-1;
         $base{$start}{$end}++;

     }

     if (scalar(@pair) == 2){ #pair mapped
        
         my @starts1 = keys %{$coverage{$fusion}{$read_root}{1}};
         my @starts2 = keys %{$coverage{$fusion}{$read_root}{2}};
         next if (scalar (@starts1) >= 2 or scalar (@starts2) >= 2);
         my $start1 = $starts1[0];
         my $start2 = $starts2[0];
         my $end1 = $start1+$read_length-1;
         my $end2 = $start2+$read_length-1;
         
         my $frag_s;
         my $frag_e;
         if ($start1 <  $start2){$frag_s=$start1;$frag_e=$end2;}
         if ($start1 >= $start2){$frag_s=$start2;$frag_e=$end1;}

         $base{$frag_s}{$frag_e}++;

     }
}


my $ptr = 0;
my @allstarts = sort { $a <=> $b } keys %base;

for (1..$length) {
    
   my $pos = $_;
   my $cov = 0;
   my $off = 0;

   while ((($ptr+$off) <= $#allstarts) and ($allstarts[$ptr+$off] <= $pos)){

      my @ends = sort {$a <=> $b} keys %{$base{$allstarts[$ptr+$off]}};

      foreach my $end (@ends){
         if ($end >= $pos){
             $cov++;
         }
      }#foreach
      
      $off++;

      if ($ends[$#ends] <= $pos){
         $ptr += $off;
         $off = 0;
      }

   }#while

   print "$pos\t$cov\n";

}


exit 22;
