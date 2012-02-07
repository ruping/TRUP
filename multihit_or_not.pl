#!/usr/bin/perl
#
#TODO:1) read in the .sam mapping file and generate uniquely mapped and multiple mapped reads
#     2) generate the numbers for basic statistics of mapping
#

use strict;
use File::Basename;
use File::Glob ':glob';
use Getopt::Long;

my $type;
my $path;
my $lanename;
my $readfile;
my $mapfile;
my $read_length;

GetOptions(      
             "type|t=s"     => \$type,     
             "path|p=s"     => \$path,
             "lanename=s"   => \$lanename,
             "readfile|r=s" => \$readfile,
             "mapfile|m=s"  => \$mapfile,
             "length|l=i"   => \$read_length, 
	     "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--type\t\tthe type of the reads, either s (single-end) or p (pair-end)\n";
                             print "\t--path\t\tthe path to which the UMR, POS files writing\n";
                             print "\t--lanename\tused as the name prefix of the output files\n";
                             print "\t--readfile\tone of the read file\n";
			     print "\t--mapfile\tthe tophat mapping hits in sam format\n";
			     print "\t--length\tlength of the reads\n";
			     print "\t--help\t\tprint this help message\n";
                             exit 0;
			    }
          );


my $total_pos = $path."$lanename\.pos\.gff";
my $umrfile = $path."$lanename\.UMR\.gff";
my $chrfile = $path."$lanename\.chr_map_stat";
my $totfile = $path."$lanename\.mapping\.stats";

#for general usage
my %mapped;
my %multi;

#for pair-end usage
my %all_aligned;
my %correct_pair;
my %wrong_pair;
my %lost_pair;
my %spliced;

#count total number of reads
my $nreads_tot = (`wc -l $readfile`)/4;
print STDERR "total number of reads = $nreads_tot\n";

#index the chr start in the mapping file;
open MAP, "$mapfile";
my $old_chr;
my %chr_start;
while ( <MAP> ) {

  next if /^@/;  #ignore comments
  chomp;

  my @cols = split /\t/;
  my $chr = $cols[2];

  $chr = 'chr'.$chr unless ($chr =~ /^chr/);
  $chr .= 'T' if ($chr eq 'chrM');

  if ($chr ne $old_chr){
     $chr_start{$chr} = tell MAP;
  }

  $old_chr = $chr;
}
close MAP;
print STDERR "mapping index loaded\n";


#
##########for single end mapping##################
#
if ($type eq 's'){

open UMR, ">$umrfile";
open POS, ">$total_pos";

my $unique = 0;
my $unique_partial = 0;
my %unique_chr;
my %pileup;
my $pos = 0;
my $pileup_pos = 0;
my %pos_chr;
my %pileup_pos_chr;


print STDERR "now output $umrfile and $total_pos files\n";
foreach my $chr (sort keys %chr_start){

  load_mapping($chr, 'single');
  print STDERR "$chr mapping loaded\n";

  %pileup = ();   #redefine the hash for pileup things

  foreach my $name (sort {my $sa = $mapped{$a}->{'start'}; my $sb = $mapped{$b}->{'start'}; $sa <=> $sb} keys %mapped){

    $unique++;  

    my $strand = $mapped{$name}->{'strand'};
    my $start  = $mapped{$name}->{'start'};
    my $end    = $mapped{$name}->{'end'};
    my $type   = $mapped{$name}->{'type'};
    my $gff    = $mapped{$name}->{'gff'};

    if ($type eq 'partial'){
      $unique_partial++;
    }

    #from here, only calculate pileup for UMR
    push (@{$pileup{$start}{$end}{$strand}}, $gff);    
    $unique_chr{$chr}++;   #uniquely mapped each chromosome
    print UMR "$gff\n";

  }
  
  foreach my $start (keys %pileup){
    foreach my $end (sort {$a <=> $b} keys %{$pileup{$start}}){
      foreach my $strand (sort keys %{$pileup{$start}{$end}}){
      
        $pos++;
        $pos_chr{$chr}++;
        
        my $phase = '.';
        my $times = scalar (@{$pileup{$start}{$end}{$strand}});
        my $tag = 'Pileup='.$times;
        
        if($times > 1){
          $pileup_pos++;
          $pileup_pos_chr{$chr}++;
        }
        
        print POS "$chr\tTopHat\tpileup\t$start\t$end\t$phase\t$strand\t$phase\t$tag\n";
        
      }#strand
    }#end
  }#start
    
}
close UMR;
close POS;

my $multi = scalar (keys %multi);
my $mapped_tot = $multi + $unique;

open CHR, ">$chrfile";
foreach my $chr (sort {$a cmp $b} keys %unique_chr){
 print CHR "$chr\t$unique_chr{$chr}\t$pos_chr{$chr}\t$pileup_pos_chr{$chr}\n";
}
close CHR;


my $unmapped_reads = $nreads_tot-$mapped_tot;
#print "mapped_total	$mapped_tot\n";
open TOT, ">$totfile";
print TOT "UMR\tUMR_partial\tMMR\tunmapped\ttotal_pos\tpileup_pos\n";
print TOT "$unique\t$unique_partial\t$multi\t$unmapped_reads\t$pos\t$pileup_pos\n";
close TOT;

}

#
##########for pair end mapping##################
#
if ($type eq 'p'){
     
  #generate file for comparison of chr and gene annotation:
  my %pileup;
  my %unique_chr;
  my $fake_correct_pair; #adjust the number of correct pair
  my $pos = 0;
  my $pileup_pos = 0;
  my %pos_chr;
  my %pileup_pos_chr;
  open UMR, ">$umrfile";
  open POS, ">$total_pos";
  
  foreach my $chr (sort keys %chr_start){
  
    load_mapping($chr, 'pair');
    print STDERR "$chr mapping loaded\n";
    
    %pileup = ();
    
    foreach my $fragment (keys %mapped){
       
         if (scalar(@{$mapped{$fragment}}) == 1){
            if (exists $correct_pair{$fragment}){
               $fake_correct_pair++;    
            }
            next;
         }
       
         my $Rname1 = $mapped{$fragment}[0]->{'Rname'};
         my $CIGAR1 = $mapped{$fragment}[0]->{'CIGAR'};
         my $code1  = $mapped{$fragment}[0]->{'code'};
         my $Pos1   = $mapped{$fragment}[0]->{'Pos'};
         my $MAPQ1  = $mapped{$fragment}[0]->{'MAPQ'};
         my $Rname2 = $mapped{$fragment}[1]->{'Rname'};
         my $CIGAR2 = $mapped{$fragment}[1]->{'CIGAR'};
         my $Pos2   = $mapped{$fragment}[1]->{'Pos'};
          
         my $fragment_length = abs($Pos2-$Pos1)+$read_length;
        
         next if ($CIGAR1 ne "$read_length\M" or $CIGAR2 ne "$read_length\M");  #skip spliced mapping
       
         my $source = 'TopHat';
         my $type;
         $type = 'right_pair' if ($code1 =~ /P/);
         $type = 'wrong_pair' if ($code1 !~ /P/);
         my $start;
         my $end; 
         my $score = $MAPQ1;
         my $strand = '+'; 
         my $phase = '.'; 
         my $tag = "Fragment=$fragment;length=$fragment_length";
       
         #generate the gff
         if (($Pos2-$Pos1) >= 0){$start = $Pos1; $end = $Pos2+$read_length-1;}
         if (($Pos2-$Pos1) <  0){$start = $Pos2; $end = $Pos1+$read_length-1;}
       
         my $gff = join ("\t", $chr,$source,$type,$start,$end,$score,$strand,$phase,$tag);
         push (@{$pileup{$start}}, $gff); 
         print UMR "$gff\n";
         $unique_chr{$chr}++;
     }
     
     foreach my $start (sort {$a <=> $b} keys %pileup){
     
        $pos++;
        $pos_chr{$chr}++;

        my @cols = split (/\t/, ${$pileup{$start}}[0]);
        my $end = $cols[4];
        my $strand = '+';
        my $phase = '.';
        my $time = scalar (@{$pileup{$start}});
        my $tag = 'Pileup='.$time;
        
        if(@{$pileup{$start}} > 1){
          $pileup_pos++;
          $pileup_pos_chr{$chr}++;
        }
        
        print POS "$chr\tTopHat\tpileup\t$start\t$end\t$phase\t$strand\t$phase\t$tag\n";
     }      
  }
  close UMR;
  close POS;

  open TOT, ">$totfile";
  my $multi = scalar (keys %multi);
  my $total_aligned = scalar (keys %all_aligned);
  my $total_correct_pair = scalar (keys %correct_pair) - $fake_correct_pair;
  my $total_wrong_pair = scalar (keys %wrong_pair) + $fake_correct_pair;
  my $total_lost_pair = scalar (keys %lost_pair);
  my $total_spliced = scalar (keys %spliced);
  printf TOT "%s\n", join("\t", 'Sequenced', 'Aligned','Proper-Pair','Wrong-Pair','Singleton','Spliced','Multi-mapping');
  printf TOT "%s\n", join("\t", $nreads_tot, $total_aligned, $total_correct_pair, $total_wrong_pair, $total_lost_pair, $total_spliced, $multi);   
  close TOT;
  
  
  open CHR, ">$chrfile";
  foreach my $chr (sort keys %unique_chr){
    print CHR "$chr\t$unique_chr{$chr}\t$pos_chr{$chr}\t$pileup_pos_chr{$chr}\n";
  }
  close CHR;
   
}

exit 22;

#sub region#########################################################################################################

sub load_mapping {

    my $chr_w = shift;
    my $style = shift; 

    %mapped = ();

    open MAPPING, "<$mapfile";
    seek(MAPPING, $chr_start{$chr_w}, 0);
    while( <MAPPING> ){

      chomp;
      my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;
      #generate the new data structure for conversion to gff format 
      my $chr = $Rname;
      $chr = 'chr'.$chr unless ($chr =~ /^chr/);
      $chr .= 'T' if ($chr eq 'chrM');
      
      last if ($chr ne $chr_w);
      
      if ($_ !~ /\tNH:i:1$/){  #ignore multiple mapping
        $multi{$Qname}++;
        next;
      }
      
      if ($style =~ /single/){       
         my $start = $Pos + 1;
         my $end = $start + $read_length - 1;
      
         my $strand_new;
         $strand_new = '+' if ($FLAG == 0);
         $strand_new = '-' if ($FLAG == 16);
      
         my $phase = '.';
         my $tag;
         my $type;
      
         if ($CIGAR eq "$read_length\M"){
           $tag  = 'ID='.$Qname;
           $type = 'whole';
         }
         else{
           $tag  = 'ID='.$Qname.';CIGAR='.$CIGAR;
           $type = 'partial';
           if ($CIGAR =~ /^\d+M(\d+)N\d+M$/){  #problematic
             $end += $1;
           }
         }
         my $gff = join("\t", $chr, 'TopHat', $type, $start, $end, '.', $strand_new, $phase, $tag);
         $mapped{$Qname} = {'chr'=>$chr, 'start'=>$start, 'end'=>$end, 'strand'=>$strand_new, 'gff'=>$gff, 'type'=>$type};
       }
       
       if ($style =~ /pair/){
          
          my $code = decode_FLAG($FLAG);
          
          my %hash_tmp = ('code'=>$code, 'Rname'=>$Rname, 'Pos'=>$Pos, 'CIGAR'=>$CIGAR, 'MAPQ'=>$MAPQ);
          
          $all_aligned{$Qname}++;
          
          if ($CIGAR ne "$read_length\M"){
            $spliced{$Qname}++;
          }
     
          if ($code =~ /P/){
            $correct_pair{$Qname}++;
            push (@{$mapped{$Qname}}, \%hash_tmp); 
          }
     
          if ($code !~ /P/ and $mateRname eq '='){
            $wrong_pair{$Qname}++;
            push (@{$mapped{$Qname}}, \%hash_tmp);
          }
     
          if ($code =~ /U/){
            $lost_pair{$Qname}++;
          }
        }
               
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
