#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

#TODO:this script is used for extract the oases assembled transcripts which can map to two different transcripts

my $refseq;
my $transcript;
my $blat_result;
my $lanename;
my $lanepath;

GetOptions (
            "refseq|r=s"     => \$refseq,
            "transcript|t=s" => \$transcript,
            "blat|b=s"       => \$blat_result,
            "lanename=s"     => \$lanename,
            "lanepath=s"     => \$lanepath,
            "help|h"         => sub{
                                    print "usage: $0 [options]\n\nOptions:\n\t--refseq\t\tthe fasta refseq mRNA file\n";
                                    print "\t--transcript\tthe transcript.fa assembled by oases\n";
                                    print "\t--blat\tthe transcript blat refseq file\n";
                                    print "\t--lanename\tthe name of the lane to be processed\n";
                                    print "\t--lanepath\tthe path root of the lane\n";
                                    print "\t--help\t\tprint this help message\n";
                                    exit 0;
                                   }

           );


#instore the refseq information into a hash
my %refseq;
open REFSEQ, "$refseq";
while ( <REFSEQ> ){
    chomp;
    if (/^>(.+\|).+\((.+?)\).+?$/){
        my $id = $1;
        my $name = $2;
        $refseq{$id} = $name;
    }
}
close REFSEQ;

#print Dumper(\%refseq);

#instore the transcript information into a hash
my %transcript;
my $trans_name;
my %trans_seq;
open TRANSCRIPT, "$transcript";
while ( <TRANSCRIPT> ){
    chomp;
    if (/^>(.+)$/){
	$trans_name = $1;
    }
    else{
	s/\n//g;
        s/\s//g;
        $transcript{$trans_name} += length($_); 
        $trans_seq{$trans_name} .= "$_\n";
    }
}
close TRANSCRIPT;

#print Dumper(\%transcript);

open BLAT, "$blat_result";
my %BLAT;
my $query;
while ( <BLAT> ){
    chomp;
  if (/^#\s+Query:\s+(.+)$/){
      $query = $1;
  }
  elsif ($_ !~ /^#/){
    my ($Qid, $Tid, $ID, $Alen, $mis, $gap, $QS, $QE, $TS, $TE, $EV, $BS) = split /\t/;
    if ($ID >= 90 and $BS >= 20) {   #threshold of idendity and bit score
      my %hash_tmp = ('id'=>$ID,'alen'=>$Alen, 'qs'=>$QS, 'qe'=>$QE, 'ts'=>$TS, 'te'=>$TE);
      push (@{$BLAT{$Qid}{$Tid}}, \%hash_tmp);
    }
  }
}
close BLAT;


my %combined;
foreach my $transcript (keys %BLAT){
    foreach my $refseq (keys %{$BLAT{$transcript}}){
        my $Alen_c = 0;  #c means combined
        my $QS_c   = 0;
        my $QE_c   = 0;
        my $TS_c   = 0;
        my $TE_c   = 0;
        my $indicator;
	foreach my $blat (@{$BLAT{$transcript}{$refseq}}){
	   my $Alen = $blat->{alen};
           my $QS   = $blat->{qs};
           my $QE   = $blat->{qe};
           my $TS   = $blat->{ts};
           my $TE   = $blat->{te};
           if (!defined $indicator){$Alen_c = $Alen; $QS_c = $QS; $QE_c = $QE; $TS_c = $TS; $TE_c = $TE; $indicator = 'SUN';}
           else {  #compare two position,there might be some risk!!!
	       if ($QE > $QE_c) {
                   $QE_c = $QE;  #move end
                   $TE_c = $TE;
                   $Alen_c += $QE-$QE_c;
               }
               if ($QS < $QS_c) {
                   $QS_c = $QS;  #move start
                   $TS_c = $TS;
                   $Alen_c += $QS_c-$QS;
               }
           }
        }
        my %hash_tmp = ('alen'=>$Alen_c, 'qs'=>$QS_c, 'qe'=>$QE_c, 'ts'=>$TS_c, 'te'=>$TE_c);
        $combined{$transcript}{$refseq} = \%hash_tmp;
    }
}


my %filtered;
my %BLAT_gn;  #restore non-overlapping gene name for each assembled transcript
foreach my $transcript (keys %combined) {

    my $transcript_length = $transcript{$transcript};

    my @tmp;
    my $indicator;
    foreach my $refseq (keys %{$combined{$transcript}}) {
        my $QS = $combined{$transcript}{$refseq}->{qs};
        my $QE = $combined{$transcript}{$refseq}->{qe};
        my $Alen = $combined{$transcript}{$refseq}->{alen};
        if (! defined $indicator){
	    push (@tmp, $refseq);
            $indicator = 'SUN';
        }
        else {  #compare with the previous refseq to generate a non-including refseq lists
            my $flag = 0;
	    for (my $i=0; $i<=$#tmp; $i++) {
               my $QS_b = $combined{$transcript}{$tmp[$i]}->{qs};
               my $QE_b = $combined{$transcript}{$tmp[$i]}->{qe};
               if ($QS >= $QS_b and $QE <= $QE_b){ #included
                   $flag = 1;
               }

               if (($QS < $QS_b and $QE >= $QE_b) or ($QS <= $QS_b and $QE > $QE_b)) { #includ the pervious
                   $tmp[$i] = $refseq; #reset the refseq
                   $flag = 4;
               }
            }
            if ($flag == 0){
               push (@tmp, $refseq);
            }
        }
    } #for each refseq

    #now in @tmp we have the non-including refseq lists
    my @tmp2;
    my $indicator2;
    foreach my $refseq (@tmp) {
        my $QS = $combined{$transcript}{$refseq}->{qs};
        my $QE = $combined{$transcript}{$refseq}->{qe};
        my $Alen = $combined{$transcript}{$refseq}->{alen};
        if (! defined $indicator2) {
            push (@tmp2, $refseq);
            $indicator2 = 'SUN';
        }
        else {  #check whether non-overlap or overlap to a limited extend
            my $flag = 0;
            for (my $i=0; $i<=$#tmp2; $i++) {
              my $QS_b = $combined{$transcript}{$tmp[$i]}->{qs};
              my $QE_b = $combined{$transcript}{$tmp[$i]}->{qe};

              if (($QS < $QE_b and $QE >= $QE_b and $QE_b - $QS > 6) or ($QS_b < $QE and $QS <= $QS_b and $QE - $QS_b > 6)) { #overlapping
                     $flag = 2;
              }

              if (($QS > $QE_b and $QS - $QE_b > 6) or ($QS_b > $QE and $QS_b - $QE > 6)) { #far apart
                     $flag = 3;
              }
            }
            unless ( $QS/$transcript_length < 0.05 or ($transcript_length-$QE)/$transcript_length < 0.05 ) {  #skip some strange in-middle blat result !(risky)
                 $flag = 6;
            }
            if ($flag == 0) {
               push (@tmp2, $refseq);
            }
        }
     }

     foreach my $refseq (@tmp2){

        $filtered{$transcript}{$refseq} = $combined{$transcript}{$refseq};
        my $gene_name = $refseq{$refseq};
        my $orientation;
        if ($combined{$transcript}{$refseq}->{te} > $combined{$transcript}{$refseq}->{ts}){
           $orientation = "->";
        }
        else{
           $orientation = "<-";
        }
        $BLAT_gn{$transcript}{$gene_name} = $orientation;

     }

}


open OUT, ">$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration.seq";
foreach my $transcript(keys %BLAT_gn) {

  my @genes = keys %{$BLAT_gn{$transcript}};
  if (scalar(@genes) == 2) {

    print "$transcript\t$transcript{$transcript}";
    print OUT ">$transcript\n$trans_seq{$transcript}";


    my @refseq = sort { $filtered{$transcript}{$a}->{qs} <=> $filtered{$transcript}{$b}->{qs} } keys %{$filtered{$transcript}};

    my $bp_s = $filtered{$transcript}{$refseq[0]}->{qe};
    my $bp_e = $filtered{$transcript}{$refseq[$#refseq]}->{qs};
    my $breakpoint = $bp_s.'..'.$bp_e;

    my @ref = map {(my $x = $_) =~ s/^gi\|.+?\|ref\|(.+?)\.\d+\|/\1/; $x;} @refseq;

    my $gene_name;
    if ($refseq{$refseq[0]} eq $genes[0]){  #genes0 is the left gene
        $gene_name = $genes[0]."\(".$ref[0]."\)\t".$genes[1]."\(".$ref[$#ref]."\)\t".$BLAT_gn{$transcript}{$genes[0]}.$BLAT_gn{$transcript}{$genes[1]};
    }
    else {
        $gene_name = $genes[1]."\(".$ref[0]."\)\t".$genes[0]."\(".$ref[$#ref]."\)\t".$BLAT_gn{$transcript}{$genes[1]}.$BLAT_gn{$transcript}{$genes[0]};
    }

    print "\t$gene_name";
    print "\t$breakpoint\n";

  }
}
close OUT;

exit;
