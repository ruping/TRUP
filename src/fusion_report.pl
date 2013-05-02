#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my %fusion;

#read data

open IN, shift;

while ( <IN> ) {

   chomp;

   my ($transcript, $idx, $blatScore, $length, $partner, $dir, $breakpoint, $rep, $selfChain, $type, $strand, $blat1, $blat2, $covSpan, $covEncom, $covAll, $covSpanAll, $spanScore) = split /\t/;

   next if ($selfChain eq 'CC');

   my $trans;
   my $fusion;

   ($trans->{'name'} = $transcript) =~ s/^#//;
   $trans->{'spanning'} = $covSpan;
   $trans->{'encompassing'} = $covEncom;
   $trans->{'supporting'} = $covAll;
   $trans->{'spanningAll'} = $covSpanAll;
   $trans->{'spanningScore'} = $spanScore;
   $trans->{'type'} = $type;
   $trans->{'rep'} = $rep;
   $trans->{'selfChain'} = $selfChain;
   $trans->{'strand'} = $strand;
   $trans->{'dir'} = $dir;

   if ($partner =~ /^((.+?)\((.+?)\))\-((.+?)\((.+?)\))$/) {
      $trans->{'gene1'} = $1;
      $trans->{'ensembl1'} = $2;
      $trans->{'symbol1'}  = $3;
      $trans->{'gene2'} = $4;
      $trans->{'ensembl2'} = $5;
      $trans->{'symbol2'}  = $6;
   } elsif ($partner =~ /^((.+?)\((.+?)\))\-IGR$/) {
      $trans->{'gene1'} = $1;
      $trans->{'ensembl1'} = $2;
      $trans->{'symbol1'}  = $3;
      $trans->{'gene2'} = 'IGR';
   } elsif ($partner =~ /^IGR\-((.+?)\((.+?)\))$/) {
      $trans->{'gene2'} = $1;
      $trans->{'ensembl2'} = $2;
      $trans->{'symbol2'}  = $3;
      $trans->{'gene1'} = 'IGR';
   } else {
      print STDERR "Two partners in the IGR...$transcript\n";
   }

   if ($breakpoint =~ /^(\d+)\.\.(\d+)$/) {
      $trans->{'bpit1'} = $1;
      $trans->{'bpit2'} = $2;
   }

   if ($blat1 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/) {
      $trans->{'blat1_qs'}  = $1;
      $trans->{'blat1_qe'}  = $2;
      $trans->{'blat1_st'}  = $3;
      $trans->{'blat1_chr'} = $4;
      $trans->{'blat1_ts'}  = $5;
      $trans->{'blat1_te'}  = $6;
      if ($trans->{'blat1_st'} eq '+'){
        $trans->{'bpig1'} = $trans->{'blat1_chr'}.':'.$trans->{'blat1_te'};
      } else {
        $trans->{'bpig1'} = $trans->{'blat1_chr'}.':'.$trans->{'blat1_ts'};
      }
   }

   if ($blat2 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/) {
      $trans->{'blat2_qs'}  = $1;
      $trans->{'blat2_qe'}  = $2;
      $trans->{'blat2_st'}  = $3;
      $trans->{'blat2_chr'} = $4;
      $trans->{'blat2_ts'}  = $5;
      $trans->{'blat2_te'}  = $6;
      if ($trans->{'blat2_st'} eq '+'){
        $trans->{'bpig2'} = $trans->{'blat2_chr'}.':'.$trans->{'blat2_ts'};
      } else {
        $trans->{'bpig2'} = $trans->{'blat2_chr'}.':'.$trans->{'blat2_te'};
      }
   }

   if ($dir eq '->->') {
      $fusion->{'identifier'} = $trans->{'gene1'}."\t".$trans->{'gene2'};
      $fusion->{'bpig1'} = $trans->{'bpig1'};
      $fusion->{'bpig2'} = $trans->{'bpig2'};
      $fusion->{'dir'} = '->->';
      $fusion->{'rep'} = $trans->{'rep'};
      $fusion->{'selfChain'} = $trans->{'selfChain'};
   } elsif ($dir eq '<-<-') {
      $fusion->{'identifier'} = $trans->{'gene2'}."\t".$trans->{'gene1'};
      $fusion->{'bpig2'} = $trans->{'bpig1'};
      $fusion->{'bpig1'} = $trans->{'bpig2'};
      $fusion->{'dir'} = '->->';
      $fusion->{'rep'} = reverse($trans->{'rep'});
      $fusion->{'selfChain'} = reverse($trans->{'selfChain'});
   } else {  # either <--> or -><- or other
      my @sorted_partenr = sort {$a cmp $b} ($trans->{'gene1'}, $trans->{'gene2'});
      $fusion->{'identifier'} = join("\t", @sorted_partenr);
      if ($sorted_partenr[0] eq $trans->{'gene1'}){
         $fusion->{'bpig1'} = $trans->{'bpig1'};
         $fusion->{'bpig2'} = $trans->{'bpig2'};
         $fusion->{'dir'} = $trans->{'dir'};
         $fusion->{'rep'} = $trans->{'rep'};
         $fusion->{'selfChain'} = $trans->{'selfChain'};
      } else {
         $fusion->{'bpig2'} = $trans->{'bpig1'};
         $fusion->{'bpig1'} = $trans->{'bpig2'};
         if ($dir eq '<-->' or $dir eq '-><-') {
           $fusion->{'dir'} = $trans->{'dir'};
         } else {
           $fusion->{'dir'} = reverse($trans->{'dir'});
           $fusion->{'dir'} =~ tr/\>\</\<\>/;
         }
         $fusion->{'rep'} = reverse($trans->{'rep'});
         $fusion->{'selfChain'} = reverse($trans->{'selfChain'});
      }
   }

   my $fusion_coors = $fusion->{'bpig1'}."\t".$fusion->{'bpig2'};
   my $fusion_genes = $fusion->{'identifier'};

   push(@{$fusion{$fusion_genes."\t".$fusion_coors."\t".$fusion->{'dir'}."\t".$fusion->{'rep'}."\t".$fusion->{'selfChain'}}}, $trans);

}

close IN;


#print Dumper(\%fusion);

#prepare data

my %report;

foreach my $fusion (keys %fusion){

   my $report;
   $report->{'id'} = $fusion;
   $report->{'spanning'} = 0;
   $report->{'encompassing'} = 0;
   $report->{'supporting'} = 0;
   $report->{'spanningAll'} = 0;
   $report->{'spanningScore'} = 0;
   $report->{'type'} = '';
   $report->{'strand'} = '';
   $report->{'transcripts'} = '';

   #every transcript
   foreach my $trans (@{$fusion{$fusion}}){

     $report->{'transcripts'} .= $trans->{'name'}.',';

     if ($trans->{'supporting'} > $report->{'supporting'}) {
       $report->{'spanning'} = $trans->{'spanning'};
       $report->{'encompassing'} = $trans->{'encompassing'};
       $report->{'supporting'} = $trans->{'supporting'};
       $report->{'spanningAll'} = $trans->{'spanningAll'};
       $report->{'spanningScore'} = $trans->{'spanningScore'};
     }

     if ($report->{'type'} eq ''){
       $report->{'type'} = $trans->{'type'};
     }
     if ($report->{'strand'} eq ''){
       $report->{'strand'} = $trans->{'strand'};
     }


   }

   $report{$fusion} = $report;

}


#write out

foreach my $fusion (sort {$report{$b}->{'supporting'} <=> $report{$a}->{'supporting'} or $report{$b}->{'spanningScore'} <=> $report{$a}->{'spanningScore'}} keys %report){
    my $span = $report{$fusion}->{'spanning'};
    my $encom = $report{$fusion}->{'encompassing'};
    my $supp = $report{$fusion}->{'supporting'};
    my $spanAll = $report{$fusion}->{'spanningAll'};
    my $SS = $report{$fusion}->{'spanningScore'};
    my $type = $report{$fusion}->{'type'};
    my $strand = $report{$fusion}->{'strand'};
    my $trans = $report{$fusion}->{'transcripts'};
    printf "%s\n", join("\t", $fusion, $type, $strand, $span, $encom, $supp, $spanAll, $SS, $trans);
}

exit 0;
