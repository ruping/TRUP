#!/usr/bin/perl

######
#TODO: 1)transforming the title of fusion candidates (addingthe 5', 3')
###### 2)get fusion coverage from the pair-end mapping and the assembly

use strict;
use Getopt::Long;
use Data::Dumper;

my $type;
my $mapping;
my $read_length;
my $accepthits;
my $encomcov;
my $gene_annotation;
my $selfChain;
my $ensemble_ref_name;
my $locnamefile;
my $ucsc_refgene;
my $repeatmasker;
my $genomeBlatPred = "SRP";
my $anchorLen = 12;

GetOptions (
              "type|t=s"              => \$type,
              "mappingfile|m=s"       => \$mapping,
	      "readlength|r=i"        => \$read_length,
	      "accepthits|ah=s"       => \$accepthits,
	      "encomcov=s"            => \$encomcov,
	      "geneanno|ga=s"         => \$gene_annotation,
              "selfChain=s"           => \$selfChain,
              "ensrefname|ern=s"      => \$ensemble_ref_name,
              "locname|loc=s"         => \$locnamefile,
              "refgene|ref=s"         => \$ucsc_refgene,
              "repeatmasker=s"        => \$repeatmasker,
              "genomeBlatPred=s"      => \$genomeBlatPred,
	      "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--type\t\tthe type of the reads, either pair or single\n";
                             print "\t--mappingfile\tthe mappingfile of the reads onto the assembled fusion candidates\n";
			     print "\t--readlength\tthe length of the read\n";
			     print "\t--geneanno\tthe ensemble gene anotation file\n";
                             print "\t--selfChain\tthe self Chain annotation downloaded from UCSC Table.\n";
                             print "\t--ensrefname\tthe ensemble, refseq, gene name list table\n";
                             print "\t--locname\tloc_gene_name2location file\n";
                             print "\t--refgene\tdownloaded from ucsc the refgene annotation file\n";
                             print "\t--genomeBlatPred\tset to the file produced by analyzing the genome blat result.\n";
			     print "\t--help\t\tprint this help message\n\n";
                             exit 0;
			     }
           );

#repeat mask
my %repeatmask;
my @rs;       #start array for each chr
my $old_chr;  #checking the chr
my $ptr;      #pointer for repeatmask sub
open REPEAT, "$repeatmasker";
while ( <REPEAT> ){
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    my $chr = $tmp[0];
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
}
close REPEAT;

#Ucsc Self Chain
my %selfChain;
my @rs_selfChain;
my $old_chr_selfChain;
my $ptr_selfChain;
open SELFCHAIN, "$selfChain";
while ( <SELFCHAIN> ){
    next if /^#/;
    chomp;
    my ($bin, $score, $tName, $tSize, $tStart, $tEnd, $qName, $qSize, $qStrand, $qStart, $qEnd, $id, $normScore) = split /\t/;

    next if ($normScore < 10);

    next if ($selfChain{$tName}{$tStart} ne '' and $selfChain{$tName}{$tStart} >= $tEnd);
    $selfChain{$tName}{$tStart} = $tEnd;

    next if ($selfChain{$qName}{$qStart} ne '' and $selfChain{$qName}{$qStart} >= $qEnd);
    $selfChain{$qName}{$qStart} = $qEnd;
}
close SELFCHAIN;


my %name2ens;
my %loc;
my %refgenes;

if ($genomeBlatPred eq 'SRP') {
    open ENS, "$ensemble_ref_name";
    while ( <ENS> ){
	chomp;
	next if /^Ensembl/;
	my ($ensemble_id, $refseq_id, $gene_name) = split /\t/;
	$name2ens{$gene_name} = $ensemble_id;
    }
    close ENS;


    open LOC, "$locnamefile";
    while (<LOC>){
	next if /^Gene/;
	chomp;
	my ($gene, $ens, $chr, $start, $end, $strand) = split /\t/;
	unless ($chr eq ''){
	    if ($chr !~ /^chr/) {
		$chr = 'chr'.$chr;
	    }
	    if ($strand == 1) {
		$strand = '+';
	    }
	    if ($strand == -1) {
		$strand = '-';
	    }
	    $loc{$gene}{'chr'}        = $chr;
	    $loc{$gene}{'start'}      = $start;
	    $loc{$gene}{'end'}        = $end;
	    $loc{$gene}{'strand'}     = $strand;
	}
    }
    close LOC;


    open REFGENE, "$ucsc_refgene";
    while ( <REFGENE> ){
	chomp;
	next if /^#/;
	my @cols = split /\t/;
	my $refid   = $cols[1];
	my $chr     = $cols[2];
	my $strand  = $cols[3];
	my $start   = $cols[4];
	my $end     = $cols[5];
	$refgenes{$refid}{'chr'}    = $chr;
	$refgenes{$refid}{'start'}  = $start;
	$refgenes{$refid}{'end'}    = $end;
	$refgenes{$refid}{'strand'} = $strand;
    }
    close REFGENE;
} #if refseq mapping

open GA, "$gene_annotation";
my %gene;
my %gene_annotation;
my %mRNA_region;  #remember the mRNA region in this gene
my %gene2exon;
my %transcript2gene;
my @rs2;          #start array for each chr
my $old_chr2;     #checking the chr
my $ptr2;         #pointer for repeatmask sub
while ( <GA> ){
   next if /^#/;
   chomp;
   my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$tag) = split /\t/;

   if ($type eq 'gene'){

     $tag =~ /^ID=(.+?);Name=(.+?);/;
     my $ensemble = $1;
     my $gene = $2;
     my $genelength= $end-$start+1;
     $gene{$ensemble}{'chr'}    = $chr;
     $gene{$ensemble}{'start'}  = $start;
     $gene{$ensemble}{'end'}    = $end;
     $gene{$ensemble}{'strand'} = $strand;
     $gene{$gene}{'chr'}        = $chr;
     $gene{$gene}{'start'}      = $start;
     $gene{$gene}{'end'}        = $end;
     $gene{$gene}{'strand'}     = $strand;
     $gene_annotation{$chr}{$start}{$end} = join(',', $ensemble.'('.$gene.')', $strand, $genelength);

   } elsif ($type eq 'mRNA') {  #remember the mRNA-part for genes

     $tag =~ /ID\=(.+?)\;/;
     my $transcriptID = $1;
     $tag =~ /Parent\=(.+?)\;/;
     my $parentGene = $1;
     $transcript2gene{$transcriptID} = $parentGene;

     my @mRNA = ($start, $end);
     if (!exists($mRNA_region{$parentGene})) {
       $mRNA_region{$parentGene} = \@mRNA;
     } else {
       if ($start <= $mRNA_region{$parentGene}->[0]) {
         $mRNA_region{$parentGene}->[0] = $start;
       }
       if ($end >= $mRNA_region{$parentGene}->[1]) {
         $mRNA_region{$parentGene}->[1] = $end;
       }
     }                          #parentGene already exists

   } elsif ($type eq 'transcript') {

     $tag =~ /ID\=(.+?)\;/;
     my $transcriptID = $1;
     $tag =~ /Parent\=(.+?)\;/;
     my $parentGene = $1;
     $transcript2gene{$transcriptID} = $parentGene;

   } elsif ($type eq 'exon') {  #remember exons for gene

     $tag =~ /Parent\=(.+?)\;/;
     my $parentTranscriptID = $1;
     my $parentGene = $transcript2gene{$parentTranscriptID};
     $gene2exon{$parentGene}{$start}{$end} = $parentTranscriptID;

   }
}
close GA;

#merge exons for gene;
my %gene2exonMerged;
foreach my $ensGene (keys %gene2exon) {

  my $last_start = 0;
  my $last_end = 0;

  foreach my $current_start (sort {$a <=> $b} keys %{$gene2exon{$ensGene}} ) {
    my @current_ends = sort {$b <=> $a} keys %{$gene2exon{$ensGene}{$current_start}};  #may be multiple ends, pick the largest one
    my $current_end = $current_ends[0];

    if ( $current_start <= $last_end ) { #overlapping
       if ($current_end > $last_end) {
          $gene2exonMerged{$ensGene}{$last_start} = $current_end;
       } else { #included in the old one
          next;
       }
    } else { #new_starts
      $gene2exonMerged{$ensGene}{$current_start} = $current_end;
      $last_start = $current_start;
    }

    $last_end = $current_end;
  }

}
%transcript2gene = (); #clear memory
%gene2exon = ();  #clear memory


my %genomeBlatPred;
if ($genomeBlatPred ne 'SRP') {
  open GBP, "$genomeBlatPred";
  while ( <GBP> ){
    chomp;
    my ($transcript, $id, $length, $bp, $blat1, $blat2, $chimeric_score) = split /\t/;
    my $blatInfo;
    $blatInfo->{'length'} = $length;
    $bp =~ /^(\d+)\.\.(\d+)$/;
    $blatInfo->{'bp_s'} = $1;
    $blatInfo->{'bp_e'} = $2;
    $blatInfo->{'blat1'} = $blat1;
    $blatInfo->{'blat2'} = $blat2;
    $blatInfo->{'cScore'} = $chimeric_score;
    $genomeBlatPred{$transcript}{$id} = $blatInfo;
  }
  close GBP;
}


open IN, "samtools view $mapping |";
my %for_rm;
my %coverage;
my %for_encompass;
while ( <IN> ){

   next if /^\@/;
   chomp;
   #my @cols = split /\t/;

   if ($type =~ /pair/) { #bowtie

      my ($read, $flag, $candidate, $start, $mapQ, $cigar, $mateR, $matePos, $insert, $read_seq, $read_qual, @tags) = split /\t/;

      #1 QNAME String [!-?A-~]{1,255} Query template NAME
      #2 FLAG Int [0,216-1] bitwise FLAG
      #3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
      #4 POS Int [0,229-1] 1-based leftmost mapping POSition
      #5 MAPQ Int [0,28-1] MAPping Quality
      #6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
      #7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next segment
      #8 PNEXT Int [0,229-1] Position of the mate/next segment
      #9 TLEN Int [-229+1,229-1] observed Template LENgth
      #10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
      #11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33

      my @transcript;
      my $transcript;

      if ($candidate =~ /^(.+?Confidence_([10]\.\d+).*?)\|(\d+)\|(.+?)\+(.+?)\|(.+?)\|(\d+)\.\.(\d+)\|(.+?)\|(.+?)\|(.+)$/){
         $transcript->{'name'} = $1;
         $transcript->{'confidence'} = $2;
         $transcript->{'length'} = $3;
         $transcript->{'part1'} = $4;
         $transcript->{'part2'} = $5;
         $transcript->{'orientation'} = $6;
         $transcript->{'bp_s'} = $7;
         $transcript->{'bp_e'} = $8;
         $transcript->{'ron'}  = $9;
         $transcript->{'blat1'} = $10;
         $transcript->{'blat2'} = $11;

         $transcript->{'part1'} =~ /^(.+?)\((.+?)\)$/;
         $transcript->{'message1'} = $1;
         $transcript->{'message_ref1'} = $2;

         $transcript->{'part2'} =~ /^(.+?)\((.+?)\)$/;
         $transcript->{'message2'} = $1;
         $transcript->{'message_ref2'} = $2;
         $transcript->{'cScore'} = 0;
         push @transcript, $transcript;

      } elsif ($genomeBlatPred ne 'SRP') {
	  $candidate =~ /^.+?Confidence_([10]\.\d+).*/;
	  my $confidence = $1;
	  my $length = $genomeBlatPred{$candidate}{1}->{'length'};
	  #gene1 gene2 orientation later
	  my $chimeric_score =  $genomeBlatPred{$candidate}{1}->{'cScore'};
	  foreach my $idx (sort {$a<=>$b} keys %{$genomeBlatPred{$candidate}}) {
	      my $trans_add_id = $candidate."\t".$idx;               #transcript name add the id
	      $transcript->{'name'} = $trans_add_id;
	      $transcript->{'confidence'} = $confidence;
	      $transcript->{'length'} = $length;
	      $transcript->{'part1'} = '';
	      $transcript->{'part2'} = '';
	      $transcript->{'orientation'} = '';
	      $transcript->{'bp_s'} = $genomeBlatPred{$candidate}{$idx}->{'bp_s'};
	      $transcript->{'bp_e'} = $genomeBlatPred{$candidate}{$idx}->{'bp_e'};
	      $transcript->{'ron'}  = '';
	      $transcript->{'blat1'} = $genomeBlatPred{$candidate}{$idx}->{'blat1'};
	      $transcript->{'blat2'} = $genomeBlatPred{$candidate}{$idx}->{'blat2'};
	      $transcript->{'cScore'} = $genomeBlatPred{$candidate}{$idx}->{'cScore'};
              push @transcript, $transcript;
	  } #for each breakpoint of the current transcript
      }

      foreach my $transcript_bp (@transcript) {

          my $transcript_name = $transcript_bp->{'name'};
	  my $bp_s = $transcript_bp->{'bp_s'};
	  my $bp_e = $transcript_bp->{'bp_e'};
	  my $confidence = $transcript_bp->{'confidence'};
          my $length = $transcript_bp->{'length'};
          my $part1 = $transcript_bp->{'part1'};
          my $part2 = $transcript_bp->{'part2'};
          my $orientation = $transcript_bp->{'orientation'};
          my $ron = $transcript_bp->{'ron'};
          my $blat1 = $transcript_bp->{'blat1'};
          my $blat2 = $transcript_bp->{'blat2'};
          my $cScore = $transcript_bp->{'cScore'};

	  if ($coverage{$transcript_name}{'info'} eq '') {

	      my ($mg1_chr, $mg1_start, $mg1_end, $mg1_strand, $mg2_chr, $mg2_start, $mg2_end, $mg2_strand);

	      if ($genomeBlatPred eq 'SRP') {   #decide messenges if not genome blat

		  #message1
		  my $message1 = $transcript_bp->{'message1'};
		  my $message_ref1 = $transcript_bp->{'message_ref1'};
		  if ($gene{$message1} ne '') {
		      $mg1_chr   = $gene{$message1}{'chr'};
		      $mg1_start = $gene{$message1}{'start'};
		      $mg1_end   = $gene{$message1}{'end'};
		      $mg1_strand= $gene{$message1}{'strand'};
		  } elsif ($refgenes{$message_ref1} ne '') {
		      $mg1_chr   = $refgenes{$message_ref1}{'chr'};
		      $mg1_start = $refgenes{$message_ref1}{'start'};
		      $mg1_end   = $refgenes{$message_ref1}{'end'};
		      $mg1_strand= $refgenes{$message_ref1}{'strand'};
		  } elsif ($loc{$message1} ne '') {
		      $mg1_chr    = $loc{$message1}{'chr'};
		      $mg1_start  = $loc{$message1}{'start'};
		      $mg1_end    = $loc{$message1}{'end'};
		      $mg1_strand = $loc{$message1}{'strand'};
		  } elsif ($name2ens{$message1} ne '') {
		      my $ensemble_id1 = $name2ens{$message1};
		      if ($gene{$ensemble_id1} ne '') {
			  $mg1_chr   = $gene{$ensemble_id1}{'chr'};
			  $mg1_start = $gene{$ensemble_id1}{'start'};
			  $mg1_end   = $gene{$ensemble_id1}{'end'};
			  $mg1_strand= $gene{$ensemble_id1}{'strand'};
		      }
		  }

		  #message2
		  my $message2 = $transcript_bp->{'message2'};
		  my $message_ref2 = $transcript_bp->{'message_ref2'};
		  if ($gene{$message2} ne '') {
		      $mg2_chr   = $gene{$message2}{'chr'};
		      $mg2_start = $gene{$message2}{'start'};
		      $mg2_end   = $gene{$message2}{'end'};
		      $mg2_strand= $gene{$message2}{'strand'};
		  } elsif ($refgenes{$message_ref2} ne '') {
		      $mg2_chr   = $refgenes{$message_ref2}{'chr'};
		      $mg2_start = $refgenes{$message_ref2}{'start'};
		      $mg2_end   = $refgenes{$message_ref2}{'end'};
		      $mg2_strand= $refgenes{$message_ref2}{'strand'};
		  } elsif ($loc{$message2} ne '') {
		      $mg2_chr    = $loc{$message2}{'chr'};
		      $mg2_start  = $loc{$message2}{'start'};
		      $mg2_end    = $loc{$message2}{'end'};
		      $mg2_strand = $loc{$message2}{'strand'};
		  } elsif (exists $name2ens{$message2}) {
		      my $ensemble_id2 = $name2ens{$message2};
		      if ($gene{$ensemble_id2} ne '') {
			  $mg2_chr   = $gene{$ensemble_id2}{'chr'};
			  $mg2_start = $gene{$ensemble_id2}{'start'};
			  $mg2_end   = $gene{$ensemble_id2}{'end'};
			  $mg2_strand= $gene{$ensemble_id2}{'strand'};
		      }
		  }

		  if ($mg1_chr ne '' and $mg2_chr ne '') {
		      $for_encompass{$transcript_name}{'chr1'}   = $mg1_chr;
		      $for_encompass{$transcript_name}{'start1'} = $mg1_start;
		      $for_encompass{$transcript_name}{'end1'}   = $mg1_end;
		      $for_encompass{$transcript_name}{'chr2'}   = $mg2_chr;
		      $for_encompass{$transcript_name}{'start2'} = $mg2_start;
		      $for_encompass{$transcript_name}{'end2'}   = $mg2_end;
		  } #if coordinates are found
	      } #refseq mapping , so check messenges


	      $blat1 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
	      my $blat1_qs  = $1;
	      my $blat1_qe  = $2;
	      my $blat1_st  = $3;
	      my $blat1_chr = $4;
	      my $blat1_ts  = $5;
	      my $blat1_te  = $6;
	      $blat2 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
	      my $blat2_qs  = $1;
	      my $blat2_qe  = $2;
	      my $blat2_st  = $3;
	      my $blat2_chr = $4;
	      my $blat2_ts  = $5;
	      my $blat2_te  = $6;


	      #for repeat masker to get the breakpoint genomic positions (also for genome Blat to get genename)
	      if ($blat1_qs < $blat2_qs) {
		  if ($blat1_st eq '+' and $blat2_st eq '+') {
		      $for_rm{$blat1_chr}{$blat1_te}{$transcript_name} = '+,1';
		      $for_rm{$blat2_chr}{$blat2_ts}{$transcript_name} = '+,2';
		  }
		  if ($blat1_st eq '+' and $blat2_st eq '-') {
		      $for_rm{$blat1_chr}{$blat1_te}{$transcript_name} = '+,1';
		      $for_rm{$blat2_chr}{$blat2_te}{$transcript_name} = '-,2';
		  }
		  if ($blat1_st eq '-' and $blat2_st eq '+') {
		      $for_rm{$blat1_chr}{$blat1_ts}{$transcript_name} = '-,1';
		      $for_rm{$blat2_chr}{$blat2_ts}{$transcript_name} = '+,2';
		  }
		  if ($blat1_st eq '-' and $blat2_st eq '-') {
		      $for_rm{$blat1_chr}{$blat1_ts}{$transcript_name} = '-,1';
		      $for_rm{$blat2_chr}{$blat2_te}{$transcript_name} = '-,2';
		  }
	      }
	      if ($blat1_qs > $blat2_qs) {
		  if ($blat1_st eq '+' and $blat2_st eq '+') {
		      $for_rm{$blat1_chr}{$blat1_ts}{$transcript_name} = '+,2';
		      $for_rm{$blat2_chr}{$blat2_te}{$transcript_name} = '+,1';
		  }
		  if ($blat1_st eq '+' and $blat2_st eq '-') {
		      $for_rm{$blat1_chr}{$blat1_ts}{$transcript_name} = '+,2';
		      $for_rm{$blat2_chr}{$blat2_ts}{$transcript_name} = '-,1';
		  }
		  if ($blat1_st eq '-' and $blat2_st eq '+') {
		      $for_rm{$blat1_chr}{$blat1_te}{$transcript_name} = '-,2';
		      $for_rm{$blat2_chr}{$blat2_te}{$transcript_name} = '+,1';
		  }
		  if ($blat1_st eq '-' and $blat2_st eq '-') {
		      $for_rm{$blat1_chr}{$blat1_te}{$transcript_name} = '-,2';
		      $for_rm{$blat2_chr}{$blat2_ts}{$transcript_name} = '-,1';
		  }
	      }

	      my ($range_s, $range_e);
	      if ($bp_s > $bp_e) {
		  $range_s = $bp_e; $range_e = $bp_s;
	      }
	      if ($bp_s <= $bp_e) {
		  $range_s = $bp_s; $range_e = $bp_e;
	      }

	      my $flag;
	      if ($blat1_chr eq $blat2_chr) {
		  if (abs($blat2_ts - $blat1_ts) < 300000) {
		      $flag = 'Read_Th';
		  } else {
		      $flag = 'Intra_C';
		  }
	      } else {
		  $flag = 'Inter_C';
	      }

	      my %newinfo = ('con'=>$confidence, 'length'=>$length, 'gene1'=>$part1, 'gene2'=>$part2, 'ron'=>$ron, 'ori'=>$orientation, 'bps'=>$range_s, 'bpe'=>$range_e, 'blat1'=>$blat1, 'blat2'=>$blat2, 'flag'=>$flag, 'cScore'=>$cScore);

	      $coverage{$transcript_name}{'info'} = \%newinfo;

	  } #if there is no information for this transcript
      } #for each transcript bp

      #store read mapping information

      #$read =~ /(.+?)[\/\s]([12])/;
      my $read_root = $read;
      my $read_end;

      if ($flag & 64) {$read_end = '1';} #the read is the first in the pair
      elsif ($flag & 128){$read_end = '2';} #the read is the second in the pair
      else {
         print STDERR "error: there is a read alignment ( $read ) seems not to be any mate end!!!!!!\n";
         exit 22;
      }

      #$start += 1;

      foreach my $transcript_bp (@transcript) {
        my $transcript_name = $transcript_bp->{'name'};
        $coverage{$transcript_name}{$read_root}{$read_end}{$start} = $_;
      }

   } #pair type
}

close IN;


#do the repeatmasker & selfChain mask for the breakpoint
my %bprepeat;
my %bpselfChain;
foreach my $chr (sort keys %for_rm) {
    foreach my $pos (sort {$a<=>$b} keys %{$for_rm{$chr}}){
	my $rmflag = &repeatmask($chr, $pos);
        my $scflag = &selfChainMask($chr, $pos);

        if ($rmflag == 1){
	    foreach my $transcript_name (keys %{$for_rm{$chr}{$pos}}){
	       if (exists $bprepeat{$transcript_name}){
		   if ($for_rm{$chr}{$pos}{$transcript_name} =~ /1$/){
                       $bprepeat{$transcript_name} = 'R'.$bprepeat{$transcript_name};
		   } else {
                       $bprepeat{$transcript_name} = $bprepeat{$transcript_name}.'R';
		   }
	       } else {
		   $bprepeat{$transcript_name} = 'R';
	       }
            }
        } #repeat
        if ($rmflag == 0){
            foreach my $transcript_name (keys %{$for_rm{$chr}{$pos}}){
               if (exists $bprepeat{$transcript_name}){
		   if ($for_rm{$chr}{$pos}{$transcript_name} =~ /1$/){
                       $bprepeat{$transcript_name} = 'N'.$bprepeat{$transcript_name};
		   } else {
                       $bprepeat{$transcript_name} = $bprepeat{$transcript_name}.'N';
		   }
	       } else {
		   $bprepeat{$transcript_name} = 'N';
	       }
            }
        } #not repeat

        if ($scflag == 1){
	    foreach my $transcript_name (keys %{$for_rm{$chr}{$pos}}){
	       if (exists $bpselfChain{$transcript_name}){
		   if ($for_rm{$chr}{$pos}{$transcript_name} =~ /1$/){
                       $bpselfChain{$transcript_name} = 'C'.$bpselfChain{$transcript_name};
		   } else {
                       $bpselfChain{$transcript_name} = $bpselfChain{$transcript_name}.'C';
		   }
	       } else {
		   $bpselfChain{$transcript_name} = 'C';
	       } 
            }
        } #selfChain
        if ($scflag == 0){
            foreach my $transcript_name (keys %{$for_rm{$chr}{$pos}}){
               if (exists $bpselfChain{$transcript_name}){
		   if ($for_rm{$chr}{$pos}{$transcript_name} =~ /1$/){
                       $bpselfChain{$transcript_name} = 'N'.$bpselfChain{$transcript_name};
		   } else {
                       $bpselfChain{$transcript_name} = $bpselfChain{$transcript_name}.'N';
		   }
	       } else {
		   $bpselfChain{$transcript_name} = 'N';
	       }
            }
        } #not selfChain

    } # for each pos
} #repeat mask & selfChain mask


#find the gene name, orientation and strands if genome Blat
if ($genomeBlatPred ne 'SRP'){
    foreach my $chr (sort keys %for_rm) {
	foreach my $pos (sort {$a<=>$b} keys %{$for_rm{$chr}}) {

	    foreach my $transcript_name (keys %{$for_rm{$chr}{$pos}}) {
                my ($strand, $oot) = split /\,/, $for_rm{$chr}{$pos}{$transcript_name};
                my ($geneName, $orient) = &breakpointInGene($chr, $pos, $strand);
                if ($oot == 1) {
                   $coverage{$transcript_name}{'info'}->{'ori'} = $orient.($coverage{$transcript_name}{'info'}->{'ori'});
                   $coverage{$transcript_name}{'info'}->{'gene1'} = $geneName;
                } elsif ($oot == 2){
                   $coverage{$transcript_name}{'info'}->{'ori'} = ($coverage{$transcript_name}{'info'}->{'ori'}).$orient;
                   $coverage{$transcript_name}{'info'}->{'gene2'} = $geneName;
                }
	    } #each transcript name

	}
    } #gene name finder

    my %name_pairs;
    foreach my $transcript_name (keys %coverage) {
       my $gene1 = $coverage{$transcript_name}{'info'}->{'gene1'};
       my $gene2 = $coverage{$transcript_name}{'info'}->{'gene2'};
       my $blat1 = $coverage{$transcript_name}{'info'}->{'blat1'};
       my $blat2 = $coverage{$transcript_name}{'info'}->{'blat2'};
       my $name_pair = $gene1.'+'.$gene2;
       $name_pairs{$name_pair} = '';

       $blat1 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
       my $blat1_st  = $3;
       my $blat1_chr = $4;
       my $blat1_ts  = $5;
       my $blat1_te  = $6;
       $blat2 =~ /^(\d+)\-(\d+)\(([+-])\)(chr\w+)\:(\d+)\-(\d+)$/;
       my $blat2_st  = $3;
       my $blat2_chr = $4;
       my $blat2_ts  = $5;
       my $blat2_te  = $6;

       my ($ensembl1, $ensembl2);
       if ( $gene1 =~ /^(.+?)\(/ ) {
         $ensembl1 = $1;
       } elsif ($gene1 == 'IGR') {
         $ensembl1 = 'IGR';
       } else {
         print STDERR "$gene1 \($transcript_name\) is not valid.\n";
       }
       if ( $gene2 =~ /^(.+?)\(/ ) {
         $ensembl2 = $1;
       } elsif ($gene2 == 'IGR') {
         $ensembl2 = 'IGR';
       } else {
         print STDERR "$gene2 \($transcript_name\) is not valid.\n";
       }

       if ($ensembl1 eq $ensembl2) { #not a fusion
          delete($coverage{$transcript_name});
          next;
       }

       #also should do for emcompassing
       if ( $ensembl1 ne 'IGR' ){
         $for_encompass{$transcript_name}{'chr1'}   = $gene{$ensembl1}{'chr'};
         if ($blat1_st eq '+'){
           $for_encompass{$transcript_name}{'start1'} = $gene{$ensembl1}{'start'};
           $for_encompass{$transcript_name}{'end1'}   = $blat1_te;
         } else {
           $for_encompass{$transcript_name}{'start1'} = $blat1_ts;
           $for_encompass{$transcript_name}{'end1'}   = $gene{$ensembl1}{'end'};
         }
       } else {
         $for_encompass{$transcript_name}{'chr1'}   = $blat1_chr;
         if ($blat1_st eq '+'){
           $for_encompass{$transcript_name}{'start1'} = $blat1_ts - 10000;
           $for_encompass{$transcript_name}{'end1'}   = $blat1_te;
         } else {
           $for_encompass{$transcript_name}{'start1'} = $blat1_ts;
           $for_encompass{$transcript_name}{'end1'}   = $blat1_te + 10000;
         }
       }
       if ( $ensembl2 ne 'IGR' ) {
         $for_encompass{$transcript_name}{'chr2'}   = $gene{$ensembl2}{'chr'};
         if ($blat2_st eq '+') {
           $for_encompass{$transcript_name}{'start2'} = $blat2_ts;
           $for_encompass{$transcript_name}{'end2'}   = $gene{$ensembl2}{'end'};
         } else {
           $for_encompass{$transcript_name}{'start2'} = $gene{$ensembl2}{'start'};
           $for_encompass{$transcript_name}{'end2'}   = $blat2_te;
         }
       } else {
          $for_encompass{$transcript_name}{'chr2'}   = $blat2_chr;
          if ($blat2_st eq '+') {
            $for_encompass{$transcript_name}{'start2'} = $blat2_ts;
            $for_encompass{$transcript_name}{'end2'}   = $blat2_te + 10000;
          } else {
            $for_encompass{$transcript_name}{'start2'} = $blat2_ts - 10000;
            $for_encompass{$transcript_name}{'end2'}   = $blat2_te;
          }
       }
    }
    foreach my $transcript_name (keys %coverage){
       my $gene1 = $coverage{$transcript_name}{'info'}->{'gene1'};
       my $gene2 = $coverage{$transcript_name}{'info'}->{'gene2'};
       my $name_pair = $gene2.'+'.$gene1;
       if (exists $name_pairs{$name_pair}){
          $coverage{$transcript_name}{'info'}->{'ron'} = 'st_both';
       } else {
          $coverage{$transcript_name}{'info'}->{'ron'} = 'st_sing';
       }
    }
} #define missing info for genome blat


#sort out the encompassing chr and positions
my %encoposA;
my %encoposB;
foreach my $transcript_name (keys %for_encompass) {

       my $chr1    = $for_encompass{$transcript_name}{'chr1'};
       my $start1  = $for_encompass{$transcript_name}{'start1'};
       my $end1    = $for_encompass{$transcript_name}{'end1'};
       my $chr2    = $for_encompass{$transcript_name}{'chr2'};
       my $start2  = $for_encompass{$transcript_name}{'start2'};
       my $end2    = $for_encompass{$transcript_name}{'end2'};

       push (@{$encoposA{$chr1}{$start1}{$end1}}, [$transcript_name, $start2, $end2]);
       push (@{$encoposB{$chr2}{$start2}{$end2}}, [$transcript_name, $start1, $end1]);

}

my @mergehits = split(/,/, $accepthits);

foreach my $accepthitsC (@mergehits) {

#generate the coverage encompassing first
my (%buffer1, %buffer2);

open AH, "samtools view $accepthitsC |";

my ($old_chrA, $old_chrB);
my (@ssA, @ssB);
my ($ptrA, $ptrB);
while ( <AH> ) {

    next if /^@/;                #ignore comments
    next if ( $_ !~ /NH\:i\:1\t/ and $_ !~ /NH\:i\:1$/ ); #ignore multiple mappable reads
    chomp;

    my ($Qname, $FLAG, $Rname, $Pos, $MAPQ, $CIGAR, $mateRname, $matePos, $ISIZE, $seq, $qual, @tag) = split /\t/;
    next if ($mateRname eq '*');
    next if ($mateRname eq '=' and abs($matePos - $Pos) < 230000);         #ignore correct pair
    if ($CIGAR =~ /^(\d+)[SH]/){
      next if ($1 > 8);
    } elsif ($CIGAR =~ /(\d+)[SH]$/) {
      next if ($1 > 8);
    }

    if ($Rname ne $old_chrA) {
       @ssA  = sort {$a <=> $b} keys %{$encoposA{$Rname}};
       $ptrA = 0;
    }
    if ($Rname ne $old_chrB) {
       @ssB  = sort {$a <=> $b} keys %{$encoposB{$Rname}};
       $ptrB = 0;
    }

    my @endsA = keys %{$encoposA{$Rname}{$ssA[$ptrA]}};
    while (($ptrA<=$#ssA) and ($endsA[0] < $Pos)){
      $ptrA++;
      @endsA = keys %{$encoposA{$Rname}{$ssA[$ptrA]}};
    }
    if ($ssA[$ptrA] <= $Pos){
       foreach my $array_ref (@{$encoposA{$Rname}{$ssA[$ptrA]}{$endsA[0]}}){
          if ($matePos >= $array_ref->[1] and $matePos <= $array_ref->[2]){
             $buffer1{$array_ref->[0]}{$Qname} = $_;
          }
       }
    }

    my @endsB = keys %{$encoposB{$Rname}{$ssB[$ptrB]}};
    while (($ptrB<=$#ssB) and ($endsB[0] < $Pos)){
       $ptrB++;
       @endsB = keys %{$encoposB{$Rname}{$ssB[$ptrB]}};
    }
    if ($ssB[$ptrB] <= $Pos){
       foreach my $array_ref (@{$encoposB{$Rname}{$ssB[$ptrB]}{$endsB[0]}}){
          if ($matePos >= $array_ref->[1] and $matePos <= $array_ref->[2]){
             $buffer2{$array_ref->[0]}{$Qname} = $_;
          }
       }
    }

    $old_chrA = $Rname;
    $old_chrB = $Rname;

}
close AH;

foreach my $transcript_name (keys %buffer1){
  foreach my $read (keys %{$buffer1{$transcript_name}}){
     if (exists $buffer2{$transcript_name}{$read}){
         push (@{$for_encompass{$transcript_name}{'hits'}}, $buffer1{$transcript_name}{$read}, $buffer2{$transcript_name}{$read});
         my @cols = split (/\t/, $buffer1{$transcript_name}{$read});
         my $mapo = join ("_", sort {$a <=> $b} ($cols[3],$cols[7]));
         $for_encompass{$transcript_name}{'cov'}{$mapo}++;
     }
  }
}

} #foreach merged bam

#end of encompassing cov

open ENCOMCOV, ">$encomcov";

foreach my $transcript_name (sort { $coverage{$b}{'info'}->{'con'} <=> $coverage{$a}{'info'}->{'con'} } keys %coverage){

   my $length = $coverage{$transcript_name}{'info'}->{'length'};
   my $gene1 = $coverage{$transcript_name}{'info'}->{'gene1'};
   my $gene2 = $coverage{$transcript_name}{'info'}->{'gene2'};
   my $bps = $coverage{$transcript_name}{'info'}->{'bps'};
   my $bpe = $coverage{$transcript_name}{'info'}->{'bpe'};
   my $blat1 = $coverage{$transcript_name}{'info'}->{'blat1'};
   my $blat2 = $coverage{$transcript_name}{'info'}->{'blat2'};
   my $flag  = $coverage{$transcript_name}{'info'}->{'flag'};
   my $ori   = $coverage{$transcript_name}{'info'}->{'ori'};
   my $ron   = $coverage{$transcript_name}{'info'}->{'ron'};
   my $cScore = $coverage{$transcript_name}{'info'}->{'cScore'};
   my $bpr = $bprepeat{$transcript_name};
   my $bpsc = $bpselfChain{$transcript_name};
   next if $bpr eq 'RR';

   foreach my $read_root (keys %{$coverage{$transcript_name}}) {

      my @pair = keys %{$coverage{$transcript_name}{$read_root}};

      if (scalar(@pair) == 1){ #only one end map
         my @starts = keys %{$coverage{$transcript_name}{$read_root}{$pair[0]}};
         next if (scalar (@starts) >= 2);
         my $start = $starts[0];
         my $end = $start+$read_length-1;
         if ($start < $bps and $end > $bpe) { #spanning
            next if ( ($bps-$start) < $anchorLen or ($end-$bpe) < $anchorLen );    #requiring minimal anchoring length
            $coverage{$transcript_name}{'spanning'}{$start}++;
            push (@{$coverage{$transcript_name}{'reads'}},  $coverage{$transcript_name}{$read_root}{$pair[0]}{$start});
         }
         elsif (($start < $bps and $end > $bps) or ($start >= $bps and $start <= $bpe)) { #other overlapping with bp, not interesting
            #push (@{$coverage{$transcript_name}{'reads'}},  $coverage{$transcript_name}{$read_root}{$pair[0]}{$start});
         }
      }

      if (scalar(@pair) == 2) { #pair mapped
         my @starts1 = keys %{$coverage{$transcript_name}{$read_root}{1}};
         my @starts2 = keys %{$coverage{$transcript_name}{$read_root}{2}};
         next if (scalar (@starts1) >= 2 or scalar (@starts2) >= 2);
         my $start1 = $starts1[0];
         my $start2 = $starts2[0];
         my $end1 = $start1+$read_length-1;
         my $end2 = $start2+$read_length-1;

         my $frag_s;
         my $frag_e;
         if ($start1 <  $start2){$frag_s=$start1;$frag_e=$end2;}
         if ($start1 >= $start2){$frag_s=$start2;$frag_e=$end1;}

         my $flag = 0;

         if ($start1 < $bps and $end1 > $bpe) { # 5'end spanning
           if ( ($bps-$start1) >= $anchorLen and ($end1-$bpe) >= $anchorLen ) {   #minimal anchoring length
             $coverage{$transcript_name}{'spanning'}{$start1}++;
             $flag = 1;
           }
         }
         elsif (($start1 < $bps and $end1 > $bps) or ($start1 >= $bps and $start1 <= $bpe)){ #5'end overlapping
            $flag = 2;
         }

         if ($flag == 0 and ($start2 < $bps and $end2 > $bpe)) { # 3'end spanning
           if ( ($bps-$start2) >= $anchorLen and ($end2-$bpe) >= $anchorLen ) { #minimal anchoring length
             $coverage{$transcript_name}{'spanning'}{$start2}++;
             $flag = 3;
           }
         }
         elsif ($flag == 0 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){ #3'end overlapping
            $flag = 4;
         }

         elsif ($flag == 1 and ($start2 < $bps and $end2 > $bpe)) { # 5' end and 3' end both spanning
            if ( ($bps-$start2) >= $anchorLen and ($end2-$bpe) >= $anchorLen ) { #minimal anchoring length
              $flag = 5;
            }
         }

         elsif ($flag == 1 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){  #3'end overlapping and 5'end spanning
            $flag = 6;
         }

         elsif ($flag == 2 and ($start2 < $bps and $end2 > $bpe)) { # 5' end overlapping and 3' end spanning
            if ( ($bps-$start2) >= $anchorLen and ($end2-$bpe) >= $anchorLen ) { #minimal anchoring length
              $coverage{$transcript_name}{'spanning'}{$start2}++;
              $flag = 7;
            }
         }

         elsif ($flag == 2 and ($start2 < $bps and $end2 > $bps) or ($start2 >= $bps and $start2 <= $bpe)){  #both 5' and 3' end overlapping
            $flag = 8;
         }


         if ($flag =~ /[0248]/ and ($frag_s < $bps and $frag_e > $bpe)){
            if ( ($bps-$frag_s) >= $anchorLen and ($frag_e-$bpe) >= $anchorLen ) { #minimal anchoring length
              $coverage{$transcript_name}{'encompass'}{$frag_s}++;
              $flag = 9;
            }
         }

         if ($flag =~ /[135679]/) {
            push (@{$coverage{$transcript_name}{'reads'}},  $coverage{$transcript_name}{$read_root}{1}{$start1});
            push (@{$coverage{$transcript_name}{'reads'}},  $coverage{$transcript_name}{$read_root}{2}{$start2});
         }
      }
   }


   my @spans= keys %{$coverage{$transcript_name}{'spanning'}};
   my $span = scalar(@spans); #pile-up reads only count once
   my $spanAll; #span redundant
   my $spanScore = $span;
   foreach my $spanStart (@spans) {
     $spanAll += $coverage{$transcript_name}{'spanning'}{$spanStart};
     my $spanEnd = $spanStart + $read_length - 1;
     my $balancePen = abs(($bps-$spanStart)-($spanEnd-$bpe))/(($bps-$spanStart)+($spanEnd-$bpe));
     $spanScore -= $balancePen;
   }
   $spanScore = sprintf("%.1f", $spanScore);
   my $enco = scalar (keys %{$coverage{$transcript_name}{'encompass'}}); #pile-up frags only count once
   my $cov  = $span.'('.$span.'+'.$enco.')';
   my $real_enco = scalar(keys (%{$for_encompass{$transcript_name}{'cov'}}));
   my $all  = $span + $real_enco;

   next if ($all <= 1);

   my $newtitle = join("\t", $transcript_name, $cScore, $length, $gene1.'-'.$gene2, $ori, $bps.'..'.$bpe, $bpr, $bpsc, $flag, $ron, $blat1, $blat2, $cov, $real_enco, $all, $spanAll, $spanScore);
   print "#$newtitle\n";
   print ENCOMCOV "#$newtitle\n";
   foreach my $razers (@{$coverage{$transcript_name}{'reads'}}){
      print "$razers\n";
   }
   foreach my $hits (@{$for_encompass{$transcript_name}{'hits'}}){
      print ENCOMCOV "$hits\n";
   }
}

close ENCOMCOV;


exit;

sub repeatmask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr){
	@rs = sort {$a <=> $b} keys %{$repeatmask{$chr}};
	$ptr = 0;
    }
    while (($ptr<=$#rs) and ($repeatmask{$chr}{$rs[$ptr]} < $coor)){
	$ptr++;
    }
    if ($rs[$ptr] <= $coor){
	$flag = 1;
    }
    $old_chr = $chr;
    return $flag;
}

sub selfChainMask {
    my ($chr, $coor) = @_;
    my $flag = 0;
    if ($chr ne $old_chr_selfChain){
	@rs_selfChain = sort {$a <=> $b} keys %{$selfChain{$chr}};
	$ptr_selfChain = 0;
    }
    while (($ptr_selfChain <= $#rs_selfChain) and ($selfChain{$chr}{$rs_selfChain[$ptr_selfChain]} < $coor)){
	$ptr_selfChain++;
    }
    if ($rs_selfChain[$ptr_selfChain] <= $coor){
	$flag = 1;
    }
    $old_chr_selfChain = $chr;
    return $flag;
}

sub breakpointInGene {
   my ($chr, $coor, $bstrand) = @_;
   my @return_info;
   my $geneName = 'IGR';
   my $orientation = '-';

   if ($chr ne $old_chr2) {
      @rs2 = sort {$a <=> $b} keys %{$gene_annotation{$chr}};
      $ptr2 = 0;
   }

   while ( $ptr2 <= $#rs2 ) {
     my @geneEnds = sort {$a <=> $b} keys(%{$gene_annotation{$chr}{$rs2[$ptr2]}});
     if ($geneEnds[$#geneEnds] < $coor) { #the farest end is beyond the coor
       $ptr2++;
       next;
     } else {
       last;
     }
   }

   my $off = 0;
   my $geneMode = "nc"; #mRNA-exon > gene-exon > nc-exon > nc
   my $max_glength = 0; #max gene length
   my $max_mlength = 0; #max rna length
   while ( ($ptr2+$off) <= $#rs2 and $rs2[$ptr2+$off] <= $coor ) { #same gene start

      my @geneEnds = sort {$a <=> $b} keys(%{$gene_annotation{$chr}{$rs2[$ptr2+$off]}});

      foreach my $geneEnd (@geneEnds) {
         my ($gname, $gstrand, $glength) = split /\,/, $gene_annotation{$chr}{$rs2[$ptr2+$off]}{$geneEnd};

         $gname =~ /^(.+?)\(.+?\)$/;
         my $ensembl = $1;
         my $mRNA_start = -1;
         my $mRNA_end   = -1;
         my $mRNA_size  = -1;
         if (exists $mRNA_region{$ensembl}) {
           $mRNA_start = $mRNA_region{$ensembl}->[0];
           $mRNA_end = $mRNA_region{$ensembl}->[1];
           $mRNA_size = $mRNA_end - $mRNA_start;
         }


         if ( $geneEnd >= $coor ) {   #seems overlapping

           my $exonInclFlag = &exonIncl(\%{$gene2exonMerged{$ensembl}}, $coor);

           if ($mRNA_start != -1) { # mRNA in this gene

             if ($mRNA_start <= $coor and $coor <= $mRNA_end) { #overlapping mRNA

               if ($exonInclFlag == 1) {            #BEST
                 if ($mRNA_size > $max_mlength) {
                   $geneName = $gname;
                   if ($bstrand eq $gstrand) {
                     $orientation = '->';
                   } else {
                     $orientation = '<-';
                   }
                   $max_mlength = $mRNA_size;
                 }
                 $max_glength = $glength if ($glength > $max_glength);
                 $geneMode = 'mRNA-exon';

               } else { #overlapping with mRNA but not with exon!?
                 if ($geneMode eq 'nc' and $glength > $max_glength) {   #treat as BAD:NC
                   $geneName = $gname;
                   if ($bstrand eq $gstrand) {
                     $orientation = '->';
                   } else {
                     $orientation = '<-';
                   }
                   $max_glength = $glength;
                 }
               }

             } else {           #not overlapping mRNA
                if ($exonInclFlag == 1) {
                  if ($geneMode eq 'nc' or $geneMode eq 'nc-exon' or ($geneMode eq 'gene-exon' and $glength > $max_glength)) {
                    $geneName = $gname;
                    if ($bstrand eq $gstrand) {
                      $orientation = '->';
                    } else {
                      $orientation = '<-';
                    }
                    $max_glength = $glength;
                    $geneMode = 'gene-exon';
                  }
                } else {
                  if ($geneMode eq 'nc' and $glength > $max_glength) {
                    $geneName = $gname;
                    if ($bstrand eq $gstrand) {
                      $orientation = '->';
                    } else {
                      $orientation = '<-';
                    }
                    $max_glength = $glength;
                  }
                }
             }

           } else { #non coding

             if ($exonInclFlag == 1) {   #non coding exon

               if ($geneMode eq 'nc' or ($geneMode eq 'nc-exon' and $glength > $max_glength)) {
                 $geneName = $gname;
                 if ($bstrand eq $gstrand) {
                   $orientation = '->';
                 } else {
                   $orientation = '<-';
                 }
                 $max_glength = $glength;
                 $geneMode = "nc-exon";
               }

             } else {       # in this nc gene but not with exon
               if ($geneMode eq 'nc' and $glength > $max_glength) {
                 $geneName = $gname;
                 if ($bstrand eq $gstrand) {
                   $orientation = '->';
                 } else {
                   $orientation = '<-';
                 }
                 $max_glength = $glength;
               }
             }

           } #non-coding

         } #seems overlapping
      } #foreach gene end

      $off++;

   } #while bp coor is greater than the gene start

   push @return_info, $geneName;
   push @return_info, $orientation;
   $old_chr2 = $chr;
   return @return_info;
}

sub exonIncl {

  my $eflag = 0;
  my ($exons, $coor) = @_;
  my @ers = sort {$a <=> $b} keys %{$exons};
  my $eptr = 0;

  while (($eptr<=$#ers) and ($exons->{$ers[$eptr]} < $coor)){
    $eptr++;
  }

  if ($ers[$eptr] <= $coor){
    $eflag = 1;
  }

  return $eflag;

}
