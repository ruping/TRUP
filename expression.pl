#!/usr/bin/perl

#deprecated script for sam format
#TODO: determine the expression level of each transcripts (the splitly mapped reads are not considered)
#      from the mapping result of TopHat
#

use Getopt::Long;
my $umr_file;
my $annot_file;
my $read_length = 80;

GetOptions(
            "umr|u=s"=>\$umr_file,
	    "anno|a=s"=>\$annot_file,
	    "length|l=i"=>\$read_length,
            "help|h" => sub{
	                     print "usage: $0 [options]\n\nOptions:\n\t--umr\t\tthe uniquely mapping file in gff format\n";
                             print "\t--anno\t\tthe gene annotation file in gff\n";
                             print "\t--length\tthe read length\n";
			     print "\t--help\t\tprint this help message\n";
                             exit 0;
                       }
          );

## prepare indices
my %umr_chr_start = ();

my %chr;
for (1..22){
  my $c = 'chr'.$_;
  $chr{$c} = '';
}
$chr{'chrX'}  = '';
$chr{'chrY'}  = '';
$chr{'chrMT'} = '';


open UMR, $umr_file;
my $old_chr = "";
my $umr_N = 0;
while(<UMR>){
    ## ignore comments;
    next if (/^[@#]/);
    chomp;
    my @cols = split(/\t/, $_);
    my $chr = $cols[0];
    $chr = 'chr'.$chr unless ($chr =~ /^chr/);
    $chr .= 'T' if ($chr eq 'chrM');
    
       if ($chr ne $old_chr){
          $umr_chr_start{$cols[0]} = tell UMR;
       }
       $umr_N++;
       $old_chr = $chr;

}
close UMR;
print STDERR "made umr index\n";
print STDERR "total number of tags = $umr_N\n";



open ANNOT, $annot_file;
$old_chr = "";
my %tags = ();
my %exon2tags = ();
my %exon2starts = ();
my %exon2transcript = ();
my %transcript2tags = ();
my %transcript2starts = ();
my %transcriptL = ();
my %transcript2gene = ();
my %gene2name = ();

while(<ANNOT>){
    next if /^#/; ## ignore comments
    chomp;
    my ($chr, $source, $feature, $start, $end, $score, $strand, $phase, @attributes) = split(/[\t\;]/, $_);
    $chr = 'chr'.$chr unless ($chr =~ /^chr/);
    $chr .= 'T' if ($chr eq 'chrM');

    next if (!exists $chr{$chr}); #skip other chromosomes

    if ($chr ne $old_chr){
	#print STDERR "loading data for $chr\n";
	load_tags($chr);
    }

    if ($feature eq 'gene'){
	my $id, $name;
	foreach my $attribute (@attributes){
	    if ($attribute =~ /ID\=(.+)$/){
		$id = $1;
	    }
	    if ($attribute =~ /Name\=(.+)$/){
		$name = $1;
	    }
	    
	    
	}
	$gene2name{$id} = $name;
    }
    if ($feature eq 'mRNA' or $feature eq 'transcript'){
	my $id, $gene;
	foreach my $attribute (@attributes){
	    if ($attribute =~ /ID\=(.+)$/){
		$id = $1;
	    }
	    if ($attribute =~ /Parent\=(.+)$/){
		$gene = $1;
	    }
	    
	}
	$transcript2gene{$id} = $gene;
    }


    if ($feature eq 'exon'){
	my $id, @transcripts;
	foreach my $attribute (@attributes){
	    if ($attribute =~ /ID\=(.+)$/){
		$id = $1;
	    }
	    if ($attribute =~ /Parent\=(.+)$/){
		@transcripts = split(",", $1);
	    }
	    
	}
	@{$exon2transcript{$id}} = @transcripts;
	$exon2tags{$id} = 0;
	my $L = $end - $start + 1;
	for my $pos ($start .. $end){
	    $exon2tags{$id} += $tags{$pos};     #for the number of reads
	    $exon2starts{$id}++ if (exists $tags{$pos});  #for the number of start positions
	}
	
	foreach my $transcript (@transcripts){
            
	    $transcript2tags{$transcript} += $exon2tags{$id}; #for the number of reads
	    $transcript2starts{$transcript} += $exon2starts{$id}; #for the number of start positions
	    $transcriptL{$transcript} += $L;
            #if ($transcript eq 'ENST00000401738'){
            #    print STDERR "$transcript2tags{$transcript}\t$transcriptL{$transcript}\n";
            #}
	}
    }

    $old_chr = $chr;
    
}
close ANNOT;

my %gene2expression = ();
my %gene2starts = ();
my %geneL = ();
foreach my $transcript (keys %transcript2tags){

    my $gene = $transcript2gene{$transcript};
    if ($gene eq ''){
        print STDERR "$transcript\n";
    }

    $gene2expression{$gene} += $transcript2tags{$transcript};
    $gene2starts{$gene} += $transcript2starts{$transcript};
    $geneL{$gene} += $transcriptL{$transcript};

}

my $total_reads;
my $total_length_tr;
my $total_g = scalar(keys %gene2expression);

print "# $umr_N\n";

foreach my $gene (sort keys %gene2expression){

	my $R = $gene2expression{$gene} * 1e9 / ($umr_N * $geneL{$gene});  #Reads per M reads per kilo base
	my $S = $gene2starts{$gene}/$geneL{$gene};
	if ($R > 0){
            print "$gene\t$gene2name{$gene}\t$geneL{$gene}\t$gene2expression{$gene}\t$R\t$S\n";
	}

    $total_reads += $gene2expression{$gene};
    $total_length_tr += $geneL{$gene};

}

my $ratio1=$total_reads/$total_g;
my $ratio2=$total_reads*1000/$total_length_tr;
print STDERR "per gene:$ratio1\nper 1000 base transcripts:$ratio2\n";

sub load_tags{
    my $chr = shift;

    %tags = ();

    open UMR, $umr_file;
    seek(UMR, $umr_chr_start{$chr}, 0);
    while(<UMR>){
	my ($c, $source, $type, $start, $end, $score, $strand, $phase, $tag) = split(/\t/, $_);
	$c = 'chr'.$c unless ($c =~ /^chr/);
        $c .= 'T' if ($c eq 'chrM');
		
	  last if ($c ne $chr);

	  $tags{$start}++;

    }
    close UMR;
}

