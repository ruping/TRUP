#!/usr/bin/perl
#$ -cwd

## Use as  perl /ifs/scratch/c2b2/ngs_lab/sz2317/runs/ALI_GHARAVI/convert_vcf_exomeannotation-all-samples.pl infilename original.vcf annovar/vcfheader_annovar.txt 

use File::Basename;

$fin = $ARGV[0];
$vcfout = $fin;
$header = $ARGV[1];	#Original VCF files
$header_ann = $ARGV[2];     #header specific describing Annovar steps

if (not $fin and not $header )
{ print "1st argument -annovar input file \n 2nd argument -original vcf file\n"; exit;}


print "\nOpen _annovar file";

open(vcf, "<" . $fin) || die("Could not open $fin file!");
open(hdr, "<" . $header) || die("Could not open original VCF file $header !");
open(hdr_ann, "<" . $header_ann) || die("Could not open file $header_ann !");
open(out, ">" . $vcfout . ".vcf" )|| die("Could not open $vcfout.vcf  file!");

print  "\nBeginning convrsion\n";
print "$fin \t $header \n ";

##Header of VCF shud look iike
###INFO=<ID=1KG,Number=1,Type=Float,Description="1KG Membership">
###INFO=<ID=dbSNP,Number=0,Type=Flag,Description="dbSNP132 Membership">
###INFO=<ID=ESP5400,Number=1,Type=Float,Description="ESP5400 Membership">
#???????????????????????????????????????##INFO=<ID=HapMapV3,Number=0,Type=Flag,Description="HapMapV3 Membership">
#???????????????????????????????????????##INFO=<ID=HEVS
#
###INFO=<ID=AB,Number=1,Type=Float,Description="Allele Balance for hets (ref/(ref+alt))">
###INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
###INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
###INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
###INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
###INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
###INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">
###INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
###INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
###INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
###INFO=<ID=GenericAnnotation,Number=1,Type=Integer,Description="For each variant in the 'variants' ROD, finds all entries in the other -B files that overlap the variant's position.">
###INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">
###INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
###INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
###INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
###INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
###INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
###INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
###INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
###INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
#
#
#
### fields of OUTPUT  vcf file has 10 fields as follows
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  136-04i ...
# 1       871146  .       C       T       923.51  PASS    AB=0.492;AC=1;AF=0.042;AN=24;BaseQRankSum=1.084;DP=1029;Dels=0.00;FS=1.545;HRun=2;HaplotypeScore=5.5593;InbreedingCoeff=-0.0435;MQ=57.01;MQ0=0;MQRankSum=-0.674;QD=15.65;ReadPosRankSum=1.003;SB=-332.61;refseq.chr=1;refseq.codingCoordStr=c.306-6;refseq.end=871146;refseq.haplotypeAlternate=T;refseq.haplotypeReference=C;refseq.inCodingRegion=false;refseq.name=NM_152486;refseq.name2=SAMD11;refseq.positionType=intron;refseq.spliceDist=-6;refseq.spliceInfo=splice-acceptor_-6;refseq.start=871146;refseq.transcriptStrand=+  GT:AD:DP:GQ:PL   0/0:117,0:117:99:0,343,4408 ... 

# header for inpu annovar file is
# 0       Func
# 1       Gene
# 2       ExonicFunc
# 3       AAChange
# 4       Conserved 	outfile.hg19_phastConsElements46way"
# 5       SegDup	_genomicSuperDups"
# 6       ESP5400_ALL
# 7       1000g2012feb_ALL
# 8       dbSNP132
# 9       AVSIFT
# 10      LJB_PhyloP
# 11      LJB_PhyloP_Pred
# 12      LJB_SIFT
# 13      LJB_SIFT_Pred
# 14      LJB_PolyPhen2
# 15      LJB_PolyPhen2_Pred
# 16      LJB_LRT
# 17      LJB_LRT_Pred
# 18      LRT_MutationTaster , now LJB_MutationTaster
# 19      LRT_MutationTaster_Pred , now LJB_MutationTaster_Pred
# 20      LJB_GERP++
# 21      Chr
# 22      Start
# 23      End
# 24      Ref
# 25      Obs
# 26,27,28,29,30,31,32,33,34,35      Otherinfo ( chr,pos,ID,ref,alt,qual,filter,info,format,data)
#

my $chrline="";
my $vcf_version=0;
while (<hdr>) { 
        if ($_ =~ m/^##/) { 
		print out $_; 
		if ($vcf_version==0){
			$vcf_version= $vcf_version + 1;
			while(<hdr_ann>){ print out $_; }
		}
	}
	else { $chrline=$_ ; last; print "Parsed header ..\n"; }
}

print out "\#\#INFO=<ID=function,Number=1,Type=String,Description=\"Function\">\n" ;
print out "\#\#INFO=<ID=functionalClass,Number=1,Type=String,Description=\"FunctionalClass\">\n";
print out "\#\#INFO=<ID=geneName,Number=1,Type=String,Description=\"GeneName\">\n";

my $sampleCount =1;
my @samplelist=();
my $buffer ="";
my $lineCount=0;
while(<vcf>) {
	if ($_ =~ m/^Func/) 
	{
# NewAnnovar:Func    Gene    ExonicFunc      AAChange        Conserved       SegDup  ESP5400_ALL     1000g2012feb_ALL        dbSNP135        AVSIFT  LJB_PhyloP      LJB_PhyloP_Pred LJB_SIFT        LJB_SIFT_Pred   LJB_PolyPhen2   LJB_PolyPhen2_Pred  LJB_LRT  LJB_LRT_Pred    LJB_MutationTaster      LJB_MutationTaster_Pred LJB_GERP++      Chr     Start   End     Ref     Obs     Otherinfo
		print out $chrline;
	}
	else 
	{
		chomp;
		$_ =~ s/"//g ;
		$buffer="";		
		my @item = split("\t");
		$buffer.="$item[26]\t$item[27]\t";	#Chr, Pos
		$buffer.=  ($item[8]) ? "$item[8]\t" : ".\t";	# ID
		$buffer.="$item[29]\t$item[30]\t$item[31]\t$item[32]\t";#Ref, Alt,Qual, Filter
		## Update INFO with annotations
		my $info="$item[33]";
		
		if ($item[6]) { $info.=";ESP5400.score=$item[6]";}
                if ($item[7]) { $info.=";1KG.score=$item[7]";}
                if ($item[8]) { $info.=";dbSNP";}
                if ($item[0]) { $info.=";function=$item[0]";}	
		if ($item[2]) { $item[2] =~ s/\ //; $info.=";functionalClass=$item[2]"; }
		if ($item[1]) { $info.=";geneName=$item[1]";}
		if ($item[3]) { $info.=";AAChange=$item[3]";}
		if ($item[4]) { $item[4] =~ s/Name=/mce.name:/ ; $item[4] =~ s/"// ;  $info.=";mce.score=$item[4]";}
		if ($item[5]) { $info.=";segdup.score=$item[5]";}
		if ($item[9]) { $info.=";avsift.score=$item[9]";}
		if ($item[10]) { $info.=";ljb_phylop.score=$item[10]";}
		if ($item[11]) { $info.=";ljb_phylop.pred=$item[11]";}
                if ($item[12]) { $info.=";ljb_sift.score=$item[12]";}
                if ($item[13]) { $info.=";ljb_sift.pred=$item[13]";}
	        if ($item[14]) { $info.=";ljb_pp2.score=$item[14]";}
                if ($item[15]) { $info.=";ljb_pp2.pred=$item[15]";}
		if ($item[16]) { $info.=";ljb_lrt.score=$item[16]";}
                if ($item[17]) { $info.=";ljb_lrt.pred=$item[17]";}
                if ($item[18]) { $info.=";ljb_mt.score=$item[18]";}
                if ($item[19]) { $info.=";ljb_mt.pred=$item[19]";}
		if ($item[20]) { $info.=";ljb_gerp++.score=$item[20]";}		
		$buffer.="$info\t";
		## Format & Data
		for (my $i=34 ; $i<scalar(@item) ; $i++)
		{ $buffer.="$item[$i]\t"; 
		}
		$buffer.="\n";
		$lineCount=$lineCount+1;
		print out $buffer;
	}
}
print "Parsed $lineCount lines . \n";
exit;

