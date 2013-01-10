#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $noexecute  = 0;
my $quiet      = 0;
my $runlevels  = 0;
my $readlen    = 0;
my $trimedlen  = 0;
my $mapper     = "gsnap";
my $seg_len    = 25;
my $ins_mean   = 0;
my $ins_mean_true = 0;
my $ins_mean_AB= 0;
my $ins_sd     = 20;
my %runlevel;
my $lanename;
my $threads    = 1;
my $help;
my $AB;      #cut reads in AB (for koeln)
my $QC;      #quality check
my $SM;      #second mapping
my $BT;      #using Blast instead of BLAT for runlevel-4
my $RA         = 1;      #regional assembly
my $idra       = 0;      #for a specific breakpoint id
my $force;   #force
my $bigWig;  #wiggle file
my $gtf_guide_assembly;  #for cufflinks
my $known_trans = 'ensembl';         #for cufflinks
my $frag_bias_correct;   #for cufflinks
my $upper_quantile_norm; #for cufflinks
my $root = "$RealBin/../PIPELINE";
my $anno = "$RealBin/../ANNOTATION";
my $bin  = "$RealBin/";
my $qual_zero = 33;
my $qual_move = 0;
my $fq_reid; #rename fastq read id (for gsnap)
my $priordf = 10;         #for edgeR
my $spaired = 1;          #for edgeR
my $patient; #the patient id for edgeR DE test
my $tissue;   #the tissue type for edgeR DE test
my $gf = "png"; #the format used in html report


if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}

GetOptions(
           "lanename=s"   => \$lanename,
           "runlevel=s"   => \$runlevels,
           "noexecute"    => \$noexecute,
           "quiet"        => \$quiet,
           "readlen=i"    => \$readlen,
           "trimedlen=i"  => \$trimedlen,
           "seglen=i"     => \$seg_len,
           "mapper=s"     => \$mapper,
           "insertmean=i" => \$ins_mean,
           "insertsd=i"   => \$ins_sd,
           "threads=i"    => \$threads,
           "AB"           => \$AB,
           "QC"           => \$QC,
           "SM"           => \$SM,
           "BT"           => \$BT,
           "RA=i"         => \$RA,
           "gf=s"         => \$gf,
           "idra=i"       => \$idra,
           "WIG"          => \$bigWig,
           "fqreid"       => \$fq_reid,
           "gtf-guide"    => \$gtf_guide_assembly,
           "known-trans"  => \$known_trans,
           "frag-bias"    => \$frag_bias_correct,
           "upper-qt"     => \$upper_quantile_norm,
           "force"        => \$force,
           "root=s"       => \$root,
           "anno=s"       => \$anno,
           "priordf=i"    => \$priordf,
           "spaired=i"    => \$spaired,
           "patient=s"    => \$patient,
           "tissue=s"     => \$tissue,
           "help|h"       => \$help,
          );

#print help
helpm() if ($help);

### Annotation paths------------------------------------------------------
my $bowtie_index = "$anno/bowtie_index/hg19/hg19";
my $genome_fasta = "$bowtie_index\.fa";
my $tophat_trans_index = "$anno/bowtie_index/hg19_trans/hg19_konw_ensemble_trans";
my $gene_annotation = "$anno/hg19\.ensembl\-for\-tophat\.gff";
my $gene_annotation_gtf = "$anno/hg19\.ensembl\-for\-tophat\.gtf";
my $ensemble_gene = "$anno/UCSC\_Ensembl\_Genes\_hg19";
my $ensemble_gene_bednew = "$anno/hg19\.ensembl\.gene\.sorted\.bed12";
my $gencode_genemap = "$anno/gencode/gencode\.v14\.genemap";
my $gencode_gene_bed = "$anno/gencode/gencode\.v14\.annotation\.gene\.bed12";
my $refseq_gene = "$anno/RefSeq\_Genes\_hg19";
my $refseq_gene_gtf = "$anno/refGene_hg19.gtf";
my $gmap_index = "$anno/gmap\_index/";
my $gmap_splicesites = "$gmap_index/hg19/hg19.splicesites.iit";
#-------------------------------------------------------------------------


if ($runlevels != 0) {
  foreach my $r (split /,/,$runlevels) {
    my $from=1;
    my $to=20;
    if ($r=~/^(\d+)/) {
      $from=$1;
    }
    if ($r=~/\-(\d+)$/) {
      $to=$1;
    } elsif ($r!~/\-/) {
      $to=$from;
    }
    for (my $i=$from;$i<=$to;$i++) {
      $runlevel{$i}=1;
    }
  }
}
else {
  print STDERR "no runlevel has been set, exit.\n";
  helpm();
}


if ($root eq "$RealBin/../PIPELINE") {
  if (-e "$RealBin/../PIPELINE") {
    print STDERR "no root dir given, analysis will be run under $root.\n";
  }
  else {
    print STDERR "no root dir given, $root does not exist, please do -h or --help to check how to set root dir.\n";
    helpm();
  }
}

if ($anno eq "$RealBin/../ANNOTATION") {
  if (-e "$RealBin/../ANNOTATION") {
    print STDERR "no annotation dir given, analysis will be run under $anno.\n";
  } else {
    print STDERR "no annotation dir given, $anno does not exist, please do -h or --help to check how to set annotation dir.\n";
    helpm();
  }
}


###
###runlevel0.5: preparation the lane and read path enviroment
###

my $lanepath;

if (defined $lanename) {

  printtime();
  print STDERR "####### lane name is set to $lanename #######\n\n";

  my @lanefile = bsd_glob("$root/$lanename*fastq\.gz");
  if (scalar(@lanefile) > 0) {
    foreach my $file (@lanefile) {
      (my $newfile = $file) =~ s/fastq\.gz/fq\.gz/;
      my $cmd = "mv $file $newfile";
      RunCommand($cmd,$noexecute,$quiet);
    }
  }
  @lanefile = ();
  @lanefile = bsd_glob("$root/$lanename*fq\.gz");

  $lanepath = "$root/$lanename";

  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$lanepath/01_READS") {
    my $cmd = "mkdir -p $lanepath/01_READS";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if (scalar(@lanefile) > 0) {
    foreach my $read_file (@lanefile) {
      my $cmd = "mv $read_file $lanepath/01_READS/";
      RunCommand($cmd,$noexecute,$quiet);
    }
  }

  if ($readlen == 0 or $trimedlen == 0) { #read length or trimed length not set
     my @original_read_files = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");
     my $first_second_line = `gzip -d -c "$original_read_files[0]" | head -2 | grep -v "^@"`;
     $readlen = length($first_second_line) - 1;
     $trimedlen = $readlen;
     print STDERR "read length and trimed length are not set, will take the original read length ($readlen bp) for both (no trimming).\n";
  }
}

###
###runlevel1: trim the reads and insert size detection using spiked in reads
###

$runlevels = 1;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my @qc_files = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");
  (my $qc_out1 = $qc_files[0]) =~ s/\.gz$/\.qc/;
  (my $qc_out2 = $qc_files[1]) =~ s/\.gz$/\.qc/;
  unless (-e "$qc_out1") {
    my $cmd = "gzip -d -c $qc_files[0] | $bin/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$noexecute,$quiet);
  }
  unless (-e "$qc_out2") {
    my $cmd = "gzip -d -c $qc_files[1] | $bin/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$noexecute,$quiet);
  }
  if ($QC) {
    print STDERR "quality check finished, please check the quality file manually.\n";
    exit;
  }

  my @read_files;
  if ( $trimedlen != $readlen ) {  #trimming
    my @trimed_read_files = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz");
    if ( scalar(@trimed_read_files) == 0 ) {
      my @ori_read_files = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");
      foreach my $read_file (@ori_read_files) {
        my $read_out = $read_file;
        $read_out =~ s/fq\.gz$/trimed\.fq\.gz/;
        if ($read_file =~ /_R?1\./) {
          $read_files[0] = $read_out;
        }
        else {
          $read_files[1] = $read_out;
        }
        my $cmd = "gzip -d -c $read_file | $bin/fastx_trimmer -l $trimedlen -z -Q33 -o $read_out";
        RunCommand($cmd,$noexecute,$quiet);
      }
    }
    else { #if trimed read file exists
      @read_files = @trimed_read_files;
      @read_files = mateorder(@read_files);
    }
  }
  else {
    @read_files = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");
    @read_files = mateorder(@read_files);
  }

  if ($fq_reid){  #renbame fastq id
    foreach my $read_file (@read_files){
       my $reid_file = $read_file;
       $reid_file =~ s/\.fq\.gz/\.reid\.fq\.gz/;
       unless (-e $reid_file){
         my $cmd = "perl $bin/fqreid.pl $read_file |gzip >$reid_file";
         RunCommand($cmd,$noexecute,$quiet);
       }
       if (-e $read_file and -e $reid_file){
         my $cmd = "mv -f $reid_file $read_file";
         RunCommand($cmd,$noexecute,$quiet);
       }
    } #foreach read file
  }

  ##### AB trimming, which is not fully tested################
  if ($AB) {  # do the AB treating
     print STDERR "AB is defined\n";
     my @AB_read_files = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq\.gz");
     my @AB_read_files_ori = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq");

     if ( scalar(@AB_read_files) == 0 and scalar(@AB_read_files_ori) == 0) {
       my ($read_1, $read_2, $AB_1, $AB_2);
       foreach my $read_file (@read_files){
         if ($read_file =~ /_R?1\./) {
           $read_1 = $read_file;
           $AB_1 = $read_1;
           $AB_1 =~ s/fq\.gz$/AB\.fq/;
           $AB_read_files_ori[0] = $AB_1;
         } else {
           $read_2 = $read_file;
           $AB_2 = $read_2;
           $AB_2 =~ s/fq\.gz$/AB\.fq/;
           $AB_read_files_ori[1] = $AB_2;
         }
       }
       my $cmd = "perl $bin/AB_reads.pl $read_1 $read_2 0 >$AB_1 2>$AB_2";
       RunCommand($cmd,$noexecute,$quiet);
     }
     if ( scalar(@AB_read_files) == 0 and scalar(@AB_read_files_ori) == 2) {
       foreach my $AB_read_files_ori (@AB_read_files_ori) {
          my $cmd = "gzip $AB_read_files_ori";
          RunCommand($cmd,$noexecute,$quiet);
       }
     }
     @AB_read_files = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq\.gz");
     if ( scalar(@AB_read_files) != 2) {
       print STDERR "AB read files were wrongly generated, please remove remnant files in 01_READS.\n";
     }
     @read_files = @AB_read_files;
  } #AB##############################################################


  unless (-e "$lanepath/00_TEST") {
    my $cmd = "mkdir -p $lanepath/00_TEST";
    RunCommand($cmd,$noexecute,$quiet);
  }

  my @spiked_in = bsd_glob("$lanepath/01_READS/$lanename*spikedin.fq");
  if (scalar(@spiked_in) == 0) {
    foreach my $read_file (@read_files) {
      my $spiked_in = $read_file;
      $spiked_in =~ s/fq\.gz$/spikedin\.fq/;
      push(@spiked_in, $spiked_in);
    }
    my $cmd = "perl $bin/select_reads_with_no_n.pl $read_files[0] $read_files[1] 200000 >$spiked_in[0] 2>$spiked_in[1]";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #do the pair-end mapping of spiked_in reads
  unless (-e "$lanepath/00_TEST/$lanename\.spikedin\.hits") {
    my $cmd= "bowtie -v 2 $tophat_trans_index -1 $spiked_in[0] -2 $spiked_in[1] $lanepath/00_TEST/$lanename\.spikedin\.hits";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/00_TEST/$lanename\.spikedin\.fragmentlength") {
    my $real_len = $trimedlen;
    $real_len = $trimedlen/2 if ($AB);
    my $cmd = "perl $bin/spike_in.pl $lanepath/00_TEST/$lanename\.spikedin\.hits $real_len >$lanepath/00_TEST/$lanename\.spikedin\.fragmentlength";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}

###
###runlevel1.5: deciding fragment/insert size
###

if (defined $lanename) {

  printtime();
  print STDERR "####### insert mean and sd calculation #######\n\n";

  my @quality_check_files = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.qc");
  if ( scalar(@quality_check_files) == 2 ) { #decide the quality shift
    open QC, "<$quality_check_files[0]";
    my $qual_min = -1;
    my $qual_max = -1;
    while ( <QC> ) {
      chomp;
      next if ($_ !~ /^\d+/);
      my @cols = split /\t/;
      $qual_min = $cols[2] if ($cols[2] < $qual_min or $qual_min = -1);
      $qual_max = $cols[3] if ($cols[3] > $qual_max or $qual_max = -1);
    }
    close QC;
    print STDERR "qual_min: $qual_min; qual_max: $qual_max; ";
    $qual_zero = $qual_min+33;
    $qual_move = -$qual_min;
    print STDERR "qual_zero: $qual_zero; qual_shift: $qual_move.\n";
  } else {
    unless ($force) {
      print STDERR "please do quality check first using option --QC.\n";
      exit;
    }
  }

  if ($ins_mean == 0 or $ins_mean_AB == 0) {

    my $fraginfo = `R --no-save --slave \'--args path=\"$lanepath/00_TEST/\" lane=\"$lanename\.spikedin\"\' < $bin/fragment_length.R`;
    $fraginfo =~ /^(.+)\s(.+)/;
    my $frag_mean = $1; $frag_mean = round($frag_mean);
    my $insert_sd   = $2; $insert_sd =~ s/\n//; $insert_sd = round($insert_sd);

    my $real_len = $trimedlen;
    $real_len = $trimedlen/2 if ($AB);

    $ins_sd = $insert_sd;
    if ($AB) {
      $ins_mean_AB = $frag_mean - 2*$real_len;
      $ins_mean = $ins_mean_AB - $real_len;
    } else {
      $ins_mean = $frag_mean - 2*$real_len;
    }

    $ins_mean_true = $ins_mean;

    print STDERR "insert mean: $ins_mean\tinsert mean AB (0 if AB is not set): $ins_mean_AB\tinsert_sd: $ins_sd\n";

    unless ($force) {
      my $ins_th = round($real_len*0.5);
      if ($ins_mean < -$ins_th) {
        print STDERR "two mates is overlapping too much, please trim more.\n";
        exit 22;
      } elsif ($ins_mean >= -$ins_th && $ins_mean < 1) {
        $ins_mean = 1;
        print STDERR "insert mean is set to 1 for processing purpose.\n";
      } else {
        print STDERR "insert mean is ok.\n";
      }
    }
  }
}

###
###runlevel2: do the mapping and generate the statistics
###

$runlevels = 2;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$lanepath/02_MAPPING") {
    my $cmd = "mkdir -p $lanepath/02_MAPPING";
    RunCommand($cmd,$noexecute,$quiet);
  }

  my @reads;
  if ($AB) {
    @reads = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq\.gz");       #AB reads
  }
  elsif ($trimedlen != $readlen) {
    @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz");   #trimmed reads
  }
  else {
    @reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");    #original reads
  }
  @reads = mateorder(@reads);

  my $real_len = $trimedlen;
  $real_len = $trimedlen/2 if ($AB);

  my $real_ins_mean = $ins_mean;
  $real_ins_mean = $ins_mean_AB if ($AB);

  my $fragment_length = 2*$real_len + $real_ins_mean;

  #do the mapping of pair - end reads
  if ($mapper eq 'tophat1') {
     unless (-e "$lanepath/02_MAPPING/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam") {
       my $cmd = "tophat --output-dir $lanepath/02_MAPPING --mate-inner-dist $real_ins_mean --mate-std-dev $ins_sd --library-type fr-unstranded -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat_trans_index $bowtie_index $reads[0] $reads[1]";
       RunCommand($cmd,$noexecute,$quiet);
     }
  } #tophat1

  elsif ($mapper eq 'tophat2') {
    unless (-e "$lanepath/02_MAPPING/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam") {
       my $cmd = "tophat --bowtie1 --output-dir $lanepath/02_MAPPING --mate-inner-dist $real_ins_mean --mate-std-dev $ins_sd --library-type fr-unstranded -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat_trans_index $bowtie_index $reads[0] $reads[1]";
       RunCommand($cmd,$noexecute,$quiet);
     }
  } #tophat2

  elsif ($mapper eq 'gsnap'){

     my $quality_options;
     if ( ($qual_zero - 33) < 10 ) {
        $quality_options = "--quality-protocol=sanger";
     } else {
        $quality_options = "--quality-zero-score=$qual_zero --quality-print-shift=$qual_move";
     }

     unless (-e "$lanepath/02_MAPPING/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.sam") {
        my $cmd = "gsnap -d hg19 -D $gmap_index -B 5 --gunzip --format=sam --nthreads=$threads -s $gmap_splicesites --npaths=5 $quality_options $reads[0] $reads[1] >$lanepath/02_MAPPING/accepted_hits\.sam";
        RunCommand($cmd,$noexecute,$quiet);
     }

     if (-e "$lanepath/02_MAPPING/accepted_hits\.sam" and (! -e "$lanepath/02_MAPPING/accepted_hits\.bam" and ! -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam")){
        my $cmd = "samtools view -Sb -@ $threads $lanepath/02_MAPPING/accepted_hits\.sam -o $lanepath/02_MAPPING/accepted_hits\.bam";
        RunCommand($cmd,$noexecute,$quiet);
     }

     if (-e "$lanepath/02_MAPPING/accepted_hits\.sam" and (-e "$lanepath/02_MAPPING/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam")){
        my $cmd  = "rm $lanepath/02_MAPPING/accepted_hits\.sam -f";
        RunCommand($cmd,$noexecute,$quiet);
     }

  } #gsnap

  else {
     print STDERR "Error: --mapper option should only be gsnap, tophat1 or tophat2. \n\n";
     exit 22;
  }


  #do the statistics
  unless (-e "$lanepath/03_STATS") {
    my $cmd = "mkdir -p $lanepath/03_STATS";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.mapping\.stats") {
    my $cmd = "$bin/Rseq_bam_stats --mapping $lanepath/02_MAPPING/accepted_hits\.bam --writer $lanepath/02_MAPPING/accepted_hits\.unique\.bam --arp $lanepath/03_STATS/$lanename\.arp --breakpoint $lanepath/03_STATS/$lanename\.breakpoints >$lanepath/03_STATS/$lanename\.mapping\.stats";
    RunCommand($cmd,$noexecute,$quiet);
  }


  my $mapping_stats_line_number = `wc -l $lanepath/03_STATS/$lanename.mapping.stats`;
  $mapping_stats_line_number =~ s/^(\d+).*$/$1/;
  chomp($mapping_stats_line_number);
  if ($mapping_stats_line_number == 12){
    my $total_reads = `gzip -d -c $reads[0] | wc -l`;
    $total_reads /= 4;
    open STATS, ">>$lanepath/03_STATS/$lanename\.mapping\.stats" || die "can not open $lanepath/03_STATS/$lanename\.mapping\.stats\n";
    print STATS "total_frag: $total_reads\n";
    close STATS;
  }

  unless (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam")  {
    my $cmd = "samtools sort -@ $threads $lanepath/02_MAPPING/accepted_hits\.unique\.bam $lanepath/02_MAPPING/accepted_hits\.unique\.sorted";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam" and -e "$lanepath/02_MAPPING/accepted_hits\.unique\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.unique\.bam -f";
    RunCommand($cmd,$noexecute,$quiet);
  }
  if (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam" and -e "$lanepath/02_MAPPING/accepted_hits\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.bam -f";
    #RunCommand($cmd,$noexecute,$quiet);
  }

  if ($bigWig) { #generate wiggle file
    unless (-e "$lanepath/03_STATS/$lanename\.bw"){
      if (-e "$lanepath/03_STATS/$lanename\.bedgraph") {
         my $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$lanename\.bedgraph $anno/chrom_sizes.hg19 $lanepath/03_STATS/$lanename\.bw";
         RunCommand($cmd,$noexecute,$quiet);
      }
      else {
         my $cmd = "genomeCoverageBed -ibam $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam -bg -split -g $anno/chrom_sizes.hg19 >$lanepath/03_STATS/$lanename\.bedgraph";
         RunCommand($cmd,$noexecute,$quiet);
         $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$lanename\.bedgraph $anno/chrom_sizes.hg19 $lanepath/03_STATS/$lanename\.bw";
         RunCommand($cmd,$noexecute,$quiet);
      }
      if (-e "$lanepath/03_STATS/$lanename\.bedgraph" and -e "$lanepath/03_STATS/$lanename\.bw") {
         my $cmd = "rm $lanepath/03_STATS/$lanename\.bedgraph -f";
         RunCommand($cmd,$noexecute,$quiet);
      }
    }
  }

  #ensembl transcripts######################################################
  unless (-e "$lanepath/03_STATS/$lanename\.expr.sorted") {
    unless (-e "$lanepath/03_STATS/$lanename\.expr") {
      my $cmd1 = "$bin/Rseq_bam_reads2expr --region $ensemble_gene --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam --posc $lanepath/03_STATS/$lanename\.pos\.gff --chrmap $lanepath/03_STATS/$lanename\.chrmap --lbias $lanepath/03_STATS/$lanename\.lbias >$lanepath/03_STATS/$lanename\.expr";
      RunCommand($cmd1,$noexecute,$quiet);
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$lanename\.expr >$lanepath/03_STATS/$lanename\.expr.sorted";
    RunCommand($cmd2,$noexecute,$quiet);
    if (-e "$lanepath/03_STATS/$lanename\.expr.sorted" and "$lanepath/03_STATS/$lanename\.expr"){
      my $cmd3 = "rm $lanepath/03_STATS/$lanename\.expr -rf";
      RunCommand($cmd3,$noexecute,$quiet);
    }
  }

  #refseq gene##############################################################
  unless (-e "$lanepath/03_STATS/$lanename\.RefSeq\.expr.sorted") {
    unless (-e "$lanepath/03_STATS/$lanename\.RefSeq\.expr") {
      my $cmd = "$bin/Rseq_bam_reads2expr --region $refseq_gene --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam >$lanepath/03_STATS/$lanename\.RefSeq\.expr";
      RunCommand($cmd,$noexecute,$quiet);
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$lanename\.RefSeq\.expr >$lanepath/03_STATS/$lanename\.RefSeq\.expr.sorted";
    RunCommand($cmd2,$noexecute,$quiet);
    if (-e "$lanepath/03_STATS/$lanename\.RefSeq\.expr.sorted" and "$lanepath/03_STATS/$lanename\.RefSeq\.expr") {
      my $cmd3 = "rm $lanepath/03_STATS/$lanename\.RefSeq\.expr -rf";
      RunCommand($cmd3,$noexecute,$quiet);
    }
  }

  #ensembl gene#############################################################
  unless (-e "$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr\.sorted") {
    unless (-e "$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr") {
      my $cmd1 = "$bin/Rseq_bam_reads2expr --region $ensemble_gene_bednew --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam >$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr";
      RunCommand($cmd1,$noexecute,$quiet);
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$lanename\.ensembl\_gene\.expr >$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr.sorted";
    RunCommand($cmd2,$noexecute,$quiet);
    if (-e "$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr.sorted" and "$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr") {
      my $cmd3 = "rm $lanepath/03_STATS/$lanename\.ensembl\_gene\.expr -rf";
      RunCommand($cmd3,$noexecute,$quiet);
    }
  }

  unless (-e "$lanepath/03_STATS/$lanename\.ensembl\_gene\.count") {
    open ENSEMBL_GENE, "$lanepath/03_STATS/$lanename\.ensembl\_gene\.expr.sorted";
    open ENSEMBL_GENE_COUNT, ">$lanepath/03_STATS/$lanename\.ensembl\_gene\.count";
    while ( <ENSEMBL_GENE> ) {
      chomp;
      my @cols = split /\t/;
      my $current_count = round($cols[7]);
      print ENSEMBL_GENE_COUNT "$cols[3]\t$current_count\n";
    }
    close ENSEMBL_GENE;
    close ENSEMBL_GENE_COUNT;
  }

  #gencode gene#############################################################
  unless (-e "$lanepath/03_STATS/$lanename\.gencode\_gene\.expr.sorted") {
    unless (-e "$lanepath/03_STATS/$lanename\.gencode\_gene\.expr") {
      my $cmd1 = "$bin/Rseq_bam_reads2expr --region $gencode_gene_bed --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam >$lanepath/03_STATS/$lanename\.gencode\_gene\.expr";
      RunCommand($cmd1,$noexecute,$quiet);
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$lanename\.gencode\_gene\.expr >$lanepath/03_STATS/$lanename\.gencode\_gene\.expr.sorted";
    RunCommand($cmd2,$noexecute,$quiet);
    if (-e "$lanepath/03_STATS/$lanename\.gencode\_gene\.expr.sorted" and "$lanepath/03_STATS/$lanename\.gencode\_gene\.expr") {
      my $cmd3 = "rm $lanepath/03_STATS/$lanename\.gencode\_gene\.expr -rf";
      RunCommand($cmd3,$noexecute,$quiet);
    }
  }


  #for RPKM normalization of gencode genes##################################
  unless (-e "$lanepath/03_STATS/$lanename\.gencode\_gene\.rpkm") {
    my $N_mapped_reads = 0;
    my $mapped = 0;
    my $singleton = 0;
    open MAPPING_STATS, "$lanepath/03_STATS/$lanename.mapping.stats" || die "can not open $lanepath/03_STATS/$lanename.mapping.stats";
    while ( <MAPPING_STATS> ) {
      chomp;
      if ($_ =~ /^Mapped\:\s+(\d+)$/) {
        $mapped = $1;
      }
      if ($_ =~ /^singletons\:\s+(\d+)$/) {
        $singleton = $1;
      }
    }
    print STDERR "mapped pairs equal to zero!!!\n" if ($mapped == 0);
    print STDERR "singletons equal to zero!!!\n" if ($singleton == 0);
    $N_mapped_reads = 2*$mapped - $singleton;
    exit if ($N_mapped_reads == 0);
    close MAPPING_STATS;

    open GENCODE_GENEMAP, "$gencode_genemap";
    my %gencode_genemap;
    while ( <GENCODE_GENEMAP> ) {
      chomp;
      my ($gene_id, $gene_type, $gene_name) = split /\t/;
      $gencode_genemap{$gene_id} = $gene_type."\t".$gene_name;
    }
    close GENCODE_GENEMAP;

    open GENCODE_GENE_EXPR, "$lanepath/03_STATS/$lanename\.gencode\_gene\.expr.sorted" || die "can not open $lanepath/03_STATS/$lanename\.gencode\_gene\.expr.sorted";
    open GENCODE_RPKM, ">$lanepath/03_STATS/$lanename\.gencode\_gene\.rpkm";
    printf GENCODE_RPKM ("%s\n", join("\t", "#gene_id","count", "RPKM", "gene_type", "gene_name"));
    while ( <GENCODE_GENE_EXPR> ) {
      chomp;
      my @cols = split /\t/;
      my $gencode_name = $cols[3];
      my $counts_dblength = $cols[4];
      my $counts = round($cols[7]);
      my $rpkm = sprintf("%.3f", $counts_dblength * 1e9/$N_mapped_reads);
      print GENCODE_RPKM "$gencode_name\t$counts\t$rpkm\t$gencode_genemap{$gencode_name}\n";
    }
    close GENCODE_GENE_EXPR;
    close GENCODE_RPKM;
  }


  unless (-e "$lanepath/03_STATS/$lanename\.cate") {
    my $cmd = "perl $bin/cate.pl $lanepath/03_STATS/$lanename\.ensembl\_gene\.expr\.sorted $gene_annotation >$lanepath/03_STATS/$lanename\.cate";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.ins") {
    my $cmd = "samtools view -f 0x2 $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam | cut -f 9 | awk \'\$1\>0 \&\& \$1\<500\' >$lanepath/03_STATS/$lanename\.ins";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.report/$lanename\.report\.html") {
    my $cmd = "R CMD BATCH --no-save --no-restore "."\'--args path=\"$lanepath\" lane=\"$lanename\" anno=\"$anno\" src=\"$bin\" readlen=$real_len gf=\"$gf\"' $bin/html_report.R $lanepath/03_STATS/R\_html\.out";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}


###
###runlevel3: select anormouls read pairs and do the assembly
###

$runlevels = 3;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  $RA = 1 if $mapper eq 'gsnap';

  if ($RA != 0) {  #regional assembly test###############################################

    my $assembly_type = "independent regional assembly";
    if ($RA == 2){$assembly_type = "columbus assembly";}

    printtime();
    print STDERR "Regional assembly process is starting \($assembly_type\)... first breakpoint processing...\n\n";

    unless (-e "$lanepath/04_ASSEMBLY") {
      my $cmd = "mkdir -p $lanepath/04_ASSEMBLY";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/03_STATS/$lanename\.breakpoints"){
      print STDERR "Error: breakpoint file does not exist, do runlevel 2 first.\n\n";
      exit 22;
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed"){
      my $cmd = "perl $bin/breakpoint\_processing.pl $lanepath/03_STATS/$lanename\.breakpoints >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted"){
      my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot"){
      my $cmd = "perl $bin/breakpoint_repeat_masker.pl $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted $anno/UCSC\_repeats\_hg19\.gff >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter"){
      my $cmd = "perl $bin/breakpoint_final_prepare.pl $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted"){
      my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p"){
      my $cmd = "grep \"\tp\t\" $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p";
      RunCommand($cmd,$noexecute,$quiet);
    }

    if ($RA == 1) { #independent regional assembly
      unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p\.reads"){
        my $mapping_bam = "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam";
        die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
        my $cmd = "$bin/reads_in_region --region $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p --mapping $mapping_bam >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p\.reads";
        RunCommand($cmd,$noexecute,$quiet);
      }

      my @RAssembly_reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.RAssembly\.fq");
      unless ( scalar(@RAssembly_reads) == 2 ) { #get raw reads
        my @reads;
        if ($trimedlen != $readlen) {
          @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
        } else {
          @reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz"); #original reads
        }
        @reads = mateorder(@reads);

        my $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p\.reads --readfile1 $reads[0] --readfile2 $reads[1] --RA";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } elsif ($RA == 2) {        #columbus assembler : get fasta sequences for regions
      unless (-e "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.regions\.fa"){
        my $cmd = "perl $bin/get_sequence_by_region.pl $lanepath/04_ASSEMBLY/$lanename\.breakpoints\.processed\.sorted\.repornot\.filter\.sorted\.p $genome_fasta >$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.regions\.fa";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } else {
      print STDERR "Error: --RA must be set to 1 or 2 if not 0 (default)...\n\n";
      exit 22;
    }


    printtime();
    print STDERR "now assembling...\n\n";

    unless (-e "$lanepath/04_ASSEMBLY/transcripts.fa") {

      if ($RA == 1) {

        my @RAssembly_reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.RAssembly\.fq");
        @RAssembly_reads = mateorder(@RAssembly_reads);
        #my $bp_regions_file = "$lanepath/04_ASSEMBLY/$lanename\.bp\_regions\.fa";
        my $ra_reads1_file = $RAssembly_reads[0];
        my $ra_reads2_file = $RAssembly_reads[1];

        #my %bp_regions_jumpers;
        my %ra_reads1_jumpers;
        my %ra_reads2_jumpers;

        #find_jumpers($bp_regions_file, \%bp_regions_jumpers);
        find_jumpers($ra_reads1_file, \%ra_reads1_jumpers);
        find_jumpers($ra_reads2_file, \%ra_reads2_jumpers);

        #my $n_ids_bp_regions = scalar(keys %bp_regions_jumpers);
        my $n_ids_ra_reads1 = scalar(keys %ra_reads1_jumpers);
        my $n_ids_ra_reads2 = scalar(keys %ra_reads2_jumpers);
        if ($n_ids_ra_reads1 != $n_ids_ra_reads2) {
          print STDERR "Error: the numbers of breakpoints are inconsistant between reads ($n_ids_ra_reads1 vs $n_ids_ra_reads2).\n";
          exit 22;
        }

        #now should do regional assembly one by one
        #open BP_REGIONS, "$bp_regions_file";
        open RA_READS1, "$ra_reads1_file";
        open RA_READS2, "$ra_reads2_file";
        my $speed_count = 0;
        my %printed_speed;
        foreach my $bp_id (sort {$a<=>$b} keys %ra_reads1_jumpers) {

          if ($idra != 0 and $bp_id != $idra) {
             next;
          }

          my $cmd = "mkdir -p $lanepath/04_ASSEMBLY/$bp_id"; #create dir
          RunCommand($cmd,$noexecute,1);

          #regional_assembly(\*BP_REGIONS, $bp_regions_jumpers{$bp_id}, $bp_id, "$lanepath/04_ASSEMBLY/$bp_id/regions.fa");
          regional_assembly(\*RA_READS1, $ra_reads1_jumpers{$bp_id}, $bp_id, "$lanepath/04_ASSEMBLY/$bp_id/reads_1.fq");
          regional_assembly(\*RA_READS2, $ra_reads2_jumpers{$bp_id}, $bp_id, "$lanepath/04_ASSEMBLY/$bp_id/reads_2.fq");

          unless (-e "$lanepath/04_ASSEMBLY/$bp_id/Roadmaps") {
            my $cmd = "velveth $lanepath/04_ASSEMBLY/$bp_id/ 21 -fastq -shortPaired -separate $lanepath/04_ASSEMBLY/$bp_id/reads_1.fq $lanepath/04_ASSEMBLY/$bp_id/reads_2.fq";
            RunCommand($cmd,$noexecute,1);
          }

          unless (-e "$lanepath/04_ASSEMBLY/$bp_id/Graph2") {
            my $cmd = "velvetg $lanepath/04_ASSEMBLY/$bp_id/ -read_trkg yes";
            RunCommand($cmd,$noexecute,1);
          }

          unless (-e "$lanepath/04_ASSEMBLY/$bp_id/transcripts.fa"){
            my $cmd = "oases $lanepath/04_ASSEMBLY/$bp_id/";
            RunCommand($cmd,$noexecute,1);
          }

          if (-e "$lanepath/04_ASSEMBLY/$bp_id/transcripts.fa") {
            open TRANSCRIPTS, "$lanepath/04_ASSEMBLY/$bp_id/transcripts.fa";
            open TRANSCRIPTSALL, ">>$lanepath/04_ASSEMBLY/transcripts.fa";
            while ( <TRANSCRIPTS> ) {
              chomp;
              if ($_ =~ /^>/) {
                $_ .= "\_id\_".$bp_id;
              }
              print TRANSCRIPTSALL "$_\n";
            }
            close TRANSCRIPTS;
            close TRANSCRIPTSALL;
          }

          if (-e "$lanepath/04_ASSEMBLY/transcripts.fa") {
            my $cmd = "rm $lanepath/04_ASSEMBLY/$bp_id/ -rf";
            RunCommand($cmd,$noexecute,1) unless $idra != 0;
          }

          $speed_count++;
          my $percentage_now = round(($speed_count/$n_ids_ra_reads1)*100);
          if ($percentage_now%5 == 0) {
            if (! exists $printed_speed{$percentage_now}) {
              print STDERR "$percentage_now\%\.\.\.";
              $printed_speed{$percentage_now} = '';
            }
          }                     #print percentage finished
        }                       #for each bp id
        print STDERR "\n";
        #close BP_REGION;
        close RA_READS1;
        close RA_READS2;
      }                         #RA == 1

      elsif ($RA == 2) { #columbus assmbly
        my $mapping_bam = "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam";
        die "Error: the name-sorted mapping bam file is not available." unless (-e $mapping_bam);
        my $regions_fa = "$lanepath/04_ASSEMBLY/$lanename\.breakpoints\.regions\.fa";
        die "Error: the region fasta file for columbus assembleris not available." unless (-e $regions_fa);

        unless (-e "$lanepath/04_ASSEMBLY/Roadmaps") {
          my $cmd = "velveth $lanepath/04_ASSEMBLY/ 21 -reference $regions_fa -shortPaired -bam $mapping_bam";
          RunCommand($cmd,$noexecute,$quiet);
        }

        unless (-e "$lanepath/04_ASSEMBLY/Graph2") {
          my $cmd = "velvetg $lanepath/04_ASSEMBLY/ -read_trkg yes";
          RunCommand($cmd,$noexecute,$quiet);
        }

        unless (-e "$lanepath/04_ASSEMBLY/transcripts.fa"){
          my $cmd = "oases $lanepath/04_ASSEMBLY/";
          RunCommand($cmd,$noexecute,$quiet);
        }
      } else {
        print STDERR "Error: --RA must be set to 1 or 2 if not 0 (default)...\n\n";
        exit 22;
      }
    }  #if transcripts.fa not exist
  }   #regional assembly test##########################################################################

  else {
    #get the ARP from unmapped by tophat ############
    unless ( -e "$lanepath/01_READS/$lanename\.ARP\.fq\.gz" ) {

      my @ARP = bsd_glob("$lanepath/01_READS/$lanename*ARP\.fq\.gz");

      if (scalar(@ARP) != 2) {
        my @reads;
        if ($trimedlen != $readlen) {
          @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
        } else {
          @reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz"); #original reads
        }
        @reads = mateorder(@reads);

        unless (-e "$lanepath/02_MAPPING/unmapped_left\.fq\.z"){
          print STDERR "Error: unmapped read file does not exist!!!\n";
          exit 22;
        }

        my $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/03_STATS/$lanename\.arp --unmap $lanepath/02_MAPPING/unmapped_left\.fq\.z --readfile1 $reads[0] --readfile2 $reads[1]";
        $cmd .= " --AB" if ($AB);
        RunCommand($cmd,$noexecute,$quiet);

        my @ori_ARP_file = bsd_glob("$lanepath/01_READS/$lanename*ARP\.fq");
        foreach my $ori_ARP_file (@ori_ARP_file){
          my $cmd = "gzip $ori_ARP_file";
          RunCommand($cmd,$noexecute,$quiet);
        }
      }

      @ARP = bsd_glob("$lanepath/01_READS/$lanename*ARP\.fq\.gz");
      @ARP = mateorder(@ARP);

      #in such case, do a second mapping based on trimed read, currently trim to 36bp
      if (($trimedlen >= 70 and $ins_mean <= 20) || $SM ) {

        #trimming
        foreach my $ARP (@ARP) {
          my $ARP_trimed36 = $ARP;
          $ARP_trimed36 =~ s/fq\.gz$/trimed36\.fq\.gz/;
          unless ( -e "$ARP_trimed36" ){
            my $cmd = "gzip -d -c $ARP | $bin/fastx_trimmer -l 36 -Q33 -z -o $ARP_trimed36";
            RunCommand($cmd,$noexecute,$quiet);
          }
        }

        #second mapping
        my @ARP_trimed36 = bsd_glob("$lanepath/01_READS/$lanename*ARP*trimed36\.fq\.gz");
        @ARP_trimed36 = mateorder(@ARP_trimed36);
        unless (-e "$lanepath/02_MAPPING/SecondMapping/") {
          my $cmd = "mkdir -p $lanepath/02_MAPPING/SecondMapping/";
          RunCommand($cmd,$noexecute,$quiet);
        }
        unless ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.unique\.bam" ) { #second mapping using gsnap
          if ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam" ) { #no bam but sam, need to compress it
            my $cmd = "samtools view -Sb $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -o $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam";
            RunCommand($cmd,$noexecute,$quiet);
          } else {

            my $quality_options;
            if ( ($qual_zero - 33) < 10 ) {
              $quality_options = "--quality-protocol=sanger";
            } else {
              $quality_options = "--quality-zero-score=$qual_zero --quality-print-shift=$qual_move";
            }
            my $cmd = "gsnap -d hg19 -D $gmap_index -B 5 --gunzip --format=sam --nthreads=$threads -s $gmap_splicesites --max-mismatches=2 --npaths=10 --trim-mismatch-score=0 --trim-indel-score=0 $quality_options $ARP_trimed36[0] $ARP_trimed36[1] >$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam";
            RunCommand($cmd,$noexecute,$quiet);
            $cmd = "samtools view -Sb $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -o $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam";
            RunCommand($cmd,$noexecute,$quiet);
          }
        }

        if ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam" and -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam") {
          my $cmd = "rm $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -f";
          RunCommand($cmd,$noexecute,$quiet);
        }

        unless ( -e "$lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.stats" )  {
          my $cmd = "$bin/Rseq_bam_stats --mapping $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam --writer $lanepath/02_MAPPING/SecondMapping/accepted_hits\.unique\.bam --arp $lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.arp >$lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.stats";
          RunCommand($cmd,$noexecute,$quiet);
        }

        #now get the real arp after second mapping
        my @ARP_sm = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq\.gz");
        if (scalar(@ARP_sm) != 2) {
          my @reads;
          if ($trimedlen != $readlen) {
            @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
          } else {
            @reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz"); #original reads
          }
          @reads = mateorder(@reads);

          my $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.arp --readfile1 $reads[0] --readfile2 $reads[1] --SM";
          $cmd .= " --AB" if ($AB);
          RunCommand($cmd,$noexecute,$quiet);
        }

        my @ori_ARP_SM_file = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq");
        foreach my $ori_ARP_SM_file (@ori_ARP_SM_file) {
          my $cmd = "gzip $ori_ARP_SM_file";
          RunCommand($cmd,$noexecute,$quiet);
        }

        @ARP = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq\.gz");
        @ARP = mateorder(@ARP);
      }

      # shuffle ARP reads to a single file
      my $cmd = "perl $bin/shuffleSequences_fastq.pl $ARP[0] $ARP[1] $lanepath/01_READS/$lanename\.ARP\.fq";
      RunCommand($cmd,$noexecute,$quiet);
      $cmd = "gzip $lanepath/01_READS/$lanename\.ARP\.fq";
      RunCommand($cmd,$noexecute,$quiet);
    }

    #do the velveth
    unless (-e "$lanepath/04_ASSEMBLY") {
      my $cmd = "mkdir -p $lanepath/04_ASSEMBLY";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/Roadmaps") {
      my $cmd = "velveth $lanepath/04_ASSEMBLY/ 21 -fastq.gz -shortPaired $lanepath/01_READS/$lanename\.ARP\.fq\.gz";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/Graph2") {
      my $frag_len = 2*$trimedlen + $ins_mean;
      my $cmd = "velvetg $lanepath/04_ASSEMBLY/ -ins_length $frag_len -ins_length_sd $ins_sd -exp_cov auto -read_trkg yes -scaffolding no";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/04_ASSEMBLY/transcripts.fa") {
      my $frag_len = 2*$trimedlen + $ins_mean;
      my $cmd = "oases $lanepath/04_ASSEMBLY/ -ins_length $frag_len -ins_length_sd $ins_sd -unused_reads yes -scaffolding no ";
      RunCommand($cmd,$noexecute,$quiet);
    }
  }
  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}

###
###runlevel4: get fusion candidates and visualize the result
###

$runlevels = 4;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$lanepath/05_FUSION") {
    my $cmd = "mkdir -p $lanepath/05_FUSION";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if ($BT) {  #using BLAST alternative
    unless (-e "$lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast") {
      my $cmd = "$bin/blastn -query $lanepath/04_ASSEMBLY/transcripts\.fa -db $anno/human\.rna\.fna\.blastdb -outfmt 7 -num_threads $threads -out $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq") {
      my $cmd = "perl $bin/pick_fusion_transcripts_from_BLAT.pl --refseq $anno/human\.rna\.fna --transcript $lanepath/04_ASSEMBLY/transcripts.fa --blat $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast --lanename $lanename --lanepath $lanepath >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap"){
      my $cmd = "gmap -D $gmap_index -d hg19 --format=psl -t $threads $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq"){
      my $cmd = "perl $bin/filter_out_FP_from_blatps.pl --fusion_bf_seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq --fusion_bf $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration --fusion_bf_blat $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq";
      RunCommand($cmd,$noexecute,$quiet);
    }

  } #BLAST

  else {  #using BLAT
    unless (-e "$lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat") {
      my $cmd = "blat $anno/human\.rna\.fna $lanepath/04_ASSEMBLY/transcripts.fa -out=blast9 $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq") {
      my $cmd = "perl $bin/pick_fusion_transcripts_from_BLAT.pl --refseq $anno/human\.rna\.fna --transcript $lanepath/04_ASSEMBLY/transcripts.fa --blat $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat --lanename $lanename --lanepath $lanepath >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.blat"){
      my $cmd = "blat $anno/hg19\_UCSC\.2bit $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration.genome\.blat";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq"){
      my $cmd = "perl $bin/filter_out_FP_from_blatps.pl --fusion_bf_seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq --fusion_bf $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration --fusion_bf_blat $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.blat >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq";
      RunCommand($cmd,$noexecute,$quiet);
    }
  } #BLAT

  #build fusion candidate index
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index\.1\.ebwt"){
    my $cmd = "bowtie-build --quiet $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #map reads to the fusion candidate index
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie"){
    my @reads;
    if ($trimedlen != $readlen) {
      @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
    } else {
      @reads = bsd_glob("$lanepath/01_READS/$lanename\_{R,}[12]\.fq\.gz");  #original reads
    }
    @reads = mateorder(@reads);

    my $cmd = "gzip -d -c $reads[0] $reads[1] | bowtie -v 1 -k 10 -m 10 -p $threads $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index \- $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #get fusion coverage
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov") {
    my $cmd = "perl $bin/get_fusion_coverage.pl --type pair --mappingfile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie --readlength $trimedlen --geneanno $gene_annotation --ensrefname $anno/Ensembl\_Ref\_Name\.tsv --locname $anno/Name2Location\.hg19 --refgene $anno/refgenes\.hg19 --repeatmasker $anno/UCSC\_repeats\_hg19\.gff --accepthits $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam --encomcov $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov.enco >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #visualize coverage
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis") {
    my $cmd = "perl $bin/further_processing_for_read_visualization.pl --fusionseqfile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq --coveragefile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov --readlength $trimedlen >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.list"){
    my $cmd = "grep \"^\#\" $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis \| sort -k 13,13nr -k 6,6d -k 8,8d >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.list";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel5: run cufflinks for gene expression or genomic guided assembly
###

$runlevels = 5;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$lanepath/06_CUFFLINKS") {
    my $cmd = "mkdir -p $lanepath/06_CUFFLINKS";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$lanepath/06_CUFFLINKS/transcripts\.gtf") {

    my $cufflinks_options = "";
    my $known_trans_file = $gene_annotation_gtf;
    if ($known_trans eq 'refseq') {
      $known_trans_file = $refseq_gene_gtf;
    }

    if ( $gtf_guide_assembly ) {
      $cufflinks_options .= "--GTF-guide $known_trans_file";
    }
    else {
      $cufflinks_options .= "--GTF $known_trans_file";
    }
    if ( $frag_bias_correct ) {
      $cufflinks_options .= "--frag-bias-correct";
    }
    if ( $upper_quantile_norm ){
      $cufflinks_options .= "--upper-quartile-norm";
    }

    my $mapping_bam = "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam";
    if (-e $mapping_bam){
      my $cmd = "cufflinks -o $lanepath/06_CUFFLINKS -p $threads $cufflinks_options --quiet $mapping_bam";
      RunCommand($cmd,$noexecute,$quiet);
    }
    else {
      print STDERR "$mapping_bam does not exist, please do the mapping first.\n";
      exit;
    }
  }

  if (-e "$lanepath/06_CUFFLINKS/transcripts\.gtf") { #collect the gtf list
     my $current_gtf = "$lanepath/06_CUFFLINKS/transcripts\.gtf";
     my $gtf_list = "$root/GTF\_list\.txt";

     if (! -e "$gtf_list") {
        open GTF_LIST, ">>$gtf_list";
        printf GTF_LIST "$current_gtf\n";
        close GTF_LIST;
     } else {
        my $writeornot = 1;
        open GTF_IN, "$gtf_list";
        while ( <GTF_IN> ) {
          chomp;
          $writeornot = 0 if ($_ eq $current_gtf);
        }
        close GTF_IN;
        if ($writeornot == 1){
          open GTF_OUT, ">>$gtf_list";
          printf GTF_OUT "$current_gtf\n";
          close GTF_OUT;
        }
     } #existing gtf list file
  } #gtf file is generated

  unless (-e "$lanepath/06_CUFFLINKS/$lanename\.transcripts\.gtf\.tmap") {
    my $cmd = "cuffcompare -o $lanepath/06_CUFFLINKS/$lanename -r $refseq_gene_gtf $lanepath/06_CUFFLINKS/transcripts\.gtf";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


DE_ANALYSIS:


###
###runlevel6: run cuffdiff for diffrential gene/isoform expression analysis
###

$runlevels = 6;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $cuffmerge_output = "$root/cuffmerge";
  my $cuffdiff_output = "$root/cuffdiff";
  my $gtf_list = "$root/GTF\_list\.txt";

  unless (-e "$gtf_list") {
    print STDERR "Error: $gtf_list file does not exist!!! Please rerun RTrace.pl for each sample for runlevel 5 (cufflinks).\n\n";
    exit 22;
  }

  #cuffmerge####################################
  unless (-e $cuffmerge_output) {
    my $cmd = "mkdir -p $cuffmerge_output";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$cuffmerge_output/merged\.gtf") {
    my $cmd = "cuffmerge -o $cuffmerge_output --ref-gtf $gene_annotation_gtf --num-threads $threads --ref-sequence $genome_fasta $gtf_list";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #cuffdiff#####################################
  unless (-e $cuffdiff_output) {
    my $cmd = "mkdir -p $cuffdiff_output";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$cuffdiff_output/splicing\.diff") {
    my %bam_files_cd;
    open TARGETS, "$root/edgeR/targets" || die "Error: could not open $root/edgeR/targets for the information of bam files for cuffdiff.\n please rerun the pipeline for each samples with --patient and --tissue arguments set.\n";
    while ( <TARGETS> ){
       chomp;
       next if /^label/;
       my ($cu_sample, $cu_expr_count, $tissue, $patient) = split /\t/;
       my $cu_bam_file = "$root/$cu_sample/02_MAPPING/accepted_hits\.unique\.sorted\.bam";
       push (@{$bam_files_cd{$tissue}}, $cu_bam_file);
    }
    close TARGETS;

    my $bam_files1;
    my $bam_files2;
    my @tissue_names = sort {$a cmp $b} keys %bam_files_cd;  #usually normal > cancer
    if (scalar(@tissue_names) != 2) {
      print STDERR "Error: the number of tissue types does not equal to 2, currently not supported.\n\n";
      exit 22;
    } else {
      $bam_files1 = join(",", @{$bam_files_cd{$tissue_names[0]}});
      $bam_files2 = join(",", @{$bam_files_cd{$tissue_names[1]}});
    }

    my $cmd = "cuffdiff --output-dir $cuffdiff_output --num-threads $threads --labels $tissue_names[0]\,$tissue_names[1] $cuffmerge_output/merged\.gtf $bam_files1 $bam_files2";
    RunCommand($cmd,$noexecute,$quiet);
  } #run cuffdiff

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}



##write target file --- for edgeR (and cuffdiff) ########################
if ($patient and $tissue){
  my $count_file = "$lanepath/03_STATS/$lanename\.ensembl\_gene\.count";
  unless (-e "$root/edgeR") {
    my $cmd = "mkdir -p $root/edgeR";
    RunCommand($cmd,$noexecute,$quiet);
  }
  if (! -e "$root/edgeR/targets") {
    open TARGETS_OUT, ">>$root/edgeR/targets";
    printf TARGETS_OUT "%s\n", join("\t", 'label', 'files', 'tissue', 'patient');
    printf TARGETS_OUT "%s\n", join("\t", $lanename, $count_file, $tissue, $patient);
    close TARGETS_OUT;
  } else {
    my $writeornot = 1;
    open TARGETS_IN, "$root/edgeR/targets";
    while ( <TARGETS_IN> ){
       chomp;
       my @cols = split /\t/;
       $writeornot = 0 if ($cols[0] eq $lanename);
    }
    close TARGETS_IN;
    if ($writeornot == 1){
      open TARGETS_OUT, ">>$root/edgeR/targets";
      printf TARGETS_OUT "%s\n", join("\t", $lanename, $count_file, $tissue, $patient);
      close TARGETS_OUT;
    }
  } #existing targets file
} elsif ($patient or $tissue) {
  print STDERR "please specify both --patient and --tissue information\n";
  exit 22;
}


###
###runlevel7: run edgeR for diffrential gene expression analysis
###

$runlevels = 7;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $edgeR_output = "$root/edgeR";

  if (! -e $edgeR_output) {
     print STDERR "Error: $edgeR_output dir does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  } elsif (! -e "$edgeR_output/targets") {
     print STDERR "Error: $edgeR_output/targets file does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  }

  unless (-e "$edgeR_output/topDE.txt") {
    my $cmd = "R2151 CMD BATCH --no-save --no-restore "."\'--args path=\"$edgeR_output\" priordf=$priordf spaired=$spaired\' $bin/edgeR_test.R $edgeR_output/R\_html\.out";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-e "$edgeR_output/DE\_table\.txt") {
    my $cmd = "perl $bin/get_DEtable.pl $edgeR_output topDE\.txt ifDE\.txt $spaired >$edgeR_output/DE\_table\.txt";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}



###
### sub-region
###

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    printtime();
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}

sub helpm {
  print STDERR "\nSynopsis: RTrace.pl --runlevel 1 --lanename <sample1> --root <dir_root> --anno <dir_anno> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 2 --lanename <sample1> --root <dir_root> --anno <dir_anno> --patient <ID> --tissue <type> --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 3-4 --lanename <sample1> --root <dir_root> --anno <dir_anno> --RA 1 --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 5 --lanename <sample1> --root <dir_root> --anno <dir_anno> --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 7 --root <dir_root> --anno <dir_anno> --priordf 1 2>>run.log\n\n";
  print STDERR "GENERAL OPTIONS (MUST SET):\n\t--runlevel\tthe steps of runlevel, from 1-7, either rl1-rl2 or rl. See below for options for each runlevel.\n";
  print STDERR "\t--lanename\tthe name of the lane needed to be processed (must set for runlevel 1-5)\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$bin/../PIPELINE/, MUST set using other dir)\n";
  print STDERR "\t--anno\t\tthe annotations directory (default is \$bin/../ANNOTATION/, MUST set using other dir)\n";
  print STDERR "\t--patient\tthe patient id, which will be written into the target file for edgeR ()\n";
  print STDERR "\t--tissue\tthe tissue type name (like \'normal\', \'cancer\'), for the target file\n\n";

  print STDERR "CONTROL OPTIONS FOR EACH RUNLEVEL:\n";
  print STDERR "runlevel 1: quality checking and insert size estimatiion using part of reads\n";
  print STDERR "\t--AB\t\tsplit reads up to generate non-overlapping paired-end reads.\n";
  print STDERR "\t--fqreid\trename the fastq id in case of \/1N, only for gsnap mapping.\n";
  print STDERR "\t--QC\t\tdo the quality check of reads, will stop the pipeline once it is finished.\n";
  print STDERR "\t--readlen\tthe sequenced read length.\n";
  print STDERR "\t--trimedlen\tthe read length after trimming (default 80). set it the same as readlen for no trimming\n";

  print STDERR "\nrunlevel 2: mapping and report of mapping statistics\n";
  print STDERR "\t--mapper\tthe mapper used in runlevel2, now support \'tophat1\', \'tophat2\' or \'gsnap\' (default).\n";
  print STDERR "\t--seglen\tthe segment length for tophat mapping (default 25)\n";
  print STDERR "\t--gf\t\tthe graphical format in mapping report, \'png\' (default) or \'pdf\' (when a x11 window is not available)\n";
  print STDERR "\t--WIG\t\tgenerate a big wiggle file in run-level 2.\n";
  print STDERR "\t--insertmean\tthe mean insert size of read mates (not required, can be decided automatically)\n";
  print STDERR "\t--insertsd\tthe SD of insert size of read mates (not required, can be decided automatically)\n";

  print STDERR "\nrunlevel 3: selecting anormalous/breakpoint-surrouding reads and preforming the assembly\n";
  print STDERR "\t--SM\t\tforce to do a second mapping of trimed initially unmapped reads reported by TopHat, for run-level 3\n";
  print STDERR "\t--RA\t\tuse regional assembly for runlevel 3. Set to 1: independent regional assembly (recommended); 2: Columbus assembly.\n\t\t\tDefault is 1. When using gsnap as mapper in runlevel 2, must set to 1 or 2.\n";

  print STDERR "\nrunlevel 4: detection of fusion candidates\n";
  print STDERR "\t--BT\t\tset if use BLAST in run-level 4 (default: use BLAT).\n";

  print STDERR "\nrunlevel 5: run cufflinks for gene/isoform quantification\n";
  print STDERR "runlevel 6: run cuffdiff for diffrential gene/isoform expression analysis\n";
  print STDERR "\t--gtf-guide\tuse gtf guided assembly method in run-level 5 (default: FALSE).\n";
  print STDERR "\t--known-trans\twhich known transcript annotation to be used for cufflinks, either \'ensembl\' (default) or 'refseq'.\n";
  print STDERR "\t--frag-bias\tcorrect fragmentation bias in run-level 5 (default: FALSE).\n";
  print STDERR "\t--upper-qt\tupper-quantile normalization in run-level 5 (default: FALSE).\n";

  print STDERR "\nrunlevel 7: run edgeR for diffrential gene expression analysis\n";
  print STDERR "\t--priordf\tthe prior.df parameter in edgeR, which determines the amount of smoothing of tagwise dispersions towards the common dispersion.\n\t\t\tThe larger the value for prior.df, the more smoothing. A prior.df of 1 gives the common likelihood the weight of one observation. \n\t\t\tDefault is 10. Set it smaller for large sample size (i.e., set to 1 for more than 20 replicates).\n";
  print STDERR "\t--spaired\twhether the experiment is a paired normal-disease design, 1 means yes (default), 0 for no.\n";

  print STDERR "\nOTHER OPTIONS\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--quiet\t\tdo not print the command line calls and time information\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n\n";

  print STDERR "Runlevel dependencies (->): 4->3->2->1, 6->5->2->1, 7->2->1\n\n";
  exit 0;
}

sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub mateorder {
  my @r = @_;
  my @tmp;

  if (scalar(@r) != 2) {
    print STDERR "there are not two mate read files, exit.\n";
    exit 1;
  }

  foreach my $r (@r){
    if ($r =~ /_R?1\./) {
      $tmp[0] = $r;
    }
    elsif ($r =~ /_R?2\./) {
      $tmp[1] = $r;
    }
    else {
      print STDERR "the mate read file dosen't contain _1 _2 information, exit.\n";
      exit 1;
    }
  }

  return @tmp;
}

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub find_jumpers {

  my $file = shift;
  my $jumpers = shift;

  my $breakpoint_id = '';
  my $jumper = 0;
  my $flag_ra = 0;

  open IN, "$file" || die "Error: could not open $file.\n";
  while ( <IN> ) {
    if ($_ =~ /^(\d+)\tchr\w+/) {
      if ($flag_ra == 0) {      #now the jumper should be set
        $breakpoint_id = $1;
        $jumpers->{$breakpoint_id} = $jumper;
        $flag_ra = 1;
      }
    } else {
      $flag_ra = 0 if ($flag_ra == 1); #old breakpoints
    }
    $jumper = tell IN;
  }
  close IN;
}

sub regional_assembly {

  my $h = shift;
  my $p = shift;
  my $id = shift;
  my $out = shift;
  die "Error: $h is not a filehandle.\n" unless defined fileno($h);
  open OUT, ">$out";

  seek($h, $p, 0);
  my $flag_ra = 0;
  while ( <$h> ) {
    if ($_ =~ /^(\d+)\tchr\w+/) {
      if ($flag_ra == 0) {
        if ($id != $1) {
          print STDERR "Error: wrong jump position.\n";
          exit 22;
        }
      } else {
        last;
      }
    } else {
      $flag_ra = 1 if ($flag_ra == 0);
      print OUT "$_";
    }
  } #while

  close OUT;
}
