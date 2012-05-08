#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my $noexecute  = 0;
my $runlevels  = 0;
my $readlen    = 0;
my $trimedlen  = 0;
my $seg_len    = 25;
my $ins_mean   = 0;
my $ins_mean_true = 0;
my $ins_mean_AB= 0;
my $ins_sd     = 20;
my %runlevel;
my $lanename;
my $threads    = 1;
my $help;
my $AB;      #cut reads in AB (it is not necessary)
my $QC;      #quality check
my $SM;      #second mapping
my $BT;      #using Blast instead of BLAT for runlevel-4
my $force;   #force
my $bigWig;  #wiggle file
my $root = "$RealBin/../PIPELINE";
my $anno = "$RealBin/../ANNOTATION";
my $bin  = "$RealBin/";
my $qual_zero = 33;
my $qual_move = 0;


if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}

GetOptions(
           "lanename=s"   => \$lanename,
           "runlevel=s"   => \$runlevels,
           "noexecute"    => \$noexecute,
           "readlen=i"    => \$readlen,
           "trimedlen=i"  => \$trimedlen,
           "seglen=i"     => \$seg_len,
           "insertmean=i" => \$ins_mean,
           "insertsd=i"   => \$ins_sd,
           "threads=i"    => \$threads,
           "AB"           => \$AB,
           "QC"           => \$QC,
           "SM"           => \$SM,
           "BT"           => \$BT,
           "WIG"          => \$bigWig,
           "force"        => \$force,
           "root=s"       => \$root,
           "anno=s"       => \$anno,
           "help|h"       => \$help,
          );

#print help
helpm() if ($help);

### Annotation paths------------------------------------------------------
my $bowtie_index = "$anno/bowtie_index/hg19/hg19";
my $tophat_trans_index = "$anno/bowtie_index/hg19_trans/hg19_konw_ensemble_trans";
my $gene_annotation = "$anno/hg19\.ensembl\-for\-tophat\.gff";
my $ensemble_gene = "$anno/UCSC\_Ensembl\_Genes\_hg19";
my $refseq_gene = "$anno/RefSeq\_Genes\_hg19";
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

my @lanefile = bsd_glob("$root/$lanename*fastq\.gz");
if (scalar(@lanefile) > 0){
  foreach my $file (@lanefile){
    (my $newfile = $file) =~ s/fastq\.gz/fq\.gz/;
    my $cmd = "mv $file $newfile";
    RunCommand($cmd,$noexecute);
  }
}
@lanefile = ();
@lanefile = bsd_glob("$root/$lanename*fq\.gz");

my $lanepath  = "$root/$lanename";

printtime();
print STDERR "####### preparing directories #######\n\n";

unless (-e "$lanepath/01_READS") {
  my $cmd = "mkdir -p $lanepath/01_READS";
  RunCommand($cmd,$noexecute);
}

if (scalar(@lanefile) > 0){
  foreach my $read_file (@lanefile) {
    my $cmd = "mv $read_file $lanepath/01_READS/";
    RunCommand($cmd,$noexecute);
  }
}

if ($readlen == 0 or $trimedlen == 0) { #read length or trimed length not set
     my @original_read_files = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");
     my $first_second_line = `gzip -d -c $original_read_files[0] | head -2 | grep -v "^@"`;
     $readlen = length($first_second_line) - 1;
     $trimedlen = $readlen;
     print STDERR "read length and trimed length are not set, will take the original read length ($readlen bp) for both (no trimming).\n";
}

###
###runlevel1: trim the reads and insert size detection using spiked in reads
###

$runlevels = 1;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my @qc_files = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");
  (my $qc_out1 = $qc_files[0]) =~ s/\.gz$/\.qc/;
  (my $qc_out2 = $qc_files[1]) =~ s/\.gz$/\.qc/;
  unless (-e "$qc_out1") {
    my $cmd = "gzip -d -c $qc_files[0] | $bin/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$noexecute);
  }
  unless (-e "$qc_out2") {
    my $cmd = "gzip -d -c $qc_files[1] | $bin/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$noexecute);
  }
  if ($QC) {
    print STDERR "quality check finished, please check the quality file manually.\n";
    exit;
  }

  my @read_files;
  if ( $trimedlen != $readlen ) {  #trimming
    my @trimed_read_files = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz");
    if ( scalar(@trimed_read_files) == 0 ) {
      my @ori_read_files = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");
      foreach my $read_file (@ori_read_files) {
        my $read_out = $read_file;
        $read_out =~ s/fq\.gz$/trimed\.fq\.gz/;
        if ($read_file =~ /_1\./) {
          $read_files[0] = $read_out;
        }
        else {
          $read_files[1] = $read_out;
        }
        my $cmd = "gzip -d -c $read_file | $bin/fastx_trimmer -l $trimedlen -z -Q33 -o $read_out";
        RunCommand($cmd,$noexecute);
      }
    }
    else { #if trimed read file exists
      @read_files = @trimed_read_files;
      @read_files = mateorder(@read_files);
    }
  }
  else {
    @read_files = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");
    @read_files = mateorder(@read_files);
  }

  ##### AB trimming, which is not fully tested################
  if ($AB) {  # do the AB treating
     print STDERR "AB is defined\n";
     my @AB_read_files = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq\.gz");
     my @AB_read_files_ori = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq");

     if ( scalar(@AB_read_files) == 0 and scalar(@AB_read_files_ori) == 0) {
       my ($read_1, $read_2, $AB_1, $AB_2);
       foreach my $read_file (@read_files){
         if ($read_file =~ /_1\./) {
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
       RunCommand($cmd,$noexecute);
     }
     if ( scalar(@AB_read_files) == 0 and scalar(@AB_read_files_ori) == 2) {
       foreach my $AB_read_files_ori (@AB_read_files_ori) {
          my $cmd = "gzip $AB_read_files_ori";
          RunCommand($cmd,$noexecute);
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
    RunCommand($cmd,$noexecute);
  }

  my @spiked_in = bsd_glob("$lanepath/01_READS/$lanename*spikedin.fq");
  if (scalar(@spiked_in) == 0) {
    foreach my $read_file (@read_files) {
      my $spiked_in = $read_file;
      $spiked_in =~ s/fq\.gz$/spikedin\.fq/;
      push(@spiked_in, $spiked_in);
    }
    my $cmd = "perl $bin/select_reads_with_no_n.pl $read_files[0] $read_files[1] 200000 >$spiked_in[0] 2>$spiked_in[1]";
    RunCommand($cmd,$noexecute);
  }

  #do the pair-end mapping of spiked_in reads
  unless (-e "$lanepath/00_TEST/$lanename\.spikedin\.hits") {
    my $cmd= "bowtie -v 2 $bowtie_index -1 $spiked_in[0] -2 $spiked_in[1] $lanepath/00_TEST/$lanename\.spikedin\.hits";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/00_TEST/$lanename\.spikedin\.fragmentlength") {
    my $real_len = $trimedlen;
    $real_len = $trimedlen/2 if ($AB);
    my $cmd = "perl $bin/spike_in.pl $lanepath/00_TEST/$lanename\.spikedin\.hits $real_len >$lanepath/00_TEST/$lanename\.spikedin\.fragmentlength";
    RunCommand($cmd,$noexecute);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}

###
###runlevel1.5: deciding fragment/insert size
###
printtime();
print STDERR "####### insert mean and sd calculation #######\n\n";

if ( -e "$lanepath/01_READS/$lanename\_1\.fq\.qc" ) { #decide the quality shift
     open QC, "<$lanepath/01_READS/$lanename\_1\.fq\.qc";
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
}
else {
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
    my $ins_th = round($real_len*0.25);
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


###
###runlevel2: do the mapping and generate the statistics
###

$runlevels = 2;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$lanepath/02_MAPPING") {
    my $cmd = "mkdir -p $lanepath/02_MAPPING";
    RunCommand($cmd,$noexecute);
  }

  my @reads;
  if ($AB) {
    @reads = bsd_glob("$lanepath/01_READS/$lanename*AB\.fq\.gz");       #AB reads
  }
  elsif ($trimedlen != $readlen) {
    @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz");   #trimmed reads
  }
  else {
    @reads = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");    #original reads
  }
  @reads = mateorder(@reads);

  my $real_len = $trimedlen;
  $real_len = $trimedlen/2 if ($AB);

  my $real_ins_mean = $ins_mean;
  $real_ins_mean = $ins_mean_AB if ($AB);

  my $fragment_length = 2*$real_len + $real_ins_mean;

  #do the mapping of pair - end reads
  unless (-e "$lanepath/02_MAPPING/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam") {
    my $cmd = "tophat --output-dir $lanepath/02_MAPPING --mate-inner-dist $real_ins_mean --mate-std-dev $ins_sd --library-type fr-unstranded -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat_trans_index $bowtie_index $reads[0] $reads[1]";
    RunCommand($cmd,$noexecute);
  }

  #do the statistics
  unless (-e "$lanepath/03_STATS") {
    my $cmd = "mkdir -p $lanepath/03_STATS";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.mapping\.stats")  {
    my $cmd = "$bin/Rseq_bam_stats --mapping $lanepath/02_MAPPING/accepted_hits\.bam --writer $lanepath/02_MAPPING/accepted_hits\.unique\.bam --arp $lanepath/03_STATS/$lanename\.arp >$lanepath/03_STATS/$lanename\.mapping\.stats";
    RunCommand($cmd,$noexecute);
  }

  my $mapping_stats_line_number = `wc -l $lanepath/03_STATS/$lanename.mapping.stats`;
  $mapping_stats_line_number =~ s/^(\d+).*$/\1/;
  chomp($mapping_stats_line_number);
  if ($mapping_stats_line_number == 12){
    my $total_reads = `gzip -d -c $lanepath/01_READS/$lanename\_1.fq.gz | wc -l`;
    $total_reads /= 4;
    open STATS, ">>$lanepath/03_STATS/$lanename\.mapping\.stats" || die "can not open $lanepath/03_STATS/$lanename\.mapping\.stats\n";
    print STATS "total_frag: $total_reads\n";
    close STATS;
  }

  unless (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam")  {
    my $cmd = "samtools sort $lanepath/02_MAPPING/accepted_hits\.unique\.bam $lanepath/02_MAPPING/accepted_hits\.unique\.sorted";
    RunCommand($cmd,$noexecute);
  }

  if (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam" and -e "$lanepath/02_MAPPING/accepted_hits\.unique\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.unique\.bam -f";
    RunCommand($cmd,$noexecute);
  }
  if (-e "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam" and -e "$lanepath/02_MAPPING/accepted_hits\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.bam -f";
    RunCommand($cmd,$noexecute);
  }

  if ($bigWig) { #generate wiggle file
    unless (-e "$lanepath/03_STATS/$lanename\.bw"){
      if (-e "$lanepath/03_STATS/$lanename\.bedgraph") {
         my $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$lanename\.bedgraph $anno/chrom_sizes.hg19 $lanepath/03_STATS/$lanename\.bw";
         RunCommand($cmd,$noexecute);
      }
      else {
         my $cmd = "genomeCoverageBed -ibam $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam -bg -split -g $anno/chrom_sizes.hg19 >$lanepath/03_STATS/$lanename\.bedgraph";
         RunCommand($cmd,$noexecute);
         $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$lanename\.bedgraph $anno/chrom_sizes.hg19 $lanepath/03_STATS/$lanename\.bw";
         RunCommand($cmd,$noexecute);
      }
      if (-e "$lanepath/03_STATS/$lanename\.bedgraph" and -e "$lanepath/03_STATS/$lanename\.bw") {
         my $cmd = "rm $lanepath/03_STATS/$lanename\.bedgraph -f";
         RunCommand($cmd,$noexecute);
      }
    }
  }

  unless (-e "$lanepath/03_STATS/$lanename\.expr") {
    my $cmd = "$bin/Rseq_bam_reads2expr --region $ensemble_gene --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam --posc $lanepath/03_STATS/$lanename\.pos\.gff --chrmap $lanepath/03_STATS/$lanename\.chrmap --lbias $lanepath/03_STATS/$lanename\.lbias >$lanepath/03_STATS/$lanename\.expr";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.RefSeq\.expr") {
    my $cmd = "$bin/Rseq_bam_reads2expr --region $refseq_gene --mapping $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam >$lanepath/03_STATS/$lanename\.RefSeq\.expr";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.cate") {
    my $cmd = "perl $bin/cate.pl $lanepath/03_STATS/$lanename\.expr $gene_annotation >$lanepath/03_STATS/$lanename\.cate";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.ins") {
    my $cmd = "samtools view -f 0x2 $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam | cut -f 9 | awk \'\$1\>0 \&\& \$1\<500\' >$lanepath/03_STATS/$lanename\.ins";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/03_STATS/$lanename\.report/$lanename\.report\.html") {
    my $cmd = "R CMD BATCH --no-save --no-restore "."\'--args path=\"$lanepath\" lane=\"$lanename\" anno=\"$anno\"\ src=\"$bin\" readlen=$real_len' $bin/html_report.R $lanepath/03_STATS/R\_html\.out";
    RunCommand($cmd,$noexecute);
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

  #get the ARP from unmapped by tophat ############
  unless ( -e "$lanepath/01_READS/$lanename\.ARP\.fq\.gz" ) {

    my @ARP = bsd_glob("$lanepath/01_READS/$lanename*ARP\.fq\.gz");

    if (scalar(@ARP) != 2) {
      my @reads;
      if ($trimedlen != $readlen) {
        @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
      } else {
        @reads = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz"); #original reads
      }
      @reads = mateorder(@reads);

      unless (-e "$lanepath/02_MAPPING/unmapped_left\.fq\.z"){
        print STDERR "unmapped read file does not exist!!!\n";
        exit 22;
      }

      my $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/03_STATS/$lanename\.arp --unmap $lanepath/02_MAPPING/unmapped_left\.fq\.z --readfile1 $reads[0] --readfile2 $reads[1]";
      $cmd .= " --AB" if ($AB);
      RunCommand($cmd,$noexecute);

      my @ori_ARP_file = bsd_glob("$lanepath/01_READS/$lanename*ARP\.fq");
      foreach my $ori_ARP_file (@ori_ARP_file){
         my $cmd = "gzip $ori_ARP_file";
         RunCommand($cmd,$noexecute);
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
          RunCommand($cmd,$noexecute);
        }
      }

      #second mapping
      my @ARP_trimed36 = bsd_glob("$lanepath/01_READS/$lanename*ARP*trimed36\.fq\.gz");
      @ARP_trimed36 = mateorder(@ARP_trimed36);
      unless (-e "$lanepath/02_MAPPING/SecondMapping/") {
        my $cmd = "mkdir -p $lanepath/02_MAPPING/SecondMapping/";
        RunCommand($cmd,$noexecute);
      }
      unless ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam" or -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.unique\.bam" ) { #second mapping using gsnap
        if ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam" ) { #no bam but sam, need to compress it
          my $cmd = "samtools view -Sb $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -o $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam";
          RunCommand($cmd,$noexecute);
        } else {

          my $cmd = "gsnap -d hg19 -D $gmap_index -B 5 --gunzip --format=sam --nthreads=$threads -s $gmap_splicesites --max-mismatches=2 --npaths=10 --trim-mismatch-score=0 --trim-indel-score=0 --quality-zero-score=$qual_zero --quality-print-shift=$qual_move $ARP_trimed36[0] $ARP_trimed36[1] >$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam";
          RunCommand($cmd,$noexecute);
          $cmd = "samtools view -Sb $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -o $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam";
          RunCommand($cmd,$noexecute);
        }
      }

      if ( -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam" and -e "$lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam"){
         my $cmd = "rm $lanepath/02_MAPPING/SecondMapping/accepted_hits\.sam -f";
         RunCommand($cmd,$noexecute);
      }

      unless ( -e "$lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.stats" )  {
        my $cmd = "$bin/Rseq_bam_stats --mapping $lanepath/02_MAPPING/SecondMapping/accepted_hits\.bam --writer $lanepath/02_MAPPING/SecondMapping/accepted_hits\.unique\.bam --arp $lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.arp >$lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.stats";
        RunCommand($cmd,$noexecute);
      }

      #now get the real arp after second mapping
      my @ARP_sm = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq\.gz");
      if (scalar(@ARP_sm) != 2) {
        my @reads;
        if ($trimedlen != $readlen) {
          @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
        } else {
          @reads = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz"); #original reads
        }
        @reads = mateorder(@reads);

        my $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/02_MAPPING/SecondMapping/$lanename\.secondmapping\.arp --readfile1 $reads[0] --readfile2 $reads[1] --SM";
        $cmd .= " --AB" if ($AB);
        RunCommand($cmd,$noexecute);
      }

      my @ori_ARP_SM_file = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq");
      foreach my $ori_ARP_SM_file (@ori_ARP_SM_file) {
         my $cmd = "gzip $ori_ARP_SM_file";
         RunCommand($cmd,$noexecute);
      }

      @ARP = bsd_glob("$lanepath/01_READS/$lanename*ARP\.secondmapping\.fq\.gz");
      @ARP = mateorder(@ARP);
    }

    # shuffle ARP reads to a single file
    my $cmd = "perl $bin/shuffleSequences_fastq.pl $ARP[0] $ARP[1] $lanepath/01_READS/$lanename\.ARP\.fq";
    RunCommand($cmd,$noexecute);
    $cmd = "gzip $lanepath/01_READS/$lanename\.ARP\.fq";
    RunCommand($cmd,$noexecute);
  }

  #do the velveth
  unless (-e "$lanepath/04_ASSEMBLY") {
    my $cmd = "mkdir -p $lanepath/04_ASSEMBLY";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/04_ASSEMBLY/Roadmaps") {
    my $cmd = "velveth $lanepath/04_ASSEMBLY/ 21 -fastq.gz -shortPaired $lanepath/01_READS/$lanename\.ARP\.fq\.gz";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/04_ASSEMBLY/Graph2") {
    my $cmd = "velvetg $lanepath/04_ASSEMBLY/ -ins_length $ins_mean -cov_cutoff auto -exp_cov auto -read_trkg yes -scaffolding no -min_contig_lgth 100";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/04_ASSEMBLY/transcripts.fa") {
    my $cmd = "oases $lanepath/04_ASSEMBLY/ -ins_length $ins_mean -ins_length_sd $ins_sd -unused_reads yes -scaffolding no ";
    RunCommand($cmd,$noexecute);
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
    RunCommand($cmd,$noexecute);
  }

  if ($BT) {  #using BLAST alternative
    unless (-e "$lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast") {
      my $cmd = "$bin/blastn -query $lanepath/04_ASSEMBLY/transcripts\.fa -db $anno/human\.rna\.fna\.blastdb -outfmt 7 -num_threads $threads -out $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq") {
      my $cmd = "perl $bin/pick_fusion_transcripts_from_BLAT.pl --refseq $anno/human\.rna\.fna --transcript $lanepath/04_ASSEMBLY/transcripts.fa --blat $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blast --lanename $lanename --lanepath $lanepath >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap"){
      my $cmd = "gmap -D $gmap_index -d hg19 --format=psl -t $threads $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq"){
      my $cmd = "perl $bin/filter_out_FP_from_blatps.pl --fusion_bf_seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq --fusion_bf $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration --fusion_bf_blat $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.gmap >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq";
      RunCommand($cmd,$noexecute);
    }

  }

  else {  #using BLAT
    unless (-e "$lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat") {
      my $cmd = "blat $anno/human\.rna\.fna $lanepath/04_ASSEMBLY/transcripts.fa -out=blast9 $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq") {
      my $cmd = "perl $bin/pick_fusion_transcripts_from_BLAT.pl --refseq $anno/human\.rna\.fna --transcript $lanepath/04_ASSEMBLY/transcripts.fa --blat $lanepath/05_FUSION/$lanename\.transcripts\.refseq\.blat --lanename $lanename --lanepath $lanepath >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.blat"){
      my $cmd = "blat $anno/hg19\_UCSC\.2bit $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration.genome\.blat";
      RunCommand($cmd,$noexecute);
    }

    unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq"){
      my $cmd = "perl $bin/filter_out_FP_from_blatps.pl --fusion_bf_seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.seq --fusion_bf $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration --fusion_bf_blat $lanepath/05_FUSION/$lanename\.fusion_transcirpts_before_filtration\.genome\.blat >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq";
      RunCommand($cmd,$noexecute);
    }
  }

  #build fusion candidate index
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index\.1\.ebwt"){
    my $cmd = "bowtie-build --quiet $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index";
    RunCommand($cmd,$noexecute);
  }

  #map reads to the fusion candidate index
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie"){
    my @reads;
    if ($trimedlen != $readlen) {
      @reads = bsd_glob("$lanepath/01_READS/$lanename*trimed\.fq\.gz"); #trimmed reads
    } else {
      @reads = bsd_glob("$lanepath/01_READS/$lanename\_[12]\.fq\.gz");  #original reads
    }
    @reads = mateorder(@reads);

    my $cmd = "gzip -d -c $reads[0] $reads[1] | bowtie -v 1 -k 10 -m 10 -p $threads $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq\.index \- $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie";
    RunCommand($cmd,$noexecute);
  }

  #get fusion coverage
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov") {
    my $cmd = "perl $bin/get_fusion_coverage.pl --type pair --mappingfile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie --readlength $trimedlen --geneanno $gene_annotation --ensrefname $anno/Ensembl\_Ref\_Name\.tsv --locname $anno/Name2Location\.hg19 --refgene $anno/refgenes\.hg19 --repeatmasker $anno/UCSC\_repeats\_hg19\.gff --accepthits $lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam --encomcov $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov.enco >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    RunCommand($cmd,$noexecute);
  }

  #visualize coverage
  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis") {
    my $cmd = "perl $bin/further_processing_for_read_visualization.pl --fusionseqfile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.seq --coveragefile $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov --readlength $trimedlen >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis";
    RunCommand($cmd,$noexecute);
  }

  unless (-e "$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.list"){
    my $cmd = "grep \"^\#\" $lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis \| sort -k 13,13nr -k 6,6d -k 8,8d >$lanepath/05_FUSION/$lanename\.fusion_transcirpts_after_filtration\.list";
    RunCommand($cmd,$noexecute);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
### sub-region
###

sub RunCommand {
  my ($command,$noexecute) = @_ ;
  print STDERR "$command\n";
  unless ($noexecute) {
    system($command);
  }
}

sub helpm {
  print STDERR "usage: $0 [options]\n\nOptions:\n\t--runlevel\tthe steps of runlevel, from 1-4, either rl1-rl2 or rl, see below\n";
  print STDERR "\t\t\t1: trim reads and decide insert size using spiked in reads\n";
  print STDERR "\t\t\t2: do the mapping and generate the statistics\n";
  print STDERR "\t\t\t3: select anormalous read pairs and do the assembly\n";
  print STDERR "\t\t\t4: get fusion candidates and visualize the result\n";
  print STDERR "\t--lanename\tthe name of the lane needed to be processed (must set for all runlevels)\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--readlen\tthe sequenced read length (default 95)\n";
  print STDERR "\t--AB\t\tsplit reads up to generate non-overlapping paired-end reads.\n";
  print STDERR "\t--QC\t\tdo the quality check of reads, will stop the pipeline once it is finished.\n";
  print STDERR "\t--SM\t\tforce to do a second mapping of trimed initially unmapped reads (using tophat)\n";
  print STDERR "\t--BT\t\tset if use BLAST in run-level 4.\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$bin/../PIPELINE/, MUST set using other dir)\n";
  print STDERR "\t--anno\t\tthe annotations directory (default is \$bin/../ANNOTATION/, MUST set using other dir)\n";
  print STDERR "\t--trimedlen\tthe read length after trimming (default 80). set it the same as readlen for no trimming\n";
  print STDERR "\t--seglen\tthe segment length for tophat mapping (default 27)\n";
  print STDERR "\t--insertmean\tthe mean insert size of read mates (not required, can be decided automatically)\n";
  print STDERR "\t--insertsd\tthe SD of insert size of read mates (not required, can be decided automatically)\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n\n\n";
  exit 0;
}

sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub mateorder {
  my @r = @_;
  my @tmp;

  if (scalar(@r) != 2){
    print STDERR "there are not two mate read files, exit.\n";
    exit 1;
  }

  foreach my $r (@r){
    if ($r =~ /_1\./) {
      $tmp[0] = $r;
    }
    elsif ($r =~ /_2\./) {
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
