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
my $ins_sd     = 20;
my %runlevel;
my $sampleName;
my $runID      = '';
my $threads    = 1;
my $help;
my $QC;      #quality check
my $qcOFF;
my $SM;      #second mapping
my $GMAP;    #using GMAP instead of BLAT for runlevel-4
my $RA         = 1;      #regional assembly
my $idra       = 0;      #for a specific breakpoint id
my $consisCount = 5;     #the threshold for the number of consistent mate pairs with discordant maping
my $force;   #force
my $bigWig;  #wiggle file
my $gtf_guide_assembly;  #for cufflinks
my $known_trans = 'ensembl';         #for cufflinks
my $frag_bias_correct;   #for cufflinks
my $upper_quantile_norm; #for cufflinks
my $root = "$RealBin/../PIPELINE";
my $anno = "$RealBin/../ANNOTATION";
my $species = 'hg19';
my $readpool = 'SRP';
my $bin  = "$RealBin/";
my $qual_zero = 33;
my $qual_move = 0;
my $fq_reid; #rename fastq read id (for gsnap)
my $priordf = 10;         #for edgeR
my $pairDE1 = 'N';        #for edgeR
my $pairDE2 = 'T';        #for edgeR
my $spaired = 1;          #for edgeR
my $patient; #the patient id for edgeR DE test
my $tissue;   #the tissue type for edgeR DE test
my $gf = "png"; #the format used in html report
my $merge = ''; #whether merge for multiple runs
my @merge;
my $bzip;  #to allow bzip compressed fastq files
my $Rbinary = 'R';
my $customMappedBam = '';
my $seqType = 'p';   #whether paired-end or single-end


if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}

GetOptions(
           "sampleName=s" => \$sampleName,
           "runID=s"      => \$runID,
           "merge=s"      => \$merge,
           "runlevel=s"   => \$runlevels,
           "seqType=s"    => \$seqType,
           "noexecute"    => \$noexecute,
           "quiet"        => \$quiet,
           "readlen=i"    => \$readlen,
           "trimedlen=i"  => \$trimedlen,
           "seglen=i"     => \$seg_len,
           "mapper=s"     => \$mapper,
           "insertmean=i" => \$ins_mean,
           "insertsd=i"   => \$ins_sd,
           "threads=i"    => \$threads,
           "consisCount=i"=> \$consisCount,
           "QC"           => \$QC,
           "qcOFF"        => \$qcOFF,
           "SM"           => \$SM,
           "GMAP"         => \$GMAP,
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
           "species=s"    => \$species,
           "readpool=s"   => \$readpool,
           "priordf=i"    => \$priordf,
           "spaired=i"    => \$spaired,
           "pairDE1=s"    => \$pairDE1,
           "pairDE2=s"    => \$pairDE2,
           "patient=s"    => \$patient,
           "tissue=s"     => \$tissue,
           "bzip"         => \$bzip,
           "Rbinary=s"    => \$Rbinary,
           "help|h"       => \$help,
           "customBam=s"  => \$customMappedBam,
          );

#print help
helpm() if ($help);

### Annotation paths------------------------------------------------------
my $bowtie_index = "$anno/$species/$species\.bowtie_index/$species/$species";
my $bowtie2_index = "$anno/$species/$species\.bowtie2_index/$species/$species";
my $genome_fasta = (-e "$bowtie_index\.fa")?"$bowtie_index\.fa":"$bowtie2_index\.fa";
my $tophat_trans_index = "$anno/$species/$species\.bowtie_index/$species\_trans/$species\_known_ensemble_trans";
my $tophat2_trans_index = "$anno/$species/$species\.bowtie2_index/$species\_trans/$species\_known_ensemble_trans";
my $gene_annotation = "$anno/$species/$species\.transcripts\_Ensembl\.gff";
my $gene_annotation_gtf = "$anno/$species/$species\.transcripts\_Ensembl\.gtf";
my $ensembl_gene_bed = "$anno/$species/$species\.genes\_Ensembl\.bed12";
my $gencode_genemap = "$anno/$species/$species\.gencode/gencode\.v14\.genemap" if ($species eq 'hg19');
my $gencode_gene_bed = "$anno/$species/$species\.gencode/gencode\.v14\.annotation\.gene\.bed12" if ($species eq 'hg19');
my $ensembl_genemap = "$anno/$species/$species\.biomart\.txt";
my $refseq_gene = "$anno/$species/$species\.genes\_RefSeq\.bed12";
my $refseq_gene_gtf = "$anno/$species/$species\.refGene\.gtf";
my $gmap_index = "$anno/$species/$species\.gmap\_index/";
my $gmap_splicesites = "$gmap_index/$species/$species\.splicesites\.iit";
my $chromosomeSize = "$anno/$species/$species\.chromosome\_size\.txt";
my $repeatMasker = "$anno/$species/$species\.repeats\_UCSC\.gff";
my $selfChain = "$anno/$species/$species\.SelfChain\_UCSC\.txt";
my $blatDatabase = "$anno/$species/$species\.genome\_UCSC\.2bit";
#-------------------------------------------------------------------------

### Frequently used names-------------------------------------------------
my $mappedBam = "accepted_hits\.mapped\.sorted\.bam";
#-------------------------------------------------------------------------

#decompression option-----------------------------------------------------
my $decompress = "gzip -d -c";
my $compress = "gzip";
my $zipSuffix = "gz";
if ($bzip) {
  $decompress = "bzip2 -d -c";
  $compress = "bzip2";
  $zipSuffix = "bz2";
}
#-------------------------------------------------------------------------

if ($runlevels) { #true runlevels
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
} else {
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
} else {
  $readpool = $root if $readpool eq 'SRP';
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

if (defined $sampleName) {

  printtime();
  print STDERR "####### lane name is set to $sampleName #######\n\n";

  my @lanefile;
  if ($bzip) {
    @lanefile = bsd_glob("$readpool/*$sampleName*fq\.bz2");
    if (scalar(@lanefile) == 0) {
      @lanefile = bsd_glob("$readpool/*$sampleName*fastq\.bz2");
    }
  } else {
    @lanefile = bsd_glob("$readpool/*$sampleName*fq\.gz");
    if (scalar(@lanefile) == 0) {
      @lanefile = bsd_glob("$readpool/*$sampleName*fastq\.gz");
    }
  }

  foreach my $lfile (@lanefile){
     print STDERR "lanefile:\t$lfile\n";
  }

  if (scalar(@lanefile) == 0) {
    print STDERR "no read files are available, please make sure that they are available under the dir of $readpool.\n";
    exit 22 unless ($force);
  }

  if ($runID ne '') {
    $lanepath = "$root/$sampleName"."\_$runID";
  } else {
    $lanepath = "$root/$sampleName";
  }

  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$lanepath/01_READS") {
    my $cmd = "mkdir -p $lanepath/01_READS";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if (scalar(@lanefile) > 0) {
    if ($seqType =~ /^p/) {
      if ($runID ne '') {
        foreach my $read_file (@lanefile) {
          if ($read_file =~ /[^a-zA-Z0-9]($sampleName)*(\_R?[123]\_$runID)\.f.+?\.([gb]z2?)$/) {
            my $prefix = $1.$2;
            my $suffix = $3;
            my $cmd = "ln -s $read_file $lanepath/01_READS/$prefix\.fq\.$suffix";
            RunCommand($cmd,$noexecute,$quiet) unless (-s "$lanepath/01_READS/$prefix\.fq\.$suffix");
          }
        }
      } else {
        foreach my $read_file (@lanefile) {
          if ($read_file =~ /[^a-zA-Z0-9]($sampleName)*(\_R?[123])\.f.+?\.([gb]z2?)$/) {
            my $prefix = $1.$2;
            my $suffix = $3;
            my $cmd = "ln -s $read_file $lanepath/01_READS/$prefix\.fq\.$suffix";
            RunCommand($cmd,$noexecute,$quiet) unless (-s "$lanepath/01_READS/$prefix\.fq\.$suffix");
          }
        }
      }
    }                           #paired-end
    elsif ($seqType =~ /^s/) {
      if (scalar(@lanefile) == 1) {
        if ($runID ne '') {
          if ($lanefile[0] =~ /$sampleName.*?$runID.*?\.([gb]z2?)$/) {
            my $prefix = $sampleName.'_'.$runID;
            my $suffix = $1;
            my $cmd = "ln -s $lanefile[0] $lanepath/01_READS/$prefix\.fq\.$suffix";
            RunCommand($cmd,$noexecute,$quiet) unless (-s "$lanepath/01_READS/$prefix\.fq\.$suffix");
          } else {
            print STDER "Error: the single-end read fastq file is wierd: $lanefile[0]\n";
            exit 22;
          }
        } else {
          if ($lanefile[0] =~ /$sampleName.*?\.([gb]z2?)$/) {
            my $prefix = $sampleName;
            my $suffix = $1;
            my $cmd = "ln -s $lanefile[0] $lanepath/01_READS/$prefix\.fq\.$suffix";
            RunCommand($cmd,$noexecute,$quiet) unless (-s "$lanepath/01_READS/$prefix\.fq\.$suffix");
          } else {
            print STDER "Error: the single-end read fastq file is wierd: $lanefile[0]\n";
            exit 22;
          }
        }
      } else {
        print STDERR "Error: multiple fastq files for single-end sequencing type.\n";
        exit 22;
      }
    }                           #single-end
    else {
      print STDERR "Error: option --seqType must be set to p or s.\n";
      exit 22;
    }
  }

  if ($readlen == 0 or $trimedlen == 0) { #read length or trimed length not set
     my @original_read_files;
     @original_read_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^p/);
     @original_read_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^s/);
     my $first_second_line = `$decompress "$original_read_files[0]" | head -2 | grep -v "^@"`;
     $readlen = length($first_second_line) - 1;
     $trimedlen = $readlen;
     print STDERR "read length and trimed length are not set, will take the original read length ($readlen bp) for both (no trimming).\n";
  }

  #prepare MERGE list
  if ($merge) { #defined merge runIDs
    @merge = split(/,/, $merge);
    my $merge_list = "$lanepath/$sampleName\.merge\_list";

    if ($runID eq '') {
       print STDERR "option \-\-merge is set, but \-\-runID is not set yet!!!\n";
       exit 22;
    }

    my %already_merged;
    if (-s $merge_list) {
       open MERGE_LIST_I, "<$merge_list";
       while ( <MERGE_LIST_I> ){
          chomp;
          $already_merged{$_} = '';
       }
       close MERGE_LIST_I;
    } #remember already merged run

    open MERGE_LIST, ">>$merge_list";
    foreach my $mergeC (@merge) {
      my $cRunPath = "$root/$sampleName"."\_$mergeC";
      my $cRunBam =  "$cRunPath/02_MAPPING/accepted_hits\.mapped\.sorted\.bam";
      if ( ($merge ne $runID) and !(-s $cRunBam) ) {
        print STDERR "runID $mergeC: the bam file $cRunBam does not exist!!!\n";
        exit 22;
      }
      next if exists($already_merged{$cRunBam}); #not print
      print MERGE_LIST "$cRunBam\n";  #print
    }
    close MERGE_LIST;
  } #merge

}

###
###runlevel1: trim the reads if necessary and insert size detection using spiked in reads
###

$runlevels = 1;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my @qc_files;
  my $qc_out1;
  my $qc_out2;
  if ($seqType =~ /^p/){
    @qc_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix");
    ($qc_out1 = $qc_files[0]) =~ s/\.$zipSuffix$/\.qc/;
    ($qc_out2 = $qc_files[1]) =~ s/\.$zipSuffix$/\.qc/;
  } elsif ($seqType =~ /^s/) {
    @qc_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix");
    ($qc_out1 = $qc_files[0]) =~ s/\.$zipSuffix$/\.qc/;
  } else {
    print STDERR "seqType must be set starting with p:paired-end or s:single-end";
    exit 22;
  }

  unless (-e "$qc_out1") {
    my $cmd = "$decompress $qc_files[0] | $bin/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$noexecute,$quiet) unless ($qcOFF);
  }
  unless (($seqType =~ /^s/) or (-e "$qc_out2")) {
    my $cmd = "$decompress $qc_files[1] | $bin/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$noexecute,$quiet) unless ($qcOFF);
  }
  if ($QC) {
    print STDERR "quality check finished, please check the quality file manually.\n";
    exit;
  }

  my @read_files;
  if ( $trimedlen != $readlen ) {  #trimming
    my @trimed_read_files = bsd_glob("$lanepath/01_READS/$sampleName*trimed\.fq\.$zipSuffix");
    if ( scalar(@trimed_read_files) == 0 ) {
      my @ori_read_files;
      @ori_read_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^p/);
      @ori_read_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^s/);
      foreach my $read_file (@ori_read_files) {
        my $read_out = $read_file;
        $read_out =~ s/fq\.$zipSuffix$/trimed\.fq\.$zipSuffix/;
        if ($seqType =~ /^p/) {
          if ($read_file =~ /\_R?1(\_$runID)?\./) {
            $read_files[0] = $read_out;
          } else {
            $read_files[1] = $read_out;
          }
        }
        my $cmd = "$decompress $read_file | $bin/fastx_trimmer -l $trimedlen -z -Q33 -o $read_out";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } else {                    #if trimed read file exists
      @read_files = @trimed_read_files;
      @read_files = &mateorder(\@read_files, $runID) if ($seqType =~ /^p/);
    }
  } else {  #no trimming
    @read_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^p/);
    @read_files = mateorder(\@read_files, $runID) if ($seqType =~ /^p/);
    @read_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^s/);
  }

  my $total_reads = `$decompress $read_files[0] | wc -l`;
  $total_reads /= 4;
  print STDERR "total number of read pairs: $total_reads\n";


  if ($fq_reid) {  #renbame fastq id
    foreach my $read_file (@read_files) {
       my $reid_file = $read_file;
       $reid_file =~ s/\.fq\.$zipSuffix/\.reid\.fq\.$zipSuffix/;
       unless (-e $reid_file) {
         my $cmd = "perl $bin/fqreid.pl $read_file | $compress >$reid_file";
         RunCommand($cmd,$noexecute,$quiet);
       }
       if (-e $read_file and -e $reid_file){
         my $cmd = "mv -f $reid_file $read_file";
         RunCommand($cmd,$noexecute,$quiet);
       }
    } #foreach read file
  }

  unless (-e "$lanepath/00_TEST") {
    my $cmd = "mkdir -p $lanepath/00_TEST";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if ($seqType =~ /^p/) {
    my @spiked_in = bsd_glob("$lanepath/01_READS/$sampleName*spikedin.fq");
    if (scalar(@spiked_in) == 0) {
      foreach my $read_file (@read_files) {
        my $spiked_in = $read_file;
        $spiked_in =~ s/fq\.$zipSuffix$/spikedin\.fq/;
        push(@spiked_in, $spiked_in);
      }
      my $cmd = "perl $bin/select_reads_with_no_n.pl $read_files[0] $read_files[1] 200000 >$spiked_in[0] 2>$spiked_in[1]";
      RunCommand($cmd,$noexecute,$quiet);
    }

    #do the pair-end mapping of spiked_in reads
    unless (-s "$lanepath/00_TEST/$sampleName\.spikedin\.hits") {
      my $cmd= "bowtie2 -p $threads --no-unal --no-hd --score-min L,-2,-0.15 -x $tophat2_trans_index -1 $spiked_in[0] -2 $spiked_in[1] >$lanepath/00_TEST/$sampleName\.spikedin\.hits";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-s "$lanepath/00_TEST/$sampleName\.spikedin\.fragmentlength") {
      my $real_len = $trimedlen;
      my $cmd = "perl $bin/spike_in.pl $lanepath/00_TEST/$sampleName\.spikedin\.hits >$lanepath/00_TEST/$sampleName\.spikedin\.fragmentlength";
      RunCommand($cmd,$noexecute,$quiet);
    }
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}

###
###runlevel1.5: deciding fragment/insert size
###

if (defined $sampleName) {

  printtime();
  print STDERR "####### insert mean and sd calculation #######\n\n";

  my @quality_check_files;
  @quality_check_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.qc") if ($seqType =~ /^p/);
  @quality_check_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.qc") if ($seqType =~ /^s/);
  if ( scalar(@quality_check_files) == 2 or (scalar(@quality_check_files) == 1 and $seqType =~ /^s/)) { #decide the quality shift
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

  if ($ins_mean == 0 and $seqType =~ /^p/) {

    my $fraginfo = `$Rbinary --no-save --slave \'--args path=\"$lanepath/00_TEST/\" lane=\"$sampleName\.spikedin\"\' < $bin/fragment_length.R`;
    $fraginfo =~ /^(.+)\s(.+)/;
    my $frag_mean = $1; $frag_mean = round($frag_mean);
    my $insert_sd   = $2; $insert_sd =~ s/\n//; $insert_sd = round($insert_sd);

    my $real_len = $trimedlen;

    $ins_sd = $insert_sd;
    $ins_mean = $frag_mean - 2*$real_len;

    $ins_mean_true = $ins_mean;

    print STDERR "insert mean: $ins_mean\tinsert_sd: $ins_sd\n";

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
  if ($trimedlen != $readlen) {
    @reads = bsd_glob("$lanepath/01_READS/$sampleName*trimed\.fq\.$zipSuffix");   #trimmed reads
  }
  else {
    @reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^p/);    #original reads
    @reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix") if ($seqType =~ /^s/);
  }
  @reads = &mateorder(\@reads, $runID) if ($seqType =~ /^p/);

  my $real_len = $trimedlen;

  my $real_ins_mean = $ins_mean;

  my $fragment_length = 2*$real_len + $real_ins_mean;

  #do the mapping of pair - end reads
  if ($mapper eq 'tophat1') {
  #   unless (-s "$lanepath/02_MAPPING/accepted_hits\.bam" or -s "$lanepath/02_MAPPING/accepted_hits\.unique\.sorted\.bam") {
  #     my $cmd = "tophat --output-dir $lanepath/02_MAPPING --mate-inner-dist $real_ins_mean --mate-std-dev $ins_sd --library-type fr-unstranded -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat_trans_index $bowtie_index $reads[0] $reads[1]";
  #     RunCommand($cmd,$noexecute,$quiet);
        print STDERR "tophat1 is deprecated, please use tophat2.\n";
        exit 22;
  } #tophat1

  elsif ($mapper eq 'tophat2') {
    unless (-s "$lanepath/02_MAPPING/accepted_hits\.bam" or -s "$lanepath/02_MAPPING/$mappedBam") {
      my $cmd;
      $cmd = "tophat2 --output-dir $lanepath/02_MAPPING --mate-inner-dist $real_ins_mean --mate-std-dev $ins_sd --library-type fr-unstranded -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat2_trans_index $bowtie2_index $reads[0] $reads[1]" if ($seqType =~ /^p/);
      $cmd = "tophat2 --output-dir $lanepath/02_MAPPING -p $threads --segment-length $seg_len --no-sort-bam --transcriptome-index $tophat2_trans_index $bowtie2_index $reads[0]" if ($seqType =~ /^s/);
      RunCommand($cmd,$noexecute,$quiet);
    }
  }                             #tophat2

  elsif ($mapper eq 'gsnap') {

     my $quality_options;
     if ( ($qual_zero - 33) < 10 ) {
        $quality_options = "--quality-protocol=sanger";
     } else {
        $quality_options = "--quality-zero-score=$qual_zero --quality-print-shift=$qual_move";
     }

     my $zipOption = "--gunzip";
     if ($bzip) {
       $zipOption = "--bunzip2";
     }
     unless (-s "$lanepath/02_MAPPING/accepted_hits\.bam" or -s "$lanepath/02_MAPPING/$mappedBam" or -s "$lanepath/02_MAPPING/accepted_hits\.sam") {
        my $cmd;
        $cmd = "gsnap -d $species -D $gmap_index $zipOption --format=sam --nthreads=$threads -s $gmap_splicesites --npaths=5 $quality_options $reads[0] $reads[1] >$lanepath/02_MAPPING/accepted_hits\.sam" if ($seqType =~ /^p/);
        $cmd = "gsnap -d $species -D $gmap_index $zipOption --format=sam --nthreads=$threads -s $gmap_splicesites --npaths=5 $quality_options $reads[0] >$lanepath/02_MAPPING/accepted_hits\.sam" if ($seqType =~ /^s/);
        RunCommand($cmd,$noexecute,$quiet);
     }

     if (-s "$lanepath/02_MAPPING/accepted_hits\.sam" and (! -s "$lanepath/02_MAPPING/accepted_hits\.bam" and ! -s "$lanepath/02_MAPPING/$mappedBam")) {
        my $cmd = "samtools view -Sb -@ $threads $lanepath/02_MAPPING/accepted_hits\.sam -o $lanepath/02_MAPPING/accepted_hits\.bam";
        RunCommand($cmd,$noexecute,$quiet);
     }

     if (-s "$lanepath/02_MAPPING/accepted_hits\.sam" and (-s "$lanepath/02_MAPPING/accepted_hits\.bam" or -s "$lanepath/02_MAPPING/$mappedBam")) {
        my $cmd  = "rm $lanepath/02_MAPPING/accepted_hits\.sam -f";
        RunCommand($cmd,$noexecute,$quiet);
     }

  } #gsnap

  else {
     print STDERR "Error: --mapper option should only be gsnap or tophat2. \n\n";
     exit 22;
  }


  #do the statistics
  unless (-e "$lanepath/03_STATS") {
    my $cmd = "mkdir -p $lanepath/03_STATS";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$lanepath/03_STATS/$sampleName\.mapping\.stats") {
    my $typeop = ($seqType =~ /^s/)? "--type s" : "--type p";
    my $arpop = ($seqType =~ /^p/)? "--arp $lanepath/03_STATS/$sampleName\.arp" : "";
    my $cmd = "$bin/Rseq_bam_stats --mapping $lanepath/02_MAPPING/accepted_hits\.bam $typeop --readlength $trimedlen --writer $lanepath/02_MAPPING/accepted_hits\.mapped\.bam --unmapped $lanepath/02_MAPPING/unmapped $arpop --breakpoint $lanepath/03_STATS/$sampleName\.breakpoints >$lanepath/03_STATS/$sampleName\.mapping\.stats";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$lanepath/03_STATS/$sampleName\.breakpoints\.gz" ) {
    if (-s "$lanepath/03_STATS/$sampleName\.breakpoints") {
       my $cmd = "gzip $lanepath/03_STATS/$sampleName\.breakpoints";
       RunCommand($cmd,$noexecute,$quiet);
    }
  }

  unless (-s "$lanepath/03_STATS/unmapped\.gz" ) {
    if (-s "$lanepath/03_STATS/unmapped") {
       my $cmd = "gzip $lanepath/03_STATS/unmapped";
       RunCommand($cmd,$noexecute,$quiet);
    }
  }

  my $mapping_stats_line_number = `wc -l $lanepath/03_STATS/$sampleName.mapping.stats`;
  $mapping_stats_line_number =~ s/^(\d+).*$/$1/;
  chomp($mapping_stats_line_number);
  if ($mapping_stats_line_number == 12) {
    my $total_reads = `$decompress $reads[0] | wc -l`;
    $total_reads /= 4;
    open STATS, ">>$lanepath/03_STATS/$sampleName\.mapping\.stats" || die "can not open $lanepath/03_STATS/$sampleName\.mapping\.stats\n";
    print STATS "total_frag: $total_reads\n";
    close STATS;
  }

  unless (-s "$lanepath/02_MAPPING/$mappedBam") {
    my $cmd = "samtools sort -@ $threads $lanepath/02_MAPPING/accepted_hits\.mapped\.bam $lanepath/02_MAPPING/accepted_hits\.mapped\.sorted";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if (-s "$lanepath/02_MAPPING/$mappedBam" and -s "$lanepath/02_MAPPING/accepted_hits\.mapped\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.mapped\.bam -f";
    RunCommand($cmd,$noexecute,$quiet);
  }
  if (-s "$lanepath/02_MAPPING/$mappedBam" and -s "$lanepath/02_MAPPING/accepted_hits\.bam") {
    my $cmd = "rm $lanepath/02_MAPPING/accepted_hits\.bam -f";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if ($bigWig) { #generate wiggle file
    unless (-s "$lanepath/03_STATS/$sampleName\.bw"){
      if (-s "$lanepath/03_STATS/$sampleName\.bedgraph") {
         my $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$sampleName\.bedgraph $chromosomeSize $lanepath/03_STATS/$sampleName\.bw";
         RunCommand($cmd,$noexecute,$quiet);
      }
      else {
         my $cmd = "genomeCoverageBed -ibam $lanepath/02_MAPPING/$mappedBam -bg -split -g $chromosomeSize >$lanepath/03_STATS/$sampleName\.bedgraph";
         RunCommand($cmd,$noexecute,$quiet);
         $cmd = "$bin/bedGraphToBigWig $lanepath/03_STATS/$sampleName\.bedgraph $chromosomeSize $lanepath/03_STATS/$sampleName\.bw";
         RunCommand($cmd,$noexecute,$quiet);
      }
      if (-s "$lanepath/03_STATS/$sampleName\.bedgraph" and -s "$lanepath/03_STATS/$sampleName\.bw") {
         my $cmd = "rm $lanepath/03_STATS/$sampleName\.bedgraph -f";
         RunCommand($cmd,$noexecute,$quiet);
      }
    }
  }

  #for expression extimation;
  my $readingBam = "$lanepath/02_MAPPING/$mappedBam";
  if ($merge) {
    if (-s "$lanepath/$sampleName\.merge\_list") {
      $readingBam = "$lanepath/$sampleName\.merge\_list";
    } else {
      print STDERR "$lanepath/$sampleName\.merge\_list does not exit or has no content!!!\n";
      exit 22;
    }
  }

  #ensembl gene#############################################################
  unless (-s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr\.sorted") {
    unless (-s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr") {
      my $cmd1;
      if ($seqType =~ /^p/){
        $cmd1 = "$bin/Rseq_bam_reads2expr --type p --region $ensembl_gene_bed --mapping $readingBam --posc $lanepath/03_STATS/$sampleName\.pos\.gff --chrmap $lanepath/03_STATS/$sampleName\.chrmap --lbias $lanepath/03_STATS/$sampleName\.lbias >$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr";
      } else {
        $cmd1 = "$bin/Rseq_bam_reads2expr --type s --region $ensembl_gene_bed --mapping $readingBam --posc $lanepath/03_STATS/$sampleName\.pos\.gff --chrmap $lanepath/03_STATS/$sampleName\.chrmap --lbias $lanepath/03_STATS/$sampleName\.lbias >$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr";
      }
      RunCommand($cmd1,$noexecute,$quiet);
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr >$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr.sorted";
    RunCommand($cmd2,$noexecute,$quiet);
    if (-s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr.sorted" and -s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr") {
      my $cmd3 = "rm $lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr -rf";
      RunCommand($cmd3,$noexecute,$quiet);
    }
  }

  unless (-s "$lanepath/03_STATS/$sampleName\.pos\.gff\.gz"){
    if (-s "$lanepath/03_STATS/$sampleName\.pos\.gff"){
      my $cmd = "gzip $lanepath/03_STATS/$sampleName\.pos\.gff";
      RunCommand($cmd,$noexecute,$quiet);
    }
  }

  unless (-s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.count") {
    open ENSEMBL_GENE, "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr.sorted";
    open ENSEMBL_GENE_COUNT, ">$lanepath/03_STATS/$sampleName\.ensembl\_gene\.count";
    while ( <ENSEMBL_GENE> ) {
      chomp;
      my @cols = split /\t/;
      my $current_count = round($cols[7]);
      print ENSEMBL_GENE_COUNT "$cols[3]\t$current_count\n";
    }
    close ENSEMBL_GENE;
    close ENSEMBL_GENE_COUNT;
  }

  #for RPKM normalization of ensembl genes##################################
  unless (-s "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.rpkm") {
    my $N_mapped_reads = 0;
    my $mapped = 0;
    my $singleton = 0;
    open MAPPING_STATS, "$lanepath/03_STATS/$sampleName.mapping.stats" || die "can not open $lanepath/03_STATS/$sampleName.mapping.stats";
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
    print STDERR "singletons equal to zero!!!\n" if ($singleton == 0 and $seqType =~ /^p/);
    $N_mapped_reads = 2*$mapped - $singleton if ($seqType =~ /^p/);
    $N_mapped_reads = $mapped if ($seqType =~ /^s/);
    exit if ($N_mapped_reads == 0);
    close MAPPING_STATS;

    open ENSEMBL_GENEMAP, "$ensembl_genemap";
    my %ensembl_genemap;
    while ( <ENSEMBL_GENEMAP> ) {
      chomp;
      my ($gene_id, $gene_name, $gene_type, $gene_desc, $entrez, $wiki_name, $wiki_desc) = split /\t/;
      $ensembl_genemap{$gene_id} = $gene_type."\t".$gene_name;
    }
    close ENSEMBL_GENEMAP;

    open ENSEMBL_GENE_EXPR, "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr.sorted" || die "can not open $lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr.sorted";
    open ENSEMBL_RPKM, ">$lanepath/03_STATS/$sampleName\.ensembl\_gene\.rpkm";
    printf ENSEMBL_RPKM ("%s\n", join("\t", "#ensembl_id","count", "RPKM", "gene_type", "gene_name"));
    while ( <ENSEMBL_GENE_EXPR> ) {
      chomp;
      my @cols = split /\t/;
      my $ensembl_name = $cols[3];
      my $counts_dblength = $cols[4];
      my $counts = round($cols[7]);
      my $rpkm = sprintf("%.3f", $counts_dblength * 1e9/$N_mapped_reads);
      if (exists($ensembl_genemap{$ensembl_name})){
        print ENSEMBL_RPKM "$ensembl_name\t$counts\t$rpkm\t$ensembl_genemap{$ensembl_name}\n";
      } else {
        print ENSEMBL_RPKM "$ensembl_name\t$counts\t$rpkm\t\n";
      }
    }
    close ENSEMBL_GENE_EXPR;
    close ENSEMBL_RPKM;
  }


  #gencode if it is human transcriptome
  if ($species eq 'hg19') {
    #gencode gene#############################################################
    unless (-s "$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr.sorted") {
      unless (-s "$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr") {
        my $cmd1;
        if ($seqType =~ /^p/){
          $cmd1 = "$bin/Rseq_bam_reads2expr --type p --region $gencode_gene_bed --mapping $readingBam >$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr";
        } else {
          $cmd1 = "$bin/Rseq_bam_reads2expr --type s --region $gencode_gene_bed --mapping $readingBam >$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr";
        }
        RunCommand($cmd1,$noexecute,$quiet);
      }
      my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $lanepath/03_STATS/$sampleName\.gencode\_gene\.expr >$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr.sorted";
      RunCommand($cmd2,$noexecute,$quiet);
      if (-s "$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr.sorted" and -s "$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr") {
        my $cmd3 = "rm $lanepath/03_STATS/$sampleName\.gencode\_gene\.expr -rf";
        RunCommand($cmd3,$noexecute,$quiet);
      }
    }

    #for RPKM normalization of gencode genes##################################
    unless (-s "$lanepath/03_STATS/$sampleName\.gencode\_gene\.rpkm") {
      my $N_mapped_reads = 0;
      my $mapped = 0;
      my $singleton = 0;
      open MAPPING_STATS, "$lanepath/03_STATS/$sampleName.mapping.stats" || die "can not open $lanepath/03_STATS/$sampleName.mapping.stats";
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
      print STDERR "singletons equal to zero!!!\n" if ($singleton == 0 and $seqType =~ /^p/);
      $N_mapped_reads = 2*$mapped - $singleton if ($seqType =~ /^p/);
      $N_mapped_reads = $mapped if ($seqType =~ /^s/);
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

      open GENCODE_GENE_EXPR, "$lanepath/03_STATS/$sampleName\.gencode\_gene\.expr.sorted" || die "can not open $lanepath/03_STATS/$sampleName\.gencode\_gene\.expr.sorted";
      open GENCODE_RPKM, ">$lanepath/03_STATS/$sampleName\.gencode\_gene\.rpkm";
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
  }  #end of the gencode processing


  unless (-s "$lanepath/03_STATS/$sampleName\.cate") {
    my $cmd = "perl $bin/cate.pl $lanepath/03_STATS/$sampleName\.ensembl\_gene\.expr\.sorted $gene_annotation >$lanepath/03_STATS/$sampleName\.cate";
    RunCommand($cmd,$noexecute,$quiet);
  }

  if ($seqType =~ /^p/) { #do the insert size only if it is paired-end
    unless (-s "$lanepath/03_STATS/$sampleName\.ins\.gz") {
      unless (-s "$lanepath/03_STATS/$sampleName\.ins") {
         my $cmd = "samtools view -f 0x2 $lanepath/02_MAPPING/$mappedBam | cut -f 9 | awk \'\$1\>0 \&\& \$1\<500\' >$lanepath/03_STATS/$sampleName\.ins";
         RunCommand($cmd,$noexecute,$quiet);
      }
      if (-s "$lanepath/03_STATS/$sampleName\.ins") {
         my $cmd = "gzip $lanepath/03_STATS/$sampleName\.ins";
         RunCommand($cmd,$noexecute,$quiet);
      }
    }
  } #insert size

  unless (-s "$lanepath/03_STATS/$sampleName\.report/$sampleName\.report\.html") {
    my @qc_files;
    if ($seqType =~ /^p/) {
       @qc_files = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.qc");
       @qc_files = mateorder(\@qc_files, $runID);
    } else {
       @qc_files = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.qc");
    }
    my $qcmatesuffix1;
    my $qcmatesuffix2;
    if ($seqType =~ /^p/){
      if ($qc_files[0] =~ /(\_R?[123](\_$runID)?\.fq\.qc)$/) {
         $qcmatesuffix1 = $1;
      }
      if ($qc_files[1] =~ /(\_R?[123](\_$runID)?\.fq\.qc)$/) {
         $qcmatesuffix2 = $1;
      }
    } else {
      if ($qc_files[0] =~ /((\_$runID)?\.fq\.qc)$/) {
         $qcmatesuffix1 = $1;
       }
    }
    my $cmd;
    if ($seqType =~ /^p/){
      $cmd  = "$Rbinary CMD BATCH --no-save --no-restore "."\'--args path=\"$lanepath\" lane=\"$sampleName\" anno=\"$anno\" species=\"$species\" src=\"$bin\" readlen=$real_len gf=\"$gf\" qcsuffix1=\"$qcmatesuffix1\" qcsuffix2=\"$qcmatesuffix2\" type=\"p\"' $bin/html_report.R $lanepath/03_STATS/R\_html\.out";
    } else {
      $cmd  = "$Rbinary CMD BATCH --no-save --no-restore "."\'--args path=\"$lanepath\" lane=\"$sampleName\" anno=\"$anno\" species=\"$species\" src=\"$bin\" readlen=$real_len gf=\"$gf\" qcsuffix1=\"$qcmatesuffix1\" type=\"s\"' $bin/html_report.R $lanepath/03_STATS/R\_html\.out";
    }
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

    printtime();
    print STDERR "Regional assembly process is starting \($assembly_type\)... first breakpoint processing...\n\n";

    unless (-e "$lanepath/04_ASSEMBLY") {
      my $cmd = "mkdir -p $lanepath/04_ASSEMBLY";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-s "$lanepath/03_STATS/$sampleName\.breakpoints.sorted\.gz") {
      unless (-s "$lanepath/03_STATS/$sampleName\.breakpoints\.gz" or -s "$lanepath/03_STATS/$sampleName\.breakpoints.sorted"){
        print STDERR "Error: breakpoint file does not exist, do runlevel 2 first.\n\n";
        exit 22;
      }
      unless (-s "$lanepath/03_STATS/$sampleName\.breakpoints.sorted"){
        my $cmd = "gzip -d -c $lanepath/03_STATS/$sampleName\.breakpoints\.gz | sort -k 1,1d -k 2,2n >$lanepath/03_STATS/$sampleName\.breakpoints.sorted";
        RunCommand($cmd,$noexecute,$quiet);
      }
      if ( -s "$lanepath/03_STATS/$sampleName\.breakpoints\.gz" and -s "$lanepath/03_STATS/$sampleName\.breakpoints.sorted" ){
         my $cmd2 = "rm $lanepath/03_STATS/$sampleName\.breakpoints\.gz -f";
         RunCommand($cmd2,$noexecute,$quiet);
      }
      if (-s "$lanepath/03_STATS/$sampleName\.breakpoints.sorted"){
         my $cmd3 = "gzip $lanepath/03_STATS/$sampleName\.breakpoints.sorted";
         RunCommand($cmd3,$noexecute,$quiet);
      }
    }

    #merge breakpoints for different runs!!
    my $breakpointSource = "$lanepath/03_STATS/$sampleName\.breakpoints\.sorted\.gz";
    if (scalar(@merge) != 0) {
      my $bpMerged = "";
      foreach my $mergeC (@merge) {
        my $cRunPath = "$root/$sampleName"."\_$mergeC";
        my $cRunBP = "$cRunPath/03_STATS/$sampleName\.breakpoints.sorted\.gz";
        my $cRunBP_us = "$cRunPath/03_STATS/$sampleName\.breakpoints\.gz";
        if (-s $cRunBP) {
          $bpMerged .= $cRunBP." ";
        } elsif (-s $cRunBP_us) {
          $bpMerged .= $cRunBP_us." ";
        } else {
          print STDERR "merge error: $cRunBP or $cRunBP_us does not exit or has no content!!!\n";
          exit 22;
        }
      }
      $breakpointSource = "$lanepath/03_STATS/$sampleName\.breakpoints\.merged\.sorted\.gz";
      unless (-s $breakpointSource) {
        my $cmd = "gzip -d -c $bpMerged |sort -k 1,1d -k 2,2n | gzip >$breakpointSource";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } #if merge,then redefine breakpointSource
    print STDERR "the breakpoint sources are: $breakpointSource\n";
    #breakpoint Source is defined

    unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed") {
      my $typeop = ($seqType =~ /^p/)? "p":"s";
      my $cmd = "perl $bin/breakpoint\_processing.pl $breakpointSource $typeop >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted"){
      my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted";
      RunCommand($cmd,$noexecute,$quiet);
    }

    unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot"){
      my $cmd = "perl $bin/breakpoint_repeat_masker.pl $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted $repeatMasker >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot";
      RunCommand($cmd,$noexecute,$quiet);
    }

    if ($seqType =~ /^p/) {
      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco") {
        my $mapping_bam = "$lanepath/02_MAPPING/$mappedBam";
        die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
        if ($merge) {
          if (-s "$lanepath/$sampleName\.merge\_list") {
            $mapping_bam = "$lanepath/$sampleName\.merge\_list";
          } else {
            print STDERR "$lanepath/$sampleName\.merge\_list does not exit or has no content!!!\n";
            exit 22;
          }
        }
        my $cmd = "$bin/discordant_consistency --region $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot --mapping $mapping_bam >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted"){
        my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-e "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter"){
        my $cmd = "perl $bin/breakpoint_final_prepare.pl $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted"){
        unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate"){
          my $mapping_bam = "$lanepath/02_MAPPING/$mappedBam";
          die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
          if ($merge) {
            if (-s "$lanepath/$sampleName\.merge\_list") {
              $mapping_bam = "$lanepath/$sampleName\.merge\_list";
            } else {
              print STDERR "$lanepath/$sampleName\.merge\_list does not exit or has no content!!!\n";
              exit 22;
            }
          }
          my $largestBPID = `tail -1 $lanepath/04_ASSEMBLY/$sampleName.breakpoints.processed | cut -f 1`;
          $largestBPID =~ s/\n//;
          my $cmd = "$bin/discordant_mate --mapping $mapping_bam --idstart $largestBPID --consisCount $consisCount >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate";
          RunCommand($cmd,$noexecute,$quiet);
        }
        if (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate"){
          my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted";
          RunCommand($cmd,$noexecute,$quiet);
        }
      }

      if (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate" and -s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted"){
        my $cmd = "rm $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate -f";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter") {
        my $cmd = "perl $bin/discordant_mate_processing.pl $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted $repeatMasker $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined") {
        my $cmd = "cat $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted") {
        my $cmd = "sort -k 3,3d -k 4,4n $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } #if it is paired-end reads

    else {  #it is single-end reads
      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted") {
        my $cmd = "awk \'\(\$2 == \"s\" && \$5 >= 5 && \$8 == \"N\") || \(\$2 == \"p\"\)' $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.sorted\.repornot >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted";
        RunCommand($cmd,$noexecute,$quiet);
      }
    } #single-end reads


    if ($RA == 1) { #independent regional assembly
      unless (-s "$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted\.reads") {
        my $mapping_bam = "$lanepath/02_MAPPING/$mappedBam";
        if ($merge) {
          if (-s "$lanepath/$sampleName\.merge\_list") {
            $mapping_bam = "$lanepath/$sampleName\.merge\_list";
          } else {
            print STDERR "$lanepath/$sampleName\.merge\_list does not exit or has no content!!!\n";
            exit 22;
          }
        }
        die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
        my $typeop = ($seqType =~ /^p/)? "--type p":"--type s";
        my $cmd = "$bin/reads_in_region --region $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted --mapping $mapping_bam $typeop >$lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted\.reads";
        RunCommand($cmd,$noexecute,$quiet);
      }

      my @RAssembly_reads;
      @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^p/);
      @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^s/);
      my @RAssembly_reads_gz;
      @RAssembly_reads_gz = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.RAssembly\.fq\.gz") if ($seqType =~ /^p/);
      @RAssembly_reads_gz = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.RAssembly\.fq\.gz") if ($seqType =~ /^s/);
      unless ( ($seqType =~ /^p/ and ($#RAssembly_reads == 1 or $#RAssembly_reads_gz == 1)) or ($seqType =~ /^s/ and ($#RAssembly_reads == 1 or $#RAssembly_reads_gz == 1)) ) { #get raw reads

        #need to do get raw reads here
        my @reads;
        if ($trimedlen != $readlen) {
          @reads = bsd_glob("$lanepath/01_READS/$sampleName*trimed\.fq\.$zipSuffix"); #trimmed reads
          @reads = mateorder(\@reads, $runID) if ($seqType =~ /^p/);
          if (scalar(@merge) != 0) {
            foreach my $mergeC (@merge) {
              next if ($mergeC eq $runID);
              my $cRunPath = "$root/$sampleName"."\_$mergeC";
              my @readsMerge = bsd_glob("$cRunPath/01_READS/$sampleName*trimed\.fq\.$zipSuffix");
              @readsMerge = mateorder(\@readsMerge, $mergeC) if ($seqType =~ /^p/);
              $reads[0] .= ",".$readsMerge[0];
              $reads[1] .= ",".$readsMerge[1] if ($seqType =~ /^p/);
            }
          } #if merge
        } else {
          if ($seqType =~ /^p/){ #paired-end
            @reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix"); #original reads
            @reads = mateorder(\@reads, $runID);
          } else { #single-end
            @reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix"); #original reads
          }
          if (scalar(@merge) != 0){
            foreach my $mergeC (@merge) {
              next if ($mergeC eq $runID);
              my $cRunPath = "$root/$sampleName"."\_$mergeC";
              my @readsMerge;
              if ($seqType =~ /^p/){
                @readsMerge = bsd_glob("$cRunPath/01_READS/$sampleName\_{R,}[123]{\_$mergeC,}\.fq\.$zipSuffix");
                @readsMerge = mateorder(\@readsMerge, $mergeC);
              } else {
                @readsMerge = bsd_glob("$cRunPath/01_READS/$sampleName*{\_$mergeC,}\.fq\.$zipSuffix");
              }
              $reads[0] .= ",".$readsMerge[0];
              $reads[1] .= ",".$readsMerge[1] if ($seqType =~ /^p/);
            }
          } #if merge
        }

        my $cmd;
        if ($seqType =~ /^p/){
          $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted\.reads --readfile1 $reads[0] --readfile2 $reads[1] --RA";
        } else {
          $cmd = "perl $bin/pick_ARP.pl --arpfile $lanepath/04_ASSEMBLY/$sampleName\.breakpoints\.processed\.filter\.combined\.sorted\.reads --readfile1 $reads[0] --RA";
        }
        RunCommand($cmd,$noexecute,$quiet);
      }
    } else {
      print STDERR "Error: --RA must be set to 1 or 0 (default)...\n\n";
      exit 22;
    }

    printtime();
    print STDERR "now assembling...\n\n";

    unless (-s "$lanepath/04_ASSEMBLY/transcripts.fa") {

      if ($RA == 1) {

        my @RAssembly_reads;
        @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^p/);
        @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^s/);
        my @RAssembly_reads_gz;
        @RAssembly_reads_gz = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.RAssembly\.fq\.gz") if ($seqType =~ /^p/);
        @RAssembly_reads_gz = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.RAssembly\.fq\.gz") if ($seqType =~ /^s/);

        if (($seqType =~ /^p/ and $#RAssembly_reads != 1) or ($seqType =~ /^s/ and $#RAssembly_reads != 0)) {
          if (($seqType =~ /^p/ and $#RAssembly_reads_gz == 1) or ($seqType =~ /^s/ and $#RAssembly_reads_gz == 0)) {
            foreach my $RA_read_file_gz (@RAssembly_reads_gz) {
               (my $RA_read_file = $RA_read_file_gz) =~ s/\.gz$//;
               my $cmd = "gzip -d -c $RA_read_file_gz >$RA_read_file";
               RunCommand($cmd,$noexecute,$quiet);
            }
          } else {
            print STDERR "Error: RA read files are not correctly generated...\n\n";
            exit 22;
          }
          @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^p/);
          @RAssembly_reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.RAssembly\.fq") if ($seqType =~ /^s/);
        }

        @RAssembly_reads = mateorder(\@RAssembly_reads, $runID) if ($seqType =~ /^p/);
        my $ra_reads1_file = $RAssembly_reads[0];
        my $ra_reads2_file = $RAssembly_reads[1] if ($seqType =~ /^p/);

        my %ra_reads1_jumpers;
        my %ra_reads2_jumpers;

        find_jumpers($ra_reads1_file, \%ra_reads1_jumpers);
        find_jumpers($ra_reads2_file, \%ra_reads2_jumpers) if ($seqType =~ /^p/);

        my $n_ids_ra_reads1 = scalar(keys %ra_reads1_jumpers);
        my $n_ids_ra_reads2 = scalar(keys %ra_reads2_jumpers) if ($seqType =~ /^p/);
        if ($seqType =~ /^p/ and ($n_ids_ra_reads1 != $n_ids_ra_reads2)) {
          print STDERR "Error: the numbers of breakpoints are inconsistant between read pairs ($n_ids_ra_reads1 vs $n_ids_ra_reads2).\n";
          exit 22;
        }

        #now should do regional assembly one by one
        open RA_READS1, "$ra_reads1_file";
        open RA_READS2, "$ra_reads2_file" if ($seqType =~ /^p/);
        my $speed_count = 0;
        my %printed_speed;
        foreach my $bp_id (sort {$a<=>$b} keys %ra_reads1_jumpers) {

          if ($idra != 0 and $bp_id != $idra) {
             next;
          }

          my $cmd = "mkdir -p $lanepath/04_ASSEMBLY/$bp_id"; #create dir
          RunCommand($cmd,$noexecute,1);

          regional_assembly(\*RA_READS1, $ra_reads1_jumpers{$bp_id}, $bp_id, "$lanepath/04_ASSEMBLY/$bp_id/reads_1.fq");
          regional_assembly(\*RA_READS2, $ra_reads2_jumpers{$bp_id}, $bp_id, "$lanepath/04_ASSEMBLY/$bp_id/reads_2.fq") if ($seqType =~ /^p/);

          unless (-s "$lanepath/04_ASSEMBLY/$bp_id/Roadmaps") {
            my $cmd;
            if ($seqType =~ /^p/){
              $cmd = "velveth $lanepath/04_ASSEMBLY/$bp_id/ 21 -fastq -shortPaired -separate $lanepath/04_ASSEMBLY/$bp_id/reads_1.fq $lanepath/04_ASSEMBLY/$bp_id/reads_2.fq";
            } else {
              $cmd = "velveth $lanepath/04_ASSEMBLY/$bp_id/ 21 -fastq -short $lanepath/04_ASSEMBLY/$bp_id/reads_1.fq";
            }
            RunCommand($cmd,$noexecute,1);
          }

          unless (-s "$lanepath/04_ASSEMBLY/$bp_id/Graph2") {
            my $cmd = "velvetg $lanepath/04_ASSEMBLY/$bp_id/ -read_trkg yes";
            RunCommand($cmd,$noexecute,1);
          }

          unless (-s "$lanepath/04_ASSEMBLY/$bp_id/transcripts.fa"){
            my $cmd = "oases $lanepath/04_ASSEMBLY/$bp_id/";
            RunCommand($cmd,$noexecute,1);
          }

          if (-s "$lanepath/04_ASSEMBLY/$bp_id/transcripts.fa") {
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

          if (-e "$lanepath/04_ASSEMBLY/$bp_id/") {
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
        close RA_READS2 if ($seqType =~ /^p/);

        if (($seqType =~ /^p/ and $#RAssembly_reads_gz != 1) or ($seqType =~ /^s/ and $#RAssembly_reads_gz != 0)) { #gz file does not exist
          if (($seqType =~ /^p/ and $#RAssembly_reads == 1) or ($seqType =~ /^s/ and $#RAssembly_reads == 0)) {
            foreach my $RA_read_file (@RAssembly_reads) {
              my $cmd = "gzip $RA_read_file";
              RunCommand($cmd,$noexecute,$quiet);
            }
          } else {
            print STDERR "Error: RA read files are not correctly generated... end of the regional assembly, very strange...\n\n";
            exit 22;
          }
        } else {
          foreach my $RA_read_file (@RAssembly_reads) {
             my $cmd = "rm $RA_read_file -f";
             RunCommand($cmd,$noexecute,$quiet);
          }
        } #if gz file exists

      } else {
        print STDERR "Error: --RA must be set to 1 if not 0 (default)...\n\n";
        exit 22;
      }
    }  #if transcripts.fa not exist
  }   #regional assembly test##########################################################################

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

  if ($GMAP) {  #using GMAP
      unless (-s "$lanepath/05_FUSION/$sampleName\.transcripts\.genome\.gmap") {
        my $speciesIndex = (-e "$gmap_index/$species\-all\/$species\-all.genomecomp")? $species."\-all":$species;
        my $cmd = "gmap -D $gmap_index -d $speciesIndex --format=psl -t $threads --intronlength=230000 $lanepath/04_ASSEMBLY/transcripts\.fa >$lanepath/05_FUSION/$sampleName\.transcripts\.genome\.gmap";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq") {
        my $cmd = "perl $bin/pick_fusion_transcripts_from_genomeBLAT.pl --transcripts $lanepath/04_ASSEMBLY/transcripts.fa $lanepath/05_FUSION/$sampleName\.transcripts\.genome\.gmap >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration 2>$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq";
        RunCommand($cmd,$noexecute,$quiet);
      }
  } #GMAP

  else {  #using BLAT
      unless (-s "$lanepath/05_FUSION/$sampleName\.transcripts\.genome\.blat") {
        my $cmd = "blat -maxIntron=230000 $blatDatabase $lanepath/04\_ASSEMBLY/transcripts\.fa $lanepath/05_FUSION/$sampleName\.transcripts\.genome\.blat";
        RunCommand($cmd,$noexecute,$quiet);
      }

      unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq") {
        my $cmd = "perl $bin/pick_fusion_transcripts_from_genomeBLAT.pl --transcripts $lanepath/04_ASSEMBLY/transcripts.fa $lanepath/05_FUSION/$sampleName\.transcripts\.genome\.blat >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration 2>$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq";
        RunCommand($cmd,$noexecute,$quiet);
      }
  } #BLAT

  #build fusion candidate index (now bowtie2)
  unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq\.index\.1\.bt2") {
    my $cmd = "bowtie2-build --quiet $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq\.index";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #map reads to the fusion candidate index
  unless ( -s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.bam" ) {

    if ( -s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam" ) { #bowtie2 is done

      my $cmd = "samtools view -Sb $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam -o $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2.bam";
      RunCommand($cmd,$noexecute,$quiet);
      if ( -s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2.bam" ) {
         my $cmd = "rm $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam -f";
         RunCommand($cmd,$noexecute,$quiet);
      }

    } else {  #bowtie2 is not done

      my @reads;
      if ($trimedlen != $readlen) {
          @reads = bsd_glob("$lanepath/01_READS/$sampleName*trimed\.fq\.$zipSuffix"); #trimmed reads
          @reads = mateorder(\@reads, $runID) if ($seqType =~ /^p/);
          if (scalar(@merge) != 0) {
            foreach my $mergeC (@merge) {
              next if ($mergeC eq $runID);
              my $cRunPath = "$root/$sampleName"."\_$mergeC";
              my @readsMerge = bsd_glob("$cRunPath/01_READS/$sampleName*trimed\.fq\.$zipSuffix");
              @readsMerge = mateorder(\@readsMerge, $mergeC) if ($seqType =~ /^p/);
              $reads[0] .= ",".$readsMerge[0];
              $reads[1] .= ",".$readsMerge[1] if ($seqType =~ /^p/);
            }
          } #if merge
       } else {
         if ($seqType =~ /^p/){
           @reads = bsd_glob("$lanepath/01_READS/$sampleName\_{R,}[123]{\_$runID,}\.fq\.$zipSuffix"); #original reads
           @reads = mateorder(\@reads, $runID);
         } else {
           @reads = bsd_glob("$lanepath/01_READS/$sampleName*{\_$runID,}\.fq\.$zipSuffix");
         }
         if (scalar(@merge) != 0) {
           foreach my $mergeC (@merge) {
             next if ($mergeC eq $runID);
             my $cRunPath = "$root/$sampleName"."\_$mergeC";
             my @readsMerge;
             if ($seqType =~ /^p/) {
               @readsMerge= bsd_glob("$cRunPath/01_READS/$sampleName\_{R,}[123]{\_$mergeC,}\.fq\.$zipSuffix");
               @readsMerge = mateorder(\@readsMerge, $mergeC);
             } else {
               @readsMerge= bsd_glob("$cRunPath/01_READS/$sampleName*{\_$mergeC,}\.fq\.$zipSuffix");
             }
             $reads[0] .= ",".$readsMerge[0];
             $reads[1] .= ",".$readsMerge[1] if ($seqType =~ /^p/);
           }
         } #if merge
       } #original read length

      my $readsop = ($seqType =~ /^p/)? "-1 $reads[0] -2 $reads[1]":"-U $reads[0]";
      my $cmd = "bowtie2 -k 22 -p $threads --no-unal --score-min L,-2,-0.15 -x $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq\.index $readsop >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam";
      RunCommand($cmd,$noexecute,$quiet);

      $cmd = "samtools view -Sb $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam -o $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2.bam";
      RunCommand($cmd,$noexecute,$quiet);

      if (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2.bam"){
         my $cmd = "rm $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.sam -f";
         RunCommand($cmd,$noexecute,$quiet);
      }
    } #bowtie2 is not done
  } #bowtie2 bam is generated

  #get fusion coverage
  unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov") {
    my $gfc_opts = '';
    $gfc_opts = "--genomeBlatPred $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration";

    my $readingBam = "$lanepath/02_MAPPING/$mappedBam";
    if ($customMappedBam ne ''){
       $readingBam = $customMappedBam;
    }
    if ($merge) {
      foreach my $mergeC (@merge) {
         next if ($mergeC eq $runID);
         my $cRunPath = "$root/$sampleName"."\_$mergeC";
         my $cRunBam = "$cRunPath/02_MAPPING/$mappedBam";
         if (-s $cRunBam){
           $readingBam .= ",".$cRunBam;
         } else {
           print STDERR "merge error: $cRunBam does not exist!!!!";
           exit 22;
         }
      }
    } #if merge

    my $cmd;
    if ($seqType =~ /^p/) { #paired-end
      $cmd = "perl $bin/get_fusion_coverage.pl $gfc_opts --type pair --mappingfile $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.bam --readlength $trimedlen --geneanno $gene_annotation --repeatmasker $repeatMasker --selfChain $selfChain --accepthits $readingBam --encomcov $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov.enco >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    } else { #single-end
      $cmd = "perl $bin/get_fusion_coverage.pl $gfc_opts --type single --mappingfile $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie2\.bam --readlength $trimedlen --geneanno $gene_annotation --repeatmasker $repeatMasker --selfChain $selfChain --accepthits $readingBam >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    }
    RunCommand($cmd,$noexecute,$quiet);
  }

  #visualize coverage
  unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis") {
    my $fpv_opts = '';
    $fpv_opts = "--genomeBlatPred";
    my $typeop = ($seqType =~ /^p/)? "--seqType p":"--seqType s";
    my $cmd = "perl $bin/further_processing_for_read_visualization.pl $fpv_opts $typeop --fusionseqfile $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.seq --coveragefile $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov --readlength $trimedlen >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.list"){
    my $sorting_opts = '';
    $sorting_opts = "-k 16,16nr -k 18,18nr -k 3,3nr -k 8,8d -k 11,11d";
    my $cmd = "grep \"^\#\" $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis \| sort $sorting_opts >$lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.list";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$lanepath/05_FUSION/$sampleName\.fusion\.report") {
    my $cmd = "perl $bin/fusion_report.pl $lanepath/05_FUSION/$sampleName\.fusion_transcirpts_after_filtration\.list >$lanepath/05_FUSION/$sampleName\.fusion\.report";
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

  unless (-s "$lanepath/06_CUFFLINKS/transcripts\.gtf") {

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

    my $mapping_bam = "$lanepath/02_MAPPING/$mappedBam";
    if (-s $mapping_bam){
      my $cmd = "cufflinks -o $lanepath/06_CUFFLINKS -p $threads $cufflinks_options --quiet $mapping_bam";
      RunCommand($cmd,$noexecute,$quiet);
    }
    else {
      print STDERR "$mapping_bam does not exist, please do the mapping first.\n";
      exit;
    }
  }

  if (-s "$lanepath/06_CUFFLINKS/transcripts\.gtf") { #collect the gtf list
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

  unless (-s "$lanepath/06_CUFFLINKS/$sampleName\.transcripts\.gtf\.tmap") {
    my $cmd = "cuffcompare -o $lanepath/06_CUFFLINKS/$sampleName -r $refseq_gene_gtf $lanepath/06_CUFFLINKS/transcripts\.gtf";
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

  unless (-s "$gtf_list") {
    print STDERR "Error: $gtf_list file does not exist!!! Please rerun RTrace.pl for each sample for runlevel 5 (cufflinks).\n\n";
    exit 22;
  }

  #cuffmerge####################################
  unless (-s $cuffmerge_output) {
    my $cmd = "mkdir -p $cuffmerge_output";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$cuffmerge_output/merged\.gtf") {
    my $cmd = "cuffmerge -o $cuffmerge_output --ref-gtf $gene_annotation_gtf --num-threads $threads --ref-sequence $genome_fasta $gtf_list";
    RunCommand($cmd,$noexecute,$quiet);
  }

  #cuffdiff#####################################
  unless (-e $cuffdiff_output) {
    my $cmd = "mkdir -p $cuffdiff_output";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$cuffdiff_output/splicing\.diff") {
    my %bam_files_cd;
    open TARGETS, "$root/edgeR/targets" || die "Error: could not open $root/edgeR/targets for the information of bam files for cuffdiff.\n please rerun the pipeline for each samples with --patient and --tissue arguments set.\n";
    while ( <TARGETS> ) {
       chomp;
       next if /^label/;
       my ($cu_sample, $cu_expr_count, $tissue, $patient) = split /\t/;
       my $cu_bam_file = "$root/$cu_sample/02_MAPPING/$mappedBam";
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
if ($patient and $tissue) {
  my $count_file = "$lanepath/03_STATS/$sampleName\.ensembl\_gene\.count";
  unless (-e "$root/edgeR") {
    my $cmd = "mkdir -p $root/edgeR";
    RunCommand($cmd,$noexecute,$quiet);
  }
  if (! -s "$root/edgeR/targets") {
    open TARGETS_OUT, ">>$root/edgeR/targets";
    printf TARGETS_OUT "%s\n", join("\t", 'label', 'files', 'tissue', 'patient');
    printf TARGETS_OUT "%s\n", join("\t", $sampleName, $count_file, $tissue, $patient);
    close TARGETS_OUT;
  } else {
    my $writeornot = 1;
    open TARGETS_IN, "$root/edgeR/targets";
    while ( <TARGETS_IN> ){
       chomp;
       my @cols = split /\t/;
       $writeornot = 0 if ($cols[0] eq $sampleName);
    }
    close TARGETS_IN;
    if ($writeornot == 1){
      open TARGETS_OUT, ">>$root/edgeR/targets";
      printf TARGETS_OUT "%s\n", join("\t", $sampleName, $count_file, $tissue, $patient);
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

  print STDERR "DE analysis parameters:\n";
  print STDERR "priordf: $priordf\; spaired: $spaired\; $pairDE1\-$pairDE2\n";

  if (! -e $edgeR_output) {
     print STDERR "Error: $edgeR_output dir does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  } elsif (! -e "$edgeR_output/targets") {
     print STDERR "Error: $edgeR_output/targets file does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  }

  unless (-s "$edgeR_output/topDE.txt") {
    my $cmd = "$Rbinary CMD BATCH --no-save --no-restore "."\'--args path=\"$edgeR_output\" priordf=$priordf spaired=$spaired pair1=\"$pairDE1\" pair2=\"$pairDE2\"\' $bin/edgeR_test.R $edgeR_output/R\_html\.out";
    RunCommand($cmd,$noexecute,$quiet);
  }

  unless (-s "$edgeR_output/DE\_table\.txt") {
    my $cmd = "perl $bin/get_DEtable.pl $edgeR_output topDE\.txt ifDE\.txt $spaired >$edgeR_output/DE\_table\.txt";
    RunCommand($cmd,$noexecute,$quiet);
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###------------###################################################################################################
##  sub-region  #################################################################################################
###------------###################################################################################################

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
  print STDERR "\nGENERAL OPTIONS (MUST SET):\n\t--runlevel\tthe steps of runlevel, from 1-7, either rl1-rl2 or rl. See below for options for each runlevel.\n";
  print STDERR "\t--sampleName\tthe name of the lane needed to be processed (must set for runlevel 1-5)\n";
  print STDERR "\t--seqType\tset to 's' if it is a single-end sequencing experiment, or 'p' for paired-end (default).\n";
  print STDERR "\t--runID\t\tthe ID of the run needed to be processed (default not set, must set if fastq files are ended with _R1_00X.fq)\n";
  print STDERR "\t--merge\t\ta comma-separated list of runIDs that are needed to be merged, the output will be written to the defined runID\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$bin/../PIPELINE/, MUST set using other dir)\n";
  print STDERR "\t--anno\t\tthe annotations directory (default is \$bin/../ANNOTATION/, MUST set using other dir)\n";
  print STDERR "\t--species\tspecify the reference version of the species, such as hg19 (default), mm10.\n";
  print STDERR "\t--patient\tthe patient id, which will be written into the target file for edgeR\n";
  print STDERR "\t--tissue\tthe tissue type name (like \'normal\', \'cancer\'), for the target file for running edgeR and cuffdiff\n";
  print STDERR "\t--Rbinary\tthe name of R executable, default is \'R\'. Set if your R binary name is different.\n\n";

  print STDERR "CONTROL OPTIONS FOR EACH RUNLEVEL:\n";
  print STDERR "runlevel 1: quality checking and insert size estimatiion using part of reads\n";
  print STDERR "\t--readpool\tthe directory where all the read files with names ending with \.fq\.gz or \.fastq\.gz located.\n";
  print STDERR "\t--fqreid\trename the fastq id in case of \/1N, only for gsnap mapping.\n";
  print STDERR "\t--QC\t\tdo the quality check of reads, will stop the pipeline once it is finished.\n";
  print STDERR "\t--qcOFF\t\tturn off quality check.\n";
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
  print STDERR "\t--RA\t\tuse regional assembly for runlevel 3. Default is 1. When using tophat as mapper in runlevel 2, set to 0.\n";
  print STDERR "\t--consisCount\tnumber of consistent read pairs with discordant mapping (default: 5). use smaller value For <70bp reads or low depth data.\n";

  print STDERR "\nrunlevel 4: detection of fusion candidates\n";
  print STDERR "\t--GMAP\t\tset if use GMAP in run-level 4 (default: use BLAT).\n";

  print STDERR "\nrunlevel 5: run cufflinks for gene/isoform quantification\n";
  print STDERR "runlevel 6: run cuffdiff for diffrential gene/isoform expression analysis\n";
  print STDERR "\t--gtf-guide\tuse gtf guided assembly method in run-level 5 (default: FALSE).\n";
  print STDERR "\t--known-trans\twhich known transcript annotation to be used for cufflinks, either \'ensembl\' (default) or 'refseq'.\n";
  print STDERR "\t--frag-bias\tcorrect fragmentation bias in run-level 5 (default: FALSE).\n";
  print STDERR "\t--upper-qt\tupper-quantile normalization in run-level 5 (default: FALSE).\n";

  print STDERR "\nrunlevel 7: run edgeR for diffrential gene expression analysis\n";
  print STDERR "\t--priordf\tthe prior.df parameter in edgeR, which determines the amount of smoothing of tagwise dispersions towards the common dispersion.\n\t\t\tThe larger the value for prior.df, the more smoothing. A prior.df of 1 gives the common likelihood the weight of one observation. \n\t\t\tDefault is 10. Set it smaller for large sample size (i.e., set to 1 for more than 20 replicates).\n";
  print STDERR "\t--spaired\twhether the experiment is a paired normal-disease design, 1 means yes (default), 0 for no.\n";
  print STDERR "\t--pairDE1\tspecify the name of the groups to be compared, e.g., Normal. (pairDE2 \-\> pairedDE1)\n";
  print STDERR "\t--pairDE2\tspecify the name of the groups to be compared, e.g., Tumour. (pairDE2 \-\> pairedDE1)\n";

  print STDERR "\nOTHER OPTIONS\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--quiet\t\tdo not print the command line calls and time information\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n";
  print STDERR "\t--bzip\t\tthe read fastq files are bziped, rather than gziped (default).\n";

  print STDERR "\nSynopsis: RTrace.pl --runlevel 1 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 2 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --patient <ID> --tissue <type> --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 3-4 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --RA 1 --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 5 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --threads <N> 2>>run.log\n";
  print STDERR "Synopsis: RTrace.pl --runlevel 7 --root <dir_root> --anno <dir_anno> --priordf 1 2>>run.log\n";
  print STDERR "remember to set --species option to a string, such as mm10, if it is not a human sample!!!\n\n";

  print STDERR "Runlevel dependencies (->): 4->3->2->1, 6->5->2->1, 7->2->1\n\n";
  exit 0;
}

sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub mateorder {
  my ($r, $laneid) = @_;
  my @tmp;

  if (scalar(@{$r}) != 2) {
    print STDERR "there are not two mate read files, exit.\n";
    exit 1;
  }

  foreach my $m (@{$r}){
    if ($m =~ /_R?1(\_$laneid)?\./) {
      $tmp[0] = $m;
    } elsif ($m =~ /_R?[23](\_$laneid)?\./) {
      $tmp[1] = $m;
    }
    else {
      print STDERR "the mate read file dosen't contain _1 _2 or _3 information, exit.\n";
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
