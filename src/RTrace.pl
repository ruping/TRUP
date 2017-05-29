#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);

my %options;
my %runlevel;
my %runTask;


$options{'configure'}  = "SRP";
$options{'noexecute'}  = 0;
$options{'quiet'}      = 0;
$options{'runlevels'}  = 0;
$options{'readlen'}    = 0;
$options{'mapper'}     = "gsnap";
$options{'seg_len'}    = 25;
$options{'ins_mean'}   = 0;
$options{'ins_mean_true'} = 0;
$options{'ins_sd'}     = 20;
$options{'sampleName'} = "SRP";
$options{'FASTQ1'}      = 'SRP';
$options{'FASTQ2'}      = 'SRP';
$options{'fastqFiles1'} = 'SRP';
$options{'fastqFiles2'} = 'SRP';
$options{'platform'}    = "ILLUMINA";
$options{'lanepath'}   = 'SRP';
$options{'threads'}    = 1;
$options{'help'}       = undef;
$options{'QC'}         = undef;          #quality check
$options{'qcOFF'}      = undef;
$options{'GMAP'}       = undef;          #using GMAP instead of BLAT for runlevel-4
$options{'RA'}         = 1;              #regional assembly
$options{'idra'}       = 0;              #for a specific breakpoint id
$options{'consisCount'} = 5;             #the threshold for the number of consistent mate pairs with discordant maping
$options{'force'}      = undef;          #force
$options{'bigWig'}     = undef;          #wiggle file
$options{'strandWig'}  = undef;          #splitstrand for wiggle
$options{'gtf_guide_assembly'} = undef;  #for cufflinks
$options{'known_trans'} = 'ensembl';     #for cufflinks
$options{'frag_bias_correct'} = undef;   #for cufflinks
$options{'upper_quantile_norm'} = undef; #for cufflinks
$options{'root'}        = "$RealBin/../PIPELINE";
$options{'species'}     = 'hg19';
$options{'readpool'}    = 'SRP';
$options{'bin'}         = "$RealBin/";
$options{'annovarbin'}  = "$options{'bin'}/../annovar/";
$options{'qual_zero'}   = 33;
$options{'qual_move'}   = 0;
$options{'fq_reid'}     = undef;         #rename fastq read id (for gsnap)
$options{'priordf'}     = 10;            #for edgeR
$options{'pairDE1'}     = 'N';           #for edgeR
$options{'pairDE2'}     = 'T';           #for edgeR
$options{'spaired'}     = 1;             #for edgeR
$options{'patient'}     = undef;         #the patient id for edgeR DE test
$options{'tissue'}      = undef;         #the tissue type for edgeR DE test
$options{'gf'}          = "png";         #the format used in html report
$options{'bzip'}        = undef;         #to allow bzip compressed fastq files
$options{'Rbinary'}     = 'R';
$options{'customMappedBam'} = '';
$options{'seqType'}     = 'p';           #whether paired-end or single-end
$options{'tmpDir'}      = '';
$options{'misPen'}      = 2;
$options{'uniqueBase'}  = 100;
$options{'maxIntron'}   = 230000;

if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}

GetOptions(
           "configure=s"  => \$options{'configure'},
           "sampleName=s" => \$options{'sampleName'},
           "FASTQ1=s"     => \$options{'FASTQ1'},
           "FASTQ2=s"     => \$options{'FASTQ2'},
           "fastqFiles1=s"=> \$options{'fastqFiles1'},
           "fastqFiles2=s"=> \$options{'fastqFiles2'},
           "platform=s"   => \$options{'platform'},
           "runlevel=s"   => \$options{'runlevels'},
           "seqType=s"    => \$options{'seqType'},
           "noexecute"    => \$options{'noexecute'},
           "quiet"        => \$options{'quiet'},
           "readlen=i"    => \$options{'readlen'},
           "seglen=i"     => \$options{'seg_len'},
           "mapper=s"     => \$options{'mapper'},
           "insertmean=i" => \$options{'ins_mean'},
           "insertsd=i"   => \$options{'ins_sd'},
           "threads=i"    => \$options{'threads'},
           "consisCount=i"=> \$options{'consisCount'},
           "QC"           => \$options{'QC'},
           "qcOFF"        => \$options{'qcOFF'},
           "GMAP"         => \$options{'GMAP'},
           "RA=i"         => \$options{'RA'},
           "gf=s"         => \$options{'gf'},
           "idra=i"       => \$options{'idra'},
           "WIG"          => \$options{'bigWig'},
           "strandWig"    => \$options{'strandWig'},
           "fqreid"       => \$options{'fq_reid'},
           "gtf-guide"    => \$options{'gtf_guide_assembly'},
           "known-trans"  => \$options{'known_trans'},
           "frag-bias"    => \$options{'frag_bias_correct'},
           "upper-qt"     => \$options{'upper_quantile_norm'},
           "force"        => \$options{'force'},
           "root=s"       => \$options{'root'},
           "species=s"    => \$options{'species'},
           "readpool=s"   => \$options{'readpool'},
           "priordf=i"    => \$options{'priordf'},
           "spaired=i"    => \$options{'spaired'},
           "pairDE1=s"    => \$options{'pairDE1'},
           "pairDE2=s"    => \$options{'pairDE2'},
           "patient=s"    => \$options{'patient'},
           "tissue=s"     => \$options{'tissue'},
           "bzip"         => \$options{'bzip'},
           "Rbinary=s"    => \$options{'Rbinary'},
           "help|h"       => \$options{'help'},
           "customBam=s"  => \$options{'customMappedBam'},
           "tmpDir=s"     => \$options{'tmpDir'},
           "misPen=f"     => \$options{'misPen'},
           "uniqueBase=i" => \$options{'uniqueBase'},
           "maxIntron=i"  => \$options{'maxIntron'},
          );

#print help
helpm() if ($options{'help'});



### Read configuration and set all paths----------------------------------
my %confs;
open IN, "$options{'configure'}";
while ( <IN> ) {
  chomp;
  next if /^#/;
  my @cols = split /\t/;
  $confs{$cols[0]} = $cols[1];
}
close IN;

#translate environment variable
foreach my $confele (keys %confs){
  while ($confs{$confele} =~ /\$([A-Za-z0-9]+)/g) {
    my $eleName = $1;
    my $eleTranslate;
    if (exists ($confs{$eleName})) {
      $eleTranslate = $confs{$eleName};
      $confs{$confele} =~ s/\$$eleName/$eleTranslate/;
    } else {
      die("can't translate eleName: $eleName\n");
    }
  }
}
print STDERR Dumper (\%confs);
#-------------------------------------------------------------------------


### Frequently used names-------------------------------------------------
my $mappedBam = "accepted_hits\.mapped\.sorted\.bam";
#-------------------------------------------------------------------------

#decompression option-----------------------------------------------------
$options{'decompress'} = "gzip -d -c";
$options{'compress'}   = "gzip";
$options{'zipSuffix'}  = "gz";
if ( $options{'bzip'} ) {
  $options{'decompress'} = "bzip2 -d -c";
  $options{'compress'}   = "bzip2";
  $options{'zipSuffix'}  = "bz2";
}
#-------------------------------------------------------------------------

### Already specified full path fastq files-------------------------------
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles1'} =~ s/\,/ /g;
}
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles2'} =~ s/\,/ /g;
}
#-------------------------------------------------------------------------

### Runlevel/Task check up------------------------------------------------
if ($options{'runlevels'}) { #true runlevels
  foreach my $r (split /\,/,$options{'runlevels'}) {
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
if ($options{'runTask'}) {
  foreach my $task (split(/\,/, $options{'runTask'})) {
    $runTask{$task} = '';
  }
}
if (! $options{'runlevels'} and ! $options{'runTask'}) {
  print STDERR "no runlevel or runTask has been set, exit.\n";
  helpm();
}
#-------------------------------------------------------------------------


if ($options{'root'} eq "$RealBin/../PIPELINE") {
  if (-e "$RealBin/../PIPELINE") {
    print STDERR "no root dir given, analysis will be run under $options{'root'}.\n";
  }
  else {
    print STDERR "no root dir given, $options{'root'} does not exist, please do -h or --help to check how to set root dir.\n";
    helpm();
  }
} else {
  $options{'readpool'} = $options{'root'} if $options{'readpool'} eq 'SRP';
}


###
###preparation the lane and read path enviroment
###

if ($options{'lanepath'} eq 'SRP' and $options{'sampleName'} ne 'SRP') {
  printtime();
  $options{'lanepath'} = "$options{'root'}/$options{'sampleName'}";   #define lane path
  print STDERR "####### lane name is set to $options{'sampleName'} #######\n\n";
  unless (-e "$options{'lanepath'}") {
    my $cmd = "mkdir -p $options{'lanepath'}";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
}

if ($options{'readpool'} ne 'SRP' and $options{'FASTQ1'} ne 'SRP' and $options{'fastqFiles1'} eq 'SRP') {

  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$options{'lanepath'}/01_READS/") {
    my $cmd = "mkdir -p $options{'lanepath'}/01_READS/";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my @fastqFile1 = split(/\,/, $options{'FASTQ1'});
  $options{'fastqFiles1'} =~ s/^SRP//;
  foreach my $fastqFile1 (@fastqFile1) {
    if ($fastqFile1 !~ /\.[bg]z2?$/){
      die "\[error\]: $fastqFile1 must be gzip or bzipped!\n";
    }
    $fastqFile1 = $options{'readpool'}.'/'.$fastqFile1;
    $options{'fastqFiles1'} .= $fastqFile1." ";
  }
  $options{'fastqFiles1'} =~ s/\s$//;

  if ($options{'FASTQ2'} ne 'SRP') {
    my @fastqFile2 = split(/\,/, $options{'FASTQ2'});
    $options{'fastqFiles2'} =~ s/^SRP//;
    foreach my $fastqFile2 (@fastqFile2) {
      $fastqFile2 = $options{'readpool'}.'/'.$fastqFile2;
      $options{'fastqFiles2'} .= $fastqFile2." ";
    }
    $options{'fastqFiles2'} =~ s/\s$//;
  }

  print STDERR "lanefile1:\t$options{'fastqFiles1'}\n";
  print STDERR "lanefile2:\t$options{'fastqFiles2'}\n"; #if paired end

}

if ($options{'fastqFiles1'} ne 'SRP') {
  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$options{'lanepath'}/01_READS/") {
    my $cmd = "mkdir -p $options{'lanepath'}/01_READS/";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  foreach my $fastqFile1 (split(" ", $options{'fastqFiles1'})) {
    my $cmd = "ln -b -s $fastqFile1 $options{'lanepath'}/01_READS/";
    my $fastqFile1Basename = basename($fastqFile1);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile1Basename");
  }
}

if ($options{'fastqFiles2'} ne 'SRP' and $options{'fastqFiles2'} ne 'interleaved') {
  foreach my $fastqFile2 (split(" ", $options{'fastqFiles2'})){
    my $cmd = "ln -b -s $fastqFile2 $options{'lanepath'}/01_READS/";
    my $fastqFile2Basename = basename($fastqFile2);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile2Basename");
  }
}


REALSTEPS:
############################ determine read length if not provided ###################################
if ( $options{'readlen'} == 0 and $options{'sampleName'} ne 'SRP') { #read length not set
  if ($options{'fastqFiles1'} ne 'SRP') {
    my @fastqFiles1Temp = split(/\s/, $options{'fastqFiles1'});
    if ( -s "$fastqFiles1Temp[0]" ) {
      my $first_second_line = `$options{'decompress'} "$fastqFiles1Temp[0]" | head -2 | grep -v "^@"`;
      $options{'readlen'} = length($first_second_line) - 1;
    }
  }
  if ( $options{'readlen'} == 0 ) {  #still zero, check bams
    my @bamTmp = bsd_glob("$options{'lanepath'}/02_MAPPING/*.bam");
    foreach my $bamTmp (@bamTmp){
      if (-s "$bamTmp") {
        my $samSix = `samtools view $bamTmp \| awk \-F\"\t\" \'\{print \$6\}\' \| awk \'\$1 \!\~ \/\[IDNHS\\\*\]\/\' \| head \-1000 \| tr \"\\n\" \"\,\"`;
        chomp($samSix);
        my @matchLen;
        foreach my $matchLen (split(/\,/, $samSix)){
          $matchLen =~ /^(\d+)M$/;
          $matchLen = $1;
          push(@matchLen, $matchLen);
        }
        $options{'readlen'} = median(\@matchLen);   #set length to the median of the first 1000 reads with complete matching
      }
      if ($options{'readlen'} > 0){
        last;
      }
    } #loop each bams to find existing bam
  } #check bam
  print STDERR "read length is not set, will take the original read length ($options{'readlen'} bp)\n";
}
#######################################################################################################



###
###runlevel1: trim the reads if necessary and insert size detection using spiked in reads
###

my $runlevels = 1;
if (exists($runlevel{$runlevels}) or exists($runTask{'QC'})) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my @qc_files;
  my $qc_out1 = "$options{'lanepath'}/01_READS/mate1.qc";    ######
  my $qc_out2;

  if ($options{'fastqFiles2'} ne 'SRP' and $options{'fastqFiles2'} ne 'interleaved') {
    $qc_out2 = "$options{'lanepath'}/01_READS/mate2.qc";
  }

  unless ((-e "$qc_out1")) {
    my $cmd = "$options{'decompress'} $options{'fastqFiles1'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }
  unless (($options{'fastqFiles2'} eq 'SRP' or $options{'fastqFiles2'} eq 'interleaved') or ($qc_out2 ne '' and -e "$qc_out2")) {
    my $cmd = "$options{'decompress'} $options{'fastqFiles2'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }

  #my $total_reads = `$decompress $read_files[0] | wc -l`;
  #$total_reads /= 4;
  #print STDERR "total number of read pairs: $total_reads\n";


  #if ( $options{'fq_reid'} ) {  #renbame fastq id
  #  foreach my $read_file (@read_files) {
  #     my $reid_file = $read_file;
  #     $reid_file =~ s/\.fq\.$options{'zipSuffix'}/\.reid\.fq\.$options{'zipSuffix'}/;
  #     unless (-e $reid_file) {
  #       my $cmd = "perl $options{'bin'}/fqreid.pl $read_file | $compress >$reid_file";
  #       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #     }
  #     if (-e $read_file and (-s $reid_file) > 1000) {
  #       my $cmd = "mv -f $reid_file $read_file";
  #       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #     } elsif (-e $reid_file and (-s $reid_file) <= 1000) {
  #       my $cmd = "rm $reid_file -f";              #the original fastq file looks good
  #       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #     }
  #  } #foreach read file
  #}

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}

###
###runlevel1.5: deciding fragment/insert size
###

if (defined $options{'sampleName'}) {

  printtime();
  print STDERR "####### insert mean and sd calculation #######\n\n";

  my @quality_check_files;
  @quality_check_files = bsd_glob("$options{'lanepath'}/01_READS/*.qc");
  if ( scalar(@quality_check_files) == 2 or (scalar(@quality_check_files) == 1 and $options{'seqType'} =~ /single-end/)) { #decide the quality shift
    open QC, "<$quality_check_files[0]";
    my $qual_min = -1;
    my $qual_max = -1;
    while ( <QC> ) {
      chomp;
      next if ($_ !~ /^\d+/);
      my @cols = split /\t/;
      $qual_min = $cols[2] if ($cols[2] < $qual_min or $qual_min = -1);
      $qual_max = $cols[3] if ($cols[3] > $qual_max or $qual_max = -1);
      print STDERR "$cols[2]\t$cols[3]\t$qual_min\t$qual_max\n";
    }
    close QC;
    print STDERR "qual_min: $qual_min; qual_max: $qual_max; ";
    $options{'qual_zero'} = $qual_min+33;
    $options{'qual_move'} = -$qual_min;
    print STDERR "qual_zero: $options{'qual_zero'}; qual_shift: $options{'qual_move'}.\n";
  } else {
    unless ( $options{'force'} ) {
      print STDERR "please do quality check first using option --QC.\n";
      exit;
    }
  }

  #if ($ins_mean == 0 and $options{'seqType'} =~ /paired-end/) {
  #  my $fraginfo = `$options{'Rbinary'} --no-save --slave \'--args path=\"$options{'lanepath'}/00_TEST/\" lane=\"$options{'sampleName'}\.spikedin\"\' < $options{'bin'}/fragment_length.R`;
  #  $fraginfo =~ /^(.+)\s(.+)/;
  #  my $frag_mean = $1; $frag_mean = round($frag_mean);
  #  my $insert_sd   = $2; $insert_sd =~ s/\n//; $insert_sd = round($insert_sd);
  #  my $real_len = $trimedlen;
  #  $ins_sd = $insert_sd;
  #  $ins_mean = $frag_mean - 2*$real_len;
  #  $ins_mean_true = $ins_mean;
  #  print STDERR "insert mean: $ins_mean\tinsert_sd: $ins_sd\n";
  #}
}


###
###runlevel2: do the mapping and generate the statistics
###

$runlevels = 2;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$options{'lanepath'}/02_MAPPING") {
    my $cmd = "mkdir -p $options{'lanepath'}/02_MAPPING";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #do the mapping of pair - end reads
  #if ($mapper eq 'star') {
  #  unless (-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam" or -s "$options{'lanepath'}/02_MAPPING/$mappedBam") {
  #    my $cmd;
  #    unless (-s "$options{'lanepath'}/02_MAPPING/starChimeric.out.bam") {
  #       my $zipcat = "zcat";
  #       $zipcat = "bzcat" if ($bzip);
  #       $cmd = "STAR --genomeDir $star_index --readFilesCommand $zipcat --readFilesIn $reads[0] $reads[1] --runThreadN $options{'threads'} --outFileNamePrefix $options{'lanepath'}/02_MAPPING/star --chimSegmentMin 18" if ($options{'seqType'} =~ /paired-end/);
  #       $cmd = "STAR --genomeDir $star_index --readFilesCommand $zipcat --readFilesIn $reads[0] --runThreadN $options{'threads'} --outFileNamePrefix $options{'lanepath'}/02_MAPPING/star --chimSegmentMin 18" if ($options{'seqType'} =~ /single-end/);
  #       unless (-s "$options{'lanepath'}/02_MAPPING/starChimeric.out.sam") {
  #         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #       }
  #       if (-s "$options{'lanepath'}/02_MAPPING/starAligned.out.sam" and !-s "$options{'lanepath'}/02_MAPPING/starAligned.out.bam"){
  #         $cmd = "samtools view -Sb $options{'lanepath'}/02_MAPPING/starAligned.out.sam -o $options{'lanepath'}/02_MAPPING/starAligned.out.bam";
  #         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #         $cmd = "rm $options{'lanepath'}/02_MAPPING/starAligned.out.sam -f";
  #         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #       }
  #       if (-s "$options{'lanepath'}/02_MAPPING/starChimeric.out.sam" and !-s "$options{'lanepath'}/02_MAPPING/starChimeric.out.bam"){
  #         $cmd = "samtools view -Sb $options{'lanepath'}/02_MAPPING/starChimeric.out.sam -o $options{'lanepath'}/02_MAPPING/starChimeric.out.bam";
  #         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #         $cmd = "rm $options{'lanepath'}/02_MAPPING/starChimeric.out.sam -f";
  #         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #       }
  #    }
  #    if (-s "$options{'lanepath'}/02_MAPPING/starChimeric.out.bam" and !-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam") {
  #       $cmd = "samtools merge -n $options{'lanepath'}/02_MAPPING/accepted_hits\.bam $options{'lanepath'}/02_MAPPING/starAligned.out.bam $options{'lanepath'}/02_MAPPING/starChimeric.out.bam";
  #       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #       $cmd = "rm $options{'lanepath'}/02_MAPPING/starAligned.out.bam -f";
  #       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #    }
  #  }
  #  if ((-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam" or -s "$options{'lanepath'}/02_MAPPING/$mappedBam") and -s "$options{'lanepath'}/02_MAPPING/starAligned.out.sam"){
  #     my $cmd = "rm $options{'lanepath'}/02_MAPPING/starAligned.out.sam -f";
  #     RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #  }
  #}    #STAR mapping

  if ($options{'mapper'} eq 'gsnap') {

     my $quality_options;
     if ( ($options{'qual_zero'} - 33) < 10 ) {
        $quality_options = "--quality-protocol=sanger";
     } else {
        $quality_options = "--quality-zero-score=$options{'qual_zero'} --quality-print-shift=$options{'qual_move'}";
     }q

     unless (-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam" or -s "$options{'lanepath'}/02_MAPPING/$mappedBam" or -s "$options{'lanepath'}/02_MAPPING/accepted_hits\.sam") {
       my $cmd;
       my $zipOption = "--gunzip";
       if ($options{'bzip'}) {
         $zipOption = "--bunzip2";
       }

       if ($options{'seqType'} =~ /paired-end/) {
         my $fastqswap = swapfastq($options{'fastqFiles1'},$options{'fastqFiles2'});
         $cmd = "gsnap -d $options{'species'} -D $confs{'gmap_index'} --format=sam --nthreads=$options{'threads'} -s $confs{'gmap_splicesites'} --npaths=5 $quality_options $zipOption $fastqswap >$options{'lanepath'}/02_MAPPING/accepted_hits\.sam";
       } elsif ($options{'seqType'} =~ /single-end/) {
         $cmd = "gsnap -d $options{'species'} -D $confs{'gmap_index'} --format=sam --nthreads=$options{'threads'} -s $confs{'gmap_splicesites'} --npaths=5 $quality_options $zipOption --force-single-end $options{'fastqFiles1'} >$options{'lanepath'}/02_MAPPING/accepted_hits\.sam";
       }
       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
     }

     if (-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.sam" and (! -s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam" and ! -s "$options{'lanepath'}/02_MAPPING/$mappedBam")) {
        my $cmd = "samtools view -Sb -@ $options{'threads'} $options{'lanepath'}/02_MAPPING/accepted_hits\.sam -o $options{'lanepath'}/02_MAPPING/accepted_hits\.bam";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
     }

     if (-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.sam" and (-s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam" or -s "$options{'lanepath'}/02_MAPPING/$mappedBam")) {
        my $cmd  = "rm $options{'lanepath'}/02_MAPPING/accepted_hits\.sam -f";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
     }

  } #gsnap

  else {
     print STDERR "Error: --mapper option should only be gsnap or STAR. \n\n";
     exit 22;
  }


  #do the statistics
  unless (-e "$options{'lanepath'}/03_STATS") {
    my $cmd = "mkdir -p $options{'lanepath'}/03_STATS";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping\.stats") {
    my $typeop = ($options{'seqType'} =~ /single-end/)? "--type s" : "--type p";
    my $arpop = ($options{'seqType'} =~ /paired-end/)? "--arp $options{'lanepath'}/03_STATS/$options{'sampleName'}\.arp" : "";
    my $cmd = "$options{'bin'}/Rseq_bam_stats --mapping $options{'lanepath'}/02_MAPPING/accepted_hits\.bam $typeop --readlength $options{'readlen'} --maxIntron $options{'maxIntron'} --writer $options{'lanepath'}/02_MAPPING/accepted_hits\.mapped\.bam --unmapped $options{'lanepath'}/02_MAPPING/unmapped $arpop --breakpoint $options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping\.stats";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.gz" ) {
    if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints") {
       my $cmd = "gzip $options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints";
       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  unless (-s "$options{'lanepath'}/03_STATS/unmapped\.gz" ) {
    if (-s "$options{'lanepath'}/03_STATS/unmapped") {
       my $cmd = "gzip $options{'lanepath'}/03_STATS/unmapped";
       RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  my $mapping_stats_line_number = `wc -l $options{'lanepath'}/03_STATS/$options{'sampleName'}.mapping.stats`;
  $mapping_stats_line_number =~ s/^(\d+).*$/$1/;
  chomp($mapping_stats_line_number);
  if ($mapping_stats_line_number == 12) {
    my $total_reads = `$options{'decompress'} $options{'fastqFiles1'} | wc -l`;
    $total_reads /= 4;
    open STATS, ">>$options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping\.stats" || die "can not open $options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping\.stats\n";
    print STATS "total_frag: $total_reads\n";
    close STATS;
  }

  unless (-s "$options{'lanepath'}/02_MAPPING/$mappedBam") {
    my $cmd = "samtools sort -@ $options{'threads'} $options{'lanepath'}/02_MAPPING/accepted_hits\.mapped\.bam $options{'lanepath'}/02_MAPPING/accepted_hits\.mapped\.sorted";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/02_MAPPING/$mappedBam\.bai") {
    my $cmd = "samtools index $options{'lanepath'}/02_MAPPING/$mappedBam";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  if (-s "$options{'lanepath'}/02_MAPPING/$mappedBam" and -s "$options{'lanepath'}/02_MAPPING/accepted_hits\.mapped\.bam") {
    my $cmd = "rm $options{'lanepath'}/02_MAPPING/accepted_hits\.mapped\.bam -f";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$options{'lanepath'}/02_MAPPING/$mappedBam" and -s "$options{'lanepath'}/02_MAPPING/accepted_hits\.bam") {
    my $cmd = "rm $options{'lanepath'}/02_MAPPING/accepted_hits\.bam -f";
    #RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  if ($options{'bigWig'}) { #generate wiggle file
    if ($options{'strandWig'}) {
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bw") {
        if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph") {
           my $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        else {
           my $cmd = "samtools view $options{'lanepath'}/02_MAPPING/$mappedBam -h \| awk -F\'\t\' \'\$1 \~ \/\^\@\/ \|\| \$2 \=\= 99 \|\| \$2 \=\= 147 \|\| \$2 \=\= 97 \|\| \$2 \=\= 145 \|\| \$2 \=\= 355 \|\| \$2 \=\= 403\' \| samtools view -Sb - >$options{'lanepath'}/02_MAPPING/plus.bam";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "samtools view $options{'lanepath'}/02_MAPPING/$mappedBam -h \| awk -F\'\t\' \'\$1 \~ \/\^\@\/ \|\| \$2 \=\= 83 \|\| \$2 \=\= 163 \|\| \$2 \=\= 81 \|\| \$2 \=\= 161 \|\| \$2 \=\= 339 \|\| \$2 \=\= 419\' \| samtools view -Sb - >$options{'lanepath'}/02_MAPPING/minus.bam";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "genomeCoverageBed -ibam $options{'lanepath'}/02_MAPPING/plus.bam -bg -split -g $confs{'chromosomeSize'} >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bedgraph";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "genomeCoverageBed -ibam $options{'lanepath'}/02_MAPPING/minus.bam -bg -split -g $confs{'chromosomeSize'} >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bedgraph";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bedgraph" and -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bw") {
           my $cmd = "rm $options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bedgraph $options{'lanepath'}/03_STATS/$options{'sampleName'}\.minus\.bedgraph -f";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.plus\.bw" and -s "$options{'lanepath'}/02_MAPPING/plus.bam"){
           my $cmd = "rm $options{'lanepath'}/02_MAPPING/plus.bam $options{'lanepath'}/02_MAPPING/minus.bam -f";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      }
    } else { #unstranded
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bw") {
        if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph") {
           my $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        else {
           my $cmd = "genomeCoverageBed -ibam $options{'lanepath'}/02_MAPPING/$mappedBam -bg -split -g $confs{'chromosomeSize'} >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
           $cmd = "$options{'bin'}/bedGraphToBigWig $options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph $confs{'chromosomeSize'} $options{'lanepath'}/03_STATS/$options{'sampleName'}\.bw";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph" and -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bw") {
           my $cmd = "rm $options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedgraph -f";
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      }
    }
  }

  #for expression extimation;
  my $readingBam = "$options{'lanepath'}/02_MAPPING/$mappedBam";

  #ensembl gene#############################################################
  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr\.sorted") {
    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr") {
      my $cmd1;
      if ($options{'seqType'} =~ /paired-end/){
        $cmd1 = "$options{'bin'}/Rseq_bam_reads2expr --type p --region $confs{'ensembl_gene_bed'} --mapping $readingBam --posc $options{'lanepath'}/03_STATS/$options{'sampleName'}\.pos\.gff --chrmap $options{'lanepath'}/03_STATS/$options{'sampleName'}\.chrmap --lbias $options{'lanepath'}/03_STATS/$options{'sampleName'}\.lbias >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr";
      } else {
        $cmd1 = "$options{'bin'}/Rseq_bam_reads2expr --type s --region $confs{'ensembl_gene_bed'} --mapping $readingBam --posc $options{'lanepath'}/03_STATS/$options{'sampleName'}\.pos\.gff --chrmap $options{'lanepath'}/03_STATS/$options{'sampleName'}\.chrmap --lbias $options{'lanepath'}/03_STATS/$options{'sampleName'}\.lbias >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr";
      }
      RunCommand($cmd1,$options{'noexecute'},$options{'quiet'});
    }
    my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr.sorted";
    RunCommand($cmd2,$options{'noexecute'},$options{'quiet'});
    if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr.sorted" and -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr") {
      my $cmd3 = "rm $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr -rf";
      RunCommand($cmd3,$options{'noexecute'},$options{'quiet'});
    }
  }

  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.pos\.gff\.gz"){
    if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.pos\.gff"){
      my $cmd = "gzip $options{'lanepath'}/03_STATS/$options{'sampleName'}\.pos\.gff";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.count") {
    open ENSEMBL_GENE, "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr.sorted";
    open ENSEMBL_GENE_COUNT, ">$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.count";
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
  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.rpkm") {
    my $N_mapped_reads = 0;
    my $mapped = 0;
    my $singleton = 0;
    open MAPPING_STATS, "$options{'lanepath'}/03_STATS/$options{'sampleName'}.mapping.stats" || die "can not open $options{'lanepath'}/03_STATS/$options{'sampleName'}.mapping.stats";
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
    print STDERR "singletons equal to zero!!!\n" if ($singleton == 0 and $options{'seqType'} =~ /paired-end/);
    $N_mapped_reads = 2*$mapped - $singleton if ($options{'seqType'} =~ /paired-end/);
    $N_mapped_reads = $mapped if ($options{'seqType'} =~ /single-end/);
    exit if ($N_mapped_reads == 0);
    close MAPPING_STATS;

    open ENSEMBL_GENEMAP, "$confs{'ensembl_genemap'}";
    my %ensembl_genemap;
    while ( <ENSEMBL_GENEMAP> ) {
      chomp;
      my ($gene_id, $gene_name, $gene_type, $gene_desc, $entrez, $wiki_name, $wiki_desc) = split /\t/;
      $ensembl_genemap{$gene_id} = $gene_type."\t".$gene_name;
    }
    close ENSEMBL_GENEMAP;

    open ENSEMBL_GENE_EXPR, "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr.sorted" || die "can not open $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr.sorted";
    open ENSEMBL_RPKM, ">$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.rpkm";
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
  if ($confs{'species'} eq 'hg19') {
    #gencode gene#############################################################
    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr.sorted") {
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr") {
        my $cmd1;
        if ($options{'seqType'} =~ /paired-end/){
          $cmd1 = "$options{'bin'}/Rseq_bam_reads2expr --type p --region $confs{'gencode_gene_bed'} --mapping $readingBam >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr";
        } else {
          $cmd1 = "$options{'bin'}/Rseq_bam_reads2expr --type s --region $confs{'gencode_gene_bed'} --mapping $readingBam >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr";
        }
        RunCommand($cmd1,$options{'noexecute'},$options{'quiet'});
      }
      my $cmd2 = "sort -k 1,1d -k 2,2n -k 3,3n $options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr.sorted";
      RunCommand($cmd2,$options{'noexecute'},$options{'quiet'});
      if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr.sorted" and -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr") {
        my $cmd3 = "rm $options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr -rf";
        RunCommand($cmd3,$options{'noexecute'},$options{'quiet'});
      }
    }

    #for RPKM normalization of gencode genes##################################
    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.rpkm") {
      my $N_mapped_reads = 0;
      my $mapped = 0;
      my $singleton = 0;
      open MAPPING_STATS, "$options{'lanepath'}/03_STATS/$options{'sampleName'}.mapping.stats" || die "can not open $options{'lanepath'}/03_STATS/$options{'sampleName'}.mapping.stats";
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
      print STDERR "singletons equal to zero!!!\n" if ($singleton == 0 and $options{'seqType'} =~ /paired-end/);
      $N_mapped_reads = 2*$mapped - $singleton if ($options{'seqType'} =~ /paired-end/);
      $N_mapped_reads = $mapped if ($options{'seqType'} =~ /single-end/);
      exit if ($N_mapped_reads == 0);
      close MAPPING_STATS;

      open GENCODE_GENEMAP, "$confs{'gencode_genemap'}";
      my %gencode_genemap;
      while ( <GENCODE_GENEMAP> ) {
        chomp;
        my ($gene_id, $gene_type, $gene_name) = split /\t/;
        $gencode_genemap{$gene_id} = $gene_type."\t".$gene_name;
      }
      close GENCODE_GENEMAP;

      open GENCODE_GENE_EXPR, "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr.sorted" || die "can not open $options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.expr.sorted";
      open GENCODE_RPKM, ">$options{'lanepath'}/03_STATS/$options{'sampleName'}\.gencode\_gene\.rpkm";
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


  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.cate") {
    my $cmd = "perl $options{'bin'}/cate.pl $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.expr\.sorted $confs{'gene_annotation'} >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.cate";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  if ($options{'seqType'} =~ /paired-end/) { #do the insert size only if it is paired-end
    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins\.gz") {
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins") {
         my $cmd = "samtools view -f 0x2 $options{'lanepath'}/02_MAPPING/$mappedBam | cut -f 9 | awk \'\$1\>0 \&\& \$1\<500\' >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins";
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins") {
         my $cmd = "gzip $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins";
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
  } #insert size

  unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.report/$options{'sampleName'}\.report\.html") {
    my $qcmatesuffix1 = '.qc';
    my $qcmatesuffix2 = '.qc';
    my $cmd;
    if ($options{'seqType'} =~ /paired-end/) {
      $cmd  = "$options{'Rbinary'} CMD BATCH --no-save --no-restore "."\'--args path=\"$options{'lanepath'}\" lane=\"$options{'sampleName'}\" anno=\"$confs{'anno'}\" species=\"$options{'species'}\" src=\"$options{'bin'}\" readlen=$options{'readlen'} gf=\"$options{'gf'}\" qcsuffix1=\"$qcmatesuffix1\" qcsuffix2=\"$qcmatesuffix2\" type=\"p\"' $options{'bin'}/html_report.R $options{'lanepath'}/03_STATS/R\_html\.out";
    } else {
      $cmd  = "$options{'Rbinary'} CMD BATCH --no-save --no-restore "."\'--args path=\"$options{'lanepath'}\" lane=\"$options{'sampleName'}\" anno=\"$confs{'anno'}\" species=\"$options{'species'}\" src=\"$options{'bin'}\" readlen=$options{'readlen'} gf=\"$options{'gf'}\" qcsuffix1=\"$qcmatesuffix1\" type=\"s\"' $options{'bin'}/html_report.R $options{'lanepath'}/03_STATS/R\_html\.out";
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
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

  $options{'RA'} = 1 if $options{'mapper'} eq 'gsnap';

  if ($options{'RA'} != 0) {  #regional assembly test###############################################

    my $assembly_type = "independent regional assembly";

    printtime();
    print STDERR "Regional assembly process is starting \($assembly_type\)... first breakpoint processing...\n\n";

    unless (-e "$options{'lanepath'}/04_ASSEMBLY") {
      my $cmd = "mkdir -p $options{'lanepath'}/04_ASSEMBLY";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted\.gz") {
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.gz" or -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted") {
        print STDERR "Error: breakpoint file does not exist, do runlevel 2 first.\n\n";
        exit 22;
      }
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted") {
        my $tmpOpt = ($options{'tmpDir'} eq '')? "":"-T $options{'tmpDir'}";
        my $cmd = "gzip -d -c $options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.gz | sort $tmpOpt -k 1,1d -k 2,2n >$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if ( -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.gz" and -s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted" ) {
         my $cmd2 = "rm $options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.gz -f";
         RunCommand($cmd2,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted") {
         my $cmd3 = "gzip $options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints.sorted";
         RunCommand($cmd3,$options{'noexecute'},$options{'quiet'});
      }
    }

    #breakpoints
    my $breakpointSource = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.breakpoints\.sorted\.gz";
    print STDERR "the breakpoint sources are: $breakpointSource\n";
    #breakpoint Source is defined

    unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed") {
      my $typeop = ($options{'seqType'} =~ /paired-end/)? "p":"s";
      my $cmd = "perl $options{'bin'}/breakpoint\_processing.pl $breakpointSource $typeop $options{'maxIntron'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted") {
      my $tmpOpt = ($options{'tmpDir'} eq '')? "":"-T $options{'tmpDir'}";
      my $cmd = "sort $tmpOpt -k 3,3d -k 4,4n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot"){
      my $cmd = "perl $options{'bin'}/breakpoint_repeat_masker.pl $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted $confs{'repeatMasker'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if ($options{'seqType'} =~ /paired-end/) {
      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco") {
        my $mapping_bam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
        die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
        my $cmd = "$options{'bin'}/discordant_consistency --region $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot --mapping $mapping_bam --maxIntron $options{'maxIntron'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted"){
        my $cmd = "sort -k 3,3d -k 4,4n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-e "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter"){
        my $cmd = "perl $options{'bin'}/breakpoint_final_prepare.pl $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted") {
        unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate") {
          my $mapping_bam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
          die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
          my $largestBPID = `tail -1 $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}.breakpoints.processed | cut -f 1`;
          $largestBPID =~ s/\n//;
          my $cmd = "$options{'bin'}/discordant_mate --mapping $mapping_bam --idstart $largestBPID --consisCount $options{'consisCount'} --maxIntron $options{'maxIntron'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        if (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate") {
          my $cmd = "sort -k 3,3d -k 4,4n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      }

      if (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate" and -s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted") {
        my $cmd = "rm $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate -f";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter") {
        my $cmd = "perl $options{'bin'}/discordant_mate_processing.pl $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted $confs{'repeatMasker'} $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined") {
        my $cmd = "cat $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.disco\.sorted\.filter $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot\.dismate\.sorted\.filter >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted") {
        my $cmd = "sort -k 3,3d -k 4,4n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      if ($options{'mapper'} eq 'star') {  #merge breakpoints for star mapper
        unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints") {
          my $largestBPID = `sort -k 1,1n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}.breakpoints.processed.sorted.repornot.dismate.sorted | tail -1 | cut -f 1`;
          $largestBPID =~ s/\n//;
          my $cmd = "perl $options{'bin'}/starjunction2bp.pl $largestBPID $options{'lanepath'}/02_MAPPING/starChimeric\.out\.junction $options{'maxIntron'} $confs{'repeatMasker'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints.masked") {
          my $cmd = "perl $options{'bin'}/intersectFiles.pl -o $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints -oichr 2 -oistart 3 -oiend 3 -m $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted -michr 2 -mistart 3 -miend 3 -t 100 -count >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints.masked";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
          $cmd = "awk \-F\"\\t\" \'\$8 \=\= \"N\" \&\& \$10 \=\= 0\' $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.star.breakpoints.masked | cut -f 1,2,3,4,5,6,7,8,9 >>$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined";    #append star breakpoints
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
          $cmd = "sort -k 3,3d -k 4,4n $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted";                 #re-sort
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      } #merge breakpoints for star mapper

    } #if it is paired-end reads

    else {  #it is single-end reads
      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted") {
        my $cmd = "awk \'\(\$2 == \"s\" && \$5 >= 5 && \$8 == \"N\") || \(\$2 == \"p\"\)' $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.sorted\.repornot >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    } #single-end reads


    if ($options{'RA'} == 1) { #independent regional assembly
      unless (-s "$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted\.reads") {
        my $mapping_bam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
        die "Error: the mapping bam file is not available." unless (-e $mapping_bam);
        my $typeop = ($options{'seqType'} =~ /paired-end/)? "--type p":"--type s";
        my $cmd = "$options{'bin'}/reads_in_region --region $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted --mapping $mapping_bam $typeop --maxIntron $options{'maxIntron'} >$options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted\.reads";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      my @RAssembly_reads;
      @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.RAssembly\.fq") if ($options{'seqType'} =~ /paired-end/);
      @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.RAssembly\.fq") if ($options{'seqType'} =~ /single-end/);
      @RAssembly_reads = uniqueArray(\@RAssembly_reads) if ($options{'seqType'} =~ /single-end/);
      my @RAssembly_reads_gz;
      @RAssembly_reads_gz = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.RAssembly\.fq\.gz") if ($options{'seqType'} =~ /paired-end/);
      @RAssembly_reads_gz = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.RAssembly\.fq\.gz") if ($options{'seqType'} =~ /single-end/);
      @RAssembly_reads_gz = uniqueArray(\@RAssembly_reads_gz) if ($options{'seqType'} =~ /single-end/);
      unless ( ($options{'seqType'} =~ /paired-end/ and ($#RAssembly_reads == 1 or $#RAssembly_reads_gz == 1)) or ($options{'seqType'} =~ /single-end/ and ($#RAssembly_reads == 0 or $#RAssembly_reads_gz == 0)) ) { #get raw reads

        #need to do get raw reads here
        my @reads;
        if ($options{'seqType'} =~ /paired-end/) {    #paired-end
          @reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.fq\.$options{'zipSuffix'}"); #original reads
          @reads = mateorder(\@reads);
        } else {                                      #single-end
          @reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.fq\.$options{'zipSuffix'}"); #original reads
          @reads = uniqueArray(\@reads);
        }

        my $cmd;
        if ($options{'seqType'} =~ /paired-end/){
          $cmd = "perl $options{'bin'}/pick_ARP.pl --arpfile $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted\.reads --readfile1 $reads[0] --readfile2 $reads[1] --RA";
        } else {
          $cmd = "perl $options{'bin'}/pick_ARP.pl --arpfile $options{'lanepath'}/04_ASSEMBLY/$options{'sampleName'}\.breakpoints\.processed\.filter\.combined\.sorted\.reads --readfile1 $reads[0] --RA";
        }
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    } else {
      print STDERR "Error: --RA must be set to 1 or 0 (default)...\n\n";
      exit 22;
    }

    printtime();
    print STDERR "now assembling...\n\n";

    unless (-s "$options{'lanepath'}/04_ASSEMBLY/transcripts.fa") {

      if ($options{'RA'} == 1) {

        my @RAssembly_reads;
        @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.RAssembly\.fq") if ($options{'seqType'} =~ /paired-end/);
        @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.RAssembly\.fq") if ($options{'seqType'} =~ /single-end/);
        @RAssembly_reads = uniqueArray(\@RAssembly_reads) if ($options{'seqType'} =~ /single-end/);
        my @RAssembly_reads_gz;
        @RAssembly_reads_gz = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.RAssembly\.fq\.gz") if ($options{'seqType'} =~ /paired-end/);
        @RAssembly_reads_gz = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.RAssembly\.fq\.gz") if ($options{'seqType'} =~ /single-end/);
        @RAssembly_reads_gz = uniqueArray(\@RAssembly_reads_gz) if ($options{'seqType'} =~ /single-end/);

        if (($options{'seqType'} =~ /paired-end/ and $#RAssembly_reads != 1) or ($options{'seqType'} =~ /single-end/ and $#RAssembly_reads != 0)) {
          if (($options{'seqType'} =~ /paired-end/ and $#RAssembly_reads_gz == 1) or ($options{'seqType'} =~ /single-end/ and $#RAssembly_reads_gz == 0)) {
            foreach my $RA_read_file_gz (@RAssembly_reads_gz) {
               (my $RA_read_file = $RA_read_file_gz) =~ s/\.gz$//;
               my $cmd = "gzip -d -c $RA_read_file_gz >$RA_read_file";
               RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            }
          } else {
            print STDERR "Error: RA read files are not correctly generated...\n\n";
            exit 22;
          }
          @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.RAssembly\.fq") if ($options{'seqType'} =~ /paired-end/);
          @RAssembly_reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.RAssembly\.fq") if ($options{'seqType'} =~ /single-end/);
          @RAssembly_reads = uniqueArray(\@RAssembly_reads) if ($options{'seqType'} =~ /single-end/);
        }

        @RAssembly_reads = mateorder(\@RAssembly_reads) if ($options{'seqType'} =~ /paired-end/);
        my $ra_reads1_file = $RAssembly_reads[0];
        my $ra_reads2_file = $RAssembly_reads[1] if ($options{'seqType'} =~ /paired-end/);

        my %ra_reads1_jumpers;
        my %ra_reads2_jumpers;

        find_jumpers($ra_reads1_file, \%ra_reads1_jumpers);
        find_jumpers($ra_reads2_file, \%ra_reads2_jumpers) if ($options{'seqType'} =~ /paired-end/);

        my $n_ids_ra_reads1 = scalar(keys %ra_reads1_jumpers);
        my $n_ids_ra_reads2 = scalar(keys %ra_reads2_jumpers) if ($options{'seqType'} =~ /paired-end/);
        if ($options{'seqType'} =~ /paired-end/ and ($n_ids_ra_reads1 != $n_ids_ra_reads2)) {
          print STDERR "Error: the numbers of breakpoints are inconsistant between read pairs ($n_ids_ra_reads1 vs $n_ids_ra_reads2).\n";
          exit 22;
        }

        #now should do regional assembly one by one
        open RA_READS1, "$ra_reads1_file";
        open RA_READS2, "$ra_reads2_file" if ($options{'seqType'} =~ /paired-end/);
        my $speed_count = 0;
        my %printed_speed;
        foreach my $bp_id (sort {$a<=>$b} keys %ra_reads1_jumpers) {

          if ($options{'idra'} != 0 and $bp_id != $options{'idra'}) {
             next;
          }

          my $cmd = "mkdir -p $options{'lanepath'}/04_ASSEMBLY/$bp_id"; #create dir
          RunCommand($cmd,$options{'noexecute'},1);

          regional_assembly(\*RA_READS1, $ra_reads1_jumpers{$bp_id}, $bp_id, "$options{'lanepath'}/04_ASSEMBLY/$bp_id/reads_1.fq");
          regional_assembly(\*RA_READS2, $ra_reads2_jumpers{$bp_id}, $bp_id, "$options{'lanepath'}/04_ASSEMBLY/$bp_id/reads_2.fq") if ($options{'seqType'} =~ /paired-end/);

          unless (-s "$options{'lanepath'}/04_ASSEMBLY/$bp_id/Roadmaps") {
            my $cmd;
            if ($options{'seqType'} =~ /paired-end/){
              $cmd = "velveth $options{'lanepath'}/04_ASSEMBLY/$bp_id/ 21 -fastq -shortPaired -separate $options{'lanepath'}/04_ASSEMBLY/$bp_id/reads_1.fq $options{'lanepath'}/04_ASSEMBLY/$bp_id/reads_2.fq";
            } else {
              $cmd = "velveth $options{'lanepath'}/04_ASSEMBLY/$bp_id/ 21 -fastq -short $options{'lanepath'}/04_ASSEMBLY/$bp_id/reads_1.fq";
            }
            RunCommand($cmd,$options{'noexecute'},1);
          }

          unless (-s "$options{'lanepath'}/04_ASSEMBLY/$bp_id/Graph2") {
            my $cmd = "velvetg $options{'lanepath'}/04_ASSEMBLY/$bp_id/ -read_trkg yes";
            RunCommand($cmd,$options{'noexecute'},1);
          }

          unless (-s "$options{'lanepath'}/04_ASSEMBLY/$bp_id/transcripts.fa"){
            my $cmd = "oases $options{'lanepath'}/04_ASSEMBLY/$bp_id/";
            RunCommand($cmd,$options{'noexecute'},1);
          }

          if (-s "$options{'lanepath'}/04_ASSEMBLY/$bp_id/transcripts.fa") {
            open TRANSCRIPTS, "$options{'lanepath'}/04_ASSEMBLY/$bp_id/transcripts.fa";
            open TRANSCRIPTSALL, ">>$options{'lanepath'}/04_ASSEMBLY/transcripts.fa";
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

          if (-e "$options{'lanepath'}/04_ASSEMBLY/$bp_id/") {
            my $cmd = "rm $options{'lanepath'}/04_ASSEMBLY/$bp_id/ -rf";
            RunCommand($cmd,$options{'noexecute'},1) unless $options{'idra'} != 0;
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
        close RA_READS2 if ($options{'seqType'} =~ /paired-end/);

        if (($options{'seqType'} =~ /paired-end/ and $#RAssembly_reads_gz != 1) or ($options{'seqType'} =~ /single-end/ and $#RAssembly_reads_gz != 0)) { #gz file does not exist
          if (($options{'seqType'} =~ /paired-end/ and $#RAssembly_reads == 1) or ($options{'seqType'} =~ /single-end/ and $#RAssembly_reads == 0)) {
            foreach my $RA_read_file (@RAssembly_reads) {
              my $cmd = "gzip $RA_read_file";
              RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            }
          } else {
            print STDERR "Error: RA read files are not correctly generated... end of the regional assembly, very strange...\n\n";
            exit 22;
          }
        } else {
          foreach my $RA_read_file (@RAssembly_reads) {
             my $cmd = "rm $RA_read_file -f";
             RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
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

  unless (-e "$options{'lanepath'}/05_FUSION") {
    my $cmd = "mkdir -p $options{'lanepath'}/05_FUSION";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  if ($options{'GMAP'}) {  #using GMAP
      unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.gmap") {
        my $speciesIndex = (-e "$confs{'gmap_index'}/$options{'species'}\-all\/$options{'species'}\-all.genomecomp")? $options{'species'}."\-all":$options{'species'};
        my $cmd = "gmap -D $confs{'gmap_index'} -d $speciesIndex --format=psl -t $options{'threads'} --intronlength=230000 $options{'lanepath'}/04_ASSEMBLY/transcripts\.fa >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.gmap";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq") {
        my $cmd = "perl $options{'bin'}/pick_fusion_transcripts_from_genomeBLAT.pl --transcripts $options{'lanepath'}/04_ASSEMBLY/transcripts.fa --uniqueBase $options{'uniqueBase'} --misPen $options{'misPen'} --maxIntron $options{'maxIntron'} $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.gmap >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration 2>$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
  } #GMAP

  else {  #using BLAT
      unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.blat") {
        my $cmd = "blat -maxIntron=230000 $confs{'blatDatabase'} $options{'lanepath'}/04\_ASSEMBLY/transcripts\.fa $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.blat";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

      unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq") {
        my $cmd = "perl $options{'bin'}/pick_fusion_transcripts_from_genomeBLAT.pl --transcripts $options{'lanepath'}/04_ASSEMBLY/transcripts.fa --uniqueBase $options{'uniqueBase'} --misPen $options{'misPen'} --maxIntron $options{'maxIntron'} $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.transcripts\.genome\.blat >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration 2>$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
  } #BLAT

  #build fusion candidate index (now bowtie2)
  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq\.index\.1\.bt2") {
    my $cmd = "bowtie2-build --quiet $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq\.index";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #map reads to the fusion candidate index
  unless ( -s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.bam" ) {

    if ( -s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam" ) { #bowtie2 is done

      my $cmd = "samtools view -Sb $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam -o $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2.bam";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      if ( -s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2.bam" ) {
         my $cmd = "rm $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam -f";
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }

    } else {  #bowtie2 is not done

      my @reads;
      if ($options{'seqType'} =~ /paired-end/) {
        @reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}\_{R,}[123]\.fq\.$options{'zipSuffix'}"); #original reads
        @reads = mateorder(\@reads);
      } else {
        @reads = bsd_glob("$options{'lanepath'}/01_READS/$options{'sampleName'}*\.fq\.$options{'zipSuffix'}");
      }

      my $readsop = ($options{'seqType'} =~ /paired-end/)? "-1 $reads[0] -2 $reads[1]":"-U $reads[0]";
      my $cmd = "bowtie2 -k 22 -p $options{'threads'} --no-unal --score-min L,-2,-0.15 -x $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq\.index $readsop >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      $cmd = "samtools view -Sb $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam -o $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2.bam";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      if (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2.bam"){
         my $cmd = "rm $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.sam -f";
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    } #bowtie2 is not done
  } #bowtie2 bam is generated

  #get fusion coverage
  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov") {
    my $gfc_opts = '';
    $gfc_opts = "--genomeBlatPred $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration";

    my $readingBam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
    if ($options{'customMappedBam'} ne '') {
       $readingBam = $options{'customMappedBam'};
    }

    my $cmd;
    if ($options{'seqType'} =~ /paired-end/) { #paired-end
      $cmd = "perl $options{'bin'}/get_fusion_coverage.pl $gfc_opts --type pair --maxIntron $options{'maxIntron'} --mappingfile $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.bam --readlength $options{'readlen'} --geneanno $confs{'gene_annotation'} --genelength $confs{'ensembl_gene_len'} --repeatmasker $confs{'repeatMasker'} --selfChain $confs{'selfChain'} --accepthits $readingBam --encomcov $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov.enco >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    } else { #single-end
      $cmd = "perl $options{'bin'}/get_fusion_coverage.pl $gfc_opts --type single --maxIntron $options{'maxIntron'} --mappingfile $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie2\.bam --readlength $options{'readlen'} --geneanno $confs{'gene_annotation'} --genelength $confs{'ensembl_gene_len'} --repeatmasker $confs{'repeatMasker'} --selfChain $confs{'selfChain'} --accepthits $readingBam >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov";
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #visualize coverage
  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis") {
    my $fpv_opts = '';
    $fpv_opts = "--genomeBlatPred";
    my $typeop = ($options{'seqType'} =~ /^paired-end/)? "--seqType p":"--seqType s";
    my $cmd = "perl $options{'bin'}/further_processing_for_read_visualization.pl $fpv_opts $typeop --fusionseqfile $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.seq --coveragefile $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov --readlength $options{'readlen'} >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.list"){
  #  my $sorting_opts = '';
  #  $sorting_opts = "-k 16,16nr -k 18,18nr -k 3,3nr -k 8,8d -k 11,11d";
  #  my $cmd = "grep \"^\#\" $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis \| sort $sorting_opts >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.list";
  #  RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  #}

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa") {
    my $cmd = "perl $options{'bin'}/vis2fa.pl $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis $options{'readlen'} >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl"){
    my $cmd = "blat -maxIntron=230000 $confs{'blatDatabase'} $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl";
    if ($options{'GMAP'}) {
       my $speciesIndex = (-e "$confs{'gmap_index'}/$options{'species'}\-all\/$options{'species'}\-all.genomecomp")? $options{'species'}."\-all":$options{'species'};
       $cmd = "gmap -D $confs{'gmap_index'} -d $speciesIndex --format=psl -t $options{'threads'} --intronlength=230000 $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl";
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl\.pass") {
    my $cmd = "perl $options{'bin'}/pick_fusion_transcripts_from_genomeBLAT.pl --identity 94 --final 1 --uniqueBase $options{'uniqueBase'} --misPen $options{'misPen'} --maxIntron $options{'maxIntron'} $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl\.pass";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.list"){
    my $sorting_opts = '';
    $sorting_opts = "-k 16,16nr -k 18,18nr -k 3,3nr -k 8,8d -k 11,11d";
    my $cmd = "perl $options{'bin'}/fusionFinalGrep.pl $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis\.fa\.psl\.pass $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.bowtie\.cov\.vis | sort $sorting_opts >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.list";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion\.report") {
    my $cmd = "perl $options{'bin'}/fusion_report.pl $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion_transcirpts_after_filtration\.list >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion\.report";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion\.report\.filtered") {
    my $cmd = "perl $options{'bin'}/fusionFinalFiltration.pl $options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion\.report >$options{'lanepath'}/05_FUSION/$options{'sampleName'}\.fusion\.report\.filtered";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
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

  unless (-e "$options{'lanepath'}/06_CUFFLINKS") {
    my $cmd = "mkdir -p $options{'lanepath'}/06_CUFFLINKS";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/06_CUFFLINKS/transcripts\.gtf") {

    my $cufflinks_options = "";
    my $known_trans_file = $confs{'gene_annotation_gtf'};
    if ($options{'known_trans'} eq 'refseq') {
      $known_trans_file = $confs{'refseq_gene_gtf'};
    }

    if ( $options{'gtf_guide_assembly'} ) {
      $cufflinks_options .= "--GTF-guide $known_trans_file ";
    }
    else {
      $cufflinks_options .= "--GTF $known_trans_file ";
    }
    if ( $options{'frag_bias_correct'} ) {
      $cufflinks_options .= "--frag-bias-correct ";
    }
    if ( $options{'upper_quantile_norm'} ){
      $cufflinks_options .= "--upper-quartile-norm ";
    }

    my $mapping_bam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
    if (-s $mapping_bam){
      my $cmd = "cufflinks -o $options{'lanepath'}/06_CUFFLINKS -p $options{'threads'} $cufflinks_options --quiet $mapping_bam";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    else {
      print STDERR "$mapping_bam does not exist, please do the mapping first.\n";
      exit;
    }
  }

  if (-s "$options{'lanepath'}/06_CUFFLINKS/transcripts\.gtf") { #collect the gtf list
     my $current_gtf = "$options{'lanepath'}/06_CUFFLINKS/transcripts\.gtf";
     my $gtf_list = "$options{'root'}/GTF\_list\.txt";

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

  unless (-s "$options{'lanepath'}/06_CUFFLINKS/$options{'sampleName'}\.transcripts\.gtf\.tmap") {
    my $cmd = "cuffcompare -o $options{'lanepath'}/06_CUFFLINKS/$options{'sampleName'} -r $confs{'refseq_gene_gtf'} $options{'lanepath'}/06_CUFFLINKS/transcripts\.gtf";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
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

  my $cuffmerge_output = "$options{'root'}/cuffmerge";
  my $cuffdiff_output = "$options{'root'}/cuffdiff";
  my $gtf_list = "$options{'root'}/GTF\_list\.txt";

  unless (-s "$gtf_list") {
    print STDERR "Error: $gtf_list file does not exist!!! Please rerun RTrace.pl for each sample for runlevel 5 (cufflinks).\n\n";
    exit 22;
  }

  #cuffmerge####################################
  unless (-s $cuffmerge_output) {
    my $cmd = "mkdir -p $cuffmerge_output";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$cuffmerge_output/merged\.gtf") {
    my $cmd = "cuffmerge -o $cuffmerge_output --ref-gtf $confs{'gene_annotation_gtf'} --num-threads $options{'threads'} --ref-sequence $confs{'GFASTA'} $gtf_list";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #cuffdiff#####################################
  unless (-e $cuffdiff_output) {
    my $cmd = "mkdir -p $cuffdiff_output";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$cuffdiff_output/splicing\.diff") {
    my %bam_files_cd;
    open TARGETS, "$options{'root'}/edgeR/targets" || die "Error: could not open $options{'root'}/edgeR/targets for the information of bam files for cuffdiff.\n please rerun the pipeline for each samples with --patient and --tissue arguments set.\n";
    while ( <TARGETS> ) {
       chomp;
       next if /^label/;
       my ($cu_sample, $cu_expr_count, $tissue, $patient) = split /\t/;
       my $cu_bam_file = "$options{'root'}/$cu_sample/02_MAPPING/$mappedBam";
       if (-e "$cu_bam_file"){
         push (@{$bam_files_cd{$tissue}}, $cu_bam_file);
       } else {
         my $dirname = dirname($cu_expr_count);
         $cu_bam_file = "$dirname/../02_MAPPING/$mappedBam";
         if (-e "$cu_bam_file"){
            push (@{$bam_files_cd{$tissue}}, $cu_bam_file);
         } else {
            print STDERR "error: $cu_bam_file does not exist!!!\n";
            exit 22;
         }
       }
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

    my $cmd = "cuffdiff --output-dir $cuffdiff_output --num-threads $options{'threads'} --labels $tissue_names[0]\,$tissue_names[1] $cuffmerge_output/merged\.gtf $bam_files1 $bam_files2";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  } #run cuffdiff

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}



##write target file --- for edgeR (and cuffdiff) ########################
if ($options{'patient'} and $options{'tissue'}) {
  my $count_file = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ensembl\_gene\.count";
  unless (-e "$options{'root'}/edgeR") {
    my $cmd = "mkdir -p $options{'root'}/edgeR";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (! -s "$options{'root'}/edgeR/targets") {
    open TARGETS_OUT, ">>$options{'root'}/edgeR/targets";
    printf TARGETS_OUT "%s\n", join("\t", 'label', 'files', 'tissue', 'patient');
    printf TARGETS_OUT "%s\n", join("\t", $options{'sampleName'}, $count_file, $options{'tissue'}, $options{'patient'});
    close TARGETS_OUT;
  } else {
    my $writeornot = 1;
    open TARGETS_IN, "$options{'root'}/edgeR/targets";
    while ( <TARGETS_IN> ){
       chomp;
       my @cols = split /\t/;
       $writeornot = 0 if ($cols[0] eq $options{'sampleName'});
    }
    close TARGETS_IN;
    if ($writeornot == 1){
      open TARGETS_OUT, ">>$options{'root'}/edgeR/targets";
      printf TARGETS_OUT "%s\n", join("\t", $options{'sampleName'}, $count_file, $options{'tissue'}, $options{'patient'});
      close TARGETS_OUT;
    }
  } #existing targets file
} elsif ($options{'patient'} or $options{'tissue'}) {
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

  my $edgeR_output = "$options{'root'}/edgeR";

  print STDERR "DE analysis parameters:\n";
  print STDERR "priordf: $options{'priordf'}\; spaired: $options{'spaired'}\; $options{'pairDE1'}\-$options{'pairDE2'}\n";

  if (! -e $edgeR_output) {
     print STDERR "Error: $edgeR_output dir does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  } elsif (! -e "$edgeR_output/targets") {
     print STDERR "Error: $edgeR_output/targets file does not exist!!! Please rerun RTrace.pl for each sample with --patient and --tissue set.\n\n";
     exit 22;
  }

  unless (-s "$edgeR_output/topDE.txt") {
    my $cmd = "$options{'Rbinary'} CMD BATCH --no-save --no-restore "."\'--args path=\"$edgeR_output\" priordf=$options{'priordf'} spaired=$options{'spaired'} pair1=\"$options{'pairDE1'}\" pair2=\"$options{'pairDE2'}\"\' $options{'bin'}/edgeR_test.R $edgeR_output/R\_html\.out";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$edgeR_output/DE\_table\.txt") {
    my $cmd = "perl $options{'bin'}/get_DEtable.pl $edgeR_output topDE\.txt ifDE\.txt $options{'spaired'} >$edgeR_output/DE\_table\.txt";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}

=pod
###
###runlevel8: SNV/INDEL calling from bam file
###

$runlevels = 8;
if (exists $runlevel{$runlevels}) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$options{'lanepath'}/08_VARIANTS") {
    my $cmd = "mkdir -p $options{'lanepath'}/08_VARIANTS";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf") {

    unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf"){
      my $mapping_bam = "$options{'lanepath'}/02_MAPPING/$mappedBam";
      if (-s $mapping_bam){
        my $cmd = "samtools mpileup -DSEugd 1000 -q 1 -C 50 -f $confs{'GFASTA'} $mapping_bam | bcftools view -p 0.9 -vcg - >$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      else {
        print STDERR "$mapping_bam does not exist, please do the mapping first.\n";
        exit;
      }
    }

    if (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf"){
      my $cmd = "cat $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf | vcf-sort >$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf" and -s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf"){
      my $cmd = "rm $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.vcf";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  #do annovar
  unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf") {

    unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar"){
      my $cmd = "perl $annovarbin/convert2annovar.pl --format vcf4 --includeinfo --allallele $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf >$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv") {
      my $cmd = "perl $annovarbin/summarize_annovar1.pl -outfile $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary -ver1000g 1000g2012feb -verdbsnp 135 -genetype=refgene --buildver hg19 $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar $annovarDB";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf") {
      my $cmd = "perl $annovarbin/convert_annovar_vcf-all-samples.pl $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf $vcfheader";
      #my $cmd = "perl $options{'bin'}/convert_annovar_vcf-all-samples.pl $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.exome_summary.csv $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf $vcfheader";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf") {
      my $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.hg19_* ";
      $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.exonic_variant_function ";
      $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.variant_function ";
      $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar ";
      $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.exome_summary.csv ";
      $needtobedeleted .= "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv ";
      my $cmd = "rm $needtobedeleted -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  } #do annovar

  unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf.snv") {
    my $cmd = "grep -v \"\^\#\" $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf \| awk -F\'\t\' \'\$8 \!\~ \/INDEL\/\' > $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf.snv";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf.indel") {
    my $cmd = "grep -v \"\^\#\" $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf \| awk -F\'\t\' \'\$8 \~ \/INDEL\/\' > $options{'lanepath'}/08_VARIANTS/$options{'sampleName'}\.transcriptome.sorted.vcf.annovar.summary.genome_summary.csv.vcf.indel";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

} #runevel 8 SNV/INDEL calling
=cut


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
  print STDERR "\nGENERAL OPTIONS (MUST SET):\n\t--runlevel\tthe steps of runlevel, from 1-8, either rl1-rl2 or rl. See below for options for each runlevel.\n";
  print STDERR "\t--sampleName\tthe name of the lane needed to be processed (must set for runlevel 1-5)\n";
  print STDERR "\t--seqType\tset to 's' if it is a single-end sequencing experiment, or 'p' for paired-end (default).\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$options{'bin'}/../PIPELINE/, MUST set using other dir)\n";
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

  print STDERR "\nrunlevel 2: mapping and report of mapping statistics\n";
  print STDERR "\t--mapper\tthe mapper used in runlevel2, now support \'gsnap\' (default), \'star\' and \'tophat2\' (no fusion detection if chosen).\n";
  print STDERR "\t--seglen\tthe segment length for tophat mapping (default 25)\n";
  print STDERR "\t--gf\t\tthe graphical format in mapping report, \'png\' (default) or \'pdf\' (when a x11 window is not available)\n";
  print STDERR "\t--WIG\t\tgenerate a big wiggle file in run-level 2.\n";
  print STDERR "\t--insertmean\tthe mean insert size of read mates (not required, can be decided automatically)\n";
  print STDERR "\t--insertsd\tthe SD of insert size of read mates (not required, can be decided automatically)\n";

  print STDERR "\nrunlevel 3: selecting anormalous/breakpoint-surrouding reads and preforming the assembly\n";
  print STDERR "\t--RA\t\tuse regional assembly for runlevel 3. Default is to set.\n";
  print STDERR "\t--consisCount\tnumber of consistent read pairs with discordant mapping (default: 5). use smaller value For <70bp reads or low depth data.\n";

  print STDERR "\nrunlevel 4: detection of fusion candidates\n";
  print STDERR "\t--GMAP\t\tset if use GMAP in run-level 4 (default: use BLAT).\n";
  print STDERR "\t--misPen\tthe penalty for mismatches for scoring each blat record (default: $options{'misPen'}).\n";
  print STDERR "\t--uniqueBase\tbase with added identity-score below this would be regarded as unique base for each blat record (default: $options{'uniqueBase'}).\n";

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

  #print STDERR "\nrunlevel 8: call SNV/INDEL from the mapping bam files using samtools\n";

  print STDERR "\nOTHER OPTIONS\n";
  print STDERR "\t--maxIntron\tthe maximum length of intron, to select for discordant mapping (fusion detection, default: $options{'maxIntron'})\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--quiet\t\tdo not print the command line calls and time information\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n";
  print STDERR "\t--tmpDir\ttmp dir for generating large tmp files\n";
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
  my $r = @_;
  my @tmp;

  if (scalar(@{$r}) != 2) {
    print STDERR "there are not two mate read files, exit.\n";
    exit 1;
  }

  foreach my $m (@{$r}){
    if ($m =~ /_R?1\./) {
      $tmp[0] = $m;
    } elsif ($m =~ /_R?[23]\./) {
      $tmp[1] = $m;
    }
    else {
      print STDERR "the mate read file dosen't contain _1 _2 or _3 information, exit.\n";
      exit 1;
    }
  }

  return @tmp;
}

sub swapfastq {
  my $fastq1 = shift;
  my $fastq2 = shift;
  my @fastq1 = split(' ', $fastq1);
  my @fastq2 = split(' ', $fastq2);
  my $outfastq = '';
  for (my $i = 0; $i <= $#fastq1; $i++){
    $outfastq .= join(' ', $fastq1[$i],$fastq2[$i]).' ';
  }
  $outfastq =~ s/\s$//;
  return $outfastq;
}

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub uniqueArray {
   my $array = shift;
   my %arraytmp;
   foreach my $item (@{$array}){
     $arraytmp{$item} = '';
   }
   my @arraytmp = keys %arraytmp;
   return @arraytmp;
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
