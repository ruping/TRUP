use strict;

my %cate;

$cate{'snRNA'} = 'snRNA';
$cate{'protein_coding'} = 'protein_coding';
$cate{'pseudogene'} = 'pseudogene';
$cate{'processed_pseudogene'} = 'pseudogene';
$cate{'processed_transcript'} = 'pseudogene';
$cate{'retained_intron'} = 'protein_coding';
$cate{'unprocessed_pseudogene'} = 'pseudogene';
$cate{'nonsense_mediated_decay'} = 'protein_coding';
$cate{'transcribed_processed_pseudogene'} = 'pseudogene';
$cate{'transcribed_unprocessed_pseudogene'} = 'pseudogene';
$cate{'miRNA'} = 'miRNA';
$cate{'lincRNA'} = 'lincRNA';
$cate{'snoRNA'} = 'snoRNA';
$cate{'sense_intronic'} = 'non-coding_other';
$cate{'scRNA_pseudogene'} = 'scRNA';
$cate{'snoRNA_pseudogene'} = 'snoRNA';
$cate{'misc_RNA'} = 'misc_RNA';
$cate{'antisense'} = 'non-coding_other';
$cate{'IG_V_gene'} = 'protein_coding';
$cate{'Mt_tRNA_pseudogene'} = 'Mt_RNA';
$cate{'non_coding'} = 'non-coding_other';
$cate{'rRNA'} = 'rRNA';
$cate{'retrotransposed'} = 'non-coding_other';
$cate{'rRNA_pseudogene'} = 'rRNA';
$cate{'TR_V_gene'} = 'protein_coding';
$cate{'TR_V_pseudogene'} = 'pseudogene';
$cate{'tRNA_pseudogene'} = 'tRNA';
$cate{'TR_J_gene'} = 'protein_coding';
$cate{'unitary_pseudogene'} = 'pseudogene';
$cate{'snRNA_pseudogene'} = 'snRNA';
$cate{'polymorphic_pseudogene'} = 'pseudogene';
$cate{'3prime_overlapping_ncrna'} = 'non-coding_other';
$cate{'TR_C_gene'} = 'protein_coding';
$cate{'miRNA_pseudogene'} = 'miRNA';
$cate{'ncrna_host'} = 'non-coding_other';
$cate{'Mt_tRNA'} = 'Mt_RNA';
$cate{'Mt_rRNA'} = 'Mt_RNA';
$cate{'misc_RNA_pseudogene'} = 'misc_RNA';
$cate{'sense_overlapping'} = 'non-coding_other';
$cate{'IG_V_pseudogene'} = 'pseudogene';
$cate{'IG_J_gene'} = 'protein_coding';
$cate{'IG_C_gene'} = 'protein_coding';
$cate{'IG_C_pseudogene'} = 'pseudogene';
$cate{'TR_J_pseudogene'} = 'pseudogene';
$cate{'TR_D_gene'} = 'protein_coding';
$cate{'IG_D_gene'} = 'protein_coding';
$cate{'IG_J_pseudogene'} = 'pseudogene';


my $expr = shift;
my $anno = shift;

my %expr;
open IN, "$expr";
while ( <IN> ) {
   chomp;
   my @cols = split /\t/;
   my $id = $cols[3];
   my $count = $cols[6];
   $expr{$id} = $count;
}
close IN;

open IN, "$anno";
while ( <IN> ){
  chomp;
  next if /^#/;

  my @cols = split /\t/;
  $cols[8] =~ /^ID=(.+?);/;
  my $id = $1;
  next if $id =~ /:/;

  if (exists $expr{$id}){
    if (exists($cate{$cols[1]})){
       print "$id\t$expr{$id}\t$cate{$cols[1]}\n";
    }
  }
}
close IN;
