use strict;
use Data::Dumper;

my @chrs = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM);
my %chrs;
foreach my $id (@chrs){
  $chrs{$id} = '';
}

my $gtf = shift;

if ($gtf !~ /\.gtf$/ or $gtf eq ''){
  print "Usage: perl exon_union_gtf2bed.pl A_GTF_file >output_bed_file\n";
  print "written by Ruping Sun, rs3412\@c2b2\.columbia\.edu\n\n";
  exit 22;
}

my %ensembl;
my $ens_old = "SRP";
my $strand_old;

open GTF, "$gtf";
while ( <GTF> ){
  chomp;
  next if /^#/;
  my ($chr, $type, $feature, $start, $end, $score, $strand, $phase, $tags) = split /\t/;

  next unless $feature eq 'exon';

  #---------chr name processing-------------------------------------------------
  my $chr_full = 'chr'.$chr;
  $chr = 'chrM' if ( ($chr eq 'chrMT') or ($chr_full eq 'chrMT') );

  if ( exists $chrs{$chr_full} ){
    $chr = $chr_full;
  }
  elsif ( exists $chrs{$chr} ){
    #do(nothing);
  }
  else {
    next;  #skip
  }
  #------------------------------------------------------------------------------

  #---------tag processing-------------------------------------------------------
  my @tags = split(";", $tags);
  my $gene_id;
  my $transcript_id;
  my $exon_number;
  my $gene_name;
  my $gene_biotype;
  foreach my $tag (@tags){
    if ($tag =~ /gene\_id\s+\"(.+?)\"/){
      $gene_id = $1;
    }
    elsif ($tag =~ /transcript\_id\s+\"(.+?)\"/){
      $transcript_id = $1;
    }
    elsif ($tag =~ /exon\_number\s+\"(.+?)\"/){
      $exon_number = $1;
    }
    elsif ($tag =~ /gene\_name\s+\"(.+?)\"/){
      $gene_name = $1;
    }
    elsif ($tag =~ /gene\_biotype\s+\"(.+?)\"/){
      $gene_biotype = $1;
    }
  }

  #print "$chr\t$start\t$end\t\t$strand\t$gene_id\t$gene_name\t$gene_biotype\t$transcript_id\t$exon_number\n";
  #------------------------------------------------------------------------------

  if ($ens_old ne "SRP" and $gene_id ne $ens_old) {
      #have to merge
      merge(\%ensembl, $ens_old, $strand_old);
      #print Dumper(\%ensembl);
      %ensembl = ();
  }

  $ensembl{$chr}{$start}{$end}  =  $transcript_id;

  $ens_old = $gene_id;
  $strand_old = $strand;

}
close GTF;

merge(\%ensembl, $ens_old, $strand_old);
exit;

sub merge {

   my ($exons, $ens_id, $ens_strand) = @_;

   my @chromosomes = keys %{$exons};
   if ( scalar(@chromosomes) != 1 ){
     print STDERR "why a gene is located in different chromosomes?\n";
       for (0..$#chromosomes) {
           print STDERR "$chromosomes[$_]\t";
       }
       print STDERR "\n";
       exit;
   }
   #ok, now $chromosome[0] is the chr
   my $chromosome = $chromosomes[0];

   my $last_start = 0;
   my $last_end = 0;
   my %merged = ();

   foreach my $current_start (sort {$a <=> $b} keys %{$exons->{$chromosome}} ) {

       my @current_ends = sort {$b <=> $a} keys %{${$exons->{$chromosome}}{$current_start}};  #may be multiple ends, pick the largest one
       my $current_end = $current_ends[0];
       #print "$current_start\t$current_end\n";

       if ( $current_start <= $last_end ){ #overlapping
          if ($current_end > $last_end) {
             $merged{$last_start}{'end'} = $current_end;
          }
          else{
             next;
          }
       }
       else { #new_starts
          $merged{$current_start}{'end'} = $current_end;
          $last_start = $current_start;
       }

       $last_end = $current_end;
   } #merged

   #print out
   my $thickstart = -1;
   my $thickend = -1;
   my $blockCount = 0;
   my @blockSizes;
   my @blockStarts;

   foreach my $merged_start (sort {$a <=> $b} keys %merged){

     if ($thickstart == -1){
        $thickstart = $merged_start - 1;  #bed format start is zero based
     }

     $blockCount++;

     my $merged_end = $merged{$merged_start}{'end'};
     my $current_size = $merged_end - $merged_start + 1;
     my $current_start = $merged_start - $thickstart - 1;
     push(@blockSizes, $current_size);
     push(@blockStarts, $current_start);

     $thickend = $merged_end if ($merged_end >= $thickend);
   }

   my $blockSizes = join(",", @blockSizes);
   my $blockStarts = join(",", @blockStarts);
   printf("%s\n", join("\t", $chromosome, $thickstart, $thickend, $ens_id, 0, $ens_strand, $thickstart, $thickend, 0, $blockCount, $blockSizes, $blockStarts));
}
