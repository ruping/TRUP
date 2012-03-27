#!/usr/bin/perl
use strict;
use Getopt::Long;

my $pileup;
my $where;
GetOptions(
             "pileup|p=i"=>\$pileup,
	     "where|w=s"=>\$where,
          );

#######
#TODO: how many reads in exon reagion?
#######

#set flag to Tophat mapping read
my $now_mapread;
my %mapread_mRNA;
my %mapread_exon;
my %mapread_CDS;
my %mapread_all;
my $flag = 0;

while ( <> ){
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $tag) = split /\t/;
  $chr =~ s/chr//;
  $chr = 'MT' if ($chr eq 'M');

  if ($pileup){
    $tag =~ /^Pileup=(\d+)$/;
    my $p = $1;
    if (($source eq 'TopHat') and ($p >= $pileup)){
      $now_mapread = $_;
      $mapread_all{$chr}++;
      $flag = 1;
    }
    if (($source eq 'TopHat') and ($p < $pileup)){
      $flag = 0;
    }
    if ($flag == 1 and $type eq $where){
       $mapread_exon{$chr}{$now_mapread}++;
    }
  }
  
  if (!$pileup){
   if ($source eq 'TopHat'){
     $now_mapread = $_;
     $mapread_all{$chr}++;
   }
   if ($type eq 'mRNA'){
     $mapread_mRNA{$chr}{$now_mapread}++;
   }
   if ($type eq 'exon'){
     $mapread_exon{$chr}{$now_mapread}++;
   }
   if ($type eq 'CDS'){
     $mapread_CDS{$chr}{$now_mapread}++;
   }
  }

}

if (!$pileup){
my $mappedread_mRNA;
my $mappedread_exon;
my $mappedread_CDS;
foreach my $chr (sort {$a cmp $b} keys %mapread_mRNA){
     my $mapread_mRNA = scalar (keys %{$mapread_mRNA{$chr}});
     $mappedread_mRNA += $mapread_mRNA;
     my $ratio1 = $mapread_mRNA/$mapread_all{$chr};

     my $mapread_exon = scalar (keys %{$mapread_exon{$chr}});
     $mappedread_exon += $mapread_exon;
     my $ratio2 = $mapread_exon/$mapread_all{$chr};

     my $mapread_CDS = scalar (keys %{$mapread_CDS{$chr}});
     $mappedread_CDS += $mapread_CDS;
     my $ratio3 = $mapread_CDS/$mapread_all{$chr};

     print "$chr\t$ratio1\t$ratio2\t$ratio3\n";
}


my $mappedread_total;
foreach my $chr (sort {$a cmp $b} keys %mapread_all){
   $mappedread_total += $mapread_all{$chr};
}

my $ratio1 = $mappedread_mRNA/$mappedread_total;
my $ratio2 = $mappedread_exon/$mappedread_total;
my $ratio3 = $mappedread_CDS/$mappedread_total;
print "all\t$ratio1\t$ratio2\t$ratio3\n";
}

if ($pileup){
my $mappedread_exon;
foreach my $chr (sort {$a cmp $b} keys %mapread_exon){

     my $mapread_exon = scalar (keys %{$mapread_exon{$chr}});
     $mappedread_exon += $mapread_exon;
     my $ratio = $mapread_exon/$mapread_all{$chr};

    # print "$chr       $mapread_exon/$mapread_all{$chr}=$ratio\n";
    # print "$chr        $ratio\n";
}

my $mappedread_total;
foreach my $chr (sort {$a cmp $b} keys %mapread_all){
   $mappedread_total += $mapread_all{$chr};
}

my $ratio = $mappedread_exon/$mappedread_total;
print "$pileup      $ratio\n";
}


exit 22;
