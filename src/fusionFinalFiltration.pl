use strict;
use Data::Dumper;

my $fusionReport = shift;

open FF, "$fusionReport" || die "could not open $fusionReport\n";

while ( <FF> ) {

    chomp;
    my ($fusion1, $fusion2, $bp1, $bp2, $dir, $rep, $sc, $type, $strand, $cov1, $cov2, $cov3, $cov4, $cov5, $transcripts) = split /\t/;

    my $gene1 = '';
    my $gene2 = '';
    if ($fusion1 =~ /^\S+\((.+?)\)$/){
      $gene1 = $1;
    } else {
      $gene1 = $fusion1;
    }

    if ($fusion2 =~ /\S+\((.+?)\)/){
      $gene2 = $1;
    } else {
      $gene2 = $fusion2;
    }

    next if ($rep =~ /R/ and $sc =~ /C/);                     #skip both rep and sc
    next if ($fusion1 =~ /KCNMB2/ and $fusion2 =~ /KCNMB2/);  #skip KCNMB2 tmp error
    next if ($fusion1 =~ /TMPRSS3/ or $fusion2 =~ /TMPRSS3/); #skip TMPRSS3 tmp error
    next if ($cov5 < 1 or $cov3 < 3);                         #skip low spanning reads
    next if ($cov5 < 2 and $sc eq 'CC');                      #skip low spanning reads
    next if ($cov4/$cov5 >= 10 and $cov5 <= 3);               #high duplication
    next if ($cov4/$cov5 >= 20 and $cov5 < 10);               #high duplication
    next if ((($gene2 ne 'IGR' and $gene1 =~ /^$gene2/) or ($gene1 ne 'IGR' and $gene2 =~ /^$gene1/)) and $sc eq 'CC');  #remaining sc
    next if (($gene1 =~ /^IGK[JV]/ or $gene2 =~ /^IGK[JV]/) and $sc eq 'CC');                                            #IGG sc
    next if (($gene1 eq 'IGR' and $gene2 =~ /^RP\d+\-\d+/) or ($gene2 eq 'IGR' and $gene1 =~ /^RP\d+\-\d+/));            #lnRNA IGR junk

    print "$_\n";

}

close FF;
exit 0;
