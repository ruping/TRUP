use strict;

open IN, shift;

my %gene;

while ( <IN> ) {
   next if /^#/;
   chomp;
   my @cols = split /\t/;
   my $tags = $cols[$#cols];

   my @tags = split(/;/, $tags);
   my $gene_id;
   my $gene_type;
   my $gene_name;

   foreach my $tag (@tags){
     if ($tag =~ /gene_id\s\"(.+?)\"/){
        $gene_id = $1;
     }
     elsif ($tag =~ /gene_type\s\"(.+?)\"/){
        $gene_type = $1;
     }
     elsif ($tag =~ /gene_name\s\"(.+?)\"/){
        $gene_name = $1;
     }
   }

   unless ( exists $gene{$gene_id} ){
       print "$gene_id\t$gene_type\t$gene_name\n";
       $gene{$gene_id} = '';
   }

}

close IN;
