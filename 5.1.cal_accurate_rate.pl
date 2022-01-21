use warnings;
use strict;
use zzIO;
use v5.10;

my $in = shift or die 'no 5.vcf.cal.pl_result input';

my $I = open_in_fh($in);

my $header=<$I>;
chomp $header;
my @header = split("\t", $header);
while(<$I>) {
    chomp;
    my @a = split("\t");
    my @say = @a[0..3];
    foreach my $i (1..5) {
        my $id = $header[5 + ($i-1)*3];
        my $mut=0;
        my $nomiss=0;
        foreach my $ii (0..2) {
            my $nid = 5 + ($i-1)*3 + $ii; 
            my $stat = $a[$nid];
            $mut++ if $stat eq '1';
            $nomiss++ if $stat ne 'NA';
        }
        
        if ($mut) {
            if ($mut == $nomiss) {
                push @say, "$id:$mut:OK";
                say join "\t", @a[0..3] , "$id:$mut:OK"; 
            } else {
                push @say, "$id:$mut:" . ($nomiss - $mut);
                say join "\t", @a[0..3] , "$id:$mut:" . ($nomiss - $mut); 
            }
        }
        #say join "\t", @say if @say > 4;
    }
}
