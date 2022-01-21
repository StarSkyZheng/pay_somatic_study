#perl

use warnings;
use strict;
use v5.10;

my $in=shift;
my $out=shift;
my $I;

$in ? open($I,"<$in") : ($I = *STDIN);
open(STDOUT,">$out") if $out;

# need
# chrom   position        ref_base        depth   q20_depth       ....
#  0         1               2              3         4             5

# base : reads : strands : avg_qual : map_qual : plus_reads : minus_reads
#   0      1        2         3          4           5             6


# source
# chr	position	reference_base	depth	...
# 0       1               2          3       4
 
# base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end  
#   0  1     2                   3               4                      5                6               7                    8                              9                          10                      11                                    12                13                         


while (<$I>) {
	chomp;
	my @a = split(/\s+/,$_);
	my @ret = (@a[0..3],$a[3]);
	foreach my $i (4..$#a) {
		my @b = split(/:/,$a[$i]);
		next if $b[1]==0;
		next if $b[0] eq 'N';
		next if $b[0] eq '=';
		my @c = ( $b[0], $b[1],'1' ,$b[3], $b[2], $b[5], $b[6] );
		push @ret, join( ':', @c);
	}
	say STDOUT join "\t", @ret;
}
