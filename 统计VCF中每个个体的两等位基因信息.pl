#!/usr/bin/env perl
#===============================================================================
#
#         FILE: 04-05.MISS_SNP_ALL.pl
#
#        USAGE: ./04-05.MISS_SNP_ALL.pl
#
#  DESCRIPTION:  剔除所有个体都是 0/1 的位点，并统计0/0、0/1、1/1的数目
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (Lanzhou University), zhengzy2014@lzu.edu.cn
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 10/24/2017 07:23:37 PM
#     REVISION: ---
#===============================================================================


use strict;
use warnings;
use utf8;
#use FileHandle;
use Getopt::Long;
#use MCE::Loop;
use v5.10;



my $rep_file="/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/v2.repeatmask.merge.gff";
my $ref='/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/genome-v2.fa';
my $GCstats_file="/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/02.map/fix.GCstat_merge.txt";


#my $rep_file='/home/share/users/zhengzeyu2018/project/Ostrya_genome/Ore.final.gff3';
#my $ref='/home/share/users/zhengzeyu2018/project/Ostrya_genome/genome/Ore.final.assembly.flt.2k.fa';
#my $GCstats_file='/home/share/users/zhengzeyu2018/project/Corylus_chinensis/GCstat_merge.txt';




my($in_vcf, $in_indel, $opt_help);
my $thread='auto';
GetOptions (
	'help|h!' => \$opt_help,
	'vcf=s' => \$in_vcf,
	'p=i' => \$thread,
);




print STDERR "Processing snp_vcf: $in_vcf, print to STDOUT\n";
open (my $INVCF,"zcat $in_vcf|") or die "Cant open $in_vcf, $!" ;

my @idline;
while ( <$INVCF> ){
	if ( /^##/ ){
		print $_;
		next;
	}
	if (/^#/){
		chomp;
		print $_."\n";
		if (/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+\w+/){
			@idline=split(/\s+/,$_);
			print join("\t",@idline)."\n";
			last;
		}else{
			print STDERR;
			last;
		}
	}
}
die "NO idline: \n" if ! @idline;


while(<$INVCF>) {
	if (/^#/) {
		print;
		next;
	}
	# my ($mce, $slurp_ref, $chunk_id) = @_;
	chomp;
	my $result=&flt($_);
	say( $result ) if $result;
} 

close $INVCF;


exit 0;





sub flt{
	my @a=split(/\s+/,$_[0]);
	#die unless @a;
	##  VCF:
	##   0    1   2  3    4   5      6     7    8      9...
	## CHROM POS ID REF  ALT QUAL FILTER INFO FORMAT   ...

	my %hash;
	my $miss=0;
	for (my $i=9;$i<@a;$i++){
		if ($a[$i]=~/^\.\/\./){
			#return undef; ###########!!!!!!!!!!!!!!!!!!!!!!!###########!!!!!!!!!
			$miss++;
			next;
		}elsif ( $a[$i] =~ /^(\d)\/(\d):/ ) {
			$hash{$1+$2}++;
		}
	}

	foreach my $key (sort keys %hash) {
		return undef if ( $hash{$key} >= $#a-8-$miss ) ;
	}

	$hash{0}=0 if ! exists $hash{0};
	$hash{1}=0 if ! exists $hash{1};
	$hash{2}=0 if ! exists $hash{2};

	$a[7] = "0:$hash{0},1:$hash{1},2:$hash{2}";

	$a[6]='PASS';
	return join("\t",@a);
}





