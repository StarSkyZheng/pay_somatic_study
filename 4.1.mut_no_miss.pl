#!/usr/bin/env perl
#===============================================================================
#
#         FILE: 04-05.MISS_SNP_ALL.pl
#
#        USAGE: ./04-05.MISS_SNP_ALL.pl
#
#  DESCRIPTION: 根据正链和反链 过滤位点，低深度位点以.代替
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

# java -jar /home/share/software/gatk/GenomeAnalysisTK-3.8-0.jar -T SelectVariants -nt 8 -R xx.fa -V xx.join.raw.vcf.gz -selectType SNP -o xx.join.raw.SNP.vcf.gz

# java -jar /home/share/software/gatk/GenomeAnalysisTK-3.8-0.jar -T SelectVariants -nt 8 -R xx.fa -V xx.join.raw.vcf.gz -selectType SNP -o xx.join.raw.SNP.vcf.gz


#use v5.18;
use strict;
use warnings;
use utf8;
use FileHandle;
use Getopt::Long;
use MCE::Map;
use v5.24;


my($in_vcf, $opt_help , $type);
my $thread = 1;
GetOptions (
	'help|h!' => \$opt_help,
	'vcf=s' => \$in_vcf,
	't=i' => \$thread,
);

MCE::Map::init {
		max_workers => $thread, chunk_size => 'auto'
};

die "usage: perl xxx.pl 
-vcf raw_snp.vcf.gz 
[-t 1] 
(| bgzip -c ) > out.vcf\n" if ($opt_help || !-e $in_vcf );


print STDERR "\n** Now in node: ";
print STDERR`hostname`;


foreach ( $in_vcf ){
	die "file not exist: $_\n" unless -e $_;
}


my $INVCF;
open ($INVCF,"zcat $in_vcf|") or die "Cant open $in_vcf, $!" ;
my @idline;
print STDERR "Processing snp_vcf: $in_vcf\n";
while ( <$INVCF> ){
	if ( /^##/ ){
		#print $_;
		next;
	}
	if (/^#/){
		chomp;
		#print $_."\n";
		if (/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+\w+/){
			@idline=split(/\s+/,$_);
			last;
		}else{
			print STDERR;
			last;
		}
	}
}

die "NO idline: \n" if ! @idline;
print STDERR "idline : @idline\n";


if ( $thread == 1) {
	while (<$INVCF>) {
		chomp;
		my $result;
		$result=&flt($_) unless /^#/;
		say $result if $result;
	}
} else {
	my @ret = grep {defined $_} mce_map_f {&flt} $INVCF;
	MCE::Map::finish;
	say STDERR 'run ok';
	say foreach @ret;
}

sub flt {
	my $in = $_;
	chomp $in;
	my @a=split(/\t/,$in);
    my $maf_nids = &cal_maf(\@a);
    return undef unless $maf_nids;
    return $in if $maf_nids eq 'all_same';
    if (scalar @$maf_nids ==1) {
        my @s = split(':', $a[$$maf_nids[0]]);
		return undef if $s[1] eq '.';
        return undef if $s[1] eq '!';
        return undef unless $s[1] == 3;
		say STDERR @s;
    }
    return $in;
}

sub cal_maf {
    my $s_line = shift or die;
    my %allecount;
    my %alle2id;
    foreach my $i (9..@$s_line-1) {
        my @a = split(':', $$s_line[$i]);
        $allecount{$a[0]} ++;
        push $alle2id{$a[0]}->@*, $i;
    }
    my ($most, @rares) = sort {$allecount{$b}<=>$allecount{$a}} keys %allecount;
    die if $allecount{$most}<0.6;
    my @ids;
    return 'all_same' unless @rares;
    push @ids, $alle2id{$_}->@* foreach @rares;
    #die join ' ',@ids, @$s_line if scalar @ids > 1;
    return undef unless @ids;
    return \@ids;
}
