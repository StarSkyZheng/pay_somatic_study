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
use v5.10;
use zzIO;


#my $GCstats_file="/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/02.map/GCstat_merge.txt";
#my $GCstats_file='/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/02.map/fix.GCstat_merge.txt';
#my $GCstats_file='/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/02.map/3.realign.unque/pay_unmerged/merge.fix.GCstat';
my $GCstats_file='/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/other/9.all/3.realign.unique/merge.merge.fix.GCstat';

my $VAF_fold_threshold=2;
my $AD_bg_fold_threadshold=10;

my $pm_min = undef;
my $pm_bg_max = 1;
my($in_vcf, $in_indel, $opt_help , $type);
my $thread = 1;
my $AD_fold_thread = 2;
GetOptions (
	'help|h!' => \$opt_help,
	'vcf=s' => \$in_vcf,
	'type=s' => \$type,
	'AD_fold_thread=s' => \$AD_fold_thread,
	't=i' => \$thread,
	'pm_min=i' => \$pm_min,
	'VAF_fold_threshold=s' => \$VAF_fold_threshold,
	'AD_bg_fold_threadshold=s' => \$AD_bg_fold_threadshold,
	'G=s' => \$GCstats_file,
	'pm_bg_max' => \$pm_bg_max,
);

my $doRC = $pm_min // undef;
my $plus_min=$pm_min;
my $mins_min=$pm_min;

if (! $type && $in_vcf) {
	$type='snp' if $in_vcf=~/\.SNP\./;
	$type='indel' if $in_vcf=~/\.INDEL\./;
	print STDERR "now type is $type ! \n";
}


MCE::Map::init {
		max_workers => $thread, chunk_size => 'auto'
};

die "usage: perl xxx.pl 
-vcf raw_snp.vcf.gz 
[-t 1] 
[-pm_min 0]
[-AD_bg_fold_threadshold 10]
[-AD_fold_thread 3]
[-VAF_fold_threshold 2]
-type (indel|snp) 
(| bgzip -c ) > out.vcf\n" if ($opt_help || !-e $in_vcf || ! $type );


print STDERR "\n** Now in node: ";
print STDERR`hostname`;


foreach ($GCstats_file , $in_vcf ){
	die "file not exist: $_\n" unless -e $_;
}

my %dp=&readdp($GCstats_file);
print STDERR "GC: ".join(' ',sort keys %dp)."\n";


my $INVCF = open_in_fh($in_vcf);
my @idline;
print STDERR "Processing snp_vcf: $in_vcf\n";
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

close $INVCF;
print STDERR "/home/share/software/samtools/htslib/tabix -p vcf xxxOUT.vcf.gz\n";
# system("/home/share/software/samtools/htslib/tabix -p vcf $outvcf");



THEEND:

print STDERR "\n** ALL Done\n";
print STDERR`date`;

exit 0;



sub readdp{
	my ($in)=@_;
	my %r;
	my $D;
	open ($D,"$in")||die"no $in\n";
	while (<$D>) {
		chomp;
		next if /^#/;
		my @a=split(/\s+/,$_);
		$r{$a[0]}=$a[5];
		# $r{$a[0]}=$a[6];
	}
	close $D;
	return %r;
}





sub flt{
	my $in = $_;
	chomp $in;
	my @a=split(/\t/,$in);
	#return undef unless $a[1] == 6354;
	##  VCF:
	##   0    1   2  3    4   5      6     7    8      9...
	## CHROM POS ID REF  ALT QUAL FILTER INFO FORMAT   ...
	#my %info= $a[7] eq '.' ? () : split(/[=;]/,$a[7]);


	my @all_alleles=split(/,/,$a[4]);
	unshift @all_alleles, $a[3];
	my $allele_count=scalar(@all_alleles);

	return undef if ( length($a[3])>1 );
	foreach ( split(/,/,$a[4]) ){
		if ( $type eq 'snp' and length($_)>1) {
			#print STDERR join("\t",@a)."\n" ;
			return undef;
		}
	}


	my $all=0;
	my ($RC_num, $GT_num, $AD_num);
	my @info=split(/:/,$a[8]);
	for (my $ii=0;$ii<@info;$ii++){
		$GT_num = $ii and next if ( (! defined $GT_num) && $info[$ii] eq 'GT');
		$AD_num = $ii and next if ( (! defined $RC_num) && $info[$ii] eq 'AD');
		$RC_num = $ii and next if ( (! defined $RC_num) && $info[$ii] eq 'RC');
	}
	my $miss=0;
	foreach my $i (9..@a-1){
		my $iddp= $dp{$idline[$i]} or die "$idline[$i] not in GCstat file";
		my @b = split(':',$a[$i]);
		my $result = &process_per_id(\@b, $GT_num, $RC_num, $AD_num, $iddp, \@all_alleles);
		if (defined $result) {
			$b[0] = './.';
			push @b, $result;
		} else {
			push @b, '.';
			$all++;
		}
		$a[$i] = join ':', @b;
	}
	#say STDERR 'no_left'. join " ", @a and return undef if $all==0;
	my (%ratio, %ratio_fix);

	$a[6]='PASS';
	#$a[7] = $a[7] eq '.' ? '.' : "$a[7];";
	$a[8].=':ERR';
	return join("\t",@a);
}


sub process_per_id($$$$$) {
	my ($b,$GT_num, $RC_num, $AD_num, $iddp, $all_alleles) = @_;
	my $allele_count = scalar @$all_alleles;
	if ($$b[$GT_num]=~/^\.\/\./){
		#return undef; ###########!!!!!!!!!!!!!!!!!!!!!!!###########!!!!!!!!!
		#$miss++;
		return undef;
	}
	$$b[$GT_num] =~ /^(\d)\/(\d)/ or die;
	my @nalleles=($1,$2);

	#$alt{$nalleles[$_]}++ foreach 0..$#nalleles;
	#$all = $all + 2;

	my (@RC, @AD);
	if (defined $doRC && !defined $RC_num) {
		die "no RC found: @$b . " 
	} elsif (defined $doRC) {
		@RC=split(/,/,$$b[$RC_num]);
	} else {}

	@AD=split(/,/,$$b[$AD_num]);
	
	my @uniq_nalleles = &uniq(@nalleles);
	my $dp_now=0;
	foreach my $now ( @uniq_nalleles ) {
		#my $now = $nalleles[$ii];
		if (
			defined $doRC &&
			($RC[$now*2] < $plus_min || $RC[$now*2+1] < $mins_min)
			) {
			return "LowPM";
		}
		$dp_now += $AD[$now];
	}
	if ($dp_now > $iddp * $AD_fold_thread || $dp_now < $iddp / $AD_fold_thread ) {
		return "dp3-$iddp-$dp_now";
	}
	my $bg_AD=0;
	my @bg_RC=(0,0);
	my $AD_add=0;
	foreach my $nalle_now (0..$allele_count-1) {
		die "!!!" if defined $doRC && ! exists $RC[$nalle_now*2] ;
		$AD_add += $AD[$nalle_now];
		next if $nalle_now ~~ @uniq_nalleles;
		if ( defined $doRC ) {
			$bg_RC[0]+=$RC[$nalle_now*2];
			$bg_RC[1]+=$RC[$nalle_now*2+1];
		}
		$bg_AD += $AD[$nalle_now];
	}
	if ( defined $doRC &&
		( $bg_RC[0] >= $pm_bg_max && $bg_RC[1] >= $pm_bg_max )
		) {
		return "HIBG_RC";
	}
	return "HiBG_AD" if ($bg_AD >= $AD_add / $AD_bg_fold_threadshold );
	if (scalar @uniq_nalleles == 1) {
		my $uniq_nallele = $uniq_nalleles[0];
		my $VAF =  $AD_add==0 ? 0 : $AD[$uniq_nallele] / $AD_add;
		if ( $VAF <= 1 / $VAF_fold_threshold || $VAF >= 1*$VAF_fold_threshold) {
			return "VAF1";
		}
	} elsif (scalar @uniq_nalleles == 2) {
		my $VAF = $AD[$uniq_nalleles[0]] / ($AD[$uniq_nalleles[0]] + $AD[$uniq_nalleles[1]]);
		if ( $VAF < 1/($VAF_fold_threshold+1) || $VAF > $VAF_fold_threshold/($VAF_fold_threshold+1) ) {
			return "VAF2";
		}
	} else { die; }
	return undef;
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}




