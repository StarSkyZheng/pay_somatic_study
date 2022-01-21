
use File::Basename;
use warnings;
use strict;
use v5.10;
use MCE::Loop;
use Getopt::Long;

#my $work_dir=shift or die;

my ($bam_list, $pre_snp_vcf, $pre_indel_vcf, $out_dir, $opt_help);

my $thread = 'auto';


GetOptions (
	'help|h!' => \$opt_help,
	'bam_list=s' => \$bam_list,
	'pre_snp_vcf=s' => \$pre_snp_vcf,
	'pre_indel_vcf=s' => \$pre_indel_vcf,
	'out_dir=s' => \$out_dir,
	'p=i' => \$thread,
);

&help() if ($opt_help or ! $bam_list or ! $out_dir);

sub help () {
	die "usage: \$0 -bam_list xx -pre_snp_vcf xx -pre_indel_vcf xx -out_dir xx -p 20\n";
}

MCE::Loop::init { max_workers => $thread,  chunk_size => 1};


my $ref='/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/genome-v2.fa';

(-e $_ || die "??? $_\n") foreach ($bam_list, $ref, "$ref.fai");

#my $work_dir=dirname($pre_flt_vcf);
#my $out_dir=$work_dir;
my @bams=&read_list($bam_list);

my %bams=&process_list(@bams);
my @samples=sort keys %bams;


my $out_readcounts_dir="$out_dir/1.readcounts";

my $BioScripts_dir='/home/share/users/zhengzeyu2018/software/BioScripts/';
my $VarScan='/home/share/users/zhengzeyu2018/software/VarScan/VarScan.v2.4.3.jar';
my $gatk='/home/share/software/gatk/GenomeAnalysisTK-3.8-0.jar';
my $vcftools_src_dir='/home/share/users/zhengzeyu2018/software/vcftools/src';
my $bam_readcount='/home/share/users/zhengzeyu2018/software/bam-readcount/bin/bam-readcount';
#my $bam_readcount2varscan='/home/share/users/zhengzeyu2018/script/99.other/bam-readcount2VarScan-readcount.pl';
my $bam_readcount2varscan = '~/script/03.call_snp/99.accurate/bam-readcount2VarScan-readcount.pl';

#&fix_header($pre_flt_vcf);

goto SNP if $pre_snp_vcf;
goto INDEL if $pre_indel_vcf;
goto EXIT;

SNP:
my $snp_bed="$pre_snp_vcf.snp.bed";
my $out_SNP_vcfgz="$out_dir/"."2.SNPwithDp.".basename($pre_snp_vcf);
my $out_SNP_flt_vcfgz="$out_SNP_vcfgz.c2t3d5m5.vcf";

{
	my $cmd1="zcat $pre_snp_vcf" . q& | perl -ne 'next if(/\\#/); my ($chrom, $pos)=(split /\\s+/)[0,1]; print "$chrom\\t",($pos-1),"\\t$pos\\n";' > & . $snp_bed ;
	#my $cmd1=qq& zcat $pre_snp_vcf | perl $vcftools_src_dir/perl/vcf-annotate & . q& --fill-type | perl -ne 'next if(/\\#/); next unless(/snp/); my ($chrom, $pos)=(split /\\s+/)[0,1]; print "$chrom\\t",($pos-1),"\\t$pos\\n";' > & . $snp_bed ;
	system $cmd1 unless -e $snp_bed;
}
mkdir $out_readcounts_dir unless -e $out_readcounts_dir;
#foreach my $sample (@samples) {
mce_loop {
	my $sample=$_;
	my $bam=$bams{$sample};
	my $out_readcounts="$out_readcounts_dir/$sample.readcounts";
	my $cmd1=qq+ samtools mpileup -d100000 -q20 -f $ref -l $snp_bed $bam | grep -vP "\t0\t" | java -jar $VarScan readcounts - --min-base-qual 20 --min-coverage 1 --output-file $out_readcounts +;
	#my $cmd1=qq+ $bam_readcount -q 20 -b 20 -l $snp_bed -f $ref $bam 2>/dev/null | perl $bam_readcount2varscan > $out_readcounts +;
	say $cmd1;
	#next;
	system $cmd1 unless -e $out_readcounts;
} @samples;

#die 'OK';

my $list="$out_readcounts_dir/readcounts.list";
open (my $L,"> $list");
say $L "$_ $_ $out_readcounts_dir/$_.readcounts" foreach @samples;
close $L;

{
	my $cmd1=qq+ perl $BioScripts_dir/fillVcfDepth.pl --vcf $pre_snp_vcf --list $list --minimum-vcf --update-AD | bgzip -c > $out_SNP_vcfgz +;
	say $cmd1;
	#&fix_header($out_SNP_vcfgz);
	system $cmd1 unless -s $out_SNP_vcfgz and system("gzip -t $out_SNP_vcfgz")==0;
	my $cmd2=qq+ tabix -p vcf $out_SNP_vcfgz +;
	my $cmd3=qq+perl $BioScripts_dir/detect_mutations.pl -v $out_SNP_vcfgz --max-cmp-depth 2 --max-cmp-total 3 --max-cmp-miss 5 --min-supp-depth 5 --min-supp-plus 1 --min-supp-minus 1 --mask-only LowDepth | perl $vcftools_src_dir/perl/vcf-annotate -f c=3,150 --fill-type > $out_SNP_flt_vcfgz +;
	#system $cmd3 unless -e $out_SNP_flt_vcfgz;
}

goto INDEL if $pre_indel_vcf;
goto EXIT;
#exit;

INDEL:
my $indel_bed="$pre_indel_vcf.indel.bed";
my $out_INDEL_vcfgz="$out_dir/"."3.INDELwithDp.".basename($pre_indel_vcf);
my $out_INDEL_flt_vcfgz="$out_INDEL_vcfgz.c2t3d5m5.vcf";
{
	my $cmd2=" bgzip -dc $pre_indel_vcf | perl $vcftools_src_dir/perl/vcf-annotate " . q& --fill-type | perl -ne 'next if(/\\#/); next unless(/ins/ || /del/); my ($chrom, $pos, $ref)=(split /\\s+/)[0,1,3]; my $start = $pos-1; my $end = $pos+length($ref); print "$chrom\\t$start\\t$end\\n";' | uniq | bedtools slop -i - -b 10 -g $ref.fai & . qq& $ref.fai > $indel_bed &;
	system $cmd2 unless -e $indel_bed;
}


{
	my $temp;
	$temp.=" -I $bams{$_}" foreach @samples;
	my $cmd1=qq+ java -jar $gatk -R $ref -T HaplotypeCaller -nct $thread -stand_call_conf 30.0 -L $indel_bed -o $out_INDEL_vcfgz $temp 2>&1 | tee $out_INDEL_vcfgz.log +;
	system $cmd1 unless -e $out_INDEL_vcfgz;
	#tabix -p vcf $out_INDEL_vcfgz
	my $cmd2=qq+ perl $BioScripts_dir/vcf_process.pl --vcf $out_INDEL_vcfgz + . q+ --quality 30 --var-type indel | awk '/\#/ || ($9 ~ /AD/)' + . " | perl $BioScripts_dir/detect_mutations.pl -v - --max-cmp-total 3 --max-cmp-depth 2 --min-supp-depth 5 --max-cmp-miss 5 --mask-only LowDepth | perl $vcftools_src_dir/perl/vcf-annotate -f c=2,150 --fill-type >  $out_INDEL_flt_vcfgz ";
	system $cmd2 unless -e $out_INDEL_flt_vcfgz;
}

# -g samples.group.txt

#merge_SNP_INDEL:
#{
#	my $cmd1=qq+ $BioScripts_dir/vcf_process.pl --vcf caller1.mut.vcf --secondary-vcf caller2.mut.vcf --combine-rows 0 1 --compare-rows 2 3 4 --primary-tag CALLER1 --secondary-tag CALLER2  --intersect-tag CALLED_BOTH > mut.combined.vcf+;
#}

EXIT:
print STDERR "\nDONE\n";
exit;

sub process_list() {
	my @t=@_;
	my %t;
	foreach my $bam (@t) {
		my $base=basename $bam;
		$base =~ /^([^.]+)/ or die;
		my $sample=$1;
		$t{$1} = $bam;
	}
	return %t;
}


sub read_list() {
	my $file=shift;
	my @t;
	open (my $F,"<$file") or die;
	while (<$F>) {
		chomp;
		push @t,$_ if $_;
	}
	return @t;
}

sub fix_header() {
	my $file=shift;
	return if ( system ("zcat $file | head -3000 | grep '^IDline: #CHROM'") != 0 );
	my $file_backup="$file.backup";
	die if -e $file_backup;
	system "mv $file $file_backup";
	system qq+zcat $file_backup | sed "s/^IDline: #CHROM/#CHROM/" | bgzip -c > $file +;
}

