use File::Basename;
use warnings;
use strict;
use v5.10;
use zzReadRefdict;
use Parallel::ForkManager;


my $out_dir=shift or &help();
my @bams=&read_list(shift);
my @type=@ARGV;

my $ref='/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/genome-v2.fa';
my %chr2len = zzReadRefdict($ref);
my $ref_split_chr_dir = '/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/Split_chr';

my %bams=&process_bam_list(@bams);
my $bams;
$bams.="$_ " foreach @bams;

my $cnv_dir="$out_dir/1.CNVnator";
my $pindel_dir="$out_dir/2.pindel";
my $breakdancer_dir="$out_dir/3.breakdancer";
my $hugeseq_bin_dir = '/home/share/users/zhengzeyu2018/software/HugeSeq/bin';

(-e $_ || mkdir $_ || die) foreach ($cnv_dir,$pindel_dir, $breakdancer_dir);

my $MAX_PROCESSES = 20;
my $cnvnator_bin_dir='/home/share/users/zhengzeyu2018/software/SV_tector/CNVnator';
my $pindel_bin_dir='/home/share/users/zhengzeyu2018/software/SV_tector/pindel';
my $BreakDancer_bin_dir='/home/share/users/zhengzeyu2018/software/SV_tector/breakdancer';

&cnvnator( "$cnv_dir/CNV" ) if 'C' ~~ @type;
&cnvnator_split_chr( "$cnv_dir/CNV" ) if 'C_S' ~~ @type;
&pindel ( "$pindel_dir/Pindel" ) if 'P'~~ @type;
&breakDancer ( "$breakdancer_dir/BreakDancer" ) if 'B' ~~ @type;

#merge gff;

my $t = qq+ sh /home/share/users/zhengzeyu2018/software/HugeSeq/bin/merge_gff.sh out in1 in2.. +;
say $t;
exit;


sub help() {
	die "usage: perl xx.pl out_dir bam_list_file [type:C/B/P]\n";
}

sub process_bam_list() {
	my %hash;
	my @bams_t=@_;
	foreach my $bam (@bams_t) {
		my $base=basename($bam);
		$base=~/^([^.]+)/ or die $base;
		$hash{$1} = $bam;
	}
	return %hash;
}


sub cnvnator() {
	my $out_prefix=shift;
	my @cmds;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.root -unique -d $ref_split_chr_dir -tree $bams +;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.root -d $ref_split_chr_dir -his 100 +;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.root -d $ref_split_chr_dir -stat 100 +;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.root -d $ref_split_chr_dir -partition 100 +;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.root -d $ref_split_chr_dir -call 100 > $out_prefix.txt+;
	push @cmds, qq+$cnvnator_bin_dir/cnvnator2VCF.pl $out_prefix.txt | bgzip -c > $out_prefix.result.vcf.gz +;
	@cmds=() if -s "$out_prefix.root";
	foreach my $temp (@cmds) {
		system $temp;
	}

	system "sh $hugeseq_bin_dir/cnvnator.sh $out_prefix.gff $bams";
}

sub cnvnator_split_chr() {
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	my $out_prefix=shift;
	my @cmds;
	foreach my $chr (sort keys %chr2len) {
		my $pid = $pm->start and next;
		say $chr;
		#die `pwd`;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.$chr.root -unique -d $ref_split_chr_dir -tree $bams -chrom $chr +;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.$chr.root -d $ref_split_chr_dir -his 100 +;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.$chr.root -d $ref_split_chr_dir -stat 100 +;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.$chr.root -d $ref_split_chr_dir -partition 100 +;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator -root $out_prefix.$chr.root -d $ref_split_chr_dir -call 100 > $out_prefix.$chr.txt+;
		push @cmds, qq+$cnvnator_bin_dir/cnvnator2VCF.pl $out_prefix.$chr.txt | bgzip -c > $out_prefix.$chr.result.vcf.gz +;
		@cmds=() if -s "$out_prefix.$chr.root";
		foreach my $temp (@cmds) {
			system "echo '$temp' >> $out_prefix.$chr.sh";
		}
		#system "sh $hugeseq_bin_dir/cnvnator.sh $out_prefix.$chr.gff $bams";
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub pindel() {
	my $out_prefix=shift;
	my $config_txt="$out_prefix.config";
	open (my $O, "> $config_txt");
	foreach my $id (sort keys %bams) {
		print $O "$bams{$id}   350   $id\n";
	}
	close $O;

	my @cmds;
	push @cmds, qq+ $pindel_bin_dir/pindel -f $ref -i $config_txt -c ALL -o $out_prefix +;
	push @cmds, qq+ $pindel_bin_dir/pindel2vcf -r $ref -R ref -d 20180000 -v $out_prefix.result.vcf.gz -P $out_prefix +;
	@cmds=() if -e "${out_prefix}_SI";
	foreach my $temp (@cmds) {
		system $temp;
	}

	system "sh $hugeseq_bin_dir/pindel.sh $out_prefix.gff $bams";
}



sub breakDancer() {
	my $out_prefix=shift;
	my @cmds;
	push @cmds, qq+ $BreakDancer_bin_dir/perl/bam2cfg.pl -g -h $bams > $out_prefix.config.txt +;
	push @cmds, qq+ $BreakDancer_bin_dir/breakdancer-max $out_prefix.config.txt > $out_prefix.config.result +;
	@cmds=() if -s "$out_prefix.config.result";
	foreach my $temp (@cmds) {
		system $temp;
	}

	system "sh $hugeseq_bin_dir/breakdancer.sh $out_prefix.config.result.gff $bams";
}


sub read_list() {
	my $f=shift;
	my @t;
	open (my $F,"<$f") or die;
	while (<$F>) {
		chomp;
		push @t, $_ if $_;
	}
	return @t;
}

