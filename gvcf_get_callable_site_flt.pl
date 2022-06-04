
use warnings;
use strict;
use zzReadRefdict;
use File::Basename;
use v5.24;
use MCE::Loop;
use zzIO;
use Data::Dumper;


my ($inlist, $gc_file, $outfile) = @ARGV;
$outfile //= '-';
#my @nodes = qw/node1 node3 node5 node6/;
my $thread = 60;
my @nodes;
my $miss_ratio=0.1;
my $dp_fold=3;


my @gvcf_dirs = (
    '/data/01/user101/project/P_pyramidalis/a1.Ptr/4.gvcf',
);

#my $gc_file = '/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/other/9.all/3.realign.unique/merge.merge.fix.GCstat';
$gc_file //= '/data/01/user101/project/P_pyramidalis/a1.Ptr/3.bams.GCstat';


my $rep_gff=undef;
my $ref_dict_file='/data/00/user/user101/project/P_pyramidalis_genome/genome-v2.dict';

my $dps = &read_gc_file($gc_file, 5);
my %ref_chr2len = zzReadRefdict($ref_dict_file);
my %pops;

my $O = open_out_fh($outfile);


my $files_raw = &read_list($inlist);
my %files;
#foreach my $gvcf_dir (@gvcf_dirs) {
#foreach my $gvcf (<$gvcf_dir/*.gvcf.gz>) {
foreach my $gvcf (@$files_raw) {
        my $base = basename $gvcf;
        $base=~m/^([^.]+)/ or die;
        my $sid = $1;
        if (exists $files{$sid}) {
            say STDERR "duplicate sample id: $sid";
            next;
        }
        if (! exists $$dps{$sid}) {
            die "no gc stat for sample id: $sid";
        }
        $files{$sid} = $gvcf;
    }
#}

say STDERR Dumper \%files;

say $O join "\t", qw/sid homo hetero all filepath/;

MCE::Loop::init {
    chunk_size => 1, max_workers => $thread
};

#foreach (sort keys %files) {
mce_loop {
    my $sid = $_;
    my $file = $files{$sid};
    my ($homo, $hetero) = &get_homo_hetero_from_gvcf($file, $sid);
    my $all = $homo+$hetero;
    say $O join "\t", $sid, $homo, $hetero, $all, $file;
} (sort keys %files);

MCE::Loop::finish;

exit;


sub get_homo_hetero_from_gvcf {
    my ($in, $sid) = @_;
    my $I = open_in_fh($in);
    my $homo = 0;
    my $hetero = 0;
    my %ret;
    my $avg_dp = $$dps{$sid} // die;
    my $min_dp = $avg_dp / $dp_fold;
    my $max_dp = $avg_dp * $dp_fold;
    while (<$I>) {
        next if /^#/;
        chomp;
        my @F = split(/\t/,$_);
        say STDERR "@F" and next unless scalar(@F) == 10;
        my @F9 = split(/:/,$F[9]);
        my $dp = &get_DP_in_line(\@F, \@F9);
        if ($dp < $min_dp or $dp > $max_dp) {
            next;
        }
        my $is_hetero = &is_hetero($F9[0]);
        next if $is_hetero == -1; # miss
        my $len = &cal_gvcf_block_len(\@F);
        $hetero += $len if $is_hetero == 1;
        $homo += $len if $is_hetero == 0;
    }
    return ($homo, $hetero);
}

sub cal_gvcf_block_len {
    my ($F) = @_;
    my $infos = $$F[7];
    if ($infos=~/END=(\d+)/) {
        my $end = $1;
        my $start = $$F[1];
        my $len = $end - $start + 1;
        return $len;
    } else {
        return 1;
    }
}

sub is_hetero {
    my ($GT) = @_;
    my @GT = split(/[\/\|]/,$GT);
    my $is_het = 0;
    # miss
    if ( $GT[0] eq '.' || $GT[1] eq '.' ) {
        $is_het = -1;
    } elsif ( $GT[0] ne $GT[1] ) {
        $is_het = 1;
    } else {
        $is_het = 0;
    }
    return $is_het;
}

sub get_DP_in_line {
    my ($in, $F9) = @_;
    my @format = split(':', $$in[8]);
    my $dp_num;
    my %ret;
    foreach my $i (0..$#format) {
        if ($format[$i] eq 'DP') {
            $dp_num=$i;
            last;
        }
    }
    #die Dumper $in if ! defined $dp_num;
    return -1 unless $dp_num;
    #die join ' ', @$in unless $dp_num;
    my $dp = $$F9[$dp_num];
    return $dp;
}

sub read_gc_file {
    my ($in, $itemi) = @_;
    $itemi //= 7;
    my $I = open_in_fh($in);
    my %ret;
    while (<$I>) {
        chomp;
        my @F = split(/\t/,$_);
        my $id = $F[0];
        my $dp = $F[$itemi];
        $ret{$id} = $dp;
    }
    return \%ret; 
}

sub read_list {
    my ($in) = @_;
    my $I = open_in_fh($in);
    my %ret;
    while (<$I>) {
        chomp;
        my @F = split(/\t/,$_);
        my $id = $F[0];
        $ret{$id} = 1;
    }
    return ([sort keys %ret]); 
}


