use zzIO;
use v5.24;
use warnings;
use strict;


my $work_dir = shift or die;

#my @files=`find $work_dir -name '*.gff'`;

my @files = ("$work_dir/1.CNVnator/CNV.gff",  
	"$work_dir/3.breakdancer/BreakDancer.config.result.gff",
	"$work_dir/2.pindel/Pindel.gff");

chomp @files;

my @uniq_feature = &get_uniq_feature(\@files);

my $mergeBed = '/home/share/users/zhengzeyu2018/software/bedtools2/bin/mergeBed';
foreach my $feature (@uniq_feature, 'ALL') {
    my $out = "$work_dir/merged.$feature.bed";
    open(my $O, qq! | $mergeBed -i -  > $out!);
    my %h;
    foreach my $file (@files) {
        my $I = open_in_fh($file);
        while (<$I>) {
            chomp;
            next unless $_;
            my @a = split(/\t/);
            next unless ($a[2] ne 'Duplication' and $feature eq 'ALL') or $a[2] eq $feature;
            next if $a[4] - $a[3] > 5e4; # 50K
            die $_ . $file unless exists $a[3];
            $h{$a[0]}{$a[3]} .= "$_\n";
        }
    }
    foreach my $chr(sort keys %h) {
        foreach my $pos (sort {$a<=>$b} keys $h{$chr}->%*) {
            print $O $h{$chr}{$pos};
        }
    }
    close $O;
}



sub get_uniq_feature($) {
    my $files = shift or die;
    my %ret;
    foreach my $file (@$files) {
        my $I = open_in_fh($file);
        while(<$I>){
            chomp;
            next if /^#/;
            next unless $_;
            my @a = split(/\t/, $_);
            $ret{$a[2]}=0 unless exists $ret{$a[2]};
        }
    }
    return sort keys %ret;
}


