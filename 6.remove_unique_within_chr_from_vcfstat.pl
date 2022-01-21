#!/usr/bin/env perl
#===============================================================================
#
#         FILE: remove_6_chr11.pl
#
#        USAGE: ./remove_6_chr11.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 2019年01月02日 11时23分21秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use zzIO;

my $in_stat = shift;
my $out_stat = shift;
my $id = shift;
my $chr = shift or die 'usage: in_stat out_stat id chr' ;

my $I = open_in_fh($in_stat);
my $O = open_out_fh($out_stat);
my @header;
while (<$I>) {
	chomp;
	next unless $_;
	if (/^#/) {
		@header = split(/\t/,$_);
		say $O $_;
		next;
	}
	die unless @header;
	say $O $_ if &need_pass($_);
}

sub need_pass {
	my $in = shift or die;
	my @a = split(/\t/, $in);
	#return 0 unless $a[0] eq 'Chr11';
	my $id_alt=0;
	my $other=0;
	foreach my $i (5..$#a) {
		$a[$i]=0 if $a[$i] eq 'NA';
		if($header[$i] eq $id) {
			$id_alt += $a[$i];
		}else {
			$other += $a[$i];
		}
	}
	return 0 if ($id_alt>0 and $other==0 and $a[0] eq $chr);
	return 1;
}

