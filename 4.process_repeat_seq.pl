# 计算merged_vcf里面，三个重复只保留基因型出现次数最多的基因型，输入应为pay01，对应pay01\d(三个重复)，输出到STDOUT

use warnings;
use strict;
use v5.24;
use zzIO;
use zzVCF;
use MCE::Map;

my $in = shift or die 'usage: in thread id1 id2 ..';
my $thread = shift or die;
my $I_VCF = open_in_fh($in);
my @ids_m = @ARGV or die;
#my $I = open_in_fh($in);
my $O = *STDOUT;

my $vcf_obj = new zzVCF({fh=>$I_VCF}) or die;
my $commit = $vcf_obj->get_commits;
my $header = $vcf_obj->get_header;
my $id_header = $vcf_obj->get_ids;


die @$id_header unless @$id_header;
my %num2id=&get_num_in_header(\@ids_m,$id_header);
say STDERR join " ", $_ , ":", $num2id{$_}->@{ sort {$a<=>$b} keys $num2id{$_}->%* } foreach @ids_m;


say $O $_ foreach $commit->@*;
say $O join "\t", $header->@[0..8], sort keys %num2id;


MCE::Map::init {
		max_workers => $thread, chunk_size => 'auto',
};

# while (<$I_VCF>) {
#     my $ret =  &process_pos($_);
#     say $ret if defined $ret;
# }


my @ret = grep {defined $_} mce_map{&process_pos} <$I_VCF>;
MCE::Map::finish;
say foreach @ret;

sub process_pos {
#POS:while (my $line_obj = $vcf_obj->next_pos) {
    my $in = $_;
    chomp $in;
    my $line_obj = new zzVCF::line(line=>$in);
    my $pos = $line_obj->get_pos;
    my $chr = $line_obj->get_chr;
    my $alts = $line_obj->get_alts;
    return undef if scalar @$alts > 2;
    length($_) == 1 or return undef foreach (@$alts);
    my $id_info = $line_obj->get_id_info;
    
    my @say = ( $line_obj->get_split->@[0..7] , 'GT:ERR' );
    foreach my $id_merge (@ids_m) {
        my %hash;
        my @nids = keys $num2id{$id_merge}->%*;
        ID:foreach my $i (@nids) {
            #my $id = $num2id{$i};
            my $alleles = $line_obj->get_alleles->[$i];
            $hash{$i}{GT}=join '/',@$alleles;
        }
        my $result = &get_most(\%hash);
        push @say, $result or die;
    }
    return join "\t", @say;
}




sub get_num_in_header {
    my ($ids,$id_header) = @_ or die;
    my %ret;
    foreach my $id_merge (@$ids) {
        foreach my $i (0..scalar @$id_header-1) {
            my $id = $$id_header[$i];
            next if $id eq $id_merge;
            $ret{$id_merge}{$i} = $id if $id=~/^$id_merge/;
        }
    }
    return %ret;
}



sub get_most {
	my ($hash)=@_ or die;
	my $t;
	#my @nums = sort {$a <=> $b} keys %num2id;
    #my @ids = @num2id{@nums};
    my %ret;
	foreach my $i (keys %$hash) {
        if ( ! exists $$hash{$i}{GT} ) {
            next;
        } else {
            $ret{ $$hash{$i}{GT} }++;
        } 
	}
    my $most_append;
    if ( exists $ret{'./.'} ) {
        my $miss = $ret{'./.'};
        delete $ret{'./.'};
        $most_append .= scalar keys %ret<=1 ? ":$miss" : ":!";
    } else {
        $most_append .= scalar keys %ret==1 ? ':.' : ':!';
    }
    my $most_t;
    if ( %ret ) {
        ($most_t) = sort {$ret{$b}<=>$ret{$a}} keys %ret;
    } else {
        $most_t = './.';
    }
	return $most_t . $most_append;
}


sub flt {
        # my $m = 999;
        #$_ eq '.' and next POS foreach (@$alleles);
        # next ID if $$id_info[$i]{DP} < 10;
        #my @AD_split = split ',',$$id_info[$i]{AD};
        #my @PL_split = split ',',$$id_info[$i]{PL};
        #$_<3 and next POS foreach @AD_split;
        #$_<100 and $_!=0 and next POS foreach @PL_split;
        #my $div = $AD_split[$$alleles[1]] eq 0 ? 1 : $AD_split[$$alleles[0]] / $AD_split[$$alleles[1]];
        #next ID if ( $div==0 or $div>$m or $div<1/$m );
        #my $allele_char = join '', map{ $$alts[$_] } @$alleles;
        #$hash{$i}{GT}=$allele_char;
        #$hash{$i}{AD}=$$id_info[$i]{AD};
        #$hash{$i}{PL}=$$id_info[$i]{PL};
}

