use warnings;
use strict;
use v5.10;
use zzIO;
use zzVCF;

my $in = shift or die;
my $id = shift or die;
my $I = open_in_fh($in);
my $O = *STDOUT;

my $vcf_obj  =  new zzVCF($in) or die;
my $header = $vcf_obj->get_ids;

die unless @$header;
my %num2id=&get_num_in_header($id,$header);
say STDERR join " ", @num2id{ sort {$a<=>$b}keys %num2id};

my $m = 999;

POS:while (my $line_obj = $vcf_obj->next_pos) {
    my $pos = $line_obj->get_pos;
    my $chr = $line_obj->get_chr;
    my $alts = $line_obj->get_alts;
    #say scalar @$alts;
    next POS if scalar @$alts <= 1;
    length($_) == 1 or next POS foreach (@$alts);
    my $id_info = $line_obj->get_id_info;
    my %hash;
    ID:foreach my $i (sort {$a<=>$b} keys %num2id) {
        my $id = $num2id{$i};
        my $alleles = $line_obj->get_alleles->[$i];
        $_ eq '.' and next POS foreach (@$alleles);
        # next ID if $$id_info[$i]{DP} < 10;
        my @AD_split = split ',',$$id_info[$i]{AD};
        my @PL_split = split ',',$$id_info[$i]{PL};
        $_<3 and next POS foreach @AD_split;
        $_<1 and $_!=0 and next POS foreach @PL_split;
        my $div = $AD_split[$$alleles[1]] eq 0 ? 1 : $AD_split[$$alleles[0]] / $AD_split[$$alleles[1]];
        next ID if ( $div==0 or $div>$m or $div<1/$m );
        my $allele_char = join '', map{ $$alts[$_] } @$alleles;
        $hash{$i}{GT}=$allele_char;
        $hash{$i}{AD}=$$id_info[$i]{AD};
        $hash{$i}{PL}=$$id_info[$i]{PL};
    }
    #die keys %hash;
    my $result = &get_del(\%hash);
    if (defined $result) {
        say join "\t", $chr, $pos, $result;
    }
}


sub get_del {
	my ($hash)=@_ or die;
	my $t;
	my @nums = sort {$a <=> $b} keys %num2id;
    my @ids = @num2id{@nums};
    my $defalut=0;
    #$defalut-- while $defalut>$#nums; 
	foreach my $i (@nums) {
        next if $nums[$defalut]==$i;
        my $id = $num2id{$i} or die $i;
        if ( ! exists $$hash{$nums[$defalut]}{GT} ) {
            $t.="$ids[$defalut] Not exists : !!!! ";
        } elsif ( ! exists $$hash{$i}{GT} ) {
            $t.="$ids[$defalut] NotIn $id : $$hash{$nums[$defalut]}{AD} !!!! ";
        } elsif ( $$hash{$i}{GT} ne $$hash{$nums[$defalut]}{GT} ) {
            $t.="$ids[$defalut] ne $id : $$hash{$i}{AD} $$hash{$nums[$defalut]}{AD} !!!! ";
        } else {

        }
		delete $$hash{$i};
	}
    #say keys %$hash;
	foreach my $i(sort {$a<=>$b} keys %$hash) {
        next if $i==$nums[$defalut];
        my $id = $num2id{$i};
		$t.="$id NotIn $ids[$defalut] : $$hash{$i}{PL} !!!! ";
	}
	return defined $t ? $t : undef;
}

sub get_num_in_header {
    my ($id,$header) = @_ or die;
    my %ret;
    foreach my $i (0..scalar @$header-1) {
        next if $$header[$i] eq $id;
        $ret{$i} = $$header[$i] if $$header[$i]=~/^$id/;
    }
    return %ret;
}

