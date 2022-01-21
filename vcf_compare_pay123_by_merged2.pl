use warnings;
use strict;
use v5.24;
use zzIO;
use zzVCF;

my $in = shift or die;
my $id = shift or die;
my $out_prefix = shift or die;
my $I = open_in_fh($in);
#my $O = *STDOUT;
my $O1 = open_out_fh("$out_prefix.ERR.txt");
my $O2 = open_out_fh("$out_prefix.OK.txt");

my $vcf_obj  =  new zzVCF($in) or die;
my $header = $vcf_obj->get_ids;

die unless @$header;
my %num2id=&get_num_in_header($id,$header);
say STDERR join " ", @num2id{ sort {$a<=>$b}keys %num2id};

my $m = 999;
my %count;

POS:while (my $line_obj = $vcf_obj->next_pos) {
    my $pos = $line_obj->get_pos;
    my $chr = $line_obj->get_chr;
    my $alts = $line_obj->get_alts;
    #say scalar @$alts;
    next POS if scalar @$alts <= 1;
    length($_) == 1 or next POS foreach (@$alts);
    my $id_info = $line_obj->get_id_info;
    my %hash;
    my %t;
    ID:foreach my $i (sort {$a<=>$b} keys %num2id) {
        my $id = $num2id{$i};
        my $alleles = $line_obj->get_alleles->[$i];
        $_ eq '.' and next ID foreach (@$alleles);
        # next ID if $$id_info[$i]{DP} < 10;
        my @AD_split = split ',',$$id_info[$i]{AD};
        $t{$AD_split[0]}{$AD_split[1]}++;
        #my $div = $AD_split[$$alleles[1]] eq 0 ? 1 : $AD_split[$$alleles[0]] / $AD_split[$$alleles[1]];
        #next ID if ( $div==0 or $div>$m or $div<1/$m );
        my $allele_char = join '', map{ $$alts[$_] } @$alleles;
        $hash{$i}{GT}=$allele_char;
        #$hash{$i}{AD}=$div;
        #$hash{$i}{PL}=$$id_info[$i]{PL};
    }
    #die keys %hash;
    my $result = &get_del(\%hash);
    #say join "\t", $chr, $pos, $result and next POS if defined $result;
    if (defined $result) {
        foreach my $i (keys %t) {
            $count{ERR}{$i}{$_}++ foreach keys $t{$i}->%*;
        }
    } else {
        foreach my $i (keys %t) {
            $count{OK}{$i}{$_}++ foreach keys $t{$i}->%*;
        }
    }
}


say $O1 join "\t",qw/i j ERR/;
foreach my $i (sort {$a<=>$b} keys $count{ERR}->%*) {
    say $O1 join "\t", $i, $_, $count{ERR}{$i}{$_} foreach sort {$a<=>$b} keys $count{ERR}{$i}->%*;
}

say $O2 join "\t",qw/i j OK/;
foreach my $i (sort {$a<=>$b} keys $count{OK}->%*) {
    say $O2 join "\t", $i, $_, $count{OK}{$i}{$_} foreach sort {$a<=>$b} keys $count{OK}{$i}->%*;
}

exit;

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
            $t.="$ids[$defalut] NotIn $id : $$hash{$nums[$defalut]}{GT} !!!! ";
        } elsif ( $$hash{$i}{GT} ne $$hash{$nums[$defalut]}{GT} ) {
            $t.="$ids[$defalut] ne $id : $$hash{$i}{GT} $$hash{$nums[$defalut]}{GT} !!!! ";
        } else {

        }
		delete $$hash{$i};
	}
    #say keys %$hash;
	foreach my $i(sort {$a<=>$b} keys %$hash) {
        next if $i==$nums[$defalut];
        my $id = $num2id{$i};
		$t.="$id NotIn $ids[$defalut] : $$hash{$i}{GT} !!!! ";
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

