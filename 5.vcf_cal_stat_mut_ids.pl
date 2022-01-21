# print $chr, $pos, $from, $to, $change_conut, @$mut_ids (0 or 1 or NA, 0 for mut, 1 for bg, NA for miss)

use warnings;
use strict;
use v5.24;
use zzIO;
use zzVCF;

my $in = shift or die 'no input vcf! ';
#my @ids_m = @ARGV or die;
#my $I = open_in_fh($in);
my $O = *STDOUT;

my $vcf_obj = new zzVCF({vcf=>$in}) or die;
my $commit = $vcf_obj->get_commits;
my $header = $vcf_obj->get_header;
my $id_header = $vcf_obj->get_ids;


die @$id_header unless @$id_header;
#say STDERR join " ", $_ , ":", $num2id{$_}->@{ sort {$a<=>$b} keys $num2id{$_}->%* } foreach @ids_m;

say $O join "\t", $header->@[0,1], qw/from to change_conut mut_ids/;
POS:while (my $line_obj = $vcf_obj->next_pos) {
    my $pos = $line_obj->get_pos;
    my $chr = $line_obj->get_chr;
    my $alts = $line_obj->get_alts;
    next POS if scalar @$alts > 2;
    #length($_) == 1 or next POS foreach (@$alts);
    my $id_info = $line_obj->get_id_info;

    my @say = ( $line_obj->get_split->@[0,1] );
    my %hash;
    #say scalar @$id_header;
    ID:foreach my $i (0..@$id_header-1) {
        #my $id = $num2id{$i};
        my $alleles = $line_obj->get_alleles->[$i];
        $hash{$i}=join '/', @$alleles;
        #die if $hash{$i}=~/\./;
    }
    my ($from, $to, $change_conut, $ids_result) = &get_most(\%hash);
	next POS unless defined $from;
    #$from = $$alts[$from];
    #$to = $$alts[$to];
    die join "\t", @say unless defined $from;
    push @say, $from, $to, $change_conut, @$ids_result;

    say $O join "\t", @say;
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
        if ( ! exists $$hash{$i} ) {
            die;
            next;
        } elsif ( $$hash{$i} eq './.' ) {
            next;
        } else {
            $ret{ $$hash{$i} }++;
        }
	}
	#die 'more than one mut allele type'. keys %ret if keys %ret > 2;
    say STDERR 'more than one mut allele type'. keys %ret and return undef if keys %ret > 2;
    # return undef if keys %ret > 2;
    my ($most, $mut) = sort{ $ret{$b}<=>$ret{$a} } keys %ret;
    die '? most < 0.6 ? ' if $ret{$most} < 0.6;
    my @t;
	foreach my $i (keys %$hash) {
        my $alles = $$hash{$i};
        if ( $most eq $alles ) {
            #$t[$i]=0;
        } elsif ($alles eq './.') {
            #$t[$i] = 'NA';
        } elsif ( $mut eq $alles ) {
            #$t[$i]=1;
            push @t, "$$id_header[$i]:$alles";
        } else {
            die;
        }
	}
    my @most = split(/\//,$most);
    my @mut = split(/\//,$mut);
    my @t2=(0,0);
    foreach my $i (@most, @mut) {
        $t2[$i]++;
    }
    my $change_conut = 1;
    $t2[$_]==2 and $change_conut=2 and last foreach 0..$#t2;
    #($most, $mut) = sort{ $t2[$b]<=>$t2[$a] } 0..$#t2;
	return ($most, $mut, $change_conut, \@t );
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

