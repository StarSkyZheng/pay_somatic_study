
# 输入一个个体的多个VCF，进行比较


#use warnings;
use strict;
use v5.10;
use zzIO;
use zzVCF;
use File::Basename;
use MCE::Loop;

#my $prefix=shift or die;
#my @vcfs=<$prefix*>;
my @vcfs=@ARGV or die;
my $max_thread = 'auto';

MCE::Loop::init {
	max_workers => $max_thread, chunk_size => 1
};
my %hash = mce_loop {
    my $vcf = $_;
	my $base = basename($vcf);
	$base=~/^([^.]+)/ or die;
	my $id = $1;
	#$hash{$id} = &read_vcf($vcf,$id);
	MCE->gather( $id, &read_vcf($vcf,$id) );
}@vcfs;


my $fh = open_in_fh($vcfs[0]);

my $t = get_del(\%hash);
# $$t{$chr}{$pos}

foreach my $chr (sort keys %$t) {
	foreach my $pos (sort {$a<=>$b} keys %{$$t{$chr}}) {
		say join "\t", $chr, $pos, $$t{$chr}{$pos};
	}
}

exit;


while(<$fh>) {
	if(/^#/) {
		print;
		next;
	}
	my @a = split(/\t/,$_);
	next if exists $$t{$a[0]}{$a[1]};
	print join"\t", @a;
}

# &compare(\%hash);

exit;

sub get_del {
	my $hash=shift;
	my %t;
	my @ids = sort keys %{$hash};
	foreach my $id (@ids[1..$#ids]) {
		foreach my $chr (sort keys %{$$hash{$id}}) {
			foreach my $pos (sort{$a<=>$b} keys %{$$hash{$ids[0]}{$chr}}) {
				unless (( exists $$hash{$id}{$chr}{$pos}{GT} and $$hash{$id}{$chr}{$pos}{GT} eq $$hash{$ids[0]}{$chr}{$pos}{GT} ) 
					#or ($$hash{$id}{$chr}{$pos}{DP}<1000 or $$hash{$ids[0]}{$chr}{$pos}{DP}<1000) 
					#or ($$hash{$id}{$chr}{$pos}{ADF}<100 or $$hash{$ids[0]}{$chr}{$pos}{ADF}<100) 
					#or ($$hash{$id}{$chr}{$pos}{ADR}<100 or $$hash{$ids[0]}{$chr}{$pos}{ADR}<100)  
					) {
					#next if exists $t{$chr}{$pos};
					$t{$chr}{$pos}.="$ids[0] ne $id! ";
					#say STDERR "$id, $ids[0],$chr,$pos  $$hash{$id}{$chr}{$pos}{GT},$$hash{$id}{$chr}{$pos}{DP} ne $$hash{$ids[0]}{$chr}{$pos}{GT},$$hash{$ids[0]}{$chr}{$pos}{DP};";
				}
				delete $$hash{$id}{$chr}{$pos};
			}
			foreach my $pos(sort {$a<=>$b}keys %{$$hash{$id}{$chr}}) {
				#next if exists $t{$chr}{$pos};
				$t{$chr}{$pos}.="$id notin $ids[0]! ";
				#say STDERR "$id, $ids[0],$chr,$pos  $$hash{$id}{$chr}{$pos}{GT},$$hash{$id}{$chr}{$pos}{DP} ne $$hash{$ids[0]}{$chr}{$pos}{GT},$$hash{$ids[0]}{$chr}{$pos}{DP};" 
			}
		}
	}
	return \%t;
}

sub compare($) {
	my $hash=shift;
	my %t;
	my @ids = sort keys %{$hash};
	foreach my $id (@ids[1..$#ids]) {
		foreach my $chr (sort keys %{$$hash{$id}}) {
			foreach my $pos (sort{$a<=>$b} keys %{$$hash{$ids[0]}{$chr}}) {
				if (( $$hash{$id}{$chr}{$pos}{GT} ne $$hash{$ids[0]}{$chr}{$pos}{GT} ) or ($$hash{$id}{$chr}{$pos}{DP}<1000 or $$hash{$ids[0]}{$chr}{$pos}{DP}<1000) ) {
					#next if exists $t{$chr}{$pos};
					$t{$chr}{$pos}++;
					say STDERR "$id, $ids[0],$chr,$pos  $$hash{$id}{$chr}{$pos}{GT},$$hash{$id}{$chr}{$pos}{DP} ne $$hash{$ids[0]}{$chr}{$pos}{GT},$$hash{$ids[0]}{$chr}{$pos}{DP};";
				}
				delete $$hash{$id}{$chr}{$pos};
			}
			foreach my $pos(sort {$a<=>$b}keys %{$$hash{$id}{$chr}}) {
				#next if exists $t{$chr}{$pos};
				$t{$chr}{$pos}++;
				say STDERR "$id, $ids[0],$chr,$pos  $$hash{$id}{$chr}{$pos}{GT},$$hash{$id}{$chr}{$pos}{DP} ne $$hash{$ids[0]}{$chr}{$pos}{GT},$$hash{$ids[0]}{$chr}{$pos}{DP};" 
			}
		}
	}
}


sub read_vcf($$) {
	my $vcf=shift;
	my $id=shift;
	say STDERR $vcf;
	my %ret;
	my $vcf_obj  =  new zzVCF($vcf) or die;
	POS:while (my $line_obj = $vcf_obj->next_pos) {
		my $pos = $line_obj->get_pos;
		my $chr = $line_obj->get_chr;
		my $alts = $line_obj->get_alts;
		next POS if scalar @$alts <= 2;
		length($_) == 1 or next POS foreach ($alts->[0..@$alts]);
		my $id_info = $line_obj->get_id_info;
		my $alleles = $line_obj->get_alleles->[0];
		next POS if $$id_info[0]{DP} < 10;
		#$_==0 and next POS foreach split ',', $$id_info[0]{SB};
		#my @SB_split = split ',', $$id_info[0]{SB};
		#next POS if $SB_split[0] + $SB_split[1] < 2;
		#next POS if $SB_split[2] + $SB_split[3] < 2;
		my @AD_split = split ',',$$id_info[0]{AD};
		my $div = $AD_split[$$alleles[1]] eq 0 ? 1 : $AD_split[$$alleles[0]] / $AD_split[$$alleles[1]];
		my $m = 2;
		next POS if ( $div==0 or $div>2 or $div<1/$m );
		#next POS if $$id_info[0]{ADF} < 100;
		#next POS if $$id_info[0]{ADR} < 100;
		#say STDERR join "!", @$alleles;
		my $allele_char = join '', map{ $$alts[$_] } @$alleles;
		$ret{$chr}{$pos}{GT}=$allele_char;
		#say "$chr $pos";
		#$ret{$chr}{$pos}{DP}=$$id_info[0]{DP};
		#$ret{$chr}{$pos}{ADF}=$$id_info[0]{ADF};
		#$ret{$chr}{$pos}{ADR}=$$id_info[0]{ADR};
	}
	return \%ret;
}


