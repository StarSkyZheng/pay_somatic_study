reference_fasta='/home/share/users/zhengzeyu2018/project/P_pyramidalis_genome/genome-v2.fa'
bam_dir='/home/share/users/zhengzeyu2018/test/t/bams' # program will use all the bam file in this dir
out_dir='/home/share/users/zhengzeyu2018/test/t/out'
perl_script_dir='/home/share/users/zhengzeyu2018/software/BioScripts'
in_vcf_gz='/home/pool_3/users/zhengzeyu2018/project/P_pyramidalis/03.callSNP/06.HDflt_by_species_all/pay/2.SNPwithDp.pay_all.HC.join.prefit.sv.vcf.gz'
samtools_019='/home/share/users/zhengzeyu2018/software/samtools/samtools-0.1.19/samtools'

for record in `zcat $in_vcf_gz | perl -ne 'next if (/^\#/); my ($chrom, $pos, $sample) = (split /\s+/)[0,1,2,7];
    if($sample =~ /\;/) {$sample = "Shared"} print "$sample;$chrom:$pos#$1\n";'`;
#$record = '.;Chr1:108311#''    Shared;Chr1:290#
do
    sample=${record/;*} # Shared
    mutation=${record/*;} # Chr1:290#
    mut_base=${mutation/*\#} # empty
    mutation=${mutation/\#*} # Chr1:290
    chrom=${mutation/:*} #Chr1
    mut_pos=${mutation/*:} # 290

    start_pos=`echo "${mut_pos}-200" | bc`
    end_pos=`echo "${mut_pos}+200" | bc`

    echo "${sample} ${chrom} ${mut_pos} ${mut_base}"
    
    if [[ -n ${out_dir}/${sample} ]]; then
        mkdir -pv ${out_dir}/${sample}
    fi
    
    #echo -e "${chrom} ${start_pos} ${end_pos}\n${chrom} ${start_pos} ${mut_pos}\n${chrom} ${mut_pos} ${end_pos}\n" | \
    #    ${perl_script_dir}/fasta_process.pl --rows 0 1 2 --subset 1 2 --query - \
    #    --fasta ${reference_fasta} > ${out_dir}/${sample}/${chrom}_${mut_pos}.fa
    
    for bam_file in `ls ${bam_dir}/*.bam`;
    do
        #`cd ${bam_dir};${samtools_019} index ${bam_file}`
        m=`basename ${bam_file} | cut -d"." -f1`
        #bam_file_base=basename
        $samtools_019 view -X -F 3844 ${bam_file} ${chrom}:${mut_pos}-${mut_pos} | \
            awk -v name=${sample} -v id=${m} 'BEGIN {OFS = FS = "\t"}; {print ">"name"|"$2"|"$1"\n"$10"|"id;}' \
            >> ${out_dir}/${sample}/${chrom}_${mut_pos}.fa
    done

    ${perl_script_dir}/reference_align.pl -i ${out_dir}/${sample}/${chrom}_${mut_pos}.fa \
        > ${out_dir}/${sample}/${chrom}_${mut_pos}.aln.fas #&& \
        #rm -v ${out_dir}/${sample}/${chrom}_${mut_pos}.fa
done
