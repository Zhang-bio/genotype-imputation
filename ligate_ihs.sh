#!/bin/bash
dp=3
sign=Scenario
#glimpse2_ligate
indir=../phase_340/d${dp}
OUT=./gwas_ligate/${sign}_depth_${dp}_chrz_ligated.vcf
LST=./gwas_ligate/${sign}_depth_${dp}_list.chrz.txt
ls -1v ${indir}/ref_d${dp}_test_chr*.vcf.gz > ${LST}
GLIMPSE2_ligate --input ${LST} --output ${OUT} --threads 16
bgzip -f ${OUT}
bcftools index -f ${OUT}.gz

bcftools annotate --rename-chrs re_chr.txt ${sign}_depth_3_chrz_ligated.vcf.gz | bgzip -c > ${sign}_depth_${dp}_chrz_rechr.vcf.gz


java -Xmx40g -jar /public/home/software/beagle.22Jul22.46e.jar nthreads=16 gt=/public/home/chenlin/1.projects/9.bam_down/22jc/22jc_mindp3/ligate_340/gwas_ligate/${sign}_depth_${dp}_chrz_rechr.vcf.gz out=${sign}_beagle_z impute=true

plink1.9 --vcf ${sign}_beagle_z.vcf.gz --allow-extra-chr --chr 1-49  --chr-set 49 --double-id --make-bed    --set-missing-var-ids @_# --out /public/home/user/1.projects/9.bam_down/22jc/22jc_mindp3/gwas/${sign}_imputed

plink1.9 --bfile ${sign}_imputed --chr-set 49  --geno 0.1   --hwe 0.000001   --maf 0.05  --make-bed  --out ./${sign}_b_clean

mkdir ihs

filechr=${sign}_b_clean
path=/public/home/user/1.projects/9.bam_down/22jc/22jc_mindp3/gwas

for i in `awk '{print $1}' $path/$filechr.bim | uniq`
do
plink1.9 --bfile $path/$filechr --chr-set 49 --keep-allele-order --chr ${i} --recode --out $filechr${i} --noweb
plink1.9 --file $filechr${i} --chr-set 49 --keep-allele-order --recode vcf --out $filechr${i} --noweb
grep -v "^#" $filechr${i}.vcf | awk '{printf("%s %s_%s %.6f %s\n", $1,$1,$2,$2/1000000,$2)}' > $filechr${i}.maptmp
#unnorm
selscan --ihs --threads 15 --vcf $filechr${i}.vcf --map $filechr${i}.maptmp --maf 0.05 --keep-low-freq --max-gap 1000000 --out ./ihs/$filechr${i}ihs
#norm
norm --ihs --files ./ihs/$filechr${i}ihs.ihs.out --bp-win --winsize 1000000

cat ./ihs/$filechr${i}ihs.ihs.out >> ./ihs/IHS.temp

#cat norm
cat ./ihs/$filechr${i}ihs.ihs.out.100bins.norm >> ./ihs/normIHS.temp

rm $filechr${i}.* ./ihs/$filechr${i}.* ./ihs/$filechr${i}ihs.* ./nsl/$filechr${i}nsl.* ./nsl/$filechr${i}.* ./ihh12/$filechr${i}ihh12.* ./ihh12/$filechr${i}.* ./pi/$filechr${i}pi.* ./pi/$filechr${i}.*
done
awk 'BEGIN{print "SNP BP 1freq ihh1 ihh0 unnormiHS"}{print $1 " ", $2 " ", $3 " ", $4 " ", $5 " ", $6 " "}' ./ihs/IHS.temp > ./ihs/IHS
awk 'BEGIN{print "SNP BP 1freq ihh1 ihh0 unnormiHS normiHS crit"}{print $1 " ", $2 " ", $3 " ", $4 " ", $5 " ", $6 " ", $7 " ", $8 " "}' ./ihs/normIHS.temp > ./ihs/normIHS_1.temp

awk '{print $1}' ./ihs/normIHS_1.temp | cut -d '_' -f1 | sed 's/SNP/CHR/g' > ./ihs/normIHSCHR.temp
paste -d ' ' ./ihs/normIHSCHR.temp ./ihs/normIHS_1.temp > ./ihs/normIHS.txt

cd ihs
R CMD BATCH ../IHS_Plot.r
mv iHS_manhattan.png ${filechr}.iHS_manhattan.png
mv significant_top001_iHS_snps.txt ${filechr}.significant_top001_iHS_snps.txt
