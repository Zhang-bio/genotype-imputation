# MAIN work pepile
# Six popelines for genotype imputation of low-coverage whole-genome sequencing
* (1) Fixed sample size to compare imputation performance at different sequencing depths (S1,STITCH)
* (2) Fixed low depth to compare imputation performance at different sample sizes (S2,STITCH)
* (3) Fixed sample size to compare imputation performance at different sequencing depths (S3,glimpse2)
* (4) Fixed low depth to compare imputation performance at different reference panel sizes (S4,glimpse2)
* (5) (S5,glimpse2) establishing a key reference panel (K_Ref) using 22 core individuals (23×) for gradient imputation of the validation set from 0.1× to 6×
* (6) (S6,glimpse2) creating a combined reference panel (C_Ref) by merging 22 high-depth individuals (23×) with 372 routine individuals (9×) for gradient imputation of the validation set from 0.1× to 6×.
* (7) This script needs to imputed all chromosomes of all samples : ligate+ihs analysis
# Requirements
* bwa
* samtools
* pandepth
* sambamba-1.0.1
* vcftools
* bcftools v1.10
* GLIMPSE2
* R v3.6.3 requires STITCH package
* beagle v4 and v5
* GATK v3.7
* plink1.9
* selscan
# Usage
```
(1) and (2)
sh ./S1-2/STITCH - diff-sample-depth-size-clean.sh 
(3) 
sh ./S3/GLIMPSE2 -N- clean.sh
(4)
sh ./S4/GLIMPSE2-diff-sample-size-clean.sh
(5)
sh ./S5/GLIMPSE2 -K- clean.sh
(6)
sh ./S6/GLIMPSE2 -C- clean.sh
(7)
sh ./ligate_ihs.sh
```
#raw_data  
All raw data has been submitted to NCBI, SUB15819465; After the article is published, it will also be published.

# id.txt  
all_id_list.txt : all tested samples id
6_id_list.txt : valid samples id

# Publication
Weiren Zhang, Lin  Chen, Yijia  Shih, Peitan  Jia, Siyi Zhou, Qionghui Qin, Yaodong Zhang, Kewei Huang, Gongyi Lin, Xiaopeng Wang, Haihui Ye. Optimization strategies for genotype imputation accuracy in low-coverage whole-genome sequencing of Scylla paramamosain. 
