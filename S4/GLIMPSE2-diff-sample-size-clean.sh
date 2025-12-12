#!/bin/bash
#1.sample-bam-file-extra

#path & file 
folder=/data1/1.Crab_lab/1.Crab/2.ReseqData/1.sp/2.Site/2.cleandata
bam_outfile_path=/public/home/user/1.projects/9.bam_down/400bam
name_list=/public/home/user/1.projects/9.bam_down/400bam/sp_inds.txt
ref_genome=/public/home/user/1.projects/9.bam_down/400bam/ref/genome.fasta


process_file() {
    local sampleID=$1
    
    echo "Samples being processed: $sampleID"
    
    # check input file
    local file1="${folder}/sp/${sampleID}/Site${sampleID}_1.clean.fq.gz"
    local file2="${folder}/sp/${sampleID}/Site${sampleID}_2.clean.fq.gz"
    
    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "wrong: input file not find: $file1 或 $file2"
        return 1
    fi
    
    # 1)mapping and analysis
    echo "blast samples: $sampleID"
    
    # bwa mapping
    bwa mem -t 4 -R "@RG\tID:2\tPL:illumina\tLB:library\tSM:${sampleID}" \
        "$ref_genome" "$file1" "$file2" > "${bam_outfile_path}/${sampleID}_aln_refgenome.sam"
    
    # check bwa running or not
    if [ $? -ne 0 ]; then
        echo "wrong: bwa failure: $sampleID"
        return 1
    fi
    
    # 2). sam file to bam file
    samtools view -Sb "${bam_outfile_path}/${sampleID}_aln_refgenome.sam" > "${bam_outfile_path}/${sampleID}_aln_refgenome.bam"
    
    # check samtools 
    if [ $? -ne 0 ]; then
        echo "wrong: samtools default: $sampleID"
        return 1
    fi
    
    echo "sample $sampleID successful"
    
    # 3). Clean up intermediate files
    rm "${bam_outfile_path}/${sampleID}_aln_refgenome.sam"
    
    # sort and index
    # samtools sort -@ 4 -o "${bam_outfile_path}/${sampleID}_aln_refgenome.sorted.bam" "${bam_outfile_path}/${sampleID}_aln_refgenome.bam"
    # samtools index "${bam_outfile_path}/${sampleID}_aln_refgenome.sorted.bam"
    
    return 0
}

# export function and parameter
export -f process_file
export folder bam_outfile_path ref_genome

# check genome
if [ ! -f "$ref_genome" ]; then
    echo "wrong: reference genome not exist: $ref_genome"
    exit 1
fi

# check samples
if [ ! -f "$name_list" ]; then
    echo "wrong: sample list does not exist : $name_list"
    exit 1
fi

echo "start to dealing..."
echo "Using reference genome: $ref_genome"
echo "output list: $bam_outfile_path"

# use xargs parallel processing
cat "$name_list" | xargs -n 1 -P 6 -I {} bash -c 'process_file "$@"' _ {}

echo "all success"

#2.Calculate downsampling coefficient

#output file path
mean_output_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/mean_depth"
mkdir -p "$mean_output_dir"

depths=("0.1" "0.5" "1" "2" "3" "4" "5" "6")

for dp in "${depths[@]}"; do
        mean_output_file="${mean_output_dir}/400bam_d${dp}.txt"
    > "$mean_output_file"  #clean ouput file
    
    cd "$bam_outfile_path"
    bams=$(ls *_refgenome.bam 2>/dev/null)
    
    if [ -z "$bams" ]; then
        echo "warning: in $bam_outfile_path not find *_refgenome.bam"
        continue
    fi
    
    # check all BAM files
    for bam_file in $bams; do
        # extra file name for PanDepth software result file's name
        pandepth_output="${bam_file%_refgenome.bam}"
        
        # PanDepth calculate depth for sample
        pandepth -i "$bam_file" -o "${mean_output_dir}/${pandepth_output}"
        
        # check PanDepth success or not
        if [ $? -eq 0 ]; then
            #echo "PanDepth success，Start processing the result file：$pandepth_output"
            
            #Enter the file name, which is the output file of PanDepth
            input_file="${mean_output_dir}/${pandepth_output}.chr.stat.gz"
            
            if [ -f "$input_file" ]; then
                #Read the last comment line starting with # # and extract the last value after:
                value=$(zcat "$input_file" 2>/dev/null | grep '##' | tail -n 1 | awk -F ':' '{print $NF}')
                
                if [ -n "$value" ] && [ "$value" != "0" ]; then
                    # downsample-targe depth
                    value=$(echo "scale=6; $dp/$value" | bc 2>/dev/null || echo "0")
                    
                    formatted_value=$(printf "%.6f" "$value" 2>/dev/null || echo "0.000000")
                    
                    echo "${pandepth_output}	$formatted_value" >> "$mean_output_file"
                    
                else
                    echo "wrong: donot extra $input_file with depth value"
                    echo "${pandepth_output}	0.000000" >> "$mean_output_file"
                fi
            else
                echo "wrong:  $input_file donot find"
                echo "${pandepth_output}	0.000000" >> "$mean_output_file"
            fi
        else
            echo "PanDepth wrong，check inputfile：$bam_file"
            echo "${pandepth_output}	0.000000" >> "$mean_output_file"
        fi
    done
done

#3.Target sample downsampling
for dp in "${depths[@]}"; do
    down_output_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/downsampling_bam_${dp}"
    mkdir -p "$down_output_dir"
    
    depth_file="${mean_output_dir}/400bam_d${dp}.txt"
    
    if [ ! -f "$depth_file" ]; then
        echo "Error: Depth file does not exist: $depth_file"
        continue
    fi
    
    # Read all ID numbers and their corresponding coefficients from the downsampling coefficient file.
    declare -A IDS
    while IFS=$'\t' read -r id depth || [ -n "$id" ]; do
        # Skip blank lines
        [ -z "$id" ] && continue
        IDS["$id"]="$depth"
    done < "$depth_file"
    
    # Process each ID number and its corresponding coefficient in a loop
    for id in "${!IDS[@]}"; do
        # Constructing the desired BAM filename
        bam_file="${id}_refgenome.bam"
        BAM_FILE="$bam_outfile_path/$bam_file"
        
        # Check whether the BAM file exists
        if [[ -f "$BAM_FILE" ]]; then
            # Obtain the decimation coefficients
            DEPTH="${IDS[$id]}"
            
            # Verify whether the depth value is valid
            if [ "$DEPTH" = "0.000000" ] || [ -z "$DEPTH" ]; then
                echo "Warning: Invalid depth value for ID $id: $DEPTH"
                continue
            fi
            
            # Constructing the output filename
            OUTPUT_BAM="$down_output_dir/${bam_file%.bam}_downsampled.bam"
            
            # Perform the downsampling operation, incorporating the parameters specified in the reference command.
            seed=123
            sambamba-1.0.1 view -h -t 4 -f bam --subsampling-seed="$seed" -s "$DEPTH" "$BAM_FILE" -o "$OUTPUT_BAM"
            
            # Verify whether the sambamba command executed successfully
            if [[ $? -eq 0 ]]; then
                echo "Downsampling completed for $id with depth $DEPTH"
            else
                echo "Error: Downsampling failed for $id"
            fi
        else
            echo "Error: BAM file not found for ID $id: $BAM_FILE"
        fi
    done
done

echo "All downsampling operations completed."

#4.Sorting and deduplication
for dp in "${depths[@]}"; do
    sort_output_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/sort_markdup_bam_${dp}"
    sort_input_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/downsampling_bam_${dp}"
    
    # Ensure the input folder exists
    if [ ! -d "$sort_input_dir" ]; then
        echo "Error: Input directory does not exist: $sort_input_dir"
        continue
    fi
    
    # Ensure the output folder exists
    mkdir -p "$sort_output_dir"
    
    # Enter the input folder

    cd "$sort_input_dir" || continue
    
    # Iterate through all BAM files in the input folder
    for bam_file in *.bam; do
        # Check whether the file exists
        if [ -f "$bam_file" ]; then
            # Constructing the output filename

            sort_output_bam="${sort_output_dir}/${bam_file%.bam}_new.bam"
            
            echo "Processing $bam_file..."
            
            # Sorting using samtools

            samtools sort -@ 16 "$bam_file" -o "${sort_output_bam%.bam}_sorted.bam"
            
            # Verify whether the samtools sorting has been successfully completed

            if [ $? -eq 0 ]; then
                # Using sambamba for deduplication

                sambamba-1.0.1 markdup -r "${sort_output_bam%.bam}_sorted.bam" "$sort_output_bam" -t 16 
                
                # Verify whether sambamba deduplication has been successfully completed.

                if [ $? -eq 0 ]; then
                    # Clear intermediate files

                    rm "${sort_output_bam%.bam}_sorted.bam"
                    echo "Successfully processed $bam_file"
                else
                    echo "sambamba-1.0.1 markdup failed for $bam_file"
                fi
            else
                echo "samtools sort failed for $bam_file"
            fi
        fi
    done
done

echo "Processing complete."

#5.Classify samples according to different depths and chromosomes
# Note: This assumes chromosomes 1 to 49 are present; adjustments may be required based on actual circumstances.
for dp in "${depths[@]}"; do
    splite_input_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/sort_markdup_bam_${dp}"
    splite_outchrdir="/public/home/user/1.projects/9.bam_down/calculate_depth/chrbam/test_d${dp}"
    
    # Ensure the input directory exists
    if [ ! -d "$splite_input_dir" ]; then
        echo "Error: Input directory does not exist: $splite_input_dir"
        continue
    fi
    
    mkdir -p "$splite_outchrdir"
    
    for chr in {1..49}; do
        mkdir -p "${splite_outchrdir}/chr${chr}"
        chrn="JAYKKS0100000$(printf "%02d" $chr).1"
        
        # Iterate through all BAM files in the input directory

        for input_bam in "$splite_input_dir"/*.bam; do
            if [ -f "$input_bam" ]; then
                # Retrieve the BAM filename (excluding the path)
                bam_file_name=$(basename "$input_bam")
                # Output filename

                output_bam="${splite_outchrdir}/chr${chr}/${bam_file_name%_aln_refgenome_downsampled_new.bam}_${chrn}_d${dp}.bam"
                
                # Using samtools to extract BAM files for specific chromosomes

                samtools view -b "$input_bam" "$chrn" > "$output_bam"
                
                # Create an index for the extracted BAM files

                samtools index "$output_bam"
            fi
        done
    done
done

echo "Chromosome splitting complete."

#6.Hard filtering of VCF files and construction of N_Ref reference panels

vcftools --gzvcf /public/home/user/Site_400inds_hardfiltered.snp.vcf.gz --max-missing 0.9 --maf 0.05 --hwe 1e-6 --minDP 3 --recode --recode-INFO-all --stdout | bgzip -f > /public/home/user/1.projects/9.bam_down/22jc/400/clean_7/400_c_mindp3.vcf.gz

samples=(50 100 150 200 250 300 350 394)

for i in "${samples[@]}"
do
    sample_id_file="/public/home/user/1.projects/9.bam_down/22jc/400/clean_7/400k394_sampeid.txt"
    if [ ! -f "$sample_id_file" ]; then
        echo "Error: Sample ID file does not exist: $sample_id_file"
        continue
    fi
    
    vcftools --gzvcf /public/home/user/1.projects/9.bam_down/22jc/400/clean_7/400_c_mindp3.vcf.gz --keep "$sample_id_file" --recode --recode-INFO-all --stdout | bgzip > "394_k${i}_mindp3.vcf.gz"
    
    # Check whether vcftools has been successfully installed
    if [ $? -ne 0 ]; then
        echo "Error: vcftools processing failed for i=$i"
        continue
    fi
    
    beagle_jar="/public/software/beagle.22Jul22.46e.jar"
    if [ ! -f "$beagle_jar" ]; then
        echo "Error: beagle jar file does not exist: $beagle_jar"
        continue
    fi
    
    java -Xmx40g -jar "$beagle_jar" nthreads=16 gt="394_k${i}_mindp3.vcf.gz" out="394_k${i}_mindp3_beagle_z" impute=true
    
    # Define input files and output folders

    n_input_vcf="394_k${i}_mindp3_beagle_z.vcf.gz"
    n_output_folder="394_k${i}_z/394_k${i}_z"
    
    # Ensure the output folder exists

    mkdir -p "$n_output_folder"
    
    # Cyclic processing of chromosome numbers

     chroms=(1 29 49)
    for chr in "${chroms[@]}"; do

        chr_name="JAYKKS0100000$(printf "%02d" "$chr").1"
        

        if [ ! -f "$n_input_vcf" ]; then
            echo "Error: Input VCF file does not exist: $n_input_vcf"
            continue
        fi
        
        # vcftools extract chromosome reference panel

        vcftools --gzvcf "$n_input_vcf" --chr "$chr_name" --recode --out "${n_output_folder}/vcf_combine_bcftools_${chr_name}.vcf"
        
        if [ $? -eq 0 ] && [ -f "${n_output_folder}/vcf_combine_bcftools_${chr_name}.vcf.recode.vcf" ]; then
            # VCF file extracted from a bgzip-compressed archive

            bgzip -f "${n_output_folder}/vcf_combine_bcftools_${chr_name}.vcf.recode.vcf"
            
            # bcftools processing ACAN values

            bcftools view -c 0 "${n_output_folder}/vcf_combine_bcftools_${chr_name}.vcf.recode.vcf.gz" -Oz -o "${n_output_folder}/vcf_combine_bcftools_${chr_name}_ACAN.vcf.recode.vcf.gz"
            
            # bcftools index file

            bcftools index -f "${n_output_folder}/vcf_combine_bcftools_${chr_name}_ACAN.vcf.recode.vcf.gz"
            
            # Clear intermediate files
            rm -f "${n_output_folder}/vcf_combine_bcftools_${chr_name}.vcf.recode.vcf.gz"
        else
            echo "Warning: vcftools failed to extract chromosome $chr_name"
        fi
    done
done

# 7.GLIMPSE2_chunk Reference Panel Setup

for sample_size in "${samples[@]}"; do
    ind=chunk
    # Define input folder and output folder

    chunk_input_dir="/public/home/user/1.projects/9.bam_down/22jc/400/clean_7_394/all_394_chr"
    chunk_output_dir="/public/home/user/1.projects/9.bam_down/22jc/400/clean_7_394/chunk"
    
    # Ensure the input folder exists

    if [ ! -d "$chunk_input_dir" ]; then
        echo "Error: Input directory does not exist: $chunk_input_dir"
        continue
    fi
    
    # Define the log file path

    log_file="./394_k${sample_size}_${ind}_z.log"
    
    # Ensure the output folder exists

    mkdir -p "$chunk_output_dir"
    
    # Process each file in turn
    for file in "$chunk_input_dir"/vcf_combine_bcftools_JAYKKS0100000*_ACAN.vcf.recode.vcf.gz; do
        if [ -f "$file" ]; then
            # Extract chromosome number

            chromosome=$(basename "$file" | awk -F'JAYKKS0100000' '{print $2}' | awk -F'.' '{printf "%02d", $1}')
            
            # Constructing the output filename

            chunk_output_file="$chunk_output_dir/JAYKKS0100000${chromosome}_ACAN.txt"
            
            # Run the GLIMPSE2_chunk command

            GLIMPSE2_chunk --input "$file" --region "JAYKKS0100000${chromosome}.1" --sequential --window-count 100 --buffer-count 50 --output "$chunk_output_file" &>> "$log_file"
            
            if [ $? -ne 0 ]; then
                echo "Error: GLIMPSE2_chunk execution failed for $file"
            fi
        fi
    done
done

# 8.Establish a link to the target sample

depths=("0.1" "0.5" "1" "2" "2.5" "3" "4" "5" "6")
base_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/norefbin/phase"

for dp in "${depths[@]}"; do
    cd "$base_dir" || { echo "Unable to access directory: $base_dir"; continue; }
    
    mkdir -p "6infile_name_chr_d${dp}"
    
    for chr in {1..49}; do
        chrn="JAYKKS0100000$(printf "%02d" "$chr").1"
        ind_input_dir="/public/home/user/1.projects/9.bam_down/calculate_depth/chrbam/test_d${dp}/chr${chr}"
        ind_output_file="6infile_name_chr_d${dp}/file_names_chr${chr}.txt"
        
        # Check whether the input directory exists

        if [ ! -d "$ind_input_dir" ]; then
            echo "Warning: Input directory does not exist: $ind_input_dir"
            continue
        fi
        
        # Clear or create output files
        > "$ind_output_file"
        
        # Verify whether the sample ID file exists

        id_file="6_id_list.txt"
	#id_file="all_id_list.txt"
        if [ ! -f "$id_file" ]; then
            echo "Error: Sample ID file does not exist: $id_file"
            continue
        fi
        
        # Iterate through all sample IDs within the directory

        while IFS= read -r k || [ -n "$k" ]; do
            # Skip blank lines

            [ -z "$k" ] && continue
            
            # Obtain the full path to the file

            ind_bam_file_path="${ind_input_dir}/${k}_${chrn}_d${dp}.bam"
            
            # Check whether the file exists

            if [ ! -f "$ind_bam_file_path" ]; then
                echo "Warning: BAM file does not exist: $ind_bam_file_path"
                continue
            fi
            
            # Remove specific suffixes from filenames
            modified_name="$k"
            
            # If the filename falls between 10_1 and 10_80, modify the filename prefix to Site.
            if [[ $modified_name =~ ^10_([1-9]|[1-7][0-9]|80)$ ]]; then
                modified_name="Site$modified_name"
            fi
            
            # Add _1 after the filename
            modified_name="${modified_name}_1"
            
            # If the filename is 10_1, prefix the individual with Site.

            if [[ $modified_name == "10_1" ]]; then
                modified_name="Site$modified_name"
            fi
            
            # Append the file path and modified filename to the output file.
            echo "$ind_bam_file_path $modified_name" >> "$ind_output_file"
        done < "$id_file"
        
        echo "The file path and modified name have been written to: $ind_output_file"
    done
done

# 9.GLIMPSE2 filling, linking and evaluation for different sample sizes and target sample depths
for i in {50,100,150,200,250,300,350,394}
do
sign=394_k${i}_z  
mkdir ./${sign}/phase
mkdir ./${sign}/ligate
mkdir ./${sign}/concordance
t_input_vcf=/public/home/user/1.projects/9.bam_down/22jc/400/clean_7/400_c_mindp3.vcf.gz  

 Input directory containing all BAM files
#range1={1..15}
#range2={16..30}
#range3={31..49}
#for i in {1..3}
#do
#{
   # eval "range=\${range${i}}"
    #eval "chr_list=($range)"
#for chr in "${chr_list[@]}"
for chr in {1,29,49}
do
for dp in {0.1,0.5,1,2,2.5,3,4,5,6}
do
{
mkdir ./${sign}/phase/d${dp}
chrn="JAYKKS0100000$(printf "%02d" $chr).1"
chrn2=JAYKKS0100000$(printf "%02d" $chr)

REF=./${sign}/${sign}/vcf_combine_bcftools_${chrn}_ACAN.vcf.recode.vcf.gz
BAMtxt=/public/home/user/1.projects/9.bam_down/calculate_depth/norefbin/phase/6infile_name_chr_d${dp}/file_names_chr${chr}.txt
while IFS="" read -r LINE || [ -n "$LINE" ];
do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	OUT=./${sign}/phase/d${dp}/ref_d${dp}_test_chr${chr}
	GLIMPSE2_phase --bam-list ${BAMtxt} --reference ${REF} --output ${OUT}_${CHR}_${REGS}_${REGE}.vcf --input-region ${IRG} --output-region ${ORG} --keep-monomorphic-ref-Sites --threads 4
                bgzip ${OUT}_${CHR}_${REGS}_${REGE}.vcf
                bcftools index -f ${OUT}_${CHR}_${REGS}_${REGE}.vcf.gz
done < ./${sign}/chunk/${chrn2}_ACAN.txt
#Module ligate
indir=./${sign}/phase/d${dp}
OUT=./${sign}/ligate/${sign}_depth_${dp}_chr${chr}_ligated.vcf
LST=./${sign}/ligate/${sign}_depth_${dp}_list.chr${chr}.txt
ls -1v ${indir}/ref_d${dp}_test_chr${chr}_*.vcf.gz > ${LST}
GLIMPSE2_ligate --input ${LST} --output ${OUT}
bgzip -f ${OUT}
#Module concordance

t_output_vcf="./${sign}/ligate/${sign}_depth_${dp}_chr${chr}_ligated.vcf.gz"
output_concordance="./${sign}/concordance/concordance_${sign}_d${dp}_chr${chr}_bin"
s_id=400k6_id_name.txt
 
# Construct and write to the concordance.txt file
echo "${chrn} ${t_input_vcf} ${t_input_vcf} ${t_output_vcf}" > ./${sign}/concordance/concordance_chr${chr}_d${dp}.txt

# Run the GLIMPSE2_concordance command
GLIMPSE2_concordance --input ./${sign}/concordance/concordance_chr${chr}_d${dp}.txt --sample ${s_id} --output ${output_concordance} --min-val-dp 8 --min-val-gl 0.999 \
  --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.02000 0.05000 0.10000 0.20000 0.30000 0.35000 0.40000 0.45000 0.50000 --thread 4 --af-tag AF 
}&
done
wait
done
done



echo "All processing completed!"
