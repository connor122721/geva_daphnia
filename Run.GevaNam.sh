#!/usr/bin/env bash
#SBATCH -J geva_nam # A single job name for the array
#SBATCH --ntasks-per-node=15 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/err/geva_nam.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/geva_nam.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load Modules
module load tabix
module load bcftools
module load vcftools

# Link programs
GEVA=/home/csm6hg/geva

# Link inout VCF
wd="/project/berglandlab/connor/for_daniel"
INPUT_VCF=${wd}/daphnia.whatshap.ann.vcf.gz

# Samples
SAMPLEIDS=${wd}/nam.samps
GUIDE=${wd}/goodChrom.txt
rename=${wd}/Scaffold_names.txt

# Performance parameters
CPU=15

# Biological Parameters
REC=1.60e-08 #Lynch et al. -- doi: https://doi.org/10.1101/2020.03.03.974485
MUT=5.69e-09 #Keith et al. -- doi: https://doi.org/10.1101/gr.191338.115
NE=700000

# Generate Internal Variables
CHR=$( cat $GUIDE | sed -n ${SLURM_ARRAY_TASK_ID}p )
chrom=$( cat $rename | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f2 )
echo "Now processing CHR:" ${CHR} ${chrom}
name=${CHR}.nam

# VCFtools portion
vcftools --gzvcf ${INPUT_VCF} \
--chr ${CHR} \
--maf 0.01 \
--keep ${SAMPLEIDS} \
--recode --recode-INFO-all \
--out ${name}

# add annotation
bcftools query \
-f '%CHROM\t%POS\t%POS\t%CHROM\_%POS\n' \
${name}.recode.vcf > ${name}.annotation.txt

#Index the annotation
bgzip ${name}.annotation.txt
tabix -s1 -b2 -e2 ${name}.annotation.txt.gz

# Add annotation
bcftools annotate \
-a ${name}.annotation.txt.gz \
-c CHROM,FROM,TO,ID \
--rename-chrs ${rename} \
${name}.recode.vcf > ${name}.SNPannotated.forGEVA.vcf

# Filter bad GT columns
bcftools view -i 'GT~"/" || GT~"|"' \
${name}.SNPannotated.forGEVA.vcf.gz > ${name}.SNPannotated.forGEVA.filt.GT.vcf
bgzip ${name}.SNPannotated.forGEVA.filt.GT.vcf
tabix ${name}.SNPannotated.forGEVA.filt.GT.vcf.gz

# Phase
module load gcc/9.2.0 shapeit4/4.1.3

# Run shapeit
shapeit4 \
--input ${name}.SNPannotated.forGEVA.filt.GT.vcf.gz \
--region ${chrom} \
-T ${CPU} \
-O ${name}.phased.vcf

# Make bin files
${GEVA}/geva_v1beta \
--vcf ${name}.phased.vcf  \
--rec ${REC} \
-t ${CPU} \
--out ${name}

### Run the TMRCA loop
for i in `cat ${name}.marker.txt \ |
sed '1d' |  \
awk '{print $3}'`
do

#i=10421608
echo ${i}

# Run GEVA
${GEVA}/geva_v1beta \
--input ${name}.bin \
--out RUN_${i} \
--treads ${CPU} \
--mut ${MUT} \
--hmm hmm_initial_probs.txt hmm_emission_probs.txt \
--Ne ${NE} \
--position ${i}

# Make TEMP file
TMP=`sed -n '7p' RUN_${i}.sites.txt`

# Conditionally output
if [ -z "$TMP" ]
then
	echo  " RUN_${i} is empty "
else
	#Print to common file
	echo -e \
	${chrom} \
	${i} \
	${TMP} \
	${SAMPLEIDS} \
	>> ${SAMPLEIDS}.TMRCA.txt
fi

# Remove Clutter inside the loop
rm RUN_${i}.err
rm RUN_${i}.log
rm RUN_${i}.pairs.txt
rm RUN_${i}.sites.txt

# Finish loop
done

# Remove Clutter outside the loop
rm ${name}.recode.vcf
rm ${name}.annotation.txt.gz
rm ${name}.annotation.txt.gz.tbi
rm ${name}.log
rm ${name}.err
rm ${name}.SNPannotated.forGEVA.vcf
rm ${name}.bin
rm ${name}.sample.txt
rm ${name}.marker.txt
