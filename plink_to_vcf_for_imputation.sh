#$ -l h_vmem=46.5G
#$ -l tmem=46.5G
#$ -l h_rt=02:00:00
#$ -j y
#$ -cwd
#$ -S /bin/bash

## Script to convert and check plink fils
filepath=/SAN/neuroscience/PEIC/projects/EEG_genetics_consortium/genetic_data/Elliot
base=MPRC_Hong.QCed

checker=/home/anjabhat/programmes/HRC-1000G-check-bim.pl 

# ### Calc SNP fre for checker script

plink --bfile $filepath/$base --freq --out $base

# ### All checking was done using the script from: http://www.well.ox.ac.uk/~wrayner/tools/
# ## HRC imputation pre checking 
$checker -b $filepath/$base.bim -f $base.frq -r /cluster/project9/PE/user/aritz/pe_imputation/fixdata/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
# ##  Not used - not our first choise of ref panel: 1000G Imputation preparation and checking - .short means only the first 20mil vars are checked
# $checker -b $filepath/$base.bim -f $vcfsteps/01.$base.frq -r /home/molpsych/data/1000G/1000GP_Phase3_combined.legend.short -g -p EUR

# ## Fix issues identified by checker
plink --bfile $filepath/$base --exclude Exclude-$base-HRC.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-$base-HRC.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-$base-HRC.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-$base-HRC.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-$base-HRC.txt --make-bed --out $base-HRC_updated
rm TEMP*

# ### Convert to VCF format -  We use --a2-allele to force the correct REF allele in the vcf file
plink --bfile $base --chr 1-22 --maf 0.01 --recode vcf-iid bgz --a2-allele Force-Allele1-$base-HRC.txt --out $base-HRC_updated

