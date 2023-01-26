#Manipulating VCF files 
#1. Describe the format of the file and the data stored 

# A VCF file (Variant Call Format) is a text file format used to store genetic variation data, such as single nucleotide polymorphisms (SNPs) and insertions/deletions (indels), as well as structural variations, such as copy number variations (CNVs) and breakpoints. The file format is widely used in bioinformatics to store and exchange genetic variation data. The file is structured in a series of rows, where each row represents a single genetic variant and the first several rows contain meta-information. The main data of a VCF file is stored in several columns including CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and sample columns. Additional optional columns can also be included for additional information such as gene annotation, functional prediction, and population frequency data.

#2. What does the header section of the file contain

#The header section of a VCF file contains metadata about the file and the data it contains. This includes information such as the file format version, the source of the data, the reference genome used, and any other relevant information about the variant calls in the file. The header section is denoted by lines starting with "##" and contains a series of key-value pairs, called INFO fields, in the format of "##KEY=VALUE" and in this case it includes.

##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">


#3. How many samples are in the file

bcftools query -l sample.vcf  | wc -l

#Answer: six samples

#4. How many variants are in the file 

bcftools query -f '%ALT\n' sample.vcf | wc -l

#Answer: 398246

# 5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited file



bcftools query -f '%CHROM\t%POS[\t%QD;%MQ]\n' sample.vcf > myoutput1.txt

# 6. Extract data that belongs to chromosomes 2,4 and MT 

 awk '$1=="2" || $1=="4" || $1=="MT"' sample.vcf > myoutput2.vcf

# 7. Print out variants that do not belong to chr20:1-30000000 

awk '$1 != "20" || ($1 == "chr20" && ($2 < 1 || $2 > 30000000)) \
        {print $1, $2, $4, $5}' sample.vcf > myvariants.vcf



# 8. Extract variants that belong to SRR13107019 

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > myvariant2.txt 

# 9. Filter out variants with a QualByDepth above 7 

bcftools filter -i 'INFO/QD>7' sample.vcf > output3.vcf

# 10. How many contigs are referred to in the file. Check the header section 

bcftools view -h sample.vcf | grep -o -w 'contig=[^;]*' | sort | uniq | wc -l

# Answer: 2211

11. Comment on the eighth and ninth columns of the file 

## The eighth and ninth columns of the file refer to the format of the genotype data for each sample. The eighth column refers to the identifier for the sample and the ninth column refers to the genotype data for that sample.

# 12. Extract data on the read depth of called variants for sample SRR13107018 

bcftools query -f '%DP\n' -s SRR13107018 sample.vcf > output4.vcf

# 13. Extract data on the allele frequency of alternate alleles. Combine this data with the chromosome and position of the alternate allele 

 bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf > output5.vcf
 
# Manipulating SAM files 
# 1. Describe the format of the file and the data stored 

#A SAM (Sequence Alignment/Map) file is a file format that stores information about how a set of sequences align to a reference genome. The file is text-based, with each line representing a single align


# 2. What does the header section of the file contain 

# @HD: This record contains general information about the alignment, such as the version of the SAM format used, the sorting order of the alignments, and the groups of read names.
# @SQ: This record describes the reference sequences used in the alignment. It contains the reference sequence name, the length of the reference sequence, and the MD5 checksum of the reference sequence.
# @RG: This record describes the read group. It contains information such as the sample name, the library name, and the sequencing center.
# @PG: This record describes the programs used to generate the alignments. It contains information such as the program name, the version, and the command-line arguments.
# @CO: This record contains any comments or additional information about the alignment

# 3. How many samples are in the file 

grep '^@RG' -c sample.sam | cut -f2
# Answer: 249


# 4. How many alignments are in the file 

grep -v '^@' sample.sam | wc -l 
# Answer: 36142

# 5. Get summary statistics for the alignments in the file 

samtools flagstat sample.sam > mysam1.txt

# 6. Count the number of fields in the file 

head -n1 sample.sam | tr '\t' '\n' | wc -l
# Answer: 4


# 7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_

grep '^@SQ.*NT_' sample.sam

# 8. Print all lines in the file that have @RG and LB tag beginning with Solexa 

grep '^@RG.*LB:Solexa' sample.sam

# 9. Extract primarily aligned sequences and save them in another file 

awk '$1 !~ /^@/ && $2 == "99" || $2 == "83"' sample.sam > my_primarily_aligned1.sam

# 10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM format
 
awk '$1 !~ /^@/ && ($3 == "1" || $3 == "3")' sample.sam | samtools view -Sb - > mysam1.bam

# 11. How would you obtain unmapped reads from the file

samtools view -f 4 sample.sam > mysam2.sam

# 12. How many reads are aligned to chromosome 4

grep -c "^4\t" sample.sam
# Answer: 0

# 13. Comment of the second and sixth column of the file 

## The second column of the file refers to the read name and the sixth column refers to the read flag. The read flag is an integer value that encodes information about the alignment of the read, such as whether it is mapped or unmapped.

# 14. Extract all optional fields of the file and save them in “optional_fields.txt

 awk '{ for (i=11; i<=NF; i++) print $i }' sample.sam > my_optional_fields.txt

