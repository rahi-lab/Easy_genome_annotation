# Easy genome annotation (EGA)
This tool consists of two modules:

1. Whole genome annotation: Given an unannotated genome in FASTA format, this module creates a GFF annotation file based on sequence similarity with labeled sequences in FASTA format.
2. Annotation of genomic modifications: Given a GFF genome annotation and a GFF file describing transforming DNA sequences, this module updates the genomic annotations by simulating homologous recombination of the transformant DNA.

# Installation

1. If you do not have conda or miniconda installed, download it from https://docs.conda.io/en/latest/miniconda.html.
2. Clone this repository (`git clone https://github.com/rahi-lab/Easy_genome_annotation`).
3. In the command line, navigate to the folder where you cloned EGA (command `cd Easy_genome_annotation`).
4. In the command line, create a virtual environment with python 3.9 with the command `conda create -n EGA python=3.8.5`.
5. Activate the environment using `conda activate EGA`. 
6. Install the Python package pysam `pip install pysam`.
7. Download Bowtie 2 from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/. This code has been tested with version 2.5.4. on Linux using bowtie2-2.5.4-sra-linux-x86_64.zip.
8. Place the Bowtie 2 zip file in the Easy_genome_annotation folder and unzip it.


# Usage

## 1. Creating a whole gene annotation file 


1. Open a terminal.
2. Activate the environment using `conda activate EGA`.
3. Navigate to the directory where you installed EGA using `cd <installation_directory>/EGA`
4. Adapt the following example with your path to Bowtie 2, the unannotated genome FASTA file, the FASTA file with sequences to label by similarity, and the name of the output file:
```
python3 whole_genome_annotation.py \
--bowtie2_path bowtie2-2.5.4-sra-linux-x86_64/ \
--unannotated_genome_FAFSA_file whole_genome_annotation_sample_files/genome_W303_Liti_2009.fsa \
--labeled_sequences_FAFSA_file whole_genome_annotation_sample_files/orf_genomic_all.fasta \
--gff_output_file_name whole_genome_annotation_sample_files/W303_Liti_2009.gff
```


### Process description:
1. Obtain the gene ORF sequences of a closely related strain in a suitable standardized format, such as FASTA. Secondly, have the full genome sequence of the reference genome of interest. <br/>
Steps 2-4 can be run by executing 'custom_annotations.loop', or by executing the next steps individually: 
3. Index the genome: 
```bash
$ bowtie2-build <genome.fasta> <genome_index>
```
3. Preform the alignment of the gene sequences to the reference genome specified by the Bowtie2 index files: 
```bash
$ bowtie2 -x <genome_index> -U <gene_sequences.fsa> -S <output.sam>
```
4. Visualize and Analyze: <br>
Using pre-installed SAMtools, you can view the resulting sorted BAM file in a genomic viewer such as IGV, or proceed with downstream analysis.
```bash
$ samtools view -bS <output.sam> | samtools sort -o output.sorted.bam 
```
5. Convert SAM format to GFF format file: <br>
Using 'sam_to_gff.py' script and precising the SAM file and the desired file name for GFF file. Script can be run with :
```bash
$ python sam_to_gff.py <output.sam> <gff_file.gff>
```
This part will generate a gff file with listed annotations, followed by FASTA formated sequqneces, both parts separated by ##FASTA. This format is used for the next steps in the pipeline. 
Please replace placeholders with actual values and descriptions relevant to your script.
#### Options

| Argument         | Description                                                | Default Value |
|------------------|------------------------------------------------------------|---------------|
| `genome_index`     | Prefix or base name for the index files                 | -             |
| `gene_sequences.fsa`     | Raw sequencing reads used as inputs            | -             |
| `output.sam`     | The output file in SAM format           | -             |

6. (OPTIONAL) Add promoter and terminator annotations: <br\>
Using 'Adding_promoters_to_gff.py', there can be added additional annotations for promoters and terminators. They are expressed as additional annotations around genes, each being 500 bp before of after the CDS as a conservative estimation of their length. This can be beneficial for the later steps as some mutations fall in these regions and therefore aren't detected if promoters and terminators are not specified. Promtoters and termiantors are distinguished by their type ("promoter" or "terminator" instead of "gene") and in  have added -pr or -ter ending in the "name" feature where gene name is specified otherwise. <br\>
8. (OPTIONAL) Evaluate the alignment: <br>
With 'compare_gff.py' script, you can evaluate the alignment. The script reports the number of gene sequences that were successfully aligned. It also reports how many are sensible. This means, obtained by extracting the sequence between the reported start and stop position of a specific gene form the reference FASTA file, how many start with the start codon, ends with a stop codon and have no in-frame stop codons.
9. (OPTIONAL) Remove non-aligned instances from annotations file: <br>
During the alignment, some genes might not align to any of the sequences in the genome. These instances are still put in the annotations file, with chromosome name denoted as 'None' and with a meaningless assign of start and stop position. They can be removed from the annotations file with the script *Removing_None_genes.py*.

## 2. Adding custom annotations to the gff file 
This part aims to add additional annotations to the already existing annotations file (*WT_genome_annotations.gff*). It sources annotations from a GFF formated file with FASTA formated sequences following the annotations part (*additional_annotations.gff*). New gff file, containing all annotations, is saved as *Synthetic_genome_annotations.gff *. 
The script can be run as: 
```bash
$ pyhon amend_gff_file.py <WT_genome_annotations.gff> <additional_annotations.gff> <Synthetic_genome_annotations.gff> 
```
#### Replacements

| Argument         | Description                                                | Default Value |
|------------------|------------------------------------------------------------|---------------|
| `WT_genome_annotations.gff`     | Existing annotations file              | -             |
| `additional_annotations.gff`     | Addtitional annptations, see description for file format            | -             |
| `Synthetic_genome_annotations.gff`     | Name of the new gff file created           | -             |


Note: If visualization in IGV or similair software is desired, the script 'split_into_genome_and_gff.py' can be used to separate the gff file containing both annotations and sequences in FASTA format into two separate files, one containing annotations and the other sequences. The script can be run with: 
```bash
$ pyhon split_into_genome_and_gff.py <gff_with_both_annotations_and_sequences.gff>
```
This script will generate two before mentioned files, in gff and .fsa format, named as '{gff_with_both_annotations_and_sequences}_annotations.gff' or '{gff_with_both_annotations_and_sequences}_sequences.fsa'. 


## 3. Alignment to the reference genome
This part can be executed by running the *process_fq_and_make_calls.loop* script from the command line. The user needs to provide sequencing data in the FASTQ format. There can be multiple file referring to the same genome, saved in a unique folder with a distinct root name (such as with a reference genome name and variant identification name). 

Script aligns sequences to the reference genome to obtain alignment files (SAM/BAM) for downstream analysis and converts it to a more readable format (VCF). 
Program structure is similar to the script under the point *1 - creating a custom gene annotations file*
The script has to be modified before execution - placeholders need to be replaced with suitable variables. 

To see genes with higher copy numbers, 'above_avg_coverage.py' can be used. It takes in a .bam format file from sequenced strains and the .gff file. Then the user can see a list of genes with the average number of reads being at least 1.5 times higher than the average of the genome. The script takes in 2 inputs: first the .bam file and second the gff file. 

#### Replacements

| Argument         | Description                                                | Default Value |
|------------------|------------------------------------------------------------|---------------|
| `reference_genome/file_root_name`     | Directory name(s) followed by the root name(s) of sequenced data files in .fq format             | -             |
| `sequence.fsa`     | Name of the file containing reference genome            | -             |
| `index`     | Index name           | -             |
| `30`     | Minimal quality score that is not filtered out from the output          | 30             |
| `10`     | Minimal depth of read that is not filtered out from the output          | 10            |

### Output
This part generates 2 vcf files containing all differences from the reference genome to which it was aligned. Filenames end with '__calls.vcf' and '__filtered_calls_no_lowDP.vcf', the later containing only filtered results from the former file.



## 4. Vcf processing
### Description 
Filters VCF files to extract mutations meeting predefined criteria, such as quality scores, read depth, IMF value.
Compares mutations between variants to identify unique and shared variants between different genomes.
Counts the occurrences of mutations across different variants.
Summarizes mutation data, considering variations in gene sequences.
<br/>

This part is executed by calling the *main.py* script as: 

```bash
$ python main.py <gff_file.gff> <ancestor.vcf> <low_qual> <low_imf>
```
### Inputs and Prerequisites:
Annotated VCF files for each variant and an annotated VCF file for the ancestral variant (wild type).
Ensure compatibility with the VCF file format (##fileformat=VCFv4.2).

Please replace placeholders with actual values and descriptions relevant to your script.
#### Options

| Argument         | Description                                                | Default Value |
|------------------|------------------------------------------------------------|---------------|
| `gff_file`     | Path to the gff file, with both annotations and FASTA parts              | -             |
| `ancesor`     | Ancestor file from part 3 in vcf format           | -             |
| `low_qual`     | Lower bound for quality of alignment at the position of mutation (filtering lower quality)       | 15           |
| `low_imf`     | Lower bound for the Inbreeding Metric Filter if imf value is defined for the mutation          | 0.15           |

### Outputs
The script generates an ./Comparison_results/ that gathers a list of mutations occurring in each lab grown variant, and a joined file called 'Comparison_results.txt' that lists mutations occurring multiple times. 
<p>
'variantname_calls_wo_ancestor_processed.txt' file is unique for each processed variant. It lists all mutations not detected in ancestral variant. The format for mutations occuring in genes is the following: <br/>
 
 *DP=x;VDB=x;SGB=x;MQSB=x;MQ0F=x;AC=x;AN=x;DP4=x,x,x,x;MQ=x* <br>
 *&gt;gene_name, :  (position_inside_the_gene) position_in_chromosome original_amino_acid -> mutated_amino_acid* <br>
  
  where the first line contains information on the quality of alignment at the position of the mutation and the second line describes the nature of mutation. Positions are expressed in bp and refer to the starting nucleotide (counting from 0) of the amino acid that was affected by the mutation. <br\>
If mutation occurs inside a promoter or terminator, then only the second line is writen. <br\>
</p>

<p>
'Comparison_results.txt' file has a specific format: <br\>
 
 *&gt;gene_name occurs x times: list of variants' names where the mutation occured on that gene separated by comma* <br/>
 
 The mutations are grouped by the gene they occur in and sorted by frequency of appearance. 
</p>



