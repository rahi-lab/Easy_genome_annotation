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
7. Download Bowtie 2 from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/. The EGA code has been tested with Bowtie 2 version 2.5.4. on Linux using bowtie2-2.5.4-sra-linux-x86_64.zip.
8. Place the Bowtie 2 zip file in the Easy_genome_annotation folder and unzip it.


# Usage

## 1. Creating a whole gene annotation file 


1. Open a terminal.
2. Activate the environment using `conda activate EGA`.
3. Navigate to the directory where you installed EGA using `cd <installation_directory>/EGA`
4. The following example runs the whole genome annotation code on the W303 strain genome sequence file by Liti et al. 2009 with the ORF genome labels for strain S288C, which are supplied with the installation. Adapt the prompt with your path to Bowtie 2, the unannotated genome FASTA file, the FASTA file with sequences to label the genome by similarity, and the name of the output annotation GFF file:
```
python3 whole_genome_annotation.py \
--bowtie2_path \
bowtie2-2.5.4-sra-linux-x86_64/ \
--unannotated_genome_FAFSA_file \
whole_genome_annotation_sample_files/genome_W303_Liti_2009.fsa \
--labeled_sequences_FAFSA_file \
whole_genome_annotation_sample_files/orf_genomic_all.fasta \
--gff_output_file_name \
whole_genome_annotation_sample_files/W303_Liti_2009.gff
```


## 2. Annotating engineered genome modifications


1. Open a terminal.
2. Navigate to the directory where you installed EGA using `cd <installation_directory>/EGA`
3. The following example runs the genome annotation code for engineered modifications for the W303 genome annotation file by Liti et al. 2009 with an annotation file for two plasmids (for tagging _HTB2_ with _mCherry_ and a _CLN2_-promoter-_yEVenus_ construct), which are supplied with the installation. Adapt the prompt with your genome GFF annotation file, your GFF file for the DNA you transformed cells with, and the name of the output annotation GFF file:
```
python3 annotate_genomic_modifications.py \
--orig_genome_gff annotate_genomic_modifications_sample_files/genome_W303_Liti_2009_formatted.gff \
--modifications_gff annotate_genomic_modifications_sample_files/additional_seqs.gff \
--final_genome_gff annotate_genomic_modifications_sample_files/genome_W303_Liti_2009_formatted_updated.gff
```

