# Cancer Genomics Assignment 1
### Prerequisities
 - unix system
    - bowtie2 
    - samtools
 - python
    - numpy
    - matplotlib

### Data
Genome data can be downloaded from https://gear.embl.de/data/.exercise/. We need to download the reference human genome of chromosome X, in this case Genome Reference Consortium Human Build 37 (GRCh37):

```
wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz
```

### Alignment to reference
Build index from the reference genome into files with name prefix human_ref:

```
bowtie2-build Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz human_ref
```

Now, we can align the tumor and normal genome to the reference into SAM files with corresponding names:

```
bowtie2 -x human_ref -1 wt.r1.fq.gz -2 wt.r2.fq.gz -S wt.sam -p 8
bowtie2 -x human_ref -1 tu.r1.fq.gz -2 tu.r2.fq.gz -S tu.sam -p 8
```

### Processing of the aligned files

Afterwards, we convert the SAM files to BAM files, sort and index the aligned files. As well as, trim the datasets to the region of interest - chromosome X from 20Mbp to 40Mbp (GRCh37/hg19 coordinates):
```
samtools view -u wt.sam | samtools sort -o wt.sorted.bam
samtools index wt.sorted.bam
samtools view -b wt.sorted.bam X:20000000-40000000 > wt.sorted.trimmed.bam
samtools index wt.sorted.trimmed.bam
```

```
samtools view -u tu.sam | samtools sort -o tu.sorted.bam
samtools index tu.sorted.bam
samtools view -b tu.sorted.bam X:20000000-40000000 > tu.sorted.trimmed.bam
samtools index tu.sorted.trimmed.bam
```

### Calculation of read-depth
Read depth can be calculated by using 'samtools', '-a' reports null depth positions and '-r' restricts the calculation to a region of interest:

```
 samtools depth -a -r X:20000000-40000000 wt.sorted.trimmed.bam -o wt.depth.data
 samtools depth -a -r X:20000000-40000000 tu.sorted.trimmed.bam -o tu.depth.data
```

### Plot the results
To plot the read-depth of the wild-type and tumor genome, and their log2 ratio use the following python script with parameters:

```
python depth_plot.py --wt_file_name wt.depth.data --tu_file_name tu.depth.data --window_size 10000
```

The resulting graphs are saved as PNG files and can be found in the folder **results**.
