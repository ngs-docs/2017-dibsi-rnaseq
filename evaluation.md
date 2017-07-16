# Evaluating your transcriptome assembly

We will be using Transrate and Busco!

## Transrate

[Transrate](http://hibberdlab.com/transrate/getting_started.html) serves two main purposes. It can compare two assemblies to see how similar they are. Or, it can give you a score which represents proportion of input reads that provide positive support for the assembly. We will use transrate to get a score for the assembly. Use the trimmed reads. For a further explanation of metrics and how to run the reference-based transrate, see the [documentation](http://hibberdlab.com/transrate/metrics.html) and the paper by [Smith-Unna et al. 2016](http://genome.cshlp.org/content/early/2016/06/01/gr.196469.115). 

Make a new directory and get the reads together:

```
cd ${PROJECT}
mkdir -p evaluation
cd evaluation

cat ${PROJECT}/quality/*R1*.qc.fq.gz > left.fq.gz
cat ${PROJECT}/quality/*R2*.qc.fq.gz > right.fq.gz
```

Transrate doesn't like pipes in sequence names. This version of Trinity doesn't output pipes into the sequence names, but others do. Let's just fix to make sure.

```
sed 's_|_-_g' ${PROJECT}/assembly/trinity_out_dir/Trinity.fasta > Trinity.fixed.fasta
```

Now, run the actual command:

```
transrate --assembly=Trinity.fixed.fasta --threads=2 \
--left=left.fq.gz \
--right=right.fq.gz \
--output=${PROJECT}/evaluation/nema
```

## BUSCO

* Eukaryota database used with 429 genes
* "Complete" lengths are within two standard deviations of the BUSCO group mean length
* Website: http://busco.ezlab.org/
* Simho et al. 2015: http://bioinformatics.oxfordjournals.org/content/31/19/3210
* http://gitlab.com/ezlab/busco/raw/master/BUSCO_v2.0_userguide.pdf

Run the actual command:

```
BUSCO.py \
-i Trinity.fixed.fasta \
-o nema_busco_metazoa -l /home/ubuntu/busco/metazoa_odb9 \
-m tran --cpu 2
```

Check the output:

```
cat run_nema_busco_metazoa/short_summary_nema_busco_metazoa.txt
```
