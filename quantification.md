# Analyzing RNAseq expression with Salmon

Tutorial based on the Eel-Pond RNAseq workflow protocol [here](http://eel-pond.readthedocs.io/en/latest/).

## Installation

If needed, start up a fresh Ubuntu 16.04 instance on Jetstream (m1.medium) using the boot instructions [here](jetstream/boot.html).

On a Jetstream instance, run the following commands to update the base
software:


We will use `salmon <http://salmon.readthedocs.org/en/latest/>`__ to
quantify expression. `Salmon
<https://github.com/COMBINE-lab/salmon>`__ is a new breed of software
for quantifying RNAseq reads that is both really fast and takes
transcript length into consideration (`Patro et al. 2015
<http://dx.doi.org/10.1038/nmeth.4197>`__).

https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/salmon.rst

http://angus.readthedocs.io/en/2016/rob_quant/tut.html

https://2016-aug-nonmodel-rnaseq.readthedocs.io/en/latest/quantification.html

Be sure you have loaded the right Python packages
::

   source ~/pondenv/bin/activate

Run Salmon
==========

First, build an index for your new transcriptome:
::

   cd ${PROJECT}
   mkdir -p quant
   cd quant
   ln -s ${PROJECT}/assembly/trinity_out_dir/Trinity.fasta .
   salmon index --index nema --transcripts Trinity.fasta --type quasi

And also link in the QC reads (produced in :doc:`1-quality`):
::

   ln -s ../quality/*R1*.qc.fq.gz .
   ln -s ../quality/*R2*.qc.fq.gz .

Then, run the salmon command:
::

  for R1 in *R1*.qc.fq.gz
  do
    sample=$(basename $R1 extract.qc.fq.gz)
    echo sample is $sample, R1 is $R1
    R2=${R1/R1/R2}
    echo R2 is $R2
    salmon quant -i nema -p 2 -l IU -1 <(gunzip -c $R1) -2 <(gunzip -c $R2) -o ${sample}quant
  done

This will create a bunch of directories named something like
``0Hour_ATCACG_L002001.quant``, containing a bunch of files. Take a
look at what files there are:
::

    find 0Hour_ATCACG_L002_R1_001* -type f

You should see::

    0Hour_ATCACG_L002_R1_001.extract.quant/quant.sf
    0Hour_ATCACG_L002_R1_001.extract.quant/aux_info/observed_bias.gz
    0Hour_ATCACG_L002_R1_001.extract.quant/aux_info/observed_bias_3p.gz
    0Hour_ATCACG_L002_R1_001.extract.quant/aux_info/fld.gz
    0Hour_ATCACG_L002_R1_001.extract.quant/aux_info/expected_bias.gz
    0Hour_ATCACG_L002_R1_001.extract.quant/aux_info/meta_info.json
    0Hour_ATCACG_L002_R1_001.extract.quant/cmd_info.json
    0Hour_ATCACG_L002_R1_001.extract.quant/logs/salmon_quant.log
    0Hour_ATCACG_L002_R1_001.extract.quant/libParams/flenDist.txt
    0Hour_ATCACG_L002_R1_001.extract.quant/lib_format_counts.json
    0Hour_ATCACG_L002_R1_001.extract.quant.counts

The two most interesting files are ``salmon_quant.log`` and
``quant.sf``. The latter contains the counts; the former contains the
log information from running things.

Working with the counts
-----------------------

The ``quant.sf`` files actually contain the relevant information about
expression -- take a look::

   head -20 0Hour_ATCACG_L002_R1_001.extract.quant/quant.sf

The first column contains the transcript names, and the
fifth column is what edgeR etc will want - the "raw counts".
However, they're not in a convenient location / format for edgeR to use;
let's fix that.

Now, grab the script...

::

   curl -L -O https://raw.githubusercontent.com/dib-lab/eel-pond/master/gather-counts.py

and run it::

   python ./gather-counts.py

This will give you a bunch of .counts files, processed from the quant.sf files
and named for the directory they are in.
