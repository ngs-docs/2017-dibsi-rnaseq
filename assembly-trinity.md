# De novo transcriptome assembly with Trinity

This tutorial will use mRNAseq reads from a small subset of data from [Nematostella vectensis](https://en.wikipedia.org/wiki/Starlet_sea_anemone) [(Tulin et al., 2013)](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16). 

Original RNAseq workflow protocol [here](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/), more updated protocol [here](http://eel-pond.readthedocs.io/en/latest/).

## Installation

On a Jetstream instance, run the following commands to update the base
software:

```
sudo apt-get update && \
sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
  default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
  python-matplotlib python-pip python-virtualenv sysstat fastqc \
  trimmomatic bowtie samtools blast2 wget bowtie2 openjdk-8-jre \
  hmmer ruby
```

Install [khmer](http://khmer.readthedocs.org) from its source code.

```
cd ~/
python2.7 -m virtualenv pondenv
source pondenv/bin/activate
cd pondenv
pip install -U setuptools
git clone --branch v2.0 https://github.com/dib-lab/khmer.git
cd khmer
make install
```

The use of `virtualenv` allows us to install Python software without having
root access. If you come back to this protocol in a different terminal session
you will need to run

```
source ~/pondenv/bin/activate
```

Install Trinity:

```
cd ${HOME}

wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz \
    -O trinity.tar.gz
tar xzf trinity.tar.gz
cd trinityrnaseq*/
make |& tee trinity-build.log
```

Assuming it succeeds, modify the path appropriately in your virtualenv
activation setup:

```
echo export PATH=$PATH:$(pwd) >> ~/pondenv/bin/activate
source ~/pondenv/bin/activate
```

You will also need to set the default Java version to 1.8

```
sudo update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java
```

## Then, let's check we still have our reads from yesterday's QC lesson
```
set -u
printf "\nMy trimmed data is in $PROJECT/quality/, and consists of $(ls -1 ${PROJECT}/quality/*.qc.fq.gz | wc -l) files\n\n"
set +u

```
where set -u should let you know if you have any unset variables, i.e. if the `$PROJECT` variable is not defined. 

If you see `-bash: PROJECT: unbound variable`, then you need to set the $PROJECT variable.  
```
export PROJECT=/mnt/work
```
and then re-run the `printf` code block.

NOTE: if you do not have files, please rerun quality trimming steps [here](quality-trimming.html)


## Let's do digital normalization with khmer

### First, interleave the sequences

Next, we need to take these R1 and R2 sequences and convert them into
interleaved form, for the next step.  To do this, we'll use scripts
from the `khmer package <http://khmer.readthedocs.org>`__, which we
installed above.

Now let's use a for loop again - you might notice this is only a minor
modification of the previous for loop...

```
cd ${PROJECT}/quality

for filename in *_R1_*.qc.fq.gz
do
# first, make the base by removing .extract.fastq.gz
  base=$(basename $filename .qc.fq.gz)
  echo $base

# now, construct the R2 filename by replacing R1 with R2
baseR2=${base/_R1/_R2}
echo $baseR2

# construct the output filename
  output=${base/_R1/}.pe.qc.fq.gz

  (interleave-reads.py ${base}.qc.fq.gz ${baseR2}.qc.fq.gz | \
  gzip > $output)
done
```

The final product of this is now a set of files named
`*.pe.qc.fq.gz` that are paired-end / interleaved and quality
filtered sequences, together with the file `orphans.qc.fq.gz` that
contains orphaned sequences.

Make the end product files read-only!

```
chmod u-w *.pe.qc.fq.gz orphans.qc.fq.gz
```

to make sure you don't accidentally delete them.

Since you linked your original data files into the ``quality`` directory, you
can now do:

```
rm *.fastq.gz
```

to remove them from this location; you don't need them for any future steps.

### Things to think about:

Note that the filenames, while ugly, are conveniently structured with the
history of what you've done to them.  This is a good strategy to keep
in mind.

## Applying Digital Normalization

In this section, we'll apply [digital normalization](http://arxiv.org/abs/1203.4802) and [variable-coverage k-mer
abundance trimming](https://peerj.com/preprints/890/) to the reads prior to assembly.  This has the effect of reducing the computational
cost of assembly [without negatively affecting the quality of the
assembly](https://peerj.com/preprints/505/).


```   
cd ${PROJECT}
mkdir -p diginorm
cd diginorm
ln -s ../quality/*.qc.fq.gz .
```

Apply digital normalization to the paired-end reads

```
normalize-by-median.py -p -k 20 -C 20 -M 4e9 \
  --savegraph normC20k20.ct -u orphans.qc.fq.gz \
  *.pe.qc.fq.gz
```

Note the `-p` in the normalize-by-median command -- when run on PE data, that ensures that no paired ends are orphaned.  The `-u` tells
noralize-by-median that the following filename is unpaired.

Also note the `-M` parameter.  This specifies how much memory diginorm should use, and should be less than the total memory on the computer
you're using. (See [choosing hash sizes for khmer](http://khmer.readthedocs.org/en/latest/choosing-hash-sizes.html)
for more information.)

### Trim off likely erroneous k-mers

Now, run through all the reads and trim off low-abundance parts of
high-coverage reads

```
filter-abund.py -V -Z 18 normC20k20.ct *.keep && \
  rm *.keep normC20k20.ct
```

This will turn some reads into orphans when their partner read is
removed by the trimming.

### Rename files

You'll have a bunch of `keep.abundfilt` files. Let's make things prettier:

First, let's break out the orphaned and still-paired reads:

```
for file in *.pe.*.abundfilt
do 
  extract-paired-reads.py ${file} && \
  rm ${file}
done
```

We can combine all of the orphaned reads into a single file

```
gzip -9c orphans.qc.fq.gz.keep.abundfilt > orphans.keep.abundfilt.fq.gz && \
  rm orphans.qc.fq.gz.keep.abundfilt
for file in *.pe.*.abundfilt.se
do
  gzip -9c ${file} >> orphans.keep.abundfilt.fq.gz && \
  rm ${file}
done
```

We can also rename the remaining PE reads & compress those files

```
for file in *.abundfilt.pe
do
  newfile=${file%%.fq.gz.keep.abundfilt.pe}.keep.abundfilt.fq
  mv ${file} ${newfile}
  gzip ${newfile}
done
```

This leaves you with a bunch of files named `*.keep.abundfilt.fq.gz`, which represent the paired-end/interleaved reads that remain after
both digital normalization and error trimming, together with `orphans.keep.abundfilt.fq.gz`.


## Running the Actual Assembly!

Let's make another working directory for the assembly

```
cd ${PROJECT}
mkdir -p assembly
cd assembly
```

For paired-end data, Trinity expects two files, 'left' and 'right'; there can be orphan sequences present, however.  So, below, we split
all of our interleaved pair files in two, and then add the single-ended seqs to one of 'em. :

```
for file in ../diginorm/*.pe.qc.keep.abundfilt.fq.gz
do
  split-paired-reads.py ${file}
done
   
cat *.1 > left.fq
cat *.2 > right.fq
   
gunzip -c ../diginorm/orphans.keep.abundfilt.fq.gz >> left.fq
```

### Assembling with Trinity

Here is the assembly command!

```
Trinity --left left.fq \
  --right right.fq --seqType fq --max_memory 14G \
  --CPU 2
```

Note that these last two parts (`--max_memory 14G --CPU 2`) configure the maximum amount of memory and CPUs to
use.  You can increase (or decrease) them based on what machines you are running on.

Once this completes, you'll have an assembled transcriptome in
`${PROJECT}/assembly/trinity_out_dir/Trinity.fasta`.

