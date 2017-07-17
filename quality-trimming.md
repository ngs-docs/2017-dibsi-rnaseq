# Short read quality and trimming

Start up a Jetstream m1.medium
[as per Jetstream startup instructions](https://angus.readthedocs.io/en/2017/jetstream/boot.html).

---

You should now be logged into your Jetstream computer!  You should see
something like this

```
titus@js-17-71:~$ 
```

## Installing some software

Run:

```
sudo apt-get -y update && \
sudo apt-get -y install trimmomatic fastqc python-pip \
   samtools zlib1g-dev ncurses-dev python-dev
pip install multiqc   
```
`apt-get install` doesn't work properly for `fastqc`. So we will update the default `fastqc` version using the following commands

```
cd ~/
wget https://launchpad.net/ubuntu/+archive/primary/+files/fastqc_0.11.5+dfsg-3_all.deb && \
sudo dpkg -i fastqc_0.11.5+dfsg-3_all.deb && \
sudo apt-get install -f
```

## Data source

We will be using mRNAseq reads from a small subset of data from [Nematostella vectensis](https://en.wikipedia.org/wiki/Starlet_sea_anemone) [(Tulin et al., 2013)](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16). 

Original RNAseq workflow protocol [here](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/), more updated protocol [here](http://eel-pond.readthedocs.io/en/latest/).

Load your data into ``/mnt/work/data``.  You may need to make the
`/mnt/` directory writeable by doing

```
sudo chmod a+rwxt /mnt
```

First, and then creating the subdirectories

```
cd /mnt
mkdir -p work work/data
cd /mnt/work/data
```

Download subset of data:

```
cd /mnt/work
curl -O https://s3.amazonaws.com/public.ged.msu.edu/mrnaseq-subset.tar
cd data
tar xvf ../mrnaseq-subset.tar
```

Define your $PROJECT variable to be the location of your work
directory; in this case, it will be ``/mnt/work``:

```
export PROJECT=/mnt/work
```

Now load your data in!


## Check that your data is where it should be

Check:

```
ls $PROJECT/data
```

If you see all the files you think you should, good!  Otherwise, debug.

These are FASTQ files -- let's take a look at them:

```
less 0Hour_ATCACG_L002_R1_001.extract.fastq.gz
```

(use the spacebar to scroll down, and type 'q' to exit 'less')

Question:

* where does the filename come from?
* why are there 1 and 2 in the file names?

Links:

* [FASTQ Format](http://en.wikipedia.org/wiki/FASTQ_format)

## Quality trimming and light quality filtering

Make sure you've got the PROJECT location defined, and your data is there:

```
set -u
printf "\nMy raw data is in $PROJECT/data/, and consists of $(ls -1 ${PROJECT}/data/*.fastq.gz | wc -l) files\n\n"
set +u
```
*Important:* If you get an error above or the count of files is wrong...STOP!! Revisit the installation instructions!

### Link your data into your working directory

Change into your project directory and make a workspace for quality trimming:

```  
cd ${PROJECT}
mkdir -p quality
cd quality
```

Now, link the data files into your new workspace

```
ln -s ../data/*.fastq.gz .
```

(Linking with `ln` avoids having to make a copy of the files, which will take up storage space.)

Check to make sure it worked

```
printf "I see $(ls -1 *.fastq.gz | wc -l) files here.\n"
```

You can also do an ``ls`` to list the files.

If you see only one entry, `*.fastq.gz`, then the ln command above didn't work properly.  One possibility is that your files aren't in your data directory; another is that their names don't end with
`.fastq.gz`.

### FastQC


We're going to use 
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
summarize the data. We already installed 'fastqc' on our computer for
you.

Now, run FastQC on two files:

```
fastqc *.fastq.gz
```

After this finishes running (has to run on each file so might take a while), type 'ls':

```
ls -d *fastqc.zip*
```

to list the files, and you should see a number of files with the extensions `.fastqc.zip`.

Inside each of the fatqc directories you will find reports from the fastqc. You can download these files using your RStudio Server console, if you like. To install and run an RStudio Server, go [here](https://angus.readthedocs.io/en/2017/visualizing-blast-scores-with-RStudio.html#installing-and-running-rstudio-on-jetstream). 

Alternatively, we can secure copy (scp) these files to our own laptops, and view them from there.
```
scp username@ip.address:/mnt/work/quality/*html ~/Desktop
```
where the first argument after `scp` is your login and location for files on your jetstream instance, and the second argument is the location to place the files on your own computer.

Questions:

* What should you pay attention to in the FastQC report?
* Which is "better", file 1 or file 2? And why?

Links:

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQC tutorial video](http://www.youtube.com/watch?v=bz93ReOv87Y)
* [Examples of fastqc after technical sequencer problems](http://angus.readthedocs.io/en/2015/_static/2015-lecture2-sequencing.pptx.pdf)(starting on slide 40)

There are several caveats about FastQC - the main one is that it only calculates certain statistics (like duplicated sequences) for subsets
of the data (e.g. duplicate sequences are only analyzed for the first 100,000 sequences in each file.

### Find the right Illumina adapters

You'll need to know which Illumina sequencing adapters were used for your library in order to trim them off. Below, we will use the TruSeq3-PE.fa adapters:

```
wget https://anonscm.debian.org/cgit/debian-med/trimmomatic.git/plain/adapters/TruSeq3-PE.fa
```

Note: If running this on your own data, make sure these are the right adapters for your data.  If they are the right adapters, you should see that some of the reads are trimmed; if they're not, you won't see anything get trimmed.
   
### Adapter trim each pair of files

See excellent paper on trimming from [MacManes 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full).

Run:

```
for filename in *_R1_*.fastq.gz
do
# first, make the base by removing fastq.gz
  base=$(basename $filename .fastq.gz)
  echo $base
        
# now, construct the R2 filename by replacing R1 with R2
  baseR2=${base/_R1_/_R2_}
  echo $baseR2
        
# finally, run Trimmomatic
  TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25
        
# save the orphans
  gzip -9c s1_se s2_se >> orphans.qc.fq.gz
  rm -f s1_se s2_se
done
```

The paired sequences output by this set of commands will be in the files ending in ``.qc.fq.gz``, with any orphaned sequences all together
in ``orphans.qc.fq.gz``.

Make these trimmed reads read-only and keep them, as we will reuse them later.

```
chmod u-w ${PROJECT}/quality/*.qc.fq.gz
```

Questions:

* How do you figure out what the parameters mean?
* How do you figure out what parameters to use?
* What adapters do you use?
* What version of Trimmomatic are we using here? (And FastQC?)
* Do you think parameters are different for RNAseq and genomic data sets?
* What's with these annoyingly long and complicated filenames?
* why are we running R1 and R2 together?

For a discussion of optimal trimming strategies, see 
[MacManes, 2014](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00013/abstract)
-- it's about RNAseq but similar arguments should apply to metagenome
assembly.

Links:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### MultiQc
[MultiQC](http://multiqc.info/) aggregates results across many samples into a single report for easy comparison.

Run Mulitqc on both the untrimmed and trimmed files

```
multiqc .
```

And now you should see output that looks like this:

```
[INFO   ]         multiqc : This is MultiQC v1.0
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
Searching 15 files..  [####################################]  100%
[INFO   ]          fastqc : Found 4 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

You can view output html file
[multiqc_report.html](_static/multiqc_report.html)

Questions:

* is the quality trimmed data "better" than before?
* Does it matter that you still have adapters!?
