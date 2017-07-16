# DIBSI Nonmodel mRNASeq Workshop (2017) 

These are the schedule and classroom materials for the
[2017 DIBSI nonmodel mRNAseq workshop at UC Davis](http://dibsi-rnaseq.readthedocs.io/en/latest/),
which will run from July 17th to July 21st, 2017.

This workshop runs under a [Code of Conduct](code-of-conduct.html). Please
respect it and be excellent to each other!

If you're not comfortable working on the command line, please work through some of this [command-line bootcamp](http://rik.smith-unna.com/command_line_bootcamp/) before the workshop.

Twitter hash tag: [#ngs2017](https://twitter.com/search?f=tweets&q=%23ngs2017&src=typd)

**Schedule and Location**:  

9am-3pm + Tu/Th evenings (7-8:30pm)  

All sessions are in **Valley Hall**, unless otherwise noted.

## Workshop materials
*in progress*

work in:
* Lecture:  RNA-Seq uses & pitfalls


### Monday, Day 1: Introduction and QC

* 11am: Introductions
* 1:15pm: All-hands gathering to introduce across all workshops 
* 3pm: Intro & setup for nonmodel mRNASeq workshop
* Hands-on:
   * [Booting a cloud computer from Jetstream](jetstream/boot.html)
   * [Quality trimming your reads](quality-trimming.html)

optional extras:  

* Hands-on: [Review and explore: Command line UNIX, and R/RStudio](command-line-and-rstudio.html)


### Tuesday, Day 2: Assembly and Evaluation

Morning: 9am-12pm

* Lecture: Kmers, de bruijn graphs, diginorm, and assembly
* Hands-on: [De novo RNAseq assembly](assembly-trinity.html) 

Afternoon: 1:15pm - 3pm.  

* Lecture: Assembly evaluation
* Hands-on: [Assembly evaluation](assembly-evaluation.html)

Evening: 7pm-9pm: student presentations and questions! (social)

### Wednesday, Day 3: Annotation and Quantification

Morning 9am-12pm
 
* Lecture: Annotation
 * Hands-on: [Transcriptome annotation](dammit_annotation.html) 


Afternoon: 1:15pm - 4pm 

* Lecture: Quantification
*  Hands-on: [Quantification with Salmon](quantification.html)


Evening: free time / social [Wed Farmers' market!](http://www.davisfarmersmarket.org/))

### Thursday, Day 4: Differential Expression and Downstream Assessment

Morning 9am-12pm
 
* Lecture: Differential Expression 
* Hands-on: [Differential Expression with DESeq2](deseq2-asthma.html) 

Afternoon: 1:15pm - 3pm.  

* Lecture: RNA-Seq Study design
* Hands-on (tbd):
  * (maybe) [Pathway Analysis](pathway_analysis.html) 

Evening 7pm-9pm

* Hands-on: own data!
* Alternative: [Reference independent analyses with k-mers; sourmash.](kmers-and-sourmash.html)


### Friday, Day 5:  Automation and repeatability

Morning 9am-12pm  

9am: All-hands wrap-up 

Since checkout is at 12, everything is optional after the hands-on meeting. Depending on interest, One or more of the instructors can stick around to teach or help you with your own data.

Some options:  

 *  [Introduction to automation](introduction-to-automation.html)
 *  [GitHub](github.html)
 *  [Jupyter Notebook, R and Python for data science.](jupyter-notebook-demo/Jupyter-Notebook-Notes.html)
 * Hands-on: own data!

  

  
  
  
  
  
  

#### Useful References:  

* [Jupyter Notebook, R and Python for data science.](jupyter-notebook-demo/Jupyter-Notebook-Notes.html)
* [GitHub](github.html)
* [where do I find the data? NCBI, ENSEMBL, ENA; how to get FASTQ out of NCBI.](database_resources.html)
* [RMarkdown](rmarkdown_rnaseq.html)
*  [Adrienne Roeder](http://roeder.wicmb.cornell.edu/), Cornell - [Reaching biological conclusions from RNA-seq: the good, the bad, and the ugly](https://osf.io/qz3m6/)
*  [Michael I Love](https://mikelove.github.io/), UNC Chapel Hill - ["Statistics and bias correction in RNAseq differential expression analysis"](https://osf.io/gbjhn/)
*  [Robert Patro](http://www.robpatro.com/redesign/), Stony Brook University - ["Don't count on it: Pragmatic and theoretical concerns and best practices for mapping and quantifying RNA-seq data"](https://osf.io/bv85u/)
*  C. Titus Brown, UC Davis - ["Effectively infinite: next steps in Data Intensive Biology."](https://osf.io/pbmeh/)
* [Assessing & assembling nanopore data](analyzing_nanopore_data.html) (Lisa Cohen and Jon Badalamenti)


*  Hands-on: [Counting](counting.html)
