==============================================
 BETA: Binding and Expression Target Analysis
==============================================


Introduction
============

This software is a tool for Transcription Factor's target prediction combine both binding data and expression data.


    
Python Version
==============

Python 2.6 or above is recommended.

Installation
============

::

    $ python setup.py install
    
Command Line
============


Help
----

::

   $ BETA --help

   BETA --- Binding Expression Target Analysis
   BETA [options]* -p <peak> -e <expression> -k <type> -b <boundary> -g <genome> 
   

Main Arguments
--------------


  -p PEAKFILE, --peakfile=PEAKFILE  Input the bed format peak file of the factor

  -e EXPREFILE, --diff_expr=EXPREFILE  Input the differential expression file get from limma for MicroArray data and cuffdiff for RNAseq data

  -k KIND, --kind=KIND  The kind of your expression file,this is required, it can be M or R. M for Microarray. R for RNAseq

  -b BOUNDARYFILE, --bound=BOUNDARYFILE  Input the conserved CTCF binding sites boundary bed format file
                             
  -g GENOME, --genome=GENOME  Select a genome file (sqlite3 file) to search refGenes.

			      
Options
-------

  --version             Show program's version number and exit
  
  -h, --help            Show this help message and exit.
                              
  --pn=PEAKNUMBER   The number of peaks you want to consider, DEFAULT=10000
                            
  -n NAME, --name=NAME  This argument is used to name the result file.If not set, the peakfile name will be used instead
                            
                             
  -d DISTANCE, --distance=DISTANCE  Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. default:100000(100kb)
                             
  --df=DIFF_FDR  Input a number 0~1 as a threshold to pick out the most significant differential expressed genes by FDR,
                 DEFAULT = 1, that is select all the genes
                            
  --da=DIFF_AMOUNT          Input a number between 0-1, so that the script will only output a percentage of most significant differential
                            expressed genes,input a number bigger than 1, for example, 2000. so that the script will only output top 2000 
                            genes DEFAULT = 0.5, that is select top 25 percentage,NOTE:If you want to use diff_fdr, please set this parameter
                            to 1, otherwose it will get the intersection of these two parameters
                            
  -c CUTOFF, --cutoff=CUTOFF  Input a number between 0~1 as a threshold to select the closer target gene list(up regulate or down regulate or both) 
                              with the p value was called by one side ks-test, DEFAULT = 0.001
                           
  --pt=PERMUTETIMES        Permutaton times,give a resonable value to get an exact FDR.Gene number and permute times decide the time it 
                           will take. DEFAULT=500    


Example
-------

::

   BETA -p 2723_peaks.bed -e gene_exp.diff -b hg19_CTCF_bound.bed -k R -g hg19.refseq

   
   
Input Files Format
==================

- ``Peak`` : BED format 

    ``chroms``  ``start``  ``end``  ``name``  ``score``  ``[strand]``
    
    If your bed don't have the name and score column, please fake one.

- ``Expression by Microarray`` : Result of Limma 

    ``ID``  ``Refseq``  ``logFC``  ``AveExpre``  ``Tscore``  ``Pvalue``  ``adj.P.Value``  ``B``

- ``Expression by RNAseq`` : Result of Cufflinks

    ``Test_id``  ``gene_id``  ``gene``  ``locus``  ``sample_1``  ``sample_2``  ``status``  ``value_1``  ``value_2``  ``Log2(foldchange)``  ``test_stat``  ``p_value``  ``q_value``  ``significant``

- ``CTCF conserved boundary`` : BED format

    ``chroms``  ``start``  ``end``  ``name``  ``score``  ``[strand]``
    
    The conserve CTCF binding sites of all the cell lines.

- ``Genome reference`` ; Downloaded from UCSC

    ``refseqID``  ``chroms``  ``strand``  ``txstart``  ``txend``  ``genesymbol``.
    
    We use that as a reference to get the gene information.


    
Output Files
============


- ``score.pdf`` : A CDF figure to test the TF's funtion, Up pr Down regulation.
- ``score.r`` : The R script to draw the ``score.pdf`` figure
- ``uptarget.txt`` : The uptarget genes, 4 column, Refseq, Gene Symbol, Rank Product, FDR
- ``downtarget.txt`` : The downregulate genes, the same format to uptarget.
    
**NOTE**: Up or Down target file depends on the test result in the PDF file, it will be not produced enless it passed the threshold you seted via -c --cutoff
    

    
