# FastTools: Python module with tools useful for reading and manipulating FASTQ/FASTA files.

## Currently only supporting Python 3.6

FastTools is a module written with the intention of making working with FASTQ/FASTA files more convenient for Python users.
This module gives the user high-level control over their NGS data, and uses a Data Science orientation to work with sequence
data.

## Setup
1. Clone or download FastTools repository.

2. Set up bioconda by following instructions at https://bioconda.github.io/

3. Create fastTools environment.
     * Use command: conda env create -f environment_<your_os>.yml
     * Tentatively, this should work for OSX and Windows.
     
     * Alternatively, you could create your own environment, and download the following libraries manually:  
        matplotlib  
        biopython  
        pandas  
        seaborn
        

4. Include fastTools in your project directory alongside your own modules or scripts.

5. Import module.
     * Place "import fastTools" at the top of the script you want to use it in.

## Module Attributes
qScoreDict: Dictionary that maps Illumina QScore symbols to their integer values.  
* Usage
  * myQualityDict = fastTools.qScoreDict
  * fastTools.qScoreDict['?']
    * Returns 30

## FastqFile class
### Usage
#### Initialization
  * myfile = fastTools.FastqFile('Sample1_S1_L001_R1_001.fastq.gz')
    * Will create interleaved FastqFile object using R1 and R2 files.
    
  * myfile = fastTools.FastqFile('Sample1_S1_L001_R2_001.fastq.gz', False)
    * Will create a FastqFile from only the file name passed.
  
#### Class attributes
self.fastq: Name of FASTQ file passed during initialization.

self.sample: Truncated name of FASTQ file, convenient for labelling.

self.paired: True if R1 and R2 files were read and combined; False if only R1 or R2 file used.

self.fastqDataFrame: Pandas DataFrame object that holds all read/calculated data for the FastqFile object.

self.numReads: Number of reads in self.fastqDataFrame.

##### Example:
* myfile.fastq
  * Returns 'Sample1_S1_L001_R2_001.fastq.gz'
* myfile.sample
  * Returns 'Sample1_S1'
* myfile.numReads
  * Returns 432281
  * Note this is the same as len(myfile.fastqDataFrame)
  
#### Class methods
**These methods create a new column in self.fastqDataFrame that contains calculated data.**  

self.averageQuality()

self.reverseComplement()

self.aminoAcid()

self.calculateGC()

**These methods create plots that can either be displayed or saved.**  

self.plotAverageQuality(outfile=False)

self.plotGCcontent(outfile=False)

**This method saves your FastqFile object as a .fastq.gz file in the current directory.**  

self.writeFASTQ(outfile)

## FastaFile class
### Usage
#### Initialization  
  * myfile = fastTools.FastaFile('Sample1.fasta')
  
#### Class attributes
self.fasta: Name of FASTQ file passed during initialization.

self.fastaDataFrame: Pandas DataFrame object that holds all read/calculated data for the FastqFile object.

self.numReads: Number of reads in self.fastqDataFrame.

##### Example:
* myfile.fastq
  * Returns 'Sample1.fasta'
* myfile.numReads
  * Returns 216140
  * Note this is the same as len(myfile.fastaDataFrame)
  
#### Class methods
**These methods create a new column in self.fastqDataFrame that contains calculated data.**

self.reverseComplement()

self.aminoAcid()

self.calculateGC()

**This method creates a plot that can either be displayed or saved.**

self.plotGCcontent(outfile=False)

**This method saves your FastaFile object as a .fasta file in the current directory.**

self.writeFASTA(outfile)
