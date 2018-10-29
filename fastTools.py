#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fastTools was created and is maintained by Dane Gellerup (https://github.com/dgellerup).

This module is intended to make reading, analyzing, manipulating, and writing FASTQ
and FASTA files in Python easier, faster, and accessible to more users.

"""
import os
import gzip
import pandas as pd
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt


class FastqFile:
    """Class that creates a FASTQ file object when given a file name. FASTQ objects
    contain a truncated sample name (self.sample) and a Pandas data frame that
    holds sequence and quality data. Methods are available for calculating basic
    informatics, such as average read quality and GC content. If paired is passed
    as False an object will be built with only reads from the file name passed;
    if paired is redundantly passed as True, FastqFile will look
    for a matching forward/reverse (R1 or R2) file and combine them into an
    interleaved FastqFile object.
    
    Args:
        fastq1 (str): Name of a .fastq or .fastq.gz file.
        fastq2 (str) (default: None): Name of a .fastq or .fastq.gz file.
        paired (bool) (default: False): Boolean indicating whether FastqFile class should look for and interleave paired read file.
                                        Coerced to True if fastq2 argument is passed.
    
    """
        
    def __init__(self, fastq1, fastq2=None, paired=False):
        
        if fastq1 in os.listdir():
            self.fastq1 = fastq1
        else:
            print(f'{fastq1} not found in current directory')
            
            return
        
        # Takes the full file name and shortens it to a readable name.
        try:
            self.sample = '_'.join(fastq1.split("_")[:2])
        except:
            self.sample = fastq1
        
        self.paired = paired
        
        if fastq2: 
            if fastq2 in os.listdir():
                self.paired = True
                self.fastq2 = fastq2
            else:
                self.paired = False
                self.fastq2 = "None"
                print((f'{fastq2} not found in current directory.\n'
                       f'Continuing with single file {self.fastq1}'))
        else:
            if self.paired:
                if "_R1_" in fastq1:
                    # Replace '_R1_' with '_R2_' in fastq1 file name.
                    if fastq1.replace("_R1_", "_R2_") in os.listdir():
                        self.fastq2 = fastq1.replace("_R1_", "_R2_")
                        fastq2 = True
                    else:
                        self.paired = False
                        self.fastq2 = "None"
                        fastq2 = False
                        print((f"Attempted to find a mate for {self.fastq1}, but none was found.\n"
                               f"Continuing with single file {self.fastq1}"))
                elif "_R2_" in fastq1:
                    # Replace '_R2_' with '_R1_' in fastq1 file name.
                    holder = fastq1.replace("_R2_", "_R1_")
                    if holder in os.listdir():
                        self.fastq2 = self.fastq1
                        self.fastq1 = holder
                        fastq2 = True
                    else:
                        self.paired = False
                        self.fastq2 = "None"
                        fastq2 = False
                        print((f"Attempted to find a mate for {self.fastq1}, but none was found.\n"
                               f"Continuing with single file {self.fastq1}"))
                else:
                    fastq2 = False
                    self.paired = False
                    self.fastq2 = "None"
                    print((f"Attempted to find a mate for {self.fastq1}, but none was found.\n"
                           f"Continuing with single file {self.fastq1}"))
            else:
                fastq2 = False
                self.fastq2 = "None"
                self.paired = False
            
        # Read fastq1 into fqlines list.
        if self.fastq1.endswith('.gz'):
            with gzip.open(self.fastq1, 'rt') as fastqFile:
                fqlines = fastqFile.readlines()
        else:
            with open(self.fastq1, 'r') as fastqFile:
                fqlines = fastqFile.readlines()
           
        if self.paired:
            # If there is a reverse file, read it and add it to fqlines.
            if self.fastq2.endswith('.gz'):
                with gzip.open(self.fastq2, 'rt') as fastq2File:
                    fq2lines = fastq2File.readlines()
            else:
                with open(self.fastq2, 'r') as fastq2File:
                    fq2lines = fastq2File.readlines()
    
            # Create a new file lines list, and add FASTQ reads in an interleaved orientation.
            pairedLinesList = []
        
            for i in range(0, len(fqlines), 4):
                pairedLinesList.append(fqlines[i])
                pairedLinesList.append(fqlines[i+1])
                pairedLinesList.append(fqlines[i+2])
                pairedLinesList.append(fqlines[i+3])
        
                pairedLinesList.append(fq2lines[i])
                pairedLinesList.append(fq2lines[i+1])
                pairedLinesList.append(fq2lines[i+2])
                pairedLinesList.append(fq2lines[i+3])
                
            # Create lists to hold FASTQ sequences and quality strings from pairedLinesList.
            fqnameList = [pairedLinesList[i] for i in range(0, len(pairedLinesList), 4)]
            fqseqList = [pairedLinesList[i] for i in range(1, len(pairedLinesList), 4)]
            fqdirectionList = [pairedLinesList[i] for i in range(2, len(pairedLinesList), 4)]
            fqqualList = [pairedLinesList[i] for i in range(3, len(pairedLinesList), 4)]
                
                       
        # If there is no reverse file, or an unpaired FastqFile object is desired.
        else:
            
            # Create lists to hold FASTQ sequences and quality strings from fqlines.
            fqnameList = [fqlines[i] for i in range(0, len(fqlines), 4)]
            fqseqList = [fqlines[i] for i in range(1, len(fqlines), 4)]
            fqdirectionList = [fqlines[i] for i in range(2, len(fqlines), 4)]
            fqqualList = [fqlines[i] for i in range(3, len(fqlines), 4)]
            
        
        # Create a data frame "fastqDataFrame" from fqseqList and fqqualList
        self.fastqDataFrame = pd.DataFrame({'Name': fqnameList, 'Seq': fqseqList, 
                                      'Direction': fqdirectionList, 'Qual': fqqualList})
        
    
    def __len__(self):
        return len(self.fastqDataFrame)
    
    def __str__(self):
        return (f'{self.sample}')
    
    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'in1={self.fastq1}, in2={self.fastq2})'
                f' Columns: {list(self.fastqDataFrame.columns)}')
                
    def numReads(self):
        return len(self.fastqDataFrame)
    
    
    def averageQuality(self):
        """Appends a column to self.fastqDataFrame that contains average quality scores for each read.
        
        Args:
            self
            
        Returns:
            None
        
        """
        
        quallist = self.fastqDataFrame['Qual']
        qualscores = []
        
        # Create a list for the quality strings and append the corresponding quality scores
        for item in quallist:
            qualstringseries = pd.Series([qScoreDict[symbol] for symbol in item.strip()])
            
            # Create a variable for the overall average quality of each line
            avequal = qualstringseries.mean()    
            qualscores.append(avequal)
        
        # Adds the quality scores in a coloumn to the dataframe
        self.fastqDataFrame['Avg Qual'] = qualscores
        
        
    def reverseComplement(self):
        """Creates a new column in self.fastqDataFrame to hold reverse complement
        DNA sequences.
        
        Args:
            self
            
        Returns:
            none
        """
        
        self.fastqDataFrame['Reverse Complement'] = self.fastqDataFrame['Seq'].apply(revComp)
        
        
    def aminoAcid(self):
        """Creates a new column in self.fastqDataFrame to hold Amino Acid sequences.
        
        Args:
            self
            
        Returns:
            none
        """

        self.fastqDataFrame['AA Sequence'] = self.fastqDataFrame['Seq'].apply(aa)
        
        
    def calculateGC(self):
        """Creates a new column in self.fastqDataFrame to hold GC content (float) for each
        sequence.
        
        Args:
            self
            
        Returns:
            none
        """
        
        self.fastqDataFrame['GC Content'] = self.fastqDataFrame['Seq'].apply(gc_content)
        
        
    def plotAverageQuality(self, outfile=False):
        """FastqFile class method that allows a user to easily plot a histogram
        of per-read average Q Score for a FastqFile object. If average quality
        has not yet been calculated, calls self.averageQuality() first.
        Requires no plotting experience.
        
        Args:
            self
            outfile (str) (optional): Name of plot output file.
            
        Returns:
            none
        """
                
        if 'Avg Qual' in self.fastqDataFrame.columns:
            
            try:
                sns
            except NameError:
                try:
                    import seaborn as sns
            
                except ModuleNotFoundError:
                        print("It appears that the Seaborn library is not installed.")
                        print("You can install Seaborn on the command line with:")
                        print("conda install seaborn \n -or- \n pip install seaborn")
                        
                        return
                
            fig = plt.figure(figsize=(9, 6))
            
            sns.distplot(self.fastqDataFrame['Avg Qual'], kde=False)
            
            plt.title("Per Sequence Average Quality")
            plt.ylabel("Count")
            plt.xlabel("Quality Score")
            
            plt.tight_layout()
            
            if outfile:
                plt.savefig(outfile)
            else:
                plt.show()
            
        else:
            print("Calculating per-read average quality data.")
            
            self.averageQuality()
            
            print("'Avg Qual' column has been added to self.fastqDataFrame")
            
            self.plotAverageQuality()
            

    def plotGCcontent(self, outfile=False):
        """FastqFile class method that allows a user to easily plot a histogram
        of per-read GC content for a FastqFile object. If GC content has not yet
        been calculated, calls self.calculateGC() first. Requires no plotting experience.
        
        Args:
            self
            outfile (str) (optional): Name of plot output file.
            
        Returns:
            none
        """        
                
        if 'GC Content' in self.fastqDataFrame.columns:
            
            try:
                sns
            except NameError:
                try:
                    import seaborn as sns
                    
                except ModuleNotFoundError:
                        print("It appears that the Seaborn library is not installed.")
                        print("You can install Seaborn on the command line with:")
                        print("conda install seaborn \n -or- \n pip install seaborn")
                        
                        return
            
            fig = plt.figure(figsize=(9, 6))
            
            sns.distplot(self.fastqDataFrame['GC Content'], kde=False)
            
            plt.title("Per Sequence GC Content")
            plt.ylabel("Count")
            plt.xlabel("GC Content (%)")
            
            plt.tight_layout()
            
            if outfile:
                plt.savefig(outfile)
            else:
                plt.show()
            
        else:
            print("Calculating per-read GC content data.")
            
            self.calculateGC()
            
            print("'GC Content' column has been added to self.fastqDataFrame")
            
            self.plotGCcontent()
            

    def writeFASTQ(self, outfile):
        """Takes an outfile string and a dataframe. Converts 2D DataFrame to a single
        column resembling a FASTQ file with each line as a row. This makes writing the
        output file much faster than using Python's built-in write function in a loop.
        
        Args:
            outfile (str): Name of the output file, cleaned and formatted by this program.
            dataframe (:obj: pd.DataFrame): Pandas dataframe containing FASTQ info, each read on one line.
        
        """
        
        try:
            # Create one-column DataFrames from the columns of dataframe.
            namedf = pd.DataFrame(self.fastqDataFrame['Name'])
            seqdf = pd.DataFrame(self.fastqDataFrame['Seq'])
            directiondf = pd.DataFrame(self.fastqDataFrame['Direction'])
            qualdf = pd.DataFrame(self.fastqDataFrame['Qual'])
        except:
            return "This doesn't appear to be a FastqFile object."
        
        # Give each row a pseudo-index, offset from 0-4, stepping 4 for each index.
        nameid = list(range(0, len(namedf)*4, 4))
        seqid = list(range(1, len(namedf)*4, 4))
        dirid = list(range(2, len(namedf)*4, 4))
        qualid = list(range(3, len(namedf)*4, 4))
        
        # Add an 'id' column to each individual DataFrame.
        namedf['id'] = nameid
        seqdf['id'] = seqid
        directiondf['id'] = dirid
        qualdf['id'] = qualid
        
        # Create new column names.
        namedf.columns = ["line", "id"]
        seqdf.columns = ["line", "id"]
        directiondf.columns = ["line", "id"]
        qualdf.columns = ["line", "id"]
        
        # Concatenate individual DataFrames into one longdf (pseudo-indeces out of order).
        longdf = pd.concat([namedf, seqdf, directiondf, qualdf])
        
        # Sort values in longdf by thier pseudo-indeces, restoring original FASTQ structure.
        longdf.sort_values('id', inplace=True)
        
        # Reset index so there are no multiple indeces.
        longdf.reset_index(drop=True, inplace=True)
        
        # Drop pseudo-indeces.
        longdf.drop('id', axis=1, inplace=True)
        
        # Remove newline characters.
        longdf['line'] = longdf['line'].apply(removeNewline)
        
        # Use pandas to write DataFrame as FASTQ file (faster than built-in Python).
        longdf.to_csv(outfile, index=False, header=False, quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\", compression='gzip')
       
        
class FastaFile:
    """Class that creates a FASTA file object when given a file name. FASTA objects
    contain a truncated sample name (self.sample) and a Pandas data frame that
    holds sequence and sequence name data. Methods are available for calculating basic
    informatics, such as GC content and reverse complement sequences.
    
    Args:
        fasta (str): Name of a .fasta file.
    
    """
    
    def __init__(self, fasta):
        self.fasta = fasta
        
        # Read fasta into falines list.
        if fasta.endswith('.gz'):
            with gzip.open(fasta, 'rt') as fastaFile:
                falines = fastaFile.readlines()
        else:
            with open(fasta, 'r') as fastaFile:
                falines = fastaFile.readlines()

        # Create lists to hold FASTQ sequences and quality strings from fqlines
        fanameList = []
        faseqList = []
        
        # Starting at index 1, add every other line to faseqList
        for i in range(0, len(falines), 2):
            fanameList.append(falines[i].strip())

        for j in range(1, len(falines), 2):
            faseqList.append(falines[j].strip())

        self.fastaDataFrame = pd.DataFrame({'Name': fanameList, 'Seq': faseqList})
        
        
    def __len__(self):
        return len(self.fastaDataFrame)
    
    def __str__(self):
        return (f'{self.fasta}')
    
    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'in1={self.fasta}'
                f' Columns: {list(self.fastaDataFrame.columns)}')
    
    
    def numReads(self):
        return len(self.fastaDataFrame)
    
    
    def reverseComplement(self):
        """Creates a new column in self.fastaDataFrame to hold reverse complement
        DNA sequences.
        
        Args:
            self
            
        Returns:
            none
        """
        
        self.fastaDataFrame['Reverse Complement'] = self.fastaDataFrame['Seq'].apply(revComp)
        
        
    def aminoAcid(self):
        """Creates a new column in self.fastaDataFrame to hold Amino Acid sequences.
        
        Args:
            self
            
        Returns:
            none
        """
        
        self.fastaDataFrame['AA Sequence'] = self.fastaDataFrame['Seq'].apply(aa)
        
        
    def calculateGC(self):
        """Creates a new column in self.fastaDataFrame to hold GC content (float) for each
        sequence.
        
        Args:
            self
            
        Returns:
            none
        """
        
        self.fastaDataFrame['GC Content'] = self.fastaDataFrame['Seq'].apply(gc_content)
        
        
    def plotGCcontent(self, outfile=False):
        """FastqFile class method that allows a user to easily plot a histogram
        of per-read GC content for a FastqFile object. If GC content has not yet
        been calculated, calls self.calculateGC() first. Requires no plotting experience.
        
        Args:
            self
            outfile (str) (optional): Name of plot output file.
            
        Returns:
            none
        """
        
        if 'GC Content' in self.fastaDataFrame.columns:
            
            try:
                sns
            except NameError:
                try:
                    import seaborn as sns
                    
                except ModuleNotFoundError:
                        print("It appears that the Seaborn library is not installed.")
                        print("You can install Seaborn on the command line with:")
                        print("conda install seaborn \n -or- \n pip install seaborn")
                        
                        return
                    
            fig = plt.figure(figsize=(9, 6))
            
            sns.distplot(self.fastaDataFrame['GC Content'], kde=False)
            
            plt.title("Per Sequence GC Content")
            plt.ylabel("Count")
            plt.xlabel("GC Content (%)")
            
            plt.tight_layout()
            
            if outfile:
                plt.savefig(outfile)
            else:
                plt.show()
            
        else:
            print("Calculating per-read GC content data.")
            
            self.calculateGC()
            
            print("'GC Content' column has been added to self.fastaDataFrame")
            
            self.plotGCcontent()
        
        
    def writeFASTA(self, outfile):
        """Takes an outfile string and a dataframe. Converts 2D DataFrame to a single
        column resembling a FASTA file with each line as a row. This makes writing the
        output file much faster than using Python's built-in write function in a loop.
        
        Args:
            self
            outfile (str): Name of the output file, cleaned and formatted by this program.
            
        Returns:
            None
        
        """
        try:
            # Create one-column DataFrames from the columns of dataframe.
            namedf = pd.DataFrame(self.fastaDataFrame['Name'])
            seqdf = pd.DataFrame(self.fastaDataFrame['Seq'])
        except:
            return "This doesn't appear to be a FastaFile object."
        
        # Give each row a pseudo-index, offset from 0-4, stepping 4 for each index.
        nameid = list(range(0, len(namedf)*4, 4))
        seqid = list(range(1, len(namedf)*4, 4))
        
        # Add an 'id' column to each individual DataFrame.
        namedf['id'] = nameid
        seqdf['id'] = seqid
        
        # Create new column names.
        namedf.columns = ["line", "id"]
        seqdf.columns = ["line", "id"]
        
        # Concatenate individual DataFrames into one longdf (pseudo-indeces out of order).
        longdf = pd.concat([namedf, seqdf])
        
        # Sort values in longdf by thier pseudo-indeces, restoring original FASTQ structure.
        longdf.sort_values('id', inplace=True)
        
        # Reset index so there are no multiple indeces.
        longdf.reset_index(drop=True, inplace=True)
        
        # Drop pseudo-indeces.
        longdf.drop('id', axis=1, inplace=True)
        
        # Remove newline characters.
        longdf['line'] = longdf['line'].apply(removeNewline)
        
        # Use pandas to write DataFrame as FASTQ file (faster than built-in Python).
        longdf.to_csv(outfile, index=False, header=False, quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\", compression='gzip')
       
        
def revComp(seqString):
    """Custom function used by FastqFile and FastaFile objects in a pandas apply() call to
    generate an reverse complement DNA sequence string from a DNA sequence string input.
    This function can also be called with fastTools.revComp(yourSequence), where it will
    return a reverse complement DNA sequence string. If you require a Bio.Seq object, 
    use biopython directly instead of this function.
    
    Args:
        seqString (str): DNA sequence string.
        
    Returns:
        revComp (str): DANE sequence string.
    """
    
    revComp = str(Seq(seqString).reverse_complement())
    
    return revComp


def gc_content(seqString):
    """Custom function used by FastqFile and FastaFile objects in a pandas apply() call to
    calculate the GC content of a DNA sequence. This function can also be called on
    a sequence manually with fastTools.gc_content(yourSequence). Same exact usage
    as biopython's Bio.SeqUtils GC() function.
    
    Args:
        seqString (str): DNA sequence string.
        
    Returns: (float): GC content of given DNA sequence in percent.
    """
    
    GCpercent = GC(seqString)
    
    return GCpercent        


def aa(seqString):
    """Custom function used by FastqFile and FastaFile objects in a pandas apply() call to
    generate an Amino Acid string translated from a DNA sequence string input.
    This function can also be called with fastTools.aa(yourSequence), where it will
    return an Amino Acid translation string. If you require a Bio.Seq object, 
    use biopython directly instead of this function.
    
    Args:
        seqString (str): DNA sequence string.
        
    Returns:
        amino (str): Amino Acid string.
    """
    
    amino = str(Seq(seqString).translate())
    
    return amino


def removeNewline(x):
    """Custom function to use with apply() in a pandas DataFrame. Simply removes
    new line characters from strings in a column.
    """
    
    x = x.strip('\n')
    
    return x
     
        

      
qScoreDict = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, '+': 10,
                ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                'J': 41, 'K': 42}