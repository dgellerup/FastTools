#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 10:43:54 2018

@author: dgelleru
"""
import os
import gzip
import pandas as pd
import csv


class FastqFile:
    """Class that creates a FASTQ file object when given a file name. FASTQ objects
    contain a truncated sample name (self.sample) and a Pandas data frame that
    holds sequence and quality data. Methods are available for calculating basic
    informatics, such as average read quality.
    
    Args:
        fastq (str): Name of a .fastq.gz file.
    
    """
        
    def __init__(self, fastq):
        
        self.fastq = fastq
        
        # Takes the full file name and  shortens it to a readable name.
        self.sample = '_'.join(fastq.split("_")[:2])
        
        self.twoFilesFlag = False
        
        # Splice R2 in place of R1 in fastq file name.
        fastq2 = fastq.split("_R1_")[0] + "_R2_" + fastq.split("_R1_")[1]
        
        # Read fastq into fqlines list.
        if fastq.endswith('.gz'):
            with gzip.open(fastq, 'rt') as fastqFile:
                fqlines = fastqFile.readlines()
        else:
            with open(fastq, 'r') as fastqFile:
                fqlines = fastqFile.readlines()
           
        # If there is a reverse file, read it and add it to fqlines.
        if fastq2 in os.listdir():
            if fastq2.endswith('.gz'):
                with gzip.open(fastq2, 'rt') as fastq2File:
                    fq2lines = fastq2File.readlines()
            else:
                with open(fastq2, 'r') as fastq2File:
                    fq2lines = fastq2File.readlines()

            # Create a new file lines list, and add fastq reads in an interleaved orientation.
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
            
            self.twoFilesFlag = True
                
            # Create lists to hold FASTQ sequences and quality strings from pairedLinesList.
            fqnameList = [pairedLinesList[i] for i in range(0, len(pairedLinesList), 4)]
            fqseqList = [pairedLinesList[i] for i in range(1, len(pairedLinesList), 4)]
            fqdirectionList = [pairedLinesList[i] for i in range(2, len(pairedLinesList), 4)]
            fqqualList = [pairedLinesList[i] for i in range(3, len(pairedLinesList), 4)]
                       
        # If there is no reverse file.
        else:
            
            # Create lists to hold FASTQ sequences and quality strings from fqlines.
            fqnameList = [fqlines[i] for i in range(0, len(fqlines), 4)]
            fqseqList = [fqlines[i] for i in range(1, len(fqlines), 4)]
            fqdirectionList = [fqlines[i] for i in range(2, len(fqlines), 4)]
            fqqualList = [fqlines[i] for i in range(3, len(fqlines), 4)]
            
        
        # Create a data frame "fqfiledf" from fqseqList and fqqualList
        self.fqfiledf = pd.DataFrame({'Name': fqnameList, 'Seq': fqseqList, 
                                      'Direction': fqdirectionList, 'Qual': fqqualList})
    
        self.numReads = len(self.fqfiledf)
            
    
    def qualAverage(self):
        """Appends a column to self.fqfiledf that contains average quality scores for each read.
        
        Args:
            self
            
        Returns:
            None
        
        """
        
        quallist = self.fqfiledf['Qual']
        qualscores = []
        
        # Create a list for the quality strings and append the corresponding quality scores
        for item in quallist:
            qualstringlist = []
            
            for symbol in item.strip():
                qualstringlist.append(qualdict[symbol])
            
            # Create a variable for the overall average quality of each line
            avequal = sum(qualstringlist)/float(len(qualstringlist))    
            qualscores.append(avequal)
        
        # Adds the quality scores in a coloumn to the dataframe
        self.fqfiledf['Avg Qual'] = qualscores
        
        
def writeFASTQ(outfile, dataframe):
    """Takes an outfile string and a dataframe. Converts 2D DataFrame to a single
    column resembling a FASTQ file with each line as a row. This makes writing the
    output file much faster than using Python's built-in write function in a loop.
    
    Args:
        outfile (str): Name of the output file, cleaned and formatted by this program.
        dataframe (:obj: pd.DataFrame): Pandas dataframe containing FASTQ info, each read on one line.
    
    """
    try:
        # Create one-column DataFrames from the columns of dataframe.
        namedf = pd.DataFrame(dataframe['Name'])
        seqdf = pd.DataFrame(dataframe['Seq'])
        directiondf = pd.DataFrame(dataframe['Direction'])
        qualdf = pd.DataFrame(dataframe['Qual'])
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
   

def removeNewline(x):
    """Custom function to use with apply() in a pandas DataFrame. Simply removes
    new line characters from strings in a column.
    """
    
    x = x.strip('\n')
    return x
     
        
class FastaFile:
    
    def __init__(self, fasta):
        self.fasta = fasta
        
        # Read fastq with gzip library.
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

        self.fafiledf = pd.DataFrame({'Name': fanameList, 'Seq': faseqList})
        
        
def writeFASTA(outfile, dataframe):
    """Takes an outfile string and a dataframe. Converts 2D DataFrame to a single
    column resembling a FASTQ file with each line as a row. This makes writing the
    output file much faster than using Python's built-in write function in a loop.
    
    Args:
        outfile (str): Name of the output file, cleaned and formatted by this program.
        dataframe (:obj: pd.DataFrame): Pandas dataframe containing FASTQ info, each read on one line.
    
    """
    try:
        # Create one-column DataFrames from the columns of dataframe.
        namedf = pd.DataFrame(dataframe['Name'])
        seqdf = pd.DataFrame(dataframe['Seq'])
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
   
      
qualdict = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, '+': 10,
                ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                'J': 41, 'K': 42}