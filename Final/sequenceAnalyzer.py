#!/usr/bin/env python3
# Name: Aurora Hofkin (ahofkin)
# Group Members: Ash Ediga, Erika Stewart, Julia Maynard
import sys

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

import subprocess

class ProtienFinder:
    """Prodigal gene finder import:
    Finds all hypothetical genes that are protein coding and gives the sequences of each.
    Credit: Prodigal Gene Finder by Doug Hyatt et al.
"""

    def runprodigal(self,inputfile = " ", outputfile = "proteincodinggenes.fa"):

        """Finds all the hypothetical genes."""

        command = f"prodigal -p meta -i {inputfile} -a {outputfile}"

        try:
            subprocess.run(command, shell=True, check=True)

            print("Prodigal execution successful.")

        except subprocess.CalledProcessError as e:

            print(f"Prodigal execution failed with error code {e.returncode}.")
            
            #print(e.output) Will Remove for Testing Purpose
            
class Analysis:

    """Analyisis each found Protein sequence to get the Name, pI, MW, then graphs the information.
"""

    def __init__(self):
        """Prepares the Protein sequence for analysis and defines dictionaries.
    """

        self.ProteinName = set()
        self.ProteinSeq = []
        self.ProteinID = set()
        self.pIs = set()
        self.molecularWeights = set()

        Reader = FastAreader("proteincodinggenes.fa")

        for head, seq in Reader.readFasta():

            sequence = seq.replace("*","")
            self.ProteinSeq.append(sequence)
            self.ProteinID.add(head)

    def BlastProtein(self):
        """Preforms a blast on each protein sequence of get the name of each
        Note: this could take some time if the virus has many protein coding genes
    """
        from Bio.Blast import NCBIWWW
        from Bio.Blast import NCBIXML

        #Iterates through protein sequence
        for protein in self.ProteinSeq:

            # Perform BLAST search
            resultHandle = NCBIWWW.qblast("blastp", "nr", protein)

            # Parse the BLAST results
            blastRecords = NCBIXML.parse(resultHandle)

            # Retrieve the scientific name from the top hit
            for blastRecord in blastRecords:
                alignment = blastRecord.alignments[0]
                scientificName = alignment.title.split("[")[-1].split("]")[0]
                
                # Saves scientific name
                self.ProteinName.add(scientificName)

    def PI (self):
        """Calculates the pI from each protein sequence.
        """

        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        #Iterates through protein sequence
        for protein in self.ProteinSeq:

            # Create a ProteinAnalysis object from the protein sequence
            proteinAnalysis = ProteinAnalysis(protein)

            # Calculate the isoelectric point (pI)
            pI = proteinAnalysis.isoelectric_point()

            # Saves the pI value
            self.pIs.add(pI)
        
    def calculateAvgpI(self):

        # Calculates average pI
        avgpI = sum(self.pIs) / len(self.pIs)

        return avgpI
    
    def MolecularWeight(self):
        """Calculates the molecular weight of each protein sequence.
        """

        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        #Iterates through protein sequence
        for protein in self.ProteinSeq:

            # Create a ProteinAnalysis object
            protein_analysis = ProteinAnalysis(protein)

            # Calculate molecular weight
            molecularweight = protein_analysis.molecular_weight() /1000
            self.molecularWeights.add(molecularweight)

    def calculateAvgMW(self):

        # Calculate average molecular weight
        avgMW = sum(self.molecularWeights) / len(self.molecularWeights)

        return avgMW


    def ProteinGraph(self,genomename):

        import matplotlib.pyplot as plt

        # Define Molecular weights and PI values
        molecularweights = list(self.molecularWeights)
        pIValues = list(self.pIs)

        # Create a scatter plot
        plt.scatter(pIValues, molecularweights)

        # Set labels and title
        plt.xlabel('Isoelectric Point (pI)')
        plt.ylabel('Molecular Weight (kDa)')
        plt.title(f'{genomename}')

        # Display the plot
        plt.show()