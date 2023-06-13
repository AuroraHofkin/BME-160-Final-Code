#!/usr/bin/env python3
# Name: Aurora Hofkin (ahofkin)
# Group Members: Ash Ediga, Erika Stewart, Julia Maynard
from sequenceAnalyzer import FastAreader,ProtienFinder,Analysis

def main():
    
    import os

    folder = "/mnt/c/Users/kaito/Desktop/Bme 160/Final/Genomes"

    finder = ProtienFinder()
    analysis = Analysis()
    fastaFiles = os.listdir(folder)

    for GenomeFasta in fastaFiles:
        if GenomeFasta.endswith(".fa"):
            finder.runprodigal(inputfile = GenomeFasta)
            #analysis.BlastProtein()
            analysis.PI()
            analysis.calculateAvgpI()
            analysis.MolecularWeight()
            analysis.calculateAvgMW()

            # Extract the genome name without the file extension
            genomename = os.path.splitext(GenomeFasta)[0]
            analysis.ProteinGraph(genomename)

            print(genomename)
            #average PI
            #average MW
            #print(proteins + Pi + MW)
main()
            
"""Reader = FastAreader()
            for head, seq in Reader.readFasta():
                finder = ProtienFinder(seq)
                finder.runprodigal()
                finder.BlastProtein()
                finder.PI()
                finder.MolecularWeight()
                finder.ProteinGraph()"""
                
  
                #gives protein and Pi
                #sort for highest value of PI
                #print figure of PI VS molecular weight KDa from proteins
                #figure of PI vs Molecule weight comapring different viruses



#output to fasta file ; May want to remove extra genes that were found but not classified. (may add for extra option)
#then after output take the length the fasta file with CDS header and save them to a list
#then iterate over that list to find the sequence and transcribe them.
#after transcription calculates PI
#graphs PI relative to the Genomes got""")