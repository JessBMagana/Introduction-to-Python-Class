from load import *
from orf import *
from dna import *
import random


# ----------------------  Functions Section  --------------------------------------------

def restOfORF(DNA):
    """Takes a DNA sequence written 5' to 3' and assumes that the sequence begins with 'ATG'.
    It then finds the next in frame stop codon, and returns the ORF from the start to that stop codon. """
    for i in range(0, len(DNA), 3):
        if DNA[i:(i + 3)] in stopList:
            seq = DNA[:i]
            return seq
    return DNA

#.........................................................................

def oneFrameV2(DNA):
    """ Ignores any ORF contained within other ORFs. """
    orfs = []
    i = 0
    for i in range(3):
        while i < len(DNA):
            if DNA[i:i + 3] == 'ATG':
                s = restOfORF(DNA[i:])
                orfs.append(s)
                i += len(s)
            i += 3
    return orfs

#.........................................................................

def longestORFV2(DNA):
    """Takes a DNA sequence, finds all ORFs within the sequence, counts the length of each ORF,
    then returns the longest ORF found."""
    orf = oneFrameV2(DNA)
    maxlen = 0
    maxseq = ""
    for i in range(len(orf)):
        if len(orf[i]) > maxlen:
            maxlen = len(orf[i])
    for j in range(len(orf)):
        if len(orf[j]) == maxlen:
            maxseq += orf[j]
    return maxseq

#......................................................................

def longestORFBothStrands(DNA):
    """Finds the longest ORF of a DNA sequence in either the forward or reverse direction."""
    orfs = []
    orfs.append(longestORFV2(DNA))
    DNA2 = reverseComplement(DNA)
    orfs.append(longestORFV2(DNA2))
    maxlen = 0
    maxseq = ''
    for i in range(len(orfs)):
        if len(orfs[i]) > maxlen:
            maxlen = len(orfs[i])
    for j in range(len(orfs)):
        if len(orfs[j]) == maxlen:
            maxseq += orfs[j]
    return maxseq

# ..................................................................

def collapse(List):  # 3rd, collapse(List)
    '''Puts Shuffled DNA and returns in a sequence form. '''
    output = ''  # initial output
    for s in List:  # for the string in the List
        output = output + s
    return output

# ..................................................................

def longestORFNoncoding(DNA, numReps):
    """
    Takes a DNA sequence in string format and first makes it into a list of nucleotides.
    It then randomly shuffles the list of nucleotides.
    The list is then turned back into a string.
    It then finds the longest ORF for both the forward and reverse sequences and appends them into a list.
    Finally, it returns the longest ORF found.
    ** The longest ORF found serves as a comparison for ORFs found in our real DNA.
    INPUT: DNA sequence, number of times we want to repeat this process.
    OUTPUT: Longest ORF found in our randomized DNA.
    """
    NonCodingORFList = []
    count = 0
    while count <= numReps:
        DNAlist = list(DNA)
        ranDNA = random.sample(DNAlist, len(DNAlist))
        ranDNAstr = collapse(ranDNA)
        orfs = longestORFBothStrands(ranDNAstr)
        NonCodingORFList.append(orfs)
        count = count + 1
    return len(max(NonCodingORFList))

#...................................................................

def findORFs(DNA):
    """
    Identifies all ORFs in a DNA sequence.
    INPUT: DNA sequence
    OUTPUT: list of ORFs
    """
    orfs = oneFrameV2(DNA)
    return orfs

#...................................................................

def findORFsBothStrands(DNA):
    """
    Searches both the forward and reverse complement strands of a DNA sequence for ORFs.
    INPUT: DNA sequence.
    OUTPUT: list of ORFs for both the forward and reverse complement strands.
    """
    DNA2 = reverseComplement(DNA)
    for i in range(3):
        forward = findORFs(DNA[i:])
        rev = findORFs(DNA2[i:])
        orfs = forward + rev
        return orfs

#....................................................................

def getCoordinates(orf,DNA):
    """
    Returns the beginning and end coordinates of an ORF in DNA.
    INPUT: an ORF and a DNA sequence.
    OUTPUT: the start and end coordinates of the ORF in the DNA in list form.
    """
    returnList = []
    startPos = 0
    endPos = 0
    revCompORF = ''
    startPos = DNA.find(orf)
    if startPos == -1:
        revCompORF = reverseComplement(orf)
        startPos = DNA.find(revCompORF)
        if startPos == -1:
            returnList = [-1, -1]
        else:
            endPos = startPos + len(revCompORF)
            returnList = [startPos, endPos]
    else:
        endPos = startPos + len(orf)
        returnList = [startPos, endPos]
    return returnList

#....................................................................

def geneFinder(DNASeq,minLen):
    '''Identifies ORFs longer that minLength and returns a list with the starting coordinate, end
    coordinate, and amino acid sequence of the longest ORFs. '''
    orfs = findORFsBothStrands(DNASeq)
    finalOutputList = []
    for i in orfs:
        if len(i) >= minLen:
            x, y = getCoordinates(i, DNASeq)
            aa = codingStrandToAA(i)
            seq = [x, y, aa]
            finalOutputList.append(seq)
    finalOutputList.sort()
    return finalOutputList

#........................................................................

def printGenes(gene):
    """
    geneL is a list of genes made by geneFinder with each sub-list showing the location of the start, end, and protein sequence.
    The output of printGenes is to go through the list and print each gene separately. In addition, print the start and end locations.
    """
    geneL = []
    i = 0
    for i, geneL in enumerate(gene):
        print(str(geneL[0]) + ',' + str(geneL[1]))
        print(geneL[2])

# .....................................................................

seq = loadSeq('X73525.fa')
minLen = longestORFNoncoding(seq,1500)
gene = geneFinder(seq,minLen)
print(printGenes(gene))

# seq2 = loadSeq('salDNA.fa')
# minLen2 = longestORFNoncoding(seq2,1500)
# gene2 = geneFinder(seq2,minLen2)
# print(printGenes(gene2))
