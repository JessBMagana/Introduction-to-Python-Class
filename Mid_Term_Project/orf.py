import random

stopList = ['TAG', 'TAA', 'TGA']


def loadSeq(filename):
    '''Takes in a FASTA file and returns the DNA sequence in the file.'''
    f = open(filename)  # open the file
    linesList = f.readlines()  # read in the file as a list of its lines
    f.close()  # close the file
    newSeq = ''
    for l in linesList:
        if l[0] != '>':
            newSeq = newSeq + l[:-1]
    return (newSeq)


def compBase(N):
    '''Returns the complimentary base'''
    if N == 'A':
        return 'T'
    elif N == 'T':
        return 'A'
    elif N == 'C':
        return 'G'
    elif N == 'G':
        return 'C'
    else:
        return "This is not a valid nucleotide"


def reverse(string):
    '''Takes a string as input and returns the reverse of the string'''
    revstring = ''
    for i in string:
        revstring = i + revstring
    return revstring


def reverseComplement(DNA):
    '''Take a DNA string as input. Creates a new string with the reverse of the input. Returns the reverse
    complementary strand.'''
    rev = reverse(DNA)
    comp = ''
    for i in rev:
        comp = comp + compBase(i)
    return comp


def amino(codon):
    '''Returns the amino acid symbol given the codon'''
    for ind in range(len(aa)):
        if codon in codons[ind]:
            return aa[ind]


aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
      '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
      'G']

codons = [['TTT', 'TTC'],
          ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
          ['ATT', 'ATC', 'ATA'],
          ['ATG'],
          ['GTT', 'GTC', 'GTA', 'GTG'],
          ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
          ['CCT', 'CCC', 'CCA', 'CCG'],
          ['ACT', 'ACC', 'ACA', 'ACG'],
          ['GCT', 'GCC', 'GCA', 'GCG'],
          ['TAT', 'TAC'],
          ['TAA', 'TAG', 'TGA'],
          ['CAT', 'CAC'],
          ['CAA', 'CAG'],
          ['AAT', 'AAC'],
          ['AAA', 'AAG'],
          ['GAT', 'GAC'],
          ['GAA', 'GAG'],
          ['TGT', 'TGC'],
          ['TGG'],
          ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
          ['GGT', 'GGC', 'GGA', 'GGG']]


def codingStrandToAA(DNA):
    '''Takes a DNA sequence and returns its corresponding sequence
    of amino acids.'''
    codonList = [DNA[i:i + 3] for i in range(0, len(DNA), 3)]
    stringAmino = ''
    for l in codonList:
        stringAmino = stringAmino + amino(l)
    return stringAmino


def restOfORF(DNA):
    """Takes a DNA sequence written 5' to 3' and assumes that the sequence begins with 'ATG'.
    It then finds the next in frame stop codon, and returns the ORF from the start to that stop codon. """
    for i in range(0, len(DNA), 3):
        if DNA[i:(i + 3)] in stopList:
            seq = DNA[:i]
            return seq
    return DNA


def oneFrame(DNA):
    '''Takes a DNA sequence and finds the ORF.'''
    orfs = []
    for k in range(3):
        for i in range(k, len(DNA), 3):
            if DNA[i:i + 3] == 'ATG':
                orfs.append(restOfORF(DNA[i:]))
    return orfs


def longestORF(DNA):
    '''Takes a DNA sequence, finds all ORFs within the sequence, counts the length of each ORF,
    then returns the longest ORF found.'''
    orfs = oneFrame(DNA)
    maxlen = 0
    maxseq = ''
    for i in range(len(orfs)):
        if len(orfs[i]) > maxlen:
            maxlen = len(orfs[i])
    for j in range(len(orfs)):
        if len(orfs[j]) == maxlen:
            maxseq += orfs[j]
    return (maxseq)


