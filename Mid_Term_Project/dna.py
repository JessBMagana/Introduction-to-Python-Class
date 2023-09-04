from aminoAcids import *
from load import *

stopList = ['TAG', 'TAA', 'TGA']

def compBase(N):
    if N=="A":
        return("T")
    elif N=="T":
        return("A")
    elif N=="C":
        return("G")
    elif N=="G":
        return("C")

def reverse(S):
    outStr=""
    for s in S:
        outStr=s+outStr
    return outStr

def reverseComplement( DNA ):
    rDNA=reverse(DNA)
    outStr=""
    for b in rDNA:
        outStr+=compBase(b)
    return(outStr)

def amino( codon ):
    for j in range(len(codons)):
        codonL=codons[j]
        if codon in codonL:
            return aa[j]   

def codingStrandToAA( DNA ):
    aaStr=""
    for i in range(0,len(DNA),3):
        thisCodon=DNA[i:(i+3)]
        aaStr+=amino(thisCodon)
    return aaStr

def restOfORF(DNA):
    """Takes a sequence starting with an ATG and finds first stop
    codon. Returns ORF. If no in-frame stop codon, return whole
    sequence."""
    for i in range(0,len(DNA),3):
        if DNA[i:(i+3)] in stopList:
            seq=DNA[:i]
            return seq
    return DNA
