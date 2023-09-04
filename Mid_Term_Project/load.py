
def loadSeq(fileName):
    """Load sequence from a fasta file with a single entry."""
    f=open(fileName,"r")
    linesL=f.readlines()
    # combine lines, skipping first
    seq=""
    for line in linesL[1:]:
        seq+=line.rstrip()
    return(seq)
