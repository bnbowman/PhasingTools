from pbcore.io.FastaIO import FastaReader

def count_fasta( fasta_file ):
    return len(list(FastaReader(fasta_file)))
