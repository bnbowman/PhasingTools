from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord

def count_fasta( fasta_file ):
    return len(list(FastaReader(fasta_file)))

def read_fasta_names( fasta_file ):
    ids = set()
    for record in FastaReader( fasta_file ):
        ids.add( record.name )

def write_fasta( filename, seq_name, sequence ):
    with FastaWriter( filename ) as writer:
        record = FastaRecord( seq_name, sequence )
        writer.writeRecord( record )
