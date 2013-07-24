from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord

def count_fasta( fasta_file ):
    """
    Count the number of records in a Fasta
    """
    return len(list(FastaReader(fasta_file)))

def read_fasta_names( fasta_file ):
    ids = set()
    for record in FastaReader( fasta_file ):
        ids.add( record.name )
    return ids

def write_fasta( filename, seq_name, sequence ):
    """
    Write a DNA sequence out to file as a Fasta
    """
    with FastaWriter( filename ) as writer:
        record = FastaRecord( seq_name, sequence )
        writer.writeRecord( record )

def is_exe( file_path ):
    """
    Check if a file exists and is executable
    """
    if file_path is None:
        return False
    return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

def which(program):
    """
    Find and return path to local executables  
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def get_file_type( filename ):
    """
    Get the filetype of a PacBio-compatible sequence data file
    """
    if filename.lower().endswith('.fofn'):
        return 'fofn'
    elif filename.lower().endswith('.bas.h5'):
        return 'bash5'
    elif filename.lower().endswith('.bax.h5'):
        return 'bash5'
    else:
        msg = 'Input file must be Bas.H5, Bax.H5, or FOFN!'
        log.error( msg )
        raise TypeError( msg )

def create_directory( directory ):
    """
    Create a directory if it doesn't already exist
    """
    if os.path.isdir( directory ):
        return 
    try:
        os.mkdir( directory )
    except:
        msg = 'Unable to create directory "%s"' % directory
        log.error( msg )
        raise IOError( msg )
