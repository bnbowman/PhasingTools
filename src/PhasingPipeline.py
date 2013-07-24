#! /usr/bin/env python
import os, logging

log = logging.getLogger()

class PhasingPipeline( object ):

    def __init__( self ):
        # Parse the options 
        parse_args()
        # Initialize the rest of the project
        self._initialize_project()
        self._initialize_logging()

    def _initialize_project( self ):
        create_directory( args.output )
        # Create the various sub-directories
        for d in ['log_files', 'HBAR', 'references', 'subreads', 
                  'alignments', 'phasing', 'phasing_results', 
                  'resequencing', 'resequencing_results', 'results', 
                  'stats', 'annotation']:
            sub_dir = os.path.join( args.output, d )
            create_directory( sub_dir )
            setattr(self, d, sub_dir)

    def _initialize_logging( self ):
        log_format = "%(asctime)s [%(levelname)s] %(funcName)s %(message)s"
        log_file = os.path.join( self.log_files, "HLA_Pipeline.log" )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format,
                             stream=sys.stdout )
        logging.basicConfig( level=logging.INFO, 
                             format=log_format, 
                             filename=log_file )
    
    def run( self ):
        pass

if __name__ == '__main__':
    PhasingPipeline().run()
