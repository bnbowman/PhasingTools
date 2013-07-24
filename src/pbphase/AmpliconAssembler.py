#! /usr/bin/env python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os, re, sys
import logging
import subprocess

from pbphase.utils import *

log = logging.getLogger()

VALID_OPTIONS = ['minLength', 'minReadScore', 'whiteList', 'noClustering']

MIN_READ_LENGTH = 3000
MIN_READ_SCORE = 0.8
NPROC = 1

class AmpliconAssembler(object):
    """
    A tool for running Amplicon Assembly 
    """

    def __init__(self, setup=None, nproc=1):
        self._setup = setup
        self._nproc = str(nproc)
        self._validate_settings()

    def _validate_settings(self):
        self._consensus_tools = which('ConsensusTools.sh')
        # If the names are provided as a filename, parse it
        if self._setup and os.path.isfile( self._setup ):
            log.info('Using supplied Setup script for SMRTanalysis')
            self._setup = os.path.abspath( self._setup )
            self._use_setup = True
            self._consensus_tools = 'ConsensusTools.sh'
        elif self._consensus_tools:
            log.info('All required SMRT Analysis tools detected')
            self._use_setup = False
        else:
            msg = 'AmpliconAssembler requires EITHER a valid copy ' + \
                  'of ConsensusTools.sh in PATH or a path ' + \
                  'to a local SMRT Analysis setup script'
            log.error( msg )
            raise Exception( msg )

    def run(self, input_file, output='Stuff', args=None):
        """
        Create and execute a shell-script that runs the Amplicon Assembler
        """
        input_file = os.path.abspath( input_file )
        filename = os.path.basename( input_file )
        log.info('Running AmplcionAssembler on "%s"' % filename)

        # Check whether the output file exists
        output = os.path.abspath( output )
        create_directory( output )
        output_file = os.path.join( output, 'consensus-all.fasta')
        if os.path.isdir( output ) and os.path.exists( output_file ):
            log.info('Existing Fasta output detected, skipping...')
            return
        log.info('No existing output detected, running Amplicon Assembly...')

        # Create the
        process_args = self.create_process_args( self, input_file, output )
        process_args = self.add_optional_args( self, process_args, args )

        # Create and run the Amplicon Assembly process
        self.run_process( process_args, output, 'AmpliconAssembler')
        log.info('Finished running Amplicon Assembler\n')

    def create_process_args( self, input_file, output ):
        """
        Create a list of args with the minimally sufficient call to AA
        """
        filetype = get_file_type( input_file )
        if filetype == 'fofn':
            process_args = [self._consensus_tools,
                            'AmpliconAssembly',
                            '--fofn', input_file,
                            '--numThreads', self._nproc,
                            '--output', output]
        elif filetype == 'bash5':
            process_args = [self._consensus_tools,
                            'AmpliconAssembly',
                            input_file,
                            '--numThreads', self._nproc,
                            '--output', output]
        return process_args

    def add_optional_args( process_args, new_args ):
        """
        Add any remaining options to the Process
        """
        for arg, value in new_args.iteritems():
            if arg in VALID_ARGS:
                arg = '--' + arg
                if value is True:
                    process_args.append( arg )
                else:
                    process_args += [arg, str(value)]
        return process_args

    def run_process(self, process_args, folder, name):
        """
        Execute a tool as a python Subprocess
        """
        log.info("Executing child '%s' process" % name)
        if self._use_setup:
            log.info('Executing subprocess indirectly via Shell Script')
            script = self.write_script( process_args, folder, name)
            log_path = os.path.join( folder, name + '.log' )
            with open( log_path, 'w' ) as log_handle:
                p = subprocess.Popen( ['source', script], 
                                       executable='/bin/bash',
                                       stderr=subprocess.STDOUT,
                                       stdout=log_handle)
                p.wait()
        else:
            log.info('Executing subprocess directly via Subprocess')
            p = subprocess.Popen( process_args )
            p.wait()
        log.info('Child process finished successfully')

    def write_script( self, process_args, folder, name ):
        """
        Write out a shell script to execute the specified process
        """
        script_path = os.path.join( folder, name + '_script.sh' )
        with open( script_path, 'w') as handle:
            handle.write('source %s\n' % self._setup)
            handle.write( ' '.join(process_args) + '\n' )
        return script_path



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    add = parser.add_argument
    add('input_file', 
        metavar='FOFN', 
        help="BasH5 or FOFN of sequence data")
    add('--output', 
        metavar='DIR',
        default='Stuff',
        help="Specify a directory for intermediate files")
    add('--setup', 
        metavar='SETUP_FILE',
        help='Path to the SMRT Analysis setup script')
    add('--nproc', 
        type=int, 
        metavar='INT',
        default=NPROC, 
        help="Number of processors to use [%s]" % NPROC)
    add('-l', '--min_read_length',
        type=int,
        metavar='INT',
        default=MIN_READ_LENGTH,
        help="Minimum length sub-reads to allow as input [%s]" % MIN_READ_LENGTH)
    add('-s', '--min_read_score',
        type=float,
        metavar='FLOAT',
        default=MIN_READ_SCORE,
        help="Minimim read score to allow from input ZMWs [%s]" % MIN_READ_SCORE)
    add('--disable_clustering',
        action='store_true',
        help="Disable the pre-phasing coarse-clustering step")
    add('-w', '--white_list',
        metavar='FILE',
        help="Text or Fasta file of white-listed reads to use")
    args = parser.parse_args()
     
    logging.basicConfig( level=logging.INFO )

    aa = AmpliconAssembler( args.setup, args.nproc )
    process_args = {'whiteList': args.white_list,
                    'minLength': args.min_read_length,
                    'minReadScore': args.min_read_score,
                    'noClustering': args.disable_clustering}
    aa.run( args.input_file, args.output, args )
