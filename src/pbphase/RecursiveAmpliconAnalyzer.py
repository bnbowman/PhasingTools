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
from pbphase.AmpliconAnalyzer import AmpliconAnalyzer
from pbhla.io.extract_subreads import extract_subreads
from pbhla.external.commandline_tools import run_blasr
from pbhla.fasta.utils import fasta_size
from pbhla.utils import check_output_file
from pbhla.io.BlasrIO import BlasrReader

VALID_ARGS = ['minLength', 'minReadScore', 'sampleName', 'whiteList', 'noClustering']
MIN_READ_LENGTH = 3300
MIN_READ_SCORE = 0.8
MIN_SIZE
NPROC = 1

log = logging.getLogger()

class RecursiveAmpliconAnalyzer(object):
    """
    A tool for running Amplicon Analysis recursively on a complex data-set
    """

    def __init__(self, input_file,
                       output,
                       white_list=None,
                       setup=None, 
                       nproc=1, 
                       min_read_length=MIN_READ_LENGTH,
                       min_read_score=MIN_READ_SCORE):
        self._input_file = os.path.abspath( input_file )
        self._output = os.path.abspath( output )
        self._white_list = white_list
        self._setup = setup
        self._nproc = str(nproc)
        self._min_read_length = min_read_length
        self._min_read_score = min_read_score
        self._output_filelist = []
        self._count = 0
        # Validate the settings and use them to create an AA wrapper
        self._validate_settings()
        self._amplicon_analyzer = AmpliconAnalyzer( setup, nproc )

    def _validate_settings(self):
        """
        Validate input settings to ensure a local valid copy of ConsensusTools.sh
        """
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

    def run( self ):
        """
        Iteratively run AmpliconAnalyzer until it stops creating new clusters
        """
        filename = os.path.basename( self._input_file )
        log.info('Running RecursiveAmplcionAnalyzer on "%s"' % filename)
        create_directory( self._output )
        self.separate_alleles( self._white_list )
        log.info('Finished running RecursiveAmpliconAnalyzer\n')
        print self._count, self._output_filelist

    def separate_alleles( self, white_list ):
        # Run the first pass, with clustering
        log.info("Beginning iteration #%s" % self._count)
        print
        print self._count, self._output_filelist
        print
        curr_output = os.path.join( self._output, 'Iteration_%s' % self._count )
        output_file = amp_assem_output_exists( curr_output )
        if output_file:
            log.info('Existing output detected, skipping...')
        else:
            log.info('No existing output detected, proceeding ...')
            if self._count == 0: # For the first pass we enable clustering
                output_file = self.run_analysis( curr_output,
                                                 white_list,
                                                 cluster=True )
            else: # For all other iterations, we disable clustering
                output_file = self.run_analysis( curr_output,
                                                 white_list,
                                                 cluster=False )
        check_output_file( output_file )
        # Outputs of a single Fasta File are returned as is:
        log.info("Finished iteration #%s" % self._count)
        self._count += 1
        fasta_count = fasta_size( output_file )
        if fasta_count == 1:
            log.info('AmpliconAnalysis generated 1 cluster, exiting...')
            self.output_filelist.append( output_file )
            return 
        log.info('Amplicon Analysis generated %s clusters, continuing splitting' % fasta_count)
        # Otherwise we partition the reads and run the process on each partition
        alignment = self.align_subreads( white_list, output_file )
        groups = group_subreads( alignment )
        output_dir = os.path.dirname( output_file )
        sub_lists = []
        for reference, group in groups.iteritems():
            group_file = '%s.ids' % reference
            group_path = os.path.join( output_dir, group_file )
            write_whitelist( group, group_path )
            white_list_seqs = self.extract_whitelist_reads( group_path )
            sub_lists.append( white_list_seqs )
            if len(group) < MIN_SIZE:
                log.info('')
                continue
        for sub_list in sub_lists:
            self.separate_alleles( sub_list )

    def run_analysis( self, output, white_list, cluster=False ):
        """
        Run Amplicon Analysis on the input file with a given white list
        """
        analyzer_args = { 'sampleName': 'Iter%s' % self._count,
                          'minLength': self._min_read_length,
                          'minReadScore': self._min_read_score,
                          'noClustering': not cluster,
                          'whiteList': white_list }
        self._amplicon_analyzer.run( self._input_file, 
                                     output, 
                                     analyzer_args )
        # Check for the process's output file
        output_file = amp_assem_output_exists( output )
        if output_file:
            return output_file
        else:
            msg = 'No Amplicon Analysis output found in "%s"' % output_folder
            log.error( msg )
            raise IOError( msg )

    def align_subreads( self, white_list, reference_file ):
        """
        Align the subreads in a Whitelist to the created reference
        """
        basename = '.'.join( reference_file.split('.')[:-1] )
        alignment_file = '%s.m1' % basename
        reference_count = fasta_size( reference_file )
        blasr_args = { 'nproc': self._nproc,
                       'out': alignment_file,
                       'bestn': 1,
                       'nCandidates': reference_count,
                       'noSplitSubreads': True }
        run_blasr( white_list,
                   reference_file,
                   blasr_args )
        check_output_file( alignment_file )
        return alignment_file

    def extract_whitelist_reads( self, white_list_ids ):
        """
        Convert a White List of Ids into a White List of Sequences
        """
        root = '.'.join( white_list_ids.split('.')[:-1] )
        white_list_seqs = '%s.fasta' % root
        extract_subreads( self._input_file,
                          white_list_seqs,
                          self._min_read_length,
                          self._min_read_score,
                          white_list_ids )
        check_output_file( white_list_seqs )
        return white_list_seqs

def write_whitelist( group, output_file ):
    """
    Write the ZMWs for a list of reads out to a new Whitelist file
    """
    with open( output_file, 'w' ) as handle:
        for read in group:
            zmw = '/'.join( read.split('/')[:2] )
            handle.write( zmw + '\n' )

def group_subreads( alignment_file ):
    """
    Stuff
    """
    groups = {}
    for record in BlasrReader( alignment_file ):
        try:
            groups[record.tname].append( record.qname )
        except:
            groups[record.tname] = [ record.qname ]
    return groups

def amp_assem_output_exists( folder ):
    """
    Check whether an AA output folder already exists
    """
    output = os.path.abspath( folder )
    if not os.path.isdir( output ):
        return False
    output_file = os.path.join( output, 'amplicon_analysis.fasta')
    if os.path.isdir( output ) and os.path.exists( output_file ):
        return output_file
    return False

if __name__ == '__main__':
    import argparse
    
    logging.basicConfig( level=logging.INFO )

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

    raa = RecursiveAmpliconAnalyzer( args.input_file,
                                     args.output,
                                     args.white_list,
                                     args.setup, 
                                     args.nproc, 
                                     args.min_read_length,
                                     args.min_read_score )
    raa.run()
