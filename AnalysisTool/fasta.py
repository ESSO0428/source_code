#!/usr/bin/env python
"""
Python Object for process sequence :

    1. `method (def)` is separate from `class`
    2. `def` use template and attribute of `class` to process data
---
Abstract:
    - class : `Dna` is all `class` core (now : only hava `Fasta`) <exclude after 23019 Andy6 update>
    - read fasta file to dict (core class is Fasta)
    - reverse complement input sequence
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Union
from pandas.core.frame import DataFrame
import pandas as pd
from .app import *

def readAsDict(fasta_file: str, force_upper: bool = True, ignore_seq_info: bool = True) -> Dict[str, str]:
    """
    read fasta into dict :
        reference : https://www.biostars.org/p/710/#1414  \n
        input : fasta file path  \n
        output : { '>seq_id' : 'seq'  } \n
        example :
        >>> dictGenome = fasta.readAsDict("..genome.fa")
        >>> dictGenome
        {'L01': 'ATCG....', ... }

        p.s. : will remove detail (default)
            ex: >L01 chr1 -> L01 \n

            set `ignore_seq_info` to `False` to keep detail \n
            ex: >L01 chr1 -> L01 chr1
    ---
        args : 
            fasta_file : fasta file path
            force_upper : true/false (force to upper case)
                atcg -> ATCG
            ignore_seq_info : true/false (remove detail)
                default : true
                    ex: >L01 chr1 -> L01
                set `false` will keep detail
                    ex: >L01 chr1 -> L01 chr1
    """

    SeqDict = {}
    with open(fasta_file, 'r') as f:
        for record in Fasta(f):
            #SeqDict[record.head] = \
            # new : remove detail ex: >L01 chr1 -> >L01
            chr=record.head
            if ignore_seq_info:
                chr=chr.split(' ')[0]
            SeqDict[chr] = \
                record.sequence.upper() if force_upper else record.sequence
    return SeqDict


def reverse_complement(dna):
    """
    reverse complement input sequence :
        input : dna sequence
        output : reverse complement sequence
        example :
        >>> reverse_complement("ATCG")
        'CGAT'
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join([complement.get(base, base) for base in dna[::-1]])

# 230907 add new def for process dict of fasta (for degradome project)
def convert_dict_fasta_sequence(dictFasta: Dict[str, str], how: str = "length") -> Dict[str, int] | Dict[str, str]:
    """
    calculate string length of each items in dict :
        input : dict of fasta (or other string)
            `key : value`
        output : 
            `key : len(value)`
        example :
        >>> dictGenome = fasta.readAsDict("..genome.fa")
        >>> dictGenome
        {'L01': 'ATCG', ... }
        >>> calculate_items_len(dictGenome)
        {'L01': 4, ...}
    """
    if how == "length":
        return {key: len(value) for key, value in dictFasta.items()}
    else:
        return dictFasta
    


class Dna:
    ''' Object representing a FASTA record. '''

    def __init__(self, header, sequence):
        self.head = header
        self.seq = sequence

    def __repr__(self):
        return '[HTML]' % (self.head)

    def __str__(self, separator=''):
        return '>%s\n%s' % (self.head, separator.join(self.seq))

    def __len__(self):
        return len(''.join(self.seq))

    @property
    def sequence(self, separator=''):
        return separator.join(self.seq)


class Fasta:
    ''' A FASTA iterator/generates DNA objects. '''

    def __init__(self, handle):
        self.handle = handle

    def __repr__(self):
        return '[HTML]' % self.handle

    def __iter__(self):
        header, sequence = '', []
        for line in self.handle:
            if line[0] == '>':
                if sequence:
                    yield Dna(header, sequence)
                header = line[1:-1]
                sequence = []
            else:
                sequence.append(line.strip())
        yield Dna(header, sequence)

