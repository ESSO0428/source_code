#!/usr/bin/env python

import collections

Position = collections.namedtuple('Position', 'start end')


class MicroRNA:
    def __init__(self, seq=None):
        self.seq = seq

    def setSeq(self, seq):
        self.seq = seq

    def setCount(self, count):
        self.count = count

    def setStrand(self, strand):
        self.strand = strand

    def setPosition(self, start, end):
        self.position = Position(start=int(start), end=int(end))

    def getArm(self, hairpin_position):
        half_point = int(sum(hairpin_position)/2)
        if hairpin_position.start < self.position.start < half_point:
            self.arm = '5p' if self.strand == '+' else '3p'
        else:
            self.arm = '3p' if self.strand == '+' else '5p'

    def setArm(self, arm):
        self.arm = arm

    def __str__(self):
        return 'MicroRNA {}'.format(self.seq)


class Hairpin:
    def __init__(self, seq=None):
        self.seq = seq

    def setSeq(self, seq):
        self.seq = seq

    def setStruct(self, struct):
        self.struct = struct

    def setStrand(self, strand):
        self.strand = strand

    def setPosition(self, start, end):
        self.position = Position(start=int(start), end=int(end))

    def __str__(self):
        return 'Hairpin {}'.format(self.seq)


class Mature(MicroRNA):
    pass


class Star(MicroRNA):
    pass
