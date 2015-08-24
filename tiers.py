# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-08-24 17:23
# modified : 2015-08-24 17:23
"""
Multitiered Sequence Representation for LingPy.
"""

__author__="Johann-Mattis List"
__date__="2015-08-24"

import lingpyd as lp
from lingpyd.algorithm import misc
from lingpyd.align.sca import Alignments
import lingpyd.sequence.sound_classes as lsc
from collections import OrderedDict



def tonal_tiers(tokens):
    """
    Assign tonal values to all members of a syllable.
    """
    out = []
    tone = '0'
    for x,y in list(zip(lp.tokens2class(tokens, lp.rc('cv')), tokens))[::-1]:
        if x == 'T':
            tone = y
            out += [tone]
        else:
            out += [tone]
    return out

def nasality_of_word(tokens):
    """
    Assign nasality values to all members of a word.
    """
    tmp = lp.tokens2class(lp.ipa2tokens(tokens, expand_nasals=True), lp.rc('dolgo'))
    if 'N' in tmp or 'M' in tmp:
        return ['1' for m in tokens]
    else:
        return ['0' for m in tokens]

def nasality_of_syllable(tokens):
    """
    Assign nasality values to all members of a syllable.
    """
    # get syllables first
    syllables = [[]]
    for k in lsc.syllabify(tokens):
        if k == lp.rc('morpheme_separator'):
            syllables += [[]]
        else:
            syllables[-1] += [k]
    # now check each syllable for nasality
    out = []
    for syllable in syllables:
        nasal = '0'
        tks = lp.tokens2class(
                lp.ipa2tokens(syllable, expand_nasals=True),
                lp.rc('dolgo')
                )
        if 'N' in tks or 'M' in tks:
            nasal = '1'
        for s in syllable:
            out += [nasal]
    return out

def get_tier(tokens, tier='sc', bigram=False):
    """
    Function produces default tiers for a given word.
    """
    # dictionary stores different functions for conversion. initially, we
    # create only sound-class systems, the decision for bigrams is done in a
    # second step
    D = dict(
            cv = lambda x:     lp.tokens2class(x, lp.rc("cv")),
            dolgo = lambda x:  lp.tokens2class(x, lp.rc("dolgo")),
            sca = lambda x:    lp.tokens2class(x, lp.rc("sca")),
            art = lambda x:    lp.tokens2class(x, lp.rc("art")),
            asjp = lambda x:   lp.tokens2class(x, lp.rc("asjp")),
            tone = lambda x: tonal_tiers(x),
            word_nasality = lambda x: nasality_of_word(x), 
            syllable_nasality = lambda x: nasality_of_syllable(x)
            )
    # convert the tiers according to the function
    tiers = D[tier](tokens)
    # check for bigram-option
    if bigram in [1, 'pre', 'preceeding']:
        tiers = ['#'] + tiers[:-1]
    elif bigram in [2, 'post', 'following']:
        tiers = tiers[1:] + ['$']
    return tiers

def check_consistency(tier_matrix):
    """
    Function checks whether a given tier matrix will uniquely produce the
    reflex sounds.
    """
    pass

def get_numeric_matrix(matrix):
    """
    Convert a value matrix into a numeric matrix.
    """
    # get a converter for all columns in the matrix
    out = [[0 for cell in row] for row in matrix]
    for i in range(len(matrix[0])):
        col = [line[i] for line in matrix]
        idx = 1
        tmp = {}
        for j,cell in enumerate(col):
            if cell in tmp:
                out[j][i] = tmp[cell]
            else:
                idx += 1
                tmp[cell] = idx
                out[j][i] = tmp[cell]
    return out

class Tiers(Alignments):

    def __init__(self, infile, proto, ref='cogid'):
        """
        Initialize the data.
        """
        Alignments.__init__(self, infile, ref="cogid")
        self.proto = proto
        # get all words from the proto-language and the daughter languages
        self.words = {}
        for taxon in self.taxa:
            self.words[taxon] = dict([
                (self[k,ref],
                    (self[k,'concept'], self[k,'alignment'])) 
                for k in self if self[k,'language'] == taxon])
        # define descendants for convenience
        self.descendants = [t for t in self.taxa if t != proto]
        # get all proto-sounds (for purpose of convenience only)
        self.sounds = []
        for key,(concept,word) in self.words[self.proto].items():
            self.sounds += word
        self.sounds = sorted(set(self.sounds))

    def make_tiers(self, tiers=None):
        """
        Create tier system for a given language.
        """
        # define the default tier system
        if not tiers:
            self.tier_system = OrderedDict(enumerate([
                    ('cv',1), 
                    ('cv',2),
                    ('dolgo',1),
                    ('dolgo',2),
                    ('sca',1),
                    ('sca',2),
                    ('art',1),
                    ('art',2),
                    ('asjp',1),
                    ('asjp',2),
                    ('word_nasality', 0),
                    ('syllable_nasality', 0),
                    ('tone', 0),
                    ]))
        else:
            self.tier_system = OrderedDict(enumerate(tiers))
        # get the tiers
        self.tiers = {}
        for taxon in self.descendants:
            self.tiers[taxon] = {}
            # get the alignment for the two strings
            for k,(concept,almA) in self.words[self.proto].items():
                if k in self.words[taxon]:
                    almB = self.words[taxon][k][1]
                    # reduce alignment automatically
                    alm = [(s,t) for s,t in zip(almA, almB) if s != '-' and t
                            != '-']
                    # get tier from reduced alm
                    ralm = [x[0] for x in alm]
                    # get the tiers for the string
                    all_tiers = [ralm]
                    for tier,bigram in self.tier_system.values():
                        tmp = get_tier(ralm, tier, bigram)
                        all_tiers += [tmp]
                    all_tiers += [[x[1] for x in alm]]
                    self.tiers[taxon][k] = all_tiers
        # now create the specific matrix representation for each taxon. we
        # start with a value matrix which contains the values
        self.value_matrices = {}
        self.matrices = {}
        for taxon in self.descendants:
            self.value_matrices[taxon] = {}
            self.matrices[taxon] = {}
            for k,tier in self.tiers[taxon].items():
                # get the columsn by transposing
                cols = misc.transpose(tier)
                for col in cols:
                    # get the three aspects of the 
                    pf,tr = col[0],col[1:]
                    try:
                        self.value_matrices[taxon][pf] += [tr]
                    except KeyError:
                        self.value_matrices[taxon][pf] = [tr]
            for pf,matrix in self.value_matrices[taxon].items():
                self.value_matrices[taxon][pf] = sorted(matrix, 
                        key=lambda x: x[-1]) 
                self.matrices[taxon][pf] = misc.transpose(get_numeric_matrix(
                        self.value_matrices[taxon][pf]))

    def remove_redundant_tiers(self):
        """
        Function checks for redundant tiers and removes them.
        """
        pass
    


if __name__ == '__main__':

    tiers = Tiers('germanic.tsv', proto="Proto-Germanic")
