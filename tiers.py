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
from lingpyd.align.sca import Alignments
import lingpyd.sequence.sound_classes as lsc



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
        return ['1' for m in tmp]
    else:
        return ['0' for m in tmp]

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
            tiers = [
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
                    ]
        # get the tiers
        self.tiers = {}
        for taxon in self.descendants:
            self.tiers[taxon] = {}
            # get the alignment for the two strings
            for k,(concept,almA) in self.words[self.proto].items():
                print(k,concept,almA)
                if k in self.words[taxon]:
                    almB = self.words[taxon][k][1]
                    # reduce alignment automatically
                    alm = [(s,t) for s,t in zip(almA, almB) if s != '-' and t
                            != '-']
                    # get tier from reduced alm
                    ralm = [x[0] for x in alm]
                    # get the tiers for the string
                    all_tiers = [ralm]
                    for tier,bigram in tiers:
                        print(tier, bigram, ' '.join(ralm))
                        tmp = get_tier(ralm, tier, bigram)
                        all_tiers += [tmp]
                    all_tiers += [[x[1] for x in alm]]
                    self.tiers[taxon][k] = all_tiers


if __name__ == '__main__':

    tiers = Tiers('germanic.tsv', proto="Proto-Germanic")
