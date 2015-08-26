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
import itertools



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
    tmp = lp.tokens2class(lp.ipa2tokens(tokens, expand_nasals=True,
        merge_geminates=False), lp.rc('dolgo'))
    if 'N' in tmp or 'M' in tmp:
        return ['1' for m in tokens]
    else:
        return ['0' for m in tokens]

def get_syllables(tokens):
    """
    Simple function to get the syllables in tokenized form and nested.
    """
    # get syllables first
    syllables = [[]]
    for k in lsc.syllabify(tokens):
        if k == lp.rc('morpheme_separator'):
            syllables += [[]]
        else:
            syllables[-1] += [k]
    return syllables

def nasality_of_syllable(tokens):
    """
    Assign nasality values to all members of a syllable.
    """
    syllables = get_syllables(tokens)
    # now check each syllable for nasality
    out = []
    for syllable in syllables:
        nasal = '0'
        tks = lp.tokens2class(
                lp.ipa2tokens(syllable, expand_nasals=True,
                    merge_geminates=False),
                lp.rc('dolgo')
                )
        if 'N' in tks or 'M' in tks:
            nasal = '1'
        for s in syllable:
            out += [nasal]
    return out

def position(tokens):
    """
    Simple algorithm to determine the position of a syllable in a word.
    """
    # get syllables first
    syllables = get_syllables(tokens)
    # get the ranges
    p_range = list(range(1,len(syllables)+1))
    n_range = [-i for i in p_range][::-1]
    # get half of the length in terms of syllables
    sylen = int(len(syllables) / 2 + 0.5)
    # combine the ranges
    c_range = p_range[:sylen] + n_range[sylen:]
    out = []
    for r,syl in zip(c_range, syllables):
        for token in syl:
            out += [r]
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
            prostring = lambda x: lp.prosodic_string(x),
            trigram = lambda x: list(zip(['$']+x[:-1], x, x[1:]+['$'])),
            asjp = lambda x:   lp.tokens2class(x, lp.rc("asjp")),
            tone = lambda x: tonal_tiers(x),
            sequence = lambda x: x,
            word_nasality = lambda x: nasality_of_word(x), 
            syllable_nasality = lambda x: nasality_of_syllable(x),
            position = lambda x: position(x),
            exception = lambda x: [''.join(x) for y in x],
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
        idx = 0
        tmp = {}
        for j,cell in enumerate(col):
            if cell in tmp:
                out[j][i] = tmp[cell]
            else:
                idx += 1
                tmp[cell] = idx
                out[j][i] = tmp[cell]
    return out

def reduce_numeric_matrix(matrix):
    """
    Reduce a given matrix by assembling identical rows into one cluster.
    """
    idxs = OrderedDict()
    for i in range(len(matrix)-1):
        line = matrix[i]
        try:
            idxs[tuple(line)] += [i]
        except KeyError:
            idxs[tuple(line)] = [i]
    active = [list(line) for line in idxs if len(set(line)) != 1]
    indices = [idxs[tuple(line)][0] for line in active]
    return active+[matrix[-1]], indices

    #return [list(line) for line in idxs], [idx[0] for idx in
    #        idxs.values()][:-1]


class Tiers(Alignments):

    def __init__(self, infile, proto, ref='cogid'):
        """
        Initialize the data.
        """
        Alignments.__init__(self, infile, ref="cogid")
        self.proto = proto
        self.ref = ref
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
                    #('cv',1), 
                    #('cv',2),
                    #('dolgo',1),
                    #('dolgo',2),
                    #('sca',1),
                    #('sca',2),
                    #('art',1),
                    #('art',2),
                    #('asjp',1),
                    #('asjp',2),
                    #('trigram',0),
                    ('sequence', 1),
                    ('sequence', 2),
                    #('word_nasality', 0),
                    ('syllable_nasality', 0),
                    #('tone', 0),
                    ('position', 0),
                    ('prostring', 0),
                    #('exception', 0),
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
        self.active_tiers = {}
        self.instances = {}
        for taxon in self.descendants:
            self.value_matrices[taxon] = {}
            self.matrices[taxon] = {}
            self.instances[taxon] = {}
            self.active_tiers[taxon] = {}
            for k,tier in self.tiers[taxon].items():
                # get the columsn by transposing
                cols = misc.transpose(tier)
                for col in cols:
                    # get the two aspects of the col (proto and rest)
                    pf,tr = col[0],col[1:]
                    try:
                        self.value_matrices[taxon][pf] += [tr]
                        self.instances[taxon][pf] += [k]
                    except KeyError:
                        self.value_matrices[taxon][pf] = [tr]
                        self.instances[taxon][pf] = [k]
            # re-iterate and create the matrix representation
            for pf,matrix in self.value_matrices[taxon].items():
                self.value_matrices[taxon][pf] = sorted(matrix, 
                        key=lambda x: x[-1]) 
                matrix = misc.transpose(get_numeric_matrix(
                        self.value_matrices[taxon][pf]))
                new_matrix, indices = reduce_numeric_matrix(matrix)
                self.matrices[taxon][pf] = new_matrix
                # assign active tiers (will be reduced later on)
                self.active_tiers[taxon][pf] = [self.tier_system[idx]
                        for idx in indices]

    def check_tiers(self, taxon, combinations=1):
        
        # make a smart search here:
        # determine regular patterns in a first run, only then start and go for
        # the irregular patterns! we want full explanation, not necessarily the
        # best of all contexts (and we can always merge regular contexts!)
        explanations = {}
        regular = 0
        irregular = 1
        for sound in self.sounds:
            if True: #lp.tokens2class(sound, 'dolgo')[0] != 'V':

                print('[i] calculate patterns for sound {0}'.format(sound))
                check = self.remove_redundant_tiers(taxon, sound,
                        combinations=combinations)
                if check:
                    explanation, reflexes = check
                    for key in explanation:
                        print("... {0} > {1}".format(sound, reflexes[key]))
                        print("... regular: {0}\n... irregular: {1}".format(
                            explanation[key][0], explanation[key][1]))
                        explanations[sound, reflexes[key]] = explanation[key][0:2]
                        regular += explanation[key][0]
                        irregular += explanation[key][1]
        total = regular + irregular
        print("Total: {0}, Regular: {1} / {2:.2f}, Irregular: {3} / {4:.2f}".format(
                    regular+irregular, regular, regular / total, irregular,
                    irregular / total))
        regular_types = len([k for k in explanations if explanations[k][1] ==
            0])
        print(regular_types, regular_types / len(explanations))
        return explanations

    def remove_redundant_tiers(self, taxon, sound, combinations=1, verbose=False):
        """
        Function checks for redundant tiers and removes them.
        
        It removes: 

        * tiers which are identical (only one is retained), and
        * tiers which are indiscriminative (they do not distinguish at all)


        # get the exceptions for a pattern, if it can be confirmed that this
        # doens't work. These is the set of cases which overlap in two
        # contexts. We still need to figure out, how they can be identified.
        
        """
        # make sure to keep track of each of the original tiers for
        # back-tracking
        if sound not in self.matrices[taxon]:
            return
        
        matrix = self.matrices[taxon][sound]
        # return if matrix is of length one
        if len(set(matrix[-1])) == 1:
            return {0 : (
                len(matrix[-1]), 0, [('contextless',0)],[i for i in
                self.instances[taxon][sound]])},[self.value_matrices[taxon][sound][0][-1]]

        atiers = self.active_tiers[taxon][sound]
        instances = self.instances[taxon][sound]
        values = self.value_matrices[taxon][sound]
        # get the basic partition of the matrix
        current = 0
        idxs = []
        for i,cell in enumerate(matrix[-1]):
            if cell != current:
                current = cell
                idxs += [i]
        idxs += [len(instances)]
        idxs = list(zip(idxs[:-1], idxs[1:]))
        # break up the matrix into the sub-parts
        reflexes = [values[a:b][0][-1] for a,b in idxs]
        # our explanations distinguish regular and irregular characters
        explanations = dict([(i,(0, len(matrix[0]),[('exception',0)],[])) for i in range(len(idxs))])
        
        # check the minimization in this context: what do we actually want to
        # minimize? Is it the number of regular tokens, is it the number of
        # good items, is it the number of items that are false?
        # maybe it would just be enough to minimize the irregularities!
        for i in range(1,int(combinations * len(matrix)-1 + 0.5)): #-1
            for indices in itertools.combinations(range(len(matrix)-1), i):
                line = list(zip(*[matrix[jj] for jj in indices]))
                for j,(a,b) in enumerate(idxs):
                    # get the current part, the rest, and the overlap
                    part = line[a:b]
                    rest = line[:a]+line[b:]
                    regular_tokens = [(p,nr) for p,nr in zip(part,range(a,b)) if p not in rest]
                    regular = len([p[0] for p in regular_tokens])
                    irregular_tokens = [(p,nr) for p,nr in zip(part, range(a,b))
                            if p in rest]
                    irregular = len([p[0] for p in irregular_tokens])
                    if verbose:
                        print(regular, irregular, reflexes[j])
                    if regular > explanations[j][0]:
                        explanations[j] = (regular, irregular, [atiers[jj] for
                            jj in indices], 
                            [self.instances[taxon][sound][n[1]] for n in irregular_tokens])
                    elif regular == explanations[j][0]:
                        if irregular < explanations[j][1]:
                            explanations[j] = (regular, irregular, [atiers[jj]
                                for jj in indices], 
                                [self.instances[taxon][sound][n[1]] for n in
                                    irregular_tokens])
        if verbose:
            tidx = self.taxa.index(taxon)
            for k in explanations:
                print(reflexes[k])
                print('---')
                print(explanations[k][0], explanations[k][1])
                print(', '.join(['.'.join(self.msa[self.ref][x]['seqs'][self.msa[self.ref][x]['taxa'].index(taxon)]) for x in explanations[k][-1]]))
        return explanations, reflexes
             

if __name__ == '__main__':
    

    tiers = Tiers('tukanox.tsv', proto="*PT")
    tiers.make_tiers()
