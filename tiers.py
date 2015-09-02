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
import networkx as nx
import json

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

def word_length(tokens):
    
    syllables = get_syllables(tokens)
    out = []
    slen = len(syllables)
    
    return [slen for x in tokens]

def nasality_of_word(tokens):
    """
    Assign nasality values to all members of a word.
    """
    tmp = lp.ipa2tokens(''.join(tokens), expand_nasals=True)
    nasals = 'nȵɳŋɴňńɲ∼mɱʍ'
    if lp.rc('nasal_placeholder') in tmp:
        return ['1' for m in tokens]
    elif [n for n in nasals if n in set(''.join(tokens))]:
        return ['1' for m in tokens]

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

def reduce_patterns(pattern, verbose=True):

    matrix = []
    reflexes = []
    for p,t,l in pattern:
        matrix += [list(p) + [t]]
        reflexes += [l]

    numatrix = get_numeric_matrix(matrix)
    rematrix,idxs = _reduce_numeric_matrix(misc.transpose(numatrix))
    rematrix = misc.transpose(numatrix)

    def compatible(pt1,pt2):
        """
        Define the compatibility of two reflex sets.
        """
        
        for a,b in zip(pt1,pt2):
            if a > 0 and b > 0:
                return False
        return True

    # check exhaustively for combinations in which tiers are missing
    # we proceed agglomeratively and stop when a good tier system is found that
    # explains everything and is compatible
    good_matrix = False
    for i in range(1,len(rematrix)):
        for idxs in itertools.combinations(range(len(rematrix)-1),i):
            print(idxs)
            bad_matrix = False
            new_matrix = [rematrix[idx] for idx in idxs] + [rematrix[-1]]

            # check for consistency
            tmp = {}
            for j in range(len(new_matrix[0])):
                col = tuple([line[j] for line in new_matrix][:-1])
                target = rematrix[-1][j]
                ptn = reflexes[j]

                if col not in tmp:
                    tmp[col] = [(target,ptn)]
                else:
                    tmp[col] += [(target,ptn)]

            # check the consistency
            for key,val in tmp.items():
                for (targetA,ptnA),(targetB,ptnB) in itertools.combinations(val,2):
                    if targetA != targetB and not compatible(ptnA,ptnB):
                        bad_matrix = True
                        break
                if bad_matrix:
                    break
            if bad_matrix:
                pass
            else:
                good_matrix = idxs
                break
        if bad_matrix:
            pass
        elif good_matrix:
            break

    # unique cols
    unique_cols = {}
    for i in range(len(new_matrix[0])):
        col = tuple([new_matrix[j][i] for j in range(len(new_matrix)-1)])
        try:
            unique_cols[col] += [i]
        except KeyError:
            unique_cols[col] = [i]
     
    # transform back now
    out = dict(
            conditions = [],
            targets = [],
            reflexes = [],
            totals = [],
            tiers = idxs
            )
    for key in sorted(unique_cols):
        
        new_refs = []
        for k in unique_cols[key]:
            new_refs += [reflexes[k]]
        ref_sum = []
        for i in range(len(new_refs[0])):
            ref_sum += [sum([row[i] for row in new_refs])]
        if verbose:
            print('  '.join([str(matrix[unique_cols[key][0]][idx]) for idx in idxs]),'\t',
                    '/'.join([matrix[i][-1] for i in unique_cols[key]]),'\t', 
                    '  '.join([str(x) for x in ref_sum+[sum(ref_sum)]]))
        
        out['conditions'] += [[matrix[unique_cols[key][0]][idx] for idx in
            idxs]]
        out['targets'] += ['/'.join([matrix[i][-1] for i in unique_cols[key]])]
        out['reflexes'] += [ref_sum]
        out['totals'] += [sum(ref_sum)]
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
            prostring = lambda x: list(lp.prosodic_string(x)),
            trigram = lambda x: list(zip(['$']+x[:-1], x, x[1:]+['$'])),
            asjp = lambda x:   lp.tokens2class(x, lp.rc("asjp")),
            tone = lambda x: tonal_tiers(x),
            sequence = lambda x: x,
            word_nasality = lambda x: nasality_of_word(x), 
            syllable_nasality = lambda x: nasality_of_syllable(x),
            position = lambda x: position(x),
            exception = lambda x: [''.join(x) for y in x],
            word_length = lambda x: word_length(x), 
            )
    # convert the tiers according to the function
    tiers = D[tier](tokens)
    # check for bigram-option
    if bigram in [1, 'pre', 'preceeding']:
        tiers = ['#'] + tiers[:-1]
    elif bigram in [2, 'post', 'following']:
        tiers = tiers[1:] + ['$']
    
    if tier in ['dolgo','sca','asjp','cv']:
        for i,(t,c) in enumerate(zip(tokens,tiers)):
            if '-' in tokens and c == '0':
                tiers[i] = 'Ø'
            elif c == '0':
                print(t,c,tiers,tier,tokens)
                input()

    if '' in tiers:
        print(tier, tiers, tokens)
        input()
    
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
    return matrix, list(range(len(matrix)-1))

def _reduce_numeric_matrix(matrix):
    idxs = OrderedDict()
    for i in range(len(matrix)-1):
        line = matrix[i]
        try:
            idxs[tuple(line)] += [i]
        except KeyError:
            idxs[tuple(line)] = [i]
    active = [list(line) for line in idxs if len(set(line)) != 1]
    if not active:
        active = [list(idxs.keys())[0]]

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
        self.descendants = [t for t in self.taxa if t != proto]
        # get all words from the proto-language and the daughter languages
        self.words = {}
        # get the cognate ids of the proto-language
        cogids = self.get_list(taxon=proto, entry=ref, flat=True)
        for taxon in self.taxa:
            self.words[taxon] = dict([
                (self[k,ref],
                    [k, self[k,'concept'], self[k,'alignment']]) 
                for k in self if self[k,'language'] == taxon and self[k,ref] in
                cogids])
        # define descendants for convenience
        self.descendants = [t for t in self.taxa if t != proto and
                len(self.words[t]) >= 100]
        # get all proto-sounds (for purpose of convenience only)
        self.sounds = []
        for key, (idx, concept, word) in self.words[self.proto].items():
            self.sounds += word
        self.sounds = sorted(set(self.sounds))


        # check sound class conversions
        warning = False
        for k in self:
            tokens = [x for x in self[k,'tokens'] if x != '-']
            classes = lp.tokens2class(tokens, 'sca')
            errors = []
            for t,c in zip(tokens,classes):
                if c == '0' and c not in errors:
                    print('[WARNING] Unrecognized symbol {0} (ID: {1})'.format(t, k))
                    warning = True
                    errors += [c]
        if warning:
            input('Press Return to continue')
        

    def make_tiers(self, tiers=None, exclude_gaps=False):
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
                    #('art',1),
                    #('art',2),
                    ('asjp',1),
                    ('asjp',2),
                    ##('trigram',0),
                    #('sequence', 1),
                    #('sequence', 2),
                    #('word_nasality', 0),
                    #('syllable_nasality', 0),
                    ('word_length', 0),
                    #('tone', 0),
                    ('position', 0),
                    ('prostring', 0),
                    #('exception', 0),
                    ]))
            self.tier_index = OrderedDict([(v,k) for k,v in
                self.tier_system.items()])
        else:
            self.tier_system = OrderedDict(enumerate(tiers))
            self.tier_index = OrderedDict([(v,k) for k,v in
                self.tier_system.items()])
        allowed_gaps = 1 if exclude_gaps else 2
        # get the tiers
        self.tiers = {}
        for taxon in self.descendants:
            self.tiers[taxon] = {}
            # get the alignment for the two strings
            for k, (idx, concept, almA) in self.words[self.proto].items():
                if k in self.words[taxon]:
                    almB = self.words[taxon][k][2]
                    # reduce alignment automatically
                    try:
                        alm = lp.read.qlc.reduce_alignment([almA,almB])
                    except IndexError:
                        raise ValueError("Miscoded alignments {0} and {1}".format(
                            ' '.join(almA), ' '.join(almB)))

                    alm = [s for s in zip(alm[0], alm[1])]
                    if alm:
                        #alm = lp.read.qlc.reduce_alignment(alm)
                        alm = [s for s in alm if s.count('-') < allowed_gaps]
                        # update the words
                        self.words[taxon][k] += [misc.transpose(alm)]
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
                for i,col in enumerate(cols):
                    # get the two aspects of the col (proto and rest)
                    pf,tr = col[0],col[1:]
                    try:
                        self.value_matrices[taxon][pf] += [tr]
                        self.instances[taxon][pf] += [(k,i)]
                    except KeyError:
                        self.value_matrices[taxon][pf] = [tr]
                        self.instances[taxon][pf] = [(k,i)]
            # re-iterate and create the matrix representation
            for pf,matrix in self.value_matrices[taxon].items():
                sorter = [x[0] for x in sorted(zip(range(len(matrix)),matrix),
                        key=lambda x: x[1][-1])]
                self.value_matrices[taxon][pf] = [matrix[i] for i in sorter]
                self.instances[taxon][pf] = [self.instances[taxon][pf][i] for i
                        in sorter]
                matrix = misc.transpose(get_numeric_matrix(
                        self.value_matrices[taxon][pf]))
                new_matrix, indices = reduce_numeric_matrix(matrix)
                if len(new_matrix) == 1:
                    print(matrix,pf)
                    input(new_matrix)
                self.matrices[taxon][pf] = new_matrix
                # assign active tiers (will be reduced later on)
                self.active_tiers[taxon][pf] = [self.tier_system[idx]
                        for idx in indices]

    def check_tiers_per_language(self, taxon, combinations=1,
            ignore_vowels=False, verbose=False):
        
        # make a smart search here:
        # determine regular patterns in a first run, only then start and go for
        # the irregular patterns! we want full explanation, not necessarily the
        # best of all contexts (and we can always merge regular contexts!)
        explanations = {}
        regular = 0
        irregular = 1
        # determine current set of sounds
        sounds = [s for s in self.sounds if lp.tokens2class(sound, 'dolgo')[0]
                != 'V'] if ignore_vowels else self.sounds
        for sound in self.sounds:
            if verbose:
                print('[i] calculate patterns for sound {0}'.format(sound))
            check = self.remove_redundant_tiers(taxon, sound,
                    combinations=combinations)
            if check:
                explanation, reflexes = check
                for key in explanation:
                    if verbose:
                        print("... {0} > {1}".format(sound, reflexes[key]))
                        print("... regular: {0}\n... irregular: {1}".format(
                            explanation[key][0], explanation[key][1]))
                    explanations[sound, reflexes[key]] = (
                            explanation[key][0],
                            explanation[key][1],
                            explanation[key][3],
                            explanation[key][4],
                            explanation[key][2],
                            )
                    regular += explanation[key][0]
                    irregular += explanation[key][1]
        if verbose:
            total = regular + irregular
            print("Total: {0}, Regular: {1} / {2:.2f}, Irregular: {3} / {4:.2f}".format(
                        regular+irregular, regular, regular / total, irregular,
                        irregular / total))
            regular_types = len([k for k in explanations if explanations[k][1] ==
                0])
            print(regular_types, regular_types / len(explanations))

        return explanations

    def check_tiers(self, combinations=1, ignore_vowels=False):
        
        # iterate over taxa, and start counting
        self.explanations = {}
        self.changes = {}
        self.regular_changes = {}
        self.irregular_changes = {}
        for taxon in self.descendants:
            self.regular_changes[taxon] = {}
            self.irregular_changes[taxon] = {}
            explanations = self.check_tiers_per_language(taxon,
                    combinations=combinations,
                    ignore_vowels=ignore_vowels)
            self.explanations[taxon] = explanations
            for s,t in explanations:
                for cid in explanations[s,t][2]:
                    # get the value
                    if explanations[s,t][4][0][0] != 'contextless':
                        idx = self.instances[taxon][s].index(cid)
                        tier = tuple([(y,self.value_matrices[taxon][s][idx][x])
                            for y,x in
                                zip(explanations[s,t][4], [self.tier_index[z]
                                    for z in
                                    explanations[s,t][4]])])
                    else:
                        tier = ((('contextless',0),''),)

                    try:
                        self.changes[s,t][tier] += 1
                    except KeyError:
                        try:
                            self.changes[s,t][tier] = 1
                        except KeyError:
                            self.changes[s,t] = {tier : 1}
                    # append changes to 
                    try:
                        self.regular_changes[taxon][s,t][tier][0] += 1
                        self.regular_changes[taxon][s,t][tier][1] += [cid]
                    except KeyError:
                        try:
                            self.regular_changes[taxon][s,t][tier] = [1,[cid]]
                        except KeyError:
                            self.regular_changes[taxon][s,t] = {tier :
                                    [1,[cid]]}

        with open('output.tsv', 'w') as f:
            for s,t in sorted(self.changes):
                for tiers in sorted(self.changes[s,t]):
                    f.write('{0} > {1} ({2}) \t/ {3}\n'.format(
                        s,t, self.changes[s,t][tiers],
                        '\t'.join([tt[0][0]+'/'+str(tt[0][1])+':'+str(tt[1]) for tt in tiers])
                        ))
    def output_changes(self):
        """
        Write data to file.
        """
        
        _translate = {
                ('cv', 1): 'PRECEDING_CV',
                ('cv', 2): 'FOLLOWING_CV',
                ('sequence', 1) : 'PRECEDING',
                ('sequence', 2) : 'FOLLOWING',
                ('word_nasality', 0): 'NASALITY',
                ('syllable_nasality',0): 'NASAL_SYLLABLE',
                ('position', 0): 'POSITION',
                ('prostring', 0): 'PROSODY',
                ('asjp', 1): 'PRECEDING_ASJP',
                ('asjp',2): 'FOLLOWING_ASJP',
                ('dolgo',1): 'PRECEDING_DOLGO',
                ('dolgo',2): 'FOLLOWING_DOLGO',
                ('art',1):'PRECEDING_SONORITY',
                ('art',2): 'FOLLOWING_SONORITY',
                }
        translate = lambda x: _translate[x] if x in _translate else str(x[0])+'/'+str(x[1])

        D = [['Source'] + [translate(tier) for tier in self.tier_system.values()] +
                ['Target'] +self.descendants + ['FREQUENCY']]

        # get the changes
        changes = []
        patterns = {}
        for taxon in self.descendants:
            for s,t in self.regular_changes[taxon]:
                for specific_conditions,frequency in self.regular_changes[taxon][s,t].items():
                    conditions = dict(specific_conditions)
                    tmp = [s]
                    for tier in self.tier_system.values():
                        if tier in conditions:
                            tmp += [conditions[tier]]
                        else:
                            tmp += ['']
                    tmp += [t]
                    tmp = tuple(tmp)
                    try:
                        patterns[tmp][taxon] = frequency[0] 
                    except KeyError:
                        patterns[tmp] = {taxon: frequency[0]}
        
        P = {}
        for pattern in sorted(patterns, key=lambda x: tuple([x[0]]+[str(y) for
            y in x[1:]])):
            if lp.tokens2class(pattern[0], 'dolgo')[0] != 'V' and pattern[0] != '-':
                tmp = list(pattern)
                score = 0
                reflexes = []
                for taxon in self.descendants:
                    try:
                        tmp += [patterns[pattern][taxon]]
                        score += patterns[pattern][taxon]
                        reflexes += [patterns[pattern][taxon]]
                    except KeyError:
                        tmp += [0]
                        reflexes += [0]
                tmp += [score]
        
                D += [tmp]

                 
                try:
                    P[pattern[0]] += [(pattern[1:-1], pattern[-1], reflexes)]
                except KeyError:
                    P[pattern[0]] = [(pattern[1:-1], pattern[-1], reflexes)]
        
        output = [D[0]]
        for p in sorted(P):
            out = reduce_patterns(P[p], verbose=False)
            sorter = [x[0] for x in
                    sorted(zip(range(len(out['reflexes'])),out['reflexes']))]

            for i in sorter:
                tmp = [p]
                conditions = out['conditions'][i]
                reflex = out['reflexes'][i]
                totals = out['totals'][i]
                targets = out['targets'][i]
                
                # assemble the conditions
                idx = 0
                for j,tier in self.tier_system.items():
                    if j in out['tiers']:
                        tmp += [conditions[idx]]
                        idx += 1
                    else:
                        tmp += ['']

                tmp += [targets]
                # now assemble the reflexes
                for taxon in reflex:
                    tmp += [taxon]
                tmp += [totals]

                output += [tmp]



        # iterate over patterns and reduce them 
        with open(self.filename.replace('.tsv', '.patterns'), 'w') as f:
            for line in output:
                f.write('\t'.join([str(x) for x in line])+'\n')

        return P
                
    def check_cognates(self, output=False, filename=''):
        
        # make a wordlist object and indicate bad and good cognates here
        D = {0:['taxon', 'concept', 'ipa', 'tokens', 'alignment','cogid']}
        goods = 0
        allwords = 0
        for taxon in self.descendants:
            all_words = self.check_cognates_per_language(taxon)
            for key,gloss,pre,target,source in all_words:
                cogid = self[key,self.ref]
                if pre == '+':
                    goods += 1
                allwords += 1
                D[key] = [taxon, gloss, self[key,'ipa'], ' '.join(self[key,'tokens']), target, self[key,'cogid']]
        for k in self.get_list(taxon=self.proto, flat=True):
            D[k] = [self.proto, self[k,'concept'], self[k,'ipa'], ' '.join(self[k,'tokens']), ' '.join(self[k,'alignment']), self[k,'cogid']]
        
        print(goods, allwords, '{0:.2f}'.format(goods/allwords))
        if output:
            wl = lp.Wordlist(D)
            wl._meta={}
            wl.output('tsv', filename=filename or self.filename)
        else:
            return D

    def check_cognates_per_language(self, taxon, output=False, filename=''):
        """
        Check the correctness of proposed cognates based on regularity
        statements.
        """
        if not hasattr(self, 'regular_changes'):
            raise ValueError('[!] You need to check the tiers first!')
        good_words = 0
        bad_word_count = 1
        bad_words = []
        all_words = []
        for key,line in sorted(self.words[taxon].items()):
            alm = line[3]
            bads = 0
            word = []
            for i,(s,t) in enumerate(misc.transpose(alm)):
                # get the tier
                actives = self.explanations[taxon][s,t][4]
                if actives[0][0] == 'contextless':
                    word += [t]
                else:
                    indices = [self.tier_index[active]+1 for active in actives]
                    tier = tuple([(y,self.tiers[taxon][key][x][i]) for x,y in
                        zip(indices,actives)])
                    # check for gaps
                    if t == '-':
                        if len(word) > 0:
                            if word[-1] == '-':
                                word += ['*'+t+'*']
                            else:
                                if (s,t) not in self.regular_changes[taxon]:
                                    bads += 1
                                    word += ['!'+t]
                                else:
                                    if tier not in self.regular_changes[taxon][s,t]:
                                        bads += 1
                                        word += ['!'+t]
                                    else:
                                        word += [t]

                    elif (s,t) not in self.regular_changes[taxon]:
                        bads += 1
                        word += ['!'+t]
                    else:
                        if tier not in self.regular_changes[taxon][s,t]:
                            bads += 1
                            word += ['!'+t]
                        else:
                            word += [t]
            if bads:
                pre = '!'
                bad_word_count += 1
                bad_words += [line[1]]
            else:
                pre = '+'
                good_words += 1

            all_words += [[line[0], line[1], pre, ' '.join(word), ' '.join(line[3][0])]]

            print('{3} [{0}] "{1}": {2}, {4}, {5}'.format(line[1], self[line[0],
                    'ipa'], bads, pre, ' '.join(word), ' '.join(line[3][0])))

        if output:
            fn = filename +'cognates.txt' or self.filename.replace('.tsv','cognates.txt')

            with open(fn,'w') as f:
                f.write('Reflexes from language {0}\n'.format(taxon))
                for key, concept, pre, target, source in all_words:
                    f.write('[{0}] ({1})\t"{2}"\t[{3}] < [*{4}]\n'.format(
                        pre,
                        key,
                        concept,
                        target,
                        source
                        ))
                

        total = good_words + bad_word_count
        print("Out of {0} words, there are {1} ({2:.2f}%) good words and {3} ({4:.2f}%) bad words.".format(
            total,
            good_words,
            100 * good_words / total,
            bad_word_count,
            100 * bad_word_count / total
                    )
                )
        return all_words

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
        if len(set(matrix[-1])) == 1 or len(matrix) == 1:
            return {0 : (
                len(matrix[-1]), 0, [('contextless',0)],[i for i in
                self.instances[taxon][sound]],[],[])},[self.value_matrices[taxon][sound][0][-1],1,0]

        atiers = self.active_tiers[taxon][sound]
        instances = self.instances[taxon][sound]
        values = misc.transpose(self.value_matrices[taxon][sound])
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
        reflexes = [values[-1][a:b][0] for a,b in idxs]
        # our explanations distinguish regular and irregular characters
        explanations = dict([(i,(0, len(matrix[0]),[('exception',0)],[],[],[],
            len(matrix), len(matrix))) for i in range(len(idxs))])

        # check the minimization in this context: what do we actually want to
        # minimize? Is it the number of regular tokens, is it the number of
        # good items, is it the number of items that are false?
        # maybe it would just be enough to minimize the irregularities!
        # note another thing, namely, that we shouldn't flag all items as
        # irregular, but only those occuring not enough, so if in the end we
        # have irregularities in two or more reflexes, we need to assemble
        # these and re-flag the one that occurs the most as regular
        #for i in range(len(matrix)-2,len(matrix)-1): #max(1,len(matrix)-2), max(2, int(len(matrix)-1))): #-1
        #    for indices in itertools.combinations(range(len(matrix)-1), i):
        line = list(zip(*[matrix[jj] for jj in range(len(matrix)-1)]))
        indices = range(len(matrix)-1)
        for j,(a,b) in enumerate(idxs):
            # get the current part, the rest, and the overlap
            part = line[a:b]
            rest = line[:a]+line[b:]
            regular_tokens = [(p,nr) for p,nr in zip(part,range(a,b)) if p not in rest]
            regular = len([p[0] for p in regular_tokens])
            irregular_tokens = [(p,nr) for p,nr in zip(part, range(a,b))
                    if p in rest]
            irregular = len([p[0] for p in irregular_tokens])   
            irregular_types = len(set([p[0] for p in
                irregular_tokens]))
            regular_types = len(set([p[0] for p in regular_tokens]))
            if verbose:
                pass #print(regular, irregular, reflexes[j])
            if regular > explanations[j][0]:
                explanations[j] = (regular, irregular, [atiers[jj] for
                    jj in indices],
                    [self.instances[taxon][sound][n[1]] for n in
                        regular_tokens],
                    [self.instances[taxon][sound][n[1]] for n in
                        irregular_tokens],
                    indices, regular_types, irregular_types)
            elif regular == explanations[j][0]:
                if regular_types < explanations[j][6]:
                    #if irregular < explanations[j][1]:
                    explanations[j] = (regular, irregular, [atiers[jj]
                        for jj in indices], 
                        [self.instances[taxon][sound][n[1]] for n in 
                            regular_tokens],
                        [self.instances[taxon][sound][n[1]] for n in
                            irregular_tokens],
                        range(len(matrix)-1), regular_types, irregular_types)
        ## we now clean the irregularities
        for j,(a,b) in enumerate(idxs):
            line = list(zip(*[matrix[x] for x in
                explanations[j][5]]))
            part = line[a:b]
            rest = line[:a] + line[b:]
            irregular_tokens = [(p,nr) for p,nr in zip(part,range(a,b)) if
                    p in rest]
            irregular_rest = len([p for p in rest if p in part])
            if len(irregular_tokens) > irregular_rest:
                print("[i] correcting irregular rest, {0} vs. {1}".format(
                    len(irregular_tokens), irregular_rest), taxon, sound, reflexes[j])
                input()
                explanations[j] = (
                        explanations[j][0],
                        explanations[j][1],
                        explanations[j][2],
                        explanations[j][3] + explanations[j][4],
                        [],
                        explanations[j][5],
                        explanations[j][6],
                        explanations[j][7]
                        )
        
        if verbose:
            tidx = self.taxa.index(taxon)
            for k in explanations:
                print(reflexes[k])
                print('---')
                print(explanations[k][0], explanations[k][1])
                #print(', '.join(['.'.join(self.msa[self.ref][x]['seqs'][self.msa[self.ref][x]['taxa'].index(taxon)]) for x in explanations[k][5]]))
        return explanations, reflexes

    def cluster_tiers(self, language, threshold=0.5, cluster_method='upgma',
            verbose=False, ignore_vowels=True, mark_singletons=True):
        
        regular_tiers = set()
        irregular_tiers = {}
        sounds = sorted(self.instances[language])
        if ignore_vowels:
            sounds = [s for s in sounds if s[0] not in lp.rc('vowels')]

        # iterate over all sounds for the language
        for sound in sounds:

            if verbose:
                print("[i] Analyzing sound {0}...".format(sound))

            tier = self.value_matrices[language][sound]
            instances = self.instances[language][sound]

            clusters, exceptions = cluster_tier(tier, threshold=threshold,
                    cluster_method=cluster_method,
                    mark_singletons=mark_singletons)

            # take all clusters which are not in exceptions and add a scorer
            # value
            for cluster in clusters:
                for idx in clusters[cluster]:
                    instance = tuple([sound] + tier[idx])
                    if instance[1:] not in exceptions:
                        regular_tiers.add(instance)
                    else:
                        irregular_tiers[instance] = exceptions[instance[1:]]

        return regular_tiers, irregular_tiers

    def check_words(self, language, threshold=0.5, cluster_method='upgma',
            verbose=True, ignore_vowels=True, mark_singletons=True,
            mark_shakyness=True):

        regular_tiers, irregular_tiers = self.cluster_tiers(language, threshold=threshold,
                cluster_method=cluster_method, verbose=verbose,
                ignore_vowels=ignore_vowels, mark_singletons=mark_singletons)

        # make a scoring dictionary to handle the output
        scorer = {}
        for tier in regular_tiers:
            scorer[tier[:-1]] = tier[-1]
            
        goods, total = 0, 0

        # evaluate function for tiers
        def evaluate_word(alignment):
            """
            Function evaluates how well a word is in concordance with a given
            system of alignments.
            """
            alms = misc.transpose(alignment)
            opt, ipt,= [], []
            gap_track = 0
            irregularity = 0
            singletons = 0
            shakyness = 0
            for i,alm in enumerate(alms):
                # we allow consecutive gaps to be excluded to allow for
                # morphological change patterns to be assumed
                if alm[-1] == '-':
                    gap_track += 1
                else:
                    gap_track = 0

                if tuple(alm) in regular_tiers:
                    opt += [alm[-1]]
                elif ignore_vowels and alm[0][0] in lp.rc('vowels'): # or alm[-1][0] in lp.rc('vowels'):
                    opt += [alm[-1]]
                elif gap_track and alm[0] != '-' and alm[-1] != '-':
                    opt += [''+alm[-1]]
                else:
                    if tuple(alm[:-1]) in scorer:
                        exp = ''+scorer[tuple(alm[:-1])]
                        irregularity += 1
                    else:
                        if tuple(alm) in irregular_tiers:
                            exp = '/'.join([a for a in
                                irregular_tiers[tuple(alm)] if a[1:] != alm[-1]])
                            shakyness += 1
                        else:
                            exp = '?'
                            singletons += 1
                    opt += ['!'+alm[-1]+'({0})'.format(exp)]
                ipt += [alm[0]]
            return ipt, opt, irregularity, shakyness, singletons

        # now we iterate over all words
        output = []
        for cogid, tier in self.tiers[language].items():
            ipt, opt, irregularity, shakyness, singletons = evaluate_word(tier)
            if mark_singletons:
                irregularity += singletons
            if mark_shakyness:
                irregularity += shakyness
                
            if irregularity == 0:
                goods += 1
            total += 1
            if verbose:
                if irregularity > 0:
                    if verbose: print(''.join(ipt), ' --> ', ''.join(opt),
                            '({0})'.format(irregularity))
                    output += [(cogid, ipt, opt, 1)]
                else:
                    output += [(cogid, ipt, opt, 0)]

        if verbose: print(goods, total, '{0:.2f}'.format(goods / total))
        
        if not hasattr(self, 'regular_tiers'):
            self.regular_tiers = {language : regular_tiers}
            self.irregular_tiers = {language: irregular_tiers}
        else:
            self.regular_tiers[language] = regular_tiers
            self.irregular_tiers[language] = irregular_tiers

        return goods, total, output
    
    def check_all_words(self, threshold=0.4, cluster_method='upgma',
            verbose=False, ignore_vowels=True, output=False,
            mark_singletons=True, mark_shakyness=True):
        """
        Check all words whether they have been correctly inferred.
        """
        all_goods, all_total = 0, 0
        out = {}
        out2 = {}
        for taxon in self.descendants:
            goods, total, words = self.check_words(taxon, threshold=threshold,
                    cluster_method=cluster_method, verbose=verbose,
                    ignore_vowels=ignore_vowels,
                    mark_singletons=mark_singletons,
                    mark_shakyness=mark_shakyness)
            if verbose:
                print('[i] intermediate scores for {0} are {1}, {2}, {3:.2f}...'.format(
                    taxon,
                    goods,
                    total,
                    goods / total))
            all_goods += goods
            all_total += total

            for cogid, ipt, opt, status in words:
                if status:
                    try:
                        out[cogid][taxon] = (ipt,opt)
                    except KeyError:
                        out[cogid] = {taxon : (ipt,opt)}
                else:
                    try:
                        out2[cogid][taxon] = (ipt,opt)
                    except KeyError:
                        out2[cogid] = {taxon : (ipt,opt)}
                
        if verbose: print('[i] summarized scores are {0}, {1}, {2:.2f}'.format(
            all_goods,
            all_total,
            all_goods / all_total))

        if output:
            scores = dict([(t,0) for t in self.descendants])
            with open(self.filename.replace('.tsv', '.errors'), 'w') as f:
                for cogid in sorted(out):
                    # get the words in the proto-language
                    pword = self.words[self.proto][cogid][2]
                    pmean = self.words[self.proto][cogid][1]

                    f.write('[{2}] Problems with item {0} «{1}»:\n'.format(
                        ' '.join(pword), pmean, cogid))
                    for t in self.descendants:
                        if t in out[cogid]:
                            f.write('{0}\t{1}\n'.format(
                                t, ' '.join(out[cogid][t][1])))
                            scores[t] += 1
                    f.write('\n')
                f.write('\nSummary for all taxa\nTaxon\tPotential Errors\tNumber of Cognates\tProportion\n')
                for t in self.descendants:
                    f.write('{0}\t{1}\t{2}\t{3:.2f}\n'.format(
                        t, scores[t],
                        len(self.words[t]),
                        scores[t] / len(self.words[t])
                        ))
            
            goodlist = []
            with open(self.filename.replace('.tsv', '.regulars'), 'w') as f:
                for cogid in sorted(out2):

                    # get the words in the proto-language
                    pword = self.words[self.proto][cogid][2]
                    pmean = self.words[self.proto][cogid][1]

                    goodlist += [self.words[self.proto][cogid][0]]

                    f.write('[{2}] Regular change for item {0} «{1}»:\n'.format(
                        ' '.join(pword), pmean, cogid))
                    for t in self.descendants:
                        if t in out2[cogid]:
                            f.write('{0}\t{1}\n'.format(
                                t, ' '.join(out2[cogid][t][1])))
                            scores[t] += 1
                            goodlist += [self.words[t][cogid][0]]
                    f.write('\n')
                f.write('\nSummary for all taxa\nTaxon\tPotential Errors\tNumber of Cognates\tProportion\n')
                for t in self.descendants:
                    f.write('{0}\t{1}\t{2}\t{3:.2f}\n'.format(
                        t, scores[t],
                        len(self.words[t]),
                        scores[t] / len(self.words[t])
                        ))
                
                if all_goods != all_total: # < 1.0:
                    print("[i] Writing exceptionless testset to file.")
                    self.output('tsv',
                            filename=self.filename.replace('.tsv','.goodies'), 
                            subset=True,
                            rows=dict(ID = 'in '+str(goodlist)),
                            )
            
        if verbose:
            print('[i] out of {0} cognate sets, {1} are irregular ({2:.2f}%).'.format(
                len(self.words[self.proto]),
                len(out),
                len(out) / len(self.words[self.proto])))

    def get_patterns(self, language, sounds = False, ignore_vowels=True,
            max_exceptions=1, verbose=True):

        # get the sounds for the language first
        if not sounds:
            sounds = sorted([s for s in self.instances[language]])
            if ignore_vowels:
                sounds = [s for s in sounds if s[0] not in lp.rc('vowels')+'-']

        patterns = []

        # now get the tiers
        for sound in sounds:
            tier = self.value_matrices[language][sound]
            instances = self.instances[language][sound]
            
            # bring order in the stuff according to the reflexes
            reflexes = {}
            for i,line in enumerate(tier):
                try:
                    reflexes[line[-1]] += [i]
                except KeyError:
                    reflexes[line[-1]] = [i]

            for reflex in reflexes:
                
                if verbose: print('---')
                
                # find optimal solution by using combinations
                all_indices = list(range(len(self.tier_system)))
                wrong = True
                best_index_set = []
                best_reflex_tiers = []
                for idx in range(1, len(all_indices)+1):
                    if not wrong:
                        break
                    for index_set in itertools.combinations(all_indices, idx):
                        if not wrong:
                            break
                        # get the reflex tiers
                        reflex_tiers = [tuple(
                                        [tier[i][j] for j in index_set]
                                        ) for i in reflexes[reflex]]


                        # get the rest of the other tiers
                        other_tiers = []
                        for k in [x for x in reflexes if x != reflex]:
                            other_tiers += [tuple(
                                [tier[i][j] for j in index_set]
                                ) for i in reflexes[k]]
                        other_tiers = sorted(set(other_tiers))
                        
                        # check for exceptions
                        exceptions = 0
                        for s in reflex_tiers:
                            if s in other_tiers:
                                exceptions += 1

                        if not exceptions or exceptions <= max_exceptions:
                            wrong = False
                            best_index_set = index_set
                            best_reflex_tiers = reflex_tiers
                        
                for i,s in zip(reflexes[reflex], reflex_tiers):
                    exp = '! ' if s in other_tiers else '  '
                    if exp.strip(): 
                        if verbose: print(i, exp+'  '.join(['{0:3}'.format(t) for t in s]),
                                exceptions, len(reflex_tiers),  
                                instances[i][0], self.words[language][instances[i][0]][1], ' '.join(self.words[language][instances[i][0]][2]))

                if exceptions:
                    if verbose: input('---pause---')

                visited = []
                for i,line in enumerate(best_reflex_tiers):
                    
                    #if line not in visited:
                    conditions = sorted(zip(['.'.join([str(y) for y in a]) for a in [self.tier_system[x] for x in
                        best_index_set]], line), 
                        key=lambda x: tuple(x[0].split('.')[::-1]))
                    if line not in visited:
                        if verbose: print(sound, '->', reflex, '/', '; '.join(['{0} : {1:3}'.format(
                            x,y) for x,y in conditions]),
                                )
                    visited += [line]
                    
                    cidx,widx = instances[reflexes[reflex][i]]
                        
                    patterns += [
                            [language, sound, reflex, dict(conditions),
                                self.tiers[language][cidx][0], 
                                self.tiers[language][cidx][-1], 
                                widx,
                                cidx
                                ]]


                if exp.strip(): 
                    if verbose: 
                        input('there were exceptions')
        return patterns

    def get_all_patterns(self, ignore_vowels=True, max_exceptions=0, output=True,
            filename=None):

        filename = filename or self.filename.replace('.tsv', '.html')

        P = []
        tiers = []
        sounds = []
        for taxon in self.descendants:
            print("Analyzing patterns for taxon {0}.".format(taxon))
            patterns = self.get_patterns(taxon, ignore_vowels=ignore_vowels,
                    max_exceptions=max_exceptions, verbose=False)  
            for pattern in patterns:
                tiers += list(pattern[3].keys())
                sounds += [pattern[1]]
            P += [p for p in patterns]
        
        template = open('tiers.html').read()
        js = open('tiers.js').read()

        with open(filename, 'w') as f:
            patterns = 'var PATTERNS = '+json.dumps(P, indent=2) + ';\n'
            
            f.write(template.format(
                DATASET=self.filename.replace('.tsv','').upper(),
                PATTERNS = patterns,
                TIERS = 'var TIERS = '+json.dumps(
                        sorted(set(tiers), key=lambda x:
                            tuple(x.split('.')[::-1])),
                            indent=2)+';\n',
                JS = js,
                SOUNDS = 'var SOUNDS = '+json.dumps(
                    sorted(set(sounds), key=lambda x:
                        (lp.tokens2class([x], 'dolgo'), x)), indent=2) +';\n'
                ))

            


def cluster_tier(value_matrix, threshold=0.5, cluster_method='upgma',
        verbose=False, mark_singletons=True):
    
    # reduce value matrix to only consider the types, not the tokens
    converter = OrderedDict()
    for i,line in enumerate(value_matrix):
        try:
            converter[tuple(line)] += [i]
        except KeyError:
            converter[tuple(line)] = [i]
    new_matrix = [k for k in converter]

    # cluster contexts 
    matrix = []
    D = {}
    for i,lineA in enumerate(new_matrix):
        try:
            D[tuple(lineA[:-1])][lineA[-1]] += converter[lineA]
        except KeyError:
            try:
                D[tuple(lineA[:-1])][lineA[-1]] = converter[lineA]
            except KeyError:
                D[tuple(lineA[:-1])] = {lineA[-1] : converter[lineA]}
        for j,lineB in enumerate(new_matrix):
            if i < j:
                d = []
                visited = []
                for a,b in zip(lineA,lineB):
                    if (a,b) not in visited:
                        if a == b:
                            d += [0]
                        else:
                            d += [1]
                    visited += [(a,b)]
                matrix += [sum(d)/len(d)]
    
    # analyse the data
    matrix = misc.squareform(matrix)
    clrs = lp.flat_cluster(cluster_method, threshold, matrix,
            list(range(len(new_matrix))))
    # reconvert the clusters according to our initial matrix
    clusters = {}
    for clr in clrs:
        clusters[clr] = []
        for idx in clrs[clr]:
            clusters[clr] += converter[new_matrix[idx]]

    # count all reflexes
    all_reflexes = {}
    for line in value_matrix:
        try:
            all_reflexes[line[-1]] += 1
        except KeyError:
            all_reflexes[line[-1]] = 1

    # identify exceptions in the data 
    singletons = {}
    ambiguous_clusters = {}
    ambiguous_contexts = {}
    reflexes = []
    for clr in sorted(clrs, key=lambda x: len(clrs[x])):
        
        indices = []
        for i in clrs[clr]:
            indices += converter[new_matrix[i]]
            reflexes += [new_matrix[i][-1]]
        
        # get basic chars
        _chars = [value_matrix[i][-1] for i in indices]
        chars = sorted(_chars, key=lambda x: _chars.count(x), reverse=True)
        char_count = dict([(char, _chars.count(char)) for char in chars])
        bchar = chars[0]
        max_val = max(char_count.values())
        bchars = [char for char in char_count if char_count[char] == max_val]

        # search for exceptions in the matrix, we define three types of
        # exceptions:
        # 1. singleton contexts which are hard to evaluate (only if they yield
        #    unexpected output, otherwise we'll 
        # 2. ambiguous contexts which lead to two or more reflexes
        # 3. the least frequent contexts in ambiguous clusters (here, we need
        #    to think of how to handle clusters with only two instances)
        
        # check for singletons
        if len(indices) == 1:
            line = value_matrix[indices[0]]
            singletons[tuple(line)] = line[-1] 
        else:
            # determine the major reflex
            main_reflex = [value_matrix[i][-1] for i in indices]
            max_ref = max([main_reflex.count(a) for a in main_reflex])
            main_reflex = [s for s in set(main_reflex) if main_reflex.count(s)
                    == max_ref]
            
            # check for major reflex and update if there's more than one
            if len(main_reflex) > 1:
                main_ref = max([all_reflexes[x] for x in main_reflex])
                main_reflex = [s for s in main_reflex if all_reflexes[s] ==
                        main_ref]
            #, key=lambda x:
            #        main_reflex.count(x))[0]
            
            for i in indices:
                line = value_matrix[i]
                key, char = tuple(line[:-1]), line[-1]

                # check for ambigous contexts
                if len(D[key]) > 1:
                    # count how often the ambiguous contexts occur
                    key_count = dict([(c, len(D[key][c])) for c in
                        D[key]])
                    max_val = max(key_count.values())
                    tkeys = [k for k in key_count if key_count[k] == max_val]
                    if char not in tkeys or len(tkeys) > 1:
                        if char not in main_reflex:
                            ambiguous_contexts[tuple(line)] = tkeys 
                
                # check for ambiguous clusters
                if line[-1] not in bchars or len(bchars) > 1:
                    if char not in main_reflex:
                        ambiguous_clusters[tuple(line)] = bchars

        # search also for exceptions by taking exceptional chars
    
    # check exceptions
    exceptions = {}

    # make count of most frequent reflexes
    if mark_singletons:
        reflexes = dict([(a,reflexes.count(a)) for a in reflexes])
        for singleton, reflex in singletons.items():
            if reflexes[reflex] == max(reflexes.values()):
                pass
            else:
                exceptions[singleton] = ['?']

    # all ambiguous contexts are exceptions per se
    for ctx, chars in ambiguous_contexts.items():
        exceptions[ctx] = ['!'+c for c in chars]

    # ambiguous clusters can be further checked 
    for clr, chars in ambiguous_clusters.items():
        if clr in exceptions:
            pass #exceptions[clr] += ['*'+c for c in chars]
        else:
            exceptions[clr] = ['*'+c for c in chars]

    return clusters, exceptions

if __name__ == '__main__':
    
    from sys import argv

    if len(argv) > 1 and argv[1] == 'tukano':
        infile = 'tukano2.tsv'
        proto = '*PT'
        lang = 'KAR'
        thr = 0.4
    else:
        infile = 'germanic.goodies.goodies.tsv'
        proto = 'Proto-Germanic'

    tiers = Tiers(infile, proto=proto)
    tiers.make_tiers(exclude_gaps=False)
    
    if proto == '*PT':
        words = []
        pwords = []
        sound = 'j'
        lang = 'BAS'
        thr = 0.25
        for tier,instance in zip(
                tiers.value_matrices[lang][sound],
                tiers.instances[lang][sound]
                ):

            words += [tiers.words[lang][instance[0]][3][1]]
            pwords += [tiers.words[lang][instance[0]][3][0]]

        tier = tiers.value_matrices[lang][sound]
        d,E = cluster_tier(tier, verbose=True, threshold=thr)
        count = 0
        for clr,char in enumerate(d):
            print('cluster {0}'.format(clr+1))
            for i in d[char]:
                count += 1
                if tuple(tier[i]) in E:
                    exc = '[!]:'+', '.join(E[tuple(tier[i])])
                else:
                    exc = ''
                print(' '.join(['{0:3}'.format(t) for t in tier[i]]),
                        ''.join(words[i]), ''.join(pwords[i]), exc)
        print(len(E),count)
        input()

    if proto == 'Proto-Germanic':
        words = []
        pwords = []
        sound = 's'
        lang = 'German'
        thr = 0.3
        #for tier,instance in zip(
        #        tiers.value_matrices[lang][sound],
        #        tiers.instances[lang][sound]
        #        ):

        #    words += [tiers.words[lang][instance[0]][3][1]]
        #    pwords += [tiers.words[lang][instance[0]][3][0]]

        #tier = tiers.value_matrices[lang][sound]
        #d,E = cluster_tier(tier, verbose=True, threshold=thr, 
        #        cluster_method='complete')
        #count = 0
        #for clr,char in enumerate(d):
        #    print('cluster {0}'.format(clr+1))
        #    for i in d[char]:
        #        count += 1
        #        if tuple(tier[i]) in E:
        #            exc = '[!]:'+', '.join(E[tuple(tier[i])])
        #        else:
        #            exc = ''
        #        print(' '.join(['{0:3}'.format(t) for t in tier[i]]),
        #                ''.join(words[i]), ''.join(pwords[i]), exc)
        #print(len(E),count)
        #input()


        #input()

    #rt = tiers.cluster_tiers('German',verbose=True)
    #tiers.check_words(lang, verbose=True, threshold=thr)
    if 'clean' in argv:
        tiers.check_all_words(threshold=thr, output=True, verbose=True,
                cluster_method='complete')
    elif 'check' in argv:
        if 'sound' in argv:
            sound = argv[argv.index('sound')+1]
        else:
            sound = 't'
        print("German")
        a = tiers.get_patterns('German', sounds=[sound], max_exceptions=0)
        print("\nEnglish")
        tiers.get_patterns('English', sounds=[sound], max_exceptions=0)
        print("\nDutch")
        tiers.get_patterns('Dutch', sounds=[sound], max_exceptions=0)

    elif 'print' in argv:
        tiers.get_all_patterns()

