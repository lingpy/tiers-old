from lingpy import *
# re-align germanic.tsv

wl = Wordlist('germanic.tsv')
blacklist = []
for k in wl:
    if wl[k,'taxon'] == 'Proto-Germanic':
        ipa = wl[k,'ipa']
        if ipa[-2:] in ['an','az','in', 'iz','uz']: #.endswith('an') or wl[k,'ipa'].endswith('az'):
            wl[k][wl.header['tokens']] = wl[k,'tokens'][:-2]
        elif ipa.endswith('janan'):
            wl[k][wl.header['tokens']] = wl[k,'tokens'][:-5]
        elif ipa.endswith('oːn'):
            wl[k][wl.header['tokens']] = wl[k,'tokens'][:-2]
    if 
    if wl[k,'taxon'] == 'English':
        tokens=wl[k,'tokens']
        tokens = ' '.join(tokens)
        tokens = tokens.replace('ḷ', 'ə l')
        tokens = tokens.replace('ṇ', 'ə n')
        wl[k][wl.header['tokens']] = tokens.split(' ')

    if k in [1125, 1199]:
        wl[k][wl.header['tokens']] = ' '.join(wl[k,'tokens']).split('/')[0][:-1].split(' ')

        
        if '+' in wl[k,'ipa'] or '/' in wl[k,'ipa'] or '_' in wl[k,'ipa']:
            blacklist += [wl[k,'concept']]
wl.output('tsv', filename='germanic.2', subset=True, rows=dict(concept = 'not in '+str(blacklist)), cols=[h for h in wl.header if h != 'alignment'])

alm = Alignments('germanic.2.tsv')

alm.align(method='library', factor=0.1,iteration=True, swap_check=True)
# make blacklist from swapped sequences
blacklist = [] 
for k,v in alm.msa['cogid'].items():
    if 'swaps' in v:
        for idx in v['ID']:
            blacklist += [idx]
alm.output('tsv', ignore='all', subset=True, filename='pgm.aligned',
        rows=dict(ID = 'not in '+str(blacklist)))

