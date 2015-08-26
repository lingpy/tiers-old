from lingpy import *
# re-align germanic.tsv
alm = Alignments('germanic.tsv')
alm.align(method='library', iteration=True, swap_check=True)
# make blacklist from swapped sequences
blacklist = [] 
for k,v in alm.msa['cogid'].items():
    if 'swaps' in v:
        for idx in v['ID']:
            blacklist += [idx]
alm.output('tsv', ignore='all', subset=True, filename='pgm.aligned',
        rows=dict(ID = 'not in '+str(blacklist)))

