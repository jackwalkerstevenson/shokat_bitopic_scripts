from modeller import *
from modeller.automodel import *

env = Environ()
a = AutoModel(env, alnfile='TvLDH-mult.ali',
              knowns=('7mswA','nsp2_137-427+rap1_m1A:B'), sequence='nsp2+rap1')
a.starting_model = 1
a.ending_model = 5
a.make()