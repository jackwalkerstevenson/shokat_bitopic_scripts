from modeller import *

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='abl1_pl1_pep', model_segment=('FIRST','LAST'))
aln.append_model(mdl, align_codes='abl1_pl1_pep', atom_files='abl1_pl1_pep.pdb')
aln.append(file='abl+pep.ali', align_codes='abl+pep')
aln.align2d(max_gap_length=50)
aln.write(file='pep_model.ali', alignment_format='PIR')
aln.write(file='pep_model.pap', alignment_format='PAP')