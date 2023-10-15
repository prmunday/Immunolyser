import sys
import os
from subprocess import call

inp = sys.argv[1]
out = sys.argv[2]
name = sys.argv[3]
nine_mers = sys.argv[4]
total_peptides = sys.argv[5]

call(['python2', 'app/tools/seq2logo-2.1/Seq2Logo.py', '-f', inp, '-o', out, '--format', '[JPEG]', '-t', '{} based on {} 9-mers'.format(name,nine_mers,total_peptides), '-S', '2', '-I', '2'])
