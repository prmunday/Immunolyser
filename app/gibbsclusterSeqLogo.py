import sys
import os
from subprocess import call
# from constants import *

task_id = sys.argv[1]
data_mount = sys.argv[2]

sample = sys.argv[3]
replicate = sys.argv[4]
replicate += '_9mer.txt'
cluster = sys.argv[5]
mhc_class = sys.argv[6]
motif_length = int(sys.argv[7])

# Set -C (perform Single sequence moves at every iteration (default is every -r iterations)) based on mhc class
arg_C = '-C' if mhc_class == "One" else ''

# Argument to specify the motif length
arg_l = f'-l {motif_length}'

print(os.popen('perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA_for_seqlogo.pl \
        -f {}/{}/{}/{} \
        -H R \
        -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
        -g {} \
        -k {} \
        -T \
        {} \
        {} \
        -R {}/app/static/images/{}/{}/gibbscluster/{} '.format(data_mount,task_id, sample, replicate, cluster, os.cpu_count(), arg_C, arg_l, os.getcwd(), task_id, sample, replicate[:-9])).read())