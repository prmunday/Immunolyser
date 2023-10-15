import sys
import os
from subprocess import call

task_id = sys.argv[1]
data_mount = sys.argv[2]

sample = sys.argv[3]
replicate = sys.argv[4]
replicate += '_9mer.txt'
cluster = sys.argv[5]

print(os.popen('perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA_for_seqlogo.pl \
        -f {}/{}/{}/{} \
        -H R \
        -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
        -g {} \
        -k 4 \
        -R {}/app/static/images/{}/{}/gibbscluster/{} '.format(data_mount,task_id, sample, replicate, cluster, os.getcwd(), task_id, sample, replicate[:-9])).read())