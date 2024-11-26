import sys
import os
from subprocess import call
# from constants import *

task_id = sys.argv[1]
data_mount = sys.argv[2]
mhc_class = sys.argv[3]

# Set -C (perform Single sequence moves at every iteration (default is every -r iterations)) based on mhc class
arg_C = '-C' if mhc_class == "One" else ''


for sample in os.listdir('{}/{}'.format(data_mount,task_id)):
    for replicate in os.listdir('{}/{}/{}'.format(data_mount,task_id,sample)):

        if replicate[-9:] == '_9mer.txt':
            print('Replicate file: ', replicate)
            print("Gibbs Command: ", 'perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl \
                    -f {}/{}/{}/{} \
                    -H R \
                    -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
                    -g 1-5 \
                    -k {} \
                    -T \
                    {} \
                    -R {}/app/static/images/{}/{}/gibbscluster/{} '.format(data_mount,task_id, sample, replicate, os.cpu_count(), arg_C, os.getcwd(), task_id, sample, replicate[:-9]))

            print(os.popen('perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl \
                    -f {}/{}/{}/{} \
                    -H R \
                    -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
                    -g 1-5 \
                    -k {} \
                    -T \
                    {} \
                    -R {}/app/static/images/{}/{}/gibbscluster/{} '.format(data_mount,task_id, sample, replicate, os.cpu_count(), arg_C, os.getcwd(), task_id, sample, replicate[:-9])).read())
