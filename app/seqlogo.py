import sys
import os
from subprocess import call

task_id = sys.argv[1]
data_mount = sys.argv[2]


for sample in os.listdir('{}/{}'.format(data_mount,task_id)):
    for replicate in os.listdir('{}/{}/{}'.format(data_mount,task_id,sample)):
        
        if replicate[-22:] == '9merwithduplicates.txt':
            call(['python2', 'app/tools/seq2logo-2.1/Seq2Logo.py', '-f', '{}/{}/{}/{}'.format(data_mount,task_id, sample, replicate), '-o', '{}/{}/{}/seqlogos/{}'.format('app/static/images', task_id, sample, replicate[:-9]), '--format', '[JPEG]', '-t', '{}'.format(replicate[:-9]), '-S', '2', '-I', '2'])
