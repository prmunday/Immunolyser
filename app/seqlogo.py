import sys
import os
from subprocess import call

task_id = sys.argv[1]

for sample in os.listdir('data/{}'.format(task_id)):
    for replicate in os.listdir('data/{}/{}'.format(task_id,sample)):
        
        if replicate[-8:] == '9mer.txt':
            call(['python2', 'app/tools/seq2logo-2.1/Seq2Logo.py', '-f', 'data/{}/{}/{}'.format(task_id, sample, replicate), '-o', '{}/{}/{}/seqlogos/{}'.format('app/static/images', task_id, sample, replicate[:-9]), '--format', '[JPEG]', '-t', '{}'.format(replicate[:-9]), '-S', '2', '-I', '2'])
