import sys
import os
from subprocess import call

task_id = sys.argv[1]

# project root path
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, '..')))
for sample in os.listdir('/immunolyser-data/data/{}'.format(task_id)):
    for replicate in os.listdir('/immunolyser-data/data/{}/{}'.format(task_id,sample)):
       
        if replicate[-8:] == '9mer.txt':
            print(os.popen('perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl \
                    -f /immunolyser-data/data/{}/{}/{} \
                    -H R \
                    -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
                    -g 1-5 \
                    -R /var/www/firstdemo/app/static/images/{}/{}/gibbscluster/{} '.format(task_id, sample, replicate, task_id, sample, replicate[:-9])).read())
