import sys
import os
from subprocess import call

task_id = sys.argv[1]

# project root path
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, '..')))
for sample in os.listdir('data\\{}'.format(task_id)):
    for replicate in os.listdir('data\\{}\\{}'.format(task_id,sample)):
       
        if replicate[-3:] == 'txt':
            print(replicate)
            # call('app\\tools\\Python2\\python.exe app\\tools\\seq2logo-2.1\\Seq2Logo.py -f data\\{}\\{}\\{} --format [JPEG] -o {}\\{}\\{}\\seqlogos\\{} -t {} -S 2 -I 2 '.format(task_id, sample, replicate, 'app\\static\\images',task_id,sample,replicate[:-4], replicate[:-4]))
            print(os.popen('powershell \
                    app\\tools\\wsl\\Ubuntu\\ubuntu2004.exe \
                    run tcsh ./app/tools/gibbscluster-2.0/gibbscluster \
                    -f ./data/{}/{}/{} \
                    -H ./app/tools/R/bin/R \
                    -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
                    -g 1-5 \
                    -R /mnt/c/Users/ASUS/Documents/sem3/thesis/firstdemo/app/static/images/{}/{}/gibbscluster/{} '.format(task_id, sample, replicate, task_id, sample, replicate[:-4])).read())