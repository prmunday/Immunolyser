import sys
import os
from subprocess import call

task_id = sys.argv[1]

for sample in os.listdir('data\\{}'.format(task_id)):
    for replicate in os.listdir('data\\{}\\{}'.format(task_id,sample)):
        
        if replicate[-3:] == 'txt':
            print(replicate)
            call('app\\tools\\Python2\\python.exe app\\tools\\seq2logo-2.1\\Seq2Logo.py \
                    -f data\\{}\\{}\\{} \
                    --format [JPEG] \
                    -o {}\\{}\\{}\\seqlogos\\{} \
                    -t {} \
                    -S 2 \
                    -I 2 '.format(task_id, sample, replicate, 'app\\static\\images',task_id,sample,replicate[:-4], replicate[:-4]))