import sys
import os
from subprocess import call

experiment_id_date = sys.argv[1]

for sample in os.listdir('data\\{}'.format(experiment_id_date)):
    for replicate in os.listdir('data\\{}\\{}'.format(experiment_id_date,sample)):
        
        try:
            os.mkdir('data\\{}\\{}\\seqlogos'.format(experiment_id_date,sample))
        except OSError:
            print ("Creation of the directory %s failed" % 'data\\{}\\{}\\seqlogos'.format(experiment_id_date,sample))
        else:
            print ("Successfully created the directory %s " % 'data\\{}\\{}\\seqlogos'.format(experiment_id_date,sample))
        if replicate[-3:] == 'txt':
            print(replicate)
            call('app\\tools\\Python2\\python.exe app\\tools\\seq2logo-2.1\\Seq2Logo.py -f data\\{}\\{}\\{} --format [JPEG] -o {}\\{}\\{}\\seqlogos\\{}'.format(experiment_id_date, sample, replicate, 'app\\static\\images',experiment_id_date,sample,replicate[:-4]))
            # call('tools\\Python2\\python.exe -m pip install numpy')

            