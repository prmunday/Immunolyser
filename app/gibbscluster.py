import sys
import os
from subprocess import call

# Get command-line arguments
task_id = sys.argv[1]
data_mount = sys.argv[2]
mhc_class = sys.argv[3]
motif_length = int(sys.argv[4])

# Determine the project root
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# Set -C argument based on MHC class
arg_C = '-C' if mhc_class == "One" else ''

# Argument to specify the motif length
arg_l = f'-l {motif_length}'

# Set seq2logo path
seq2logo_path = os.path.join(project_root, "app/tools/seq2logo-2.1/Seq2Logo.py")

task_path = os.path.join(data_mount, task_id)
for sample in os.listdir(task_path):
    sample_path = os.path.join(task_path, sample)
    
    for replicate in os.listdir(sample_path):
        if replicate.endswith(f'_{motif_length}mer.txt'):
            print('Replicate file:', replicate)
            gibbs_command = (
                f'perl {project_root}/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl '
                f'-f {sample_path}/{replicate} '
                f'-H R '
                f'-G {seq2logo_path} '
                f'-g 1-5 '
                f'-k {os.cpu_count()} '
                f'-T {arg_C} {arg_l} '
                f'-R {project_root}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-9]}'
            )
            
            print("Gibbs Command:", gibbs_command)
            print(os.popen(gibbs_command).read())
