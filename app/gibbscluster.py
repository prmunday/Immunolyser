import sys
import os

task_id = sys.argv[1]
data_mount = sys.argv[2]
mhc_class = sys.argv[3]
motif_length = int(sys.argv[4])

# Set arguments based on MHC class
if mhc_class == "One":
    arg_C = '-C'  # Perform Single sequence moves at every iteration
    arg_s = ''  # No shift moves
    arg_r = ''  # No single peptide move interval
else:
    arg_C = ''  # No -C flag
    arg_s = '-s 100'  # Interval between Phase Shift moves
    arg_r = '-r 20'  # Interval between Single Peptide moves

# Common arguments
arg_l = f'-l {motif_length}'
arg_I = '-I 0'  # Max insertion length
arg_D = '-D 0'  # Max deletion length
arg_S = '-S 5'  # Number of initial seeds
arg_b = '-b 0.8'  # Penalty lambda
arg_q = '-q 5'  # Weight on small clusters
arg_c = '-c 0'  # Sequence weighting type
arg_z = '-z 1'  # Background model (Uniprot pre-calculated)
arg_j = '-j 2'  # Threshold for trash cluster

for sample in os.listdir(f'{data_mount}/{task_id}'):
    for replicate in os.listdir(f'{data_mount}/{task_id}/{sample}'):
        if replicate.endswith(f'_{motif_length}mer.txt'):
            print('Replicate file: ', replicate)
            gibbs_command = f'''
                perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl \
                -f {data_mount}/{task_id}/{sample}/{replicate} \
                -H R \
                -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
                -g 1-5 \
                -k {os.cpu_count()} \
                -T \
                {arg_C} \
                {arg_l} \
                {arg_I} \
                {arg_D} \
                {arg_S} \
                {arg_b} \
                {arg_q} \
                {arg_c} \
                {arg_z} \
                {arg_j} \
                {arg_s} \
                {arg_r} \
                -R {os.getcwd()}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-9]}
            '''.replace("\n", " ")

            print("Gibbs Command: ", gibbs_command)
            print(os.popen(gibbs_command).read())
