import sys
import os
import subprocess

# Get command-line arguments
task_id = sys.argv[1]
data_mount = sys.argv[2]
sample = sys.argv[3]
replicate = sys.argv[4] + '_9mer.txt'
cluster = sys.argv[5]
mhc_class = sys.argv[6]
motif_length = int(sys.argv[7])

cmd = [
    'perl', './app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA_for_seqlogo.pl',
    '-f', f'{data_mount}/{task_id}/{sample}/{replicate}',
    '-H', 'R',
    '-G', './../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py',
    '-g', cluster,
    '-k', str(os.cpu_count()),
    '-T',
]

if mhc_class == "I":
    cmd.append('-C')
    cmd += ['-I', '0', '-D', '0', '-S', '5', '-b', '0.8', '-q', '5', '-c', '0', '-z', '1', '-j', '2']
else:
    cmd += ['-I', '0', '-D', '0', '-S', '5', '-b', '0.8', '-q', '5', '-c', '0', '-z', '1', '-j', '2',
            '-s', '100', '-r', '20']

cmd += [
    '-l', str(motif_length),
    '-R', f'{os.getcwd()}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-9]}',
]

print("Gibbs Command:", ' '.join(cmd))
result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print(result.stderr, file=sys.stderr)
