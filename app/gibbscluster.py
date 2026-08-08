import sys
import os
import subprocess

# Get command-line arguments
task_id = sys.argv[1]
data_mount = sys.argv[2]
mhc_class = sys.argv[3]
motif_length = int(sys.argv[4])

# Determine the project root
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# Set seq2logo path
seq2logo_path = os.path.join(project_root, "app/tools/seq2logo-2.1/Seq2Logo.py")

task_path = os.path.join(data_mount, task_id)

# Default settings for both classes
num_clusters = '1-6'
mc_temperature = '1.5'
num_temp_steps = '20'
num_iterations = '10'
num_seeds = '5'
penalty_lambda = '0.8'
small_cluster_weight = '5'
sequence_weighting_type = '0'
background_model = '1'  # Pre-calculated Uniprot
use_trash_cluster = '1'
trash_cluster_threshold = '2'

# Command for MHC class I
if mhc_class == "I":
    max_insertion_length = '1'
    max_deletion_length = '4'
    indel_move_interval = '10'
    shift_move_interval = '20'  # Class I specific
    phase_shift_move_interval = '100'  # Class I specific
    hydrophobic_p1_preference = '0'

    input_file_ends_with = '_8to14mer.txt'

    for sample in os.listdir(task_path):
        sample_path = os.path.join(task_path, sample)
        if not os.path.isdir(sample_path):
            continue
        for replicate in os.listdir(sample_path):
            if replicate.endswith(input_file_ends_with):
                print('Replicate file:', replicate)
                cmd = [
                    'perl', f'{project_root}/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl',
                    '-f', f'{sample_path}/{replicate}',
                    '-H', 'R',
                    '-G', seq2logo_path,
                    '-g', num_clusters,
                    '-k', str(os.cpu_count()),
                    '-l', str(motif_length),
                    '-R', f'{project_root}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-13]}',
                    '-i', num_iterations,
                    '-t', mc_temperature,
                    '-n', num_temp_steps,
                    '-b', penalty_lambda,
                    '-q', small_cluster_weight,
                    '-S', num_seeds,
                    '-c', sequence_weighting_type,
                    '-z', background_model,
                    '-j', trash_cluster_threshold,
                    '-D', max_deletion_length,
                    '-I', max_insertion_length,
                    '-u', indel_move_interval,
                    '-r', shift_move_interval,
                    '-s', phase_shift_move_interval,
                ]
                # Boolean-only flags (no argument in GibbsCluster's getopts spec) must come
                # last: Getopt::Std stops parsing at the first token that doesn't start with
                # '-', so a bare flag followed by a value silently drops every flag after it.
                if use_trash_cluster == '1':
                    cmd.append('-T')
                cmd.append('-C')
                if hydrophobic_p1_preference == '1':
                    cmd.append('-p')
                print("Gibbs Command for MHC Class I:", ' '.join(cmd))
                result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
                print(result.stdout)
                if result.returncode != 0:
                    print(result.stderr, file=sys.stderr)

# Command for MHC class II
elif mhc_class == "II":
    max_insertion_length = '0'  # No insertion for Class II
    max_deletion_length = '0'   # No deletion for Class II
    indel_move_interval = '20'  # For Class II, different move interval
    shift_move_interval = '20'  # Shift moves are activated
    phase_shift_move_interval = '100'  # Phase shift interval for Class II
    hydrophobic_p1_preference = '1'  # Class II specific

    input_file_ends_with = '_12to20mer.txt'

    for sample in os.listdir(task_path):
        sample_path = os.path.join(task_path, sample)
        if not os.path.isdir(sample_path):
            continue
        for replicate in os.listdir(sample_path):
            if replicate.endswith(input_file_ends_with):
                print('Replicate file:', replicate)
                cmd = [
                    'perl', f'{project_root}/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl',
                    '-f', f'{sample_path}/{replicate}',
                    '-H', 'R',
                    '-G', seq2logo_path,
                    '-g', num_clusters,
                    '-k', str(os.cpu_count()),
                    '-l', str(motif_length),
                    '-R', f'{project_root}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-14]}',
                    '-i', num_iterations,
                    '-t', mc_temperature,
                    '-n', num_temp_steps,
                    '-b', penalty_lambda,
                    '-q', small_cluster_weight,
                    '-S', num_seeds,
                    '-c', sequence_weighting_type,
                    '-z', background_model,
                    '-j', trash_cluster_threshold,
                    '-D', max_deletion_length,
                    '-I', max_insertion_length,
                    '-u', indel_move_interval,
                    '-r', shift_move_interval,
                    '-s', phase_shift_move_interval,
                ]
                # Boolean-only flags (no argument in GibbsCluster's getopts spec) must come
                # last: Getopt::Std stops parsing at the first token that doesn't start with
                # '-', so a bare flag followed by a value silently drops every flag after it.
                if use_trash_cluster == '1':
                    cmd.append('-T')
                if hydrophobic_p1_preference == '1':
                    cmd.append('-p')
                print("Gibbs Command for MHC Class II:", ' '.join(cmd))
                result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
                print(result.stdout)
                if result.returncode != 0:
                    print(result.stderr, file=sys.stderr)
