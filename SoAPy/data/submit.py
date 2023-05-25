import os

# Get current working directory.
cwd = os.getcwd()

# List of directories.
relative_dir_list = DIRECTORY_LIST
dir_parameters = DIRECTORY_PARAMETERS

# Change to test specific directory.
for a in range(len(relative_dir_list)):
    #print(f'Spectroscopy: {relative_dir_parameters[a][0]} || Solvent Shell Type: {relative_dir_parameters[a][1]} || Functional: {relative_dir_parameters[a][2]} || Basis: {relative_dir_parameters[a][3]} || Distance: {relative_dir_parameters[a][4]} || Snapshots: {relative_dir_parameters[a][5]}')
    print ('Spectroscopy: ', dir_parameters[a][0], '|| Solvent Shell Type: ', dir_parameters[a][1], '|| Functional: ', dir_parameters[a][2], '|| Basis: ', dir_parameters[a][3], '|| Distance: ', dir_parameters[a][4], '|| Snapshots: ', dir_parameters[a][5])
    print(" ")
    conformer_count = 1
    while conformer_count <= int(dir_parameters[a][5]):
            
        # Submit jobs to ARC.
        os.chdir('{}/{}/cmpd_{}'.format(cwd,relative_dir_list[a],conformer_count))
        os.system("dos2unix G09_sub_SLURM.sh")
        os.system("sbatch G09_sub_SLURM.sh")
        print('Job submitted for compound', conformer_count, '.')
        
        conformer_count += 1

