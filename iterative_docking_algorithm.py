#!/usr/bin/env python3

import os
import re
from os.path import isfile, splitext
from numpy import array
import pandas as pd
from math import ceil
import subprocess

# Get working directory for data
def get_working_directory():
    while True:
        dir_path = input('Please input working directory (copy full path): ')

        if os.path.isdir(dir_path):
            return dir_path
        else:
            print(f'Error: Directory {dir_path} not found...')

#get the ligand directory
def get_ligand_directory():
   while True:
        dir_path = input('Please input the ligand directory:\n' + os.getcwd() + '\\')
        if dir_path in os.listdir():
            return dir_path
        else:
            print(f'Error: Directory {dir_path} not found...')

#get a list of ligands found in the folder (with extension)
def identify_ligands(lig_dir):
    while True:
        ligand_files = [
            f
            for f in os.listdir(lig_dir)
            if isfile(os.path.join(lig_dir,f)) and splitext(f)[1] == '.pdbqt']
        if len(ligand_files) > 0:
            print('The following ligands have been identified:\n')
            print(*ligand_files, sep = '\n')
            return ligand_files
        else:
            input('No ligand files found, check the folder and press Enter to try again.')

# Collect receptor file name for log files
def get_receptor():
    while True:
        receptor_file = input('Receptor file (with extension): ')

        # Check that the file exists and returns filename if True
        if os.path.isfile(receptor_file):
            print('Receptor file accepted...')
            return receptor_file
        # If file not found - re-request file name
        else:
            print(f'Error: {receptor_file} not found')
            print('Please check filenames and re-enter')

#not tested nor used
def get_flex_receptor():
    enable_flex = input('Do you want to use flexible docking?(y/n) ')

    if enable_flex == 'y':

        # Check that the file exists and returns filename if True
        flex_receptor_file = input('Receptor file: ')
        if os.path.isfile(flex_receptor_file):
            print('Flexible Receptor file accepted...')
            return flex_receptor_file

        # If file not found - re-request file name
        else:
            print(f'Error: {flex_receptor_file} not found')
            print('Please check filenames and re-enter')
            retry = get_flex_receptor()
            return retry

    else:
        return None

# Collect coordinates for log files
def get_coordinates():
    # Function to determine if a value is a float, returns True or False accordingly
    def is_float(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    while True:
        min_coordinates_input = input('Min limit co-ordinates (x y z): ')  # Receive coordinates in format x y z
        max_coordinates_input = input('Max limit co-ordinates (x y z): ')  # Receive coordinates in format x y z
        min_coord = min_coordinates_input.split()  # Generate list of lowest coordinate point
        max_coord = max_coordinates_input.split()  # Generate list of highest coordinate point
        all_coord = min_coord + max_coord

        # Check all 3 coordinates have been provided
        if len(min_coord) != 3 or len(max_coord) != 3:
            print('Error: Please enter exactly 3 coordinates for both min and max separated by a space')
        # Check all coordinates are numeric values
        elif not all(map(is_float, all_coord)):
            print('Error: not all coordinates are numbers')
        # check if values are min and max
        elif not all(x <= y for x, y in zip(min_coord, max_coord)):
            print('At least one minimum coordinate exceeds the respective maximum coordinate.')
        # accept and return input if all values are present and correct
        else:
            print('Coordinates accepted...')
            return min_coord, max_coord

# Collect box size parameters for log files
def get_box_size():
    # Check if a given value is a positive integer
    def check_pos(value):
        n = int(value)
        if n > 0:
            return True
        else:
            return False

    while True:
        box_size_input = input('Box size (x y z) integers: ')
        box = box_size_input.split()

        # Check all coordinates are integers
        try:
            all(map(int, box))
            # Check all 3 coordinates have been provided
            if len(box) != 3:
                print('Error: Please enter exactly 3 coordinates separated by a space.')
            # Check all integers are positive
            elif not all([check_pos(n) for n in box]):
                print('Error: Box values must be positive integers.')
            # accept and return input if all values are present and correct
            else:
                print('Box parameters accepted...')
                return box
        except ValueError:
            print('Error: Box values must be positive integers')

# Collect seeds if any for log files not tested for multiple
def get_seeds():
    while True:
        seed = input('Seed(s) (optional, default is 0) [a b c ...] integer: ')  # Receive seeds in format 1 2 3
        seed_list_input = seed.split()  # separate seeds into list

        # Check seeds are integers (can be positive or negative)
        try:
            seed_map = map(int, seed_list_input)  # generate a list of seeds as integers
            # If no input give - default seed is 0
            if len(seed_list_input) < 1:
                return [0]
            # accept and return input if all values are present and correct
            else:
                print('Seed(s) accepted...')
                return list(seed_map)
        # If not integer value is received ask for input again
        except ValueError:
            print('Error: seed(s) must be an integer')

# Write a basic  config file for AutoDock Vina - ligands, receptors, outputs, seed (0), grid size and coordinates.
def config_writer(ligand, receptor_input, flex_recept, coord, box, iteration, lig_folder, seeding=0):
    file_name = os.path.join(ligand, f'{iteration}-config-{seeding}.txt')
    output_file = f'{ligand}/conformations-{seeding}/{iteration}-{seeding}.pdbqt'

    # Generate the config file using inputs
    with open(file_name, 'w') as config:
        if flex_recept is not None:
            config.write(f'flex = {flex_recept}\n')

        config.write(f'receptor = {receptor_input}\n'  # receptor file
                     f'ligand = {lig_folder}/{ligand}.pdbqt\n\n'  # ligand file
                     f'center_x = {coord[0]}\n'  # grid coordinates for box
                     f'center_y = {coord[1]}\n'
                     f'center_z = {coord[2]}\n\n'
                     f'size_x = {box[0]}\n'  # box size parameters 
                     f'size_y = {box[1]}\n'
                     f'size_z = {box[2]}\n\n'
                     f'out = {output_file}\n'  # Docking conformation output
                     'exhaustiveness = 8\n' #number of attempts
                     f'seed = {seeding}'  # Seed for experiment
                     )

        file_names.append(file_name)  # append file name to list for automatic docking

def get_binding_data_csv(ligand):
    # Extract data from log files
    def extract_data(file):
        file_name = os.path.basename(file)
        file_name_split = file_name.split('-')  # extract iteration number and seed number from file name as a list
        iteration = file_name_split[0]  # store iteration number
        seed = file_name_split[-1].replace('.txt', '')  # store seed number

        # Open log file to collect binding data
        with open(file, 'r') as log:
            try:
                logs = log.readlines()  # convert each line into a string
                for i in range(38, 47, 1): #38-47 in the case of exhaustiveness 8 producing 9 lines, adjust accordingly
                    raw_data = logs[i].split()  # split values into list
                    map_data = map(float, raw_data)  # convert strings to floats
                    list_data = list(map_data)  # convert map to list
                    list_data.insert(0, iteration)  # Add ligand name as 1st element
                    list_data.insert(1, seed)  # Add seed number as 2nd element
                    log_list.append(list_data)  # Append data to overall data list
            except IndexError:
                print(f'Error could not read file: {file}')
                pass
            except UnicodeDecodeError:
                print(f'Error could not read file: {file}')
                pass
            except ValueError:
                print(f'Error could not read file: {file}')
                pass

    if "results" not in os.listdir(ligand):
        os.mkdir(os.path.join(ligand, 'results'))
        print('results directory created')

    # this is only needed to support multiple seeds
    log_dirs = [os.path.join(ligand, folder) for folder in os.listdir(ligand) if re.match(r"^logs-\d+$", folder)]  # Identify logs directories within pwd

    log_list = []  # store rows of data extracted from log files

    # Access and collect log data from log files in each directory
    for directory in log_dirs:
        print(f'Extracting Log data from {directory}...')
        [extract_data(os.path.join(directory, file)) for file in os.listdir(directory)]
        print(f'Log data extracted from {directory}...')

    log_array = array(log_list)  # convert collected data into a 2D NumPy array
    # convert 2D NumPy array to DataFrame
    log_df = pd.DataFrame(log_array, columns=['iteration', 'seed', 'binding mode', 'affinity (kcal/mol)',
                                              'distance from best mode rmsd l.b', 'distance from best mode rmsd u.b'])

    path = os.path.join(ligand, 'results', 'summary-across-iterations.csv')
    log_df.to_csv(path, index=False, header=True)  # Save DataFrame to CSV file
    print(f'Log data retrieved and saved as {path}')  # Tell user location of CSV file

def get_covering_matrix(minim, maxim, two_step):
    x_min = float(minim[0]) # convert string to float for math operators
    x_max = float(maxim[0])
    y_min = float(minim[1])
    y_max = float(maxim[1])
    z_min = float(minim[2])
    z_max = float(maxim[2])

    x_range = x_max - x_min # get the max range in each dimension
    y_range = y_max - y_min
    z_range = z_max - z_min
    x_step = float(two_step[0]) / 2 # get the step size that is equal to half box size in each dimension
    y_step = float(two_step[1]) / 2
    z_step = float(two_step[2]) / 2
    x_counter = ceil(x_range / x_step) + 1 # calculate how many iteration is required to exhaustively cover the range
    y_counter = ceil(y_range / y_step) + 1
    z_counter = ceil(z_range / z_step) + 1
    vector_list = []

    # generate a list of coordinate vectors to cover the entire cubic space
    for a in range(0,x_counter, 1):
        for b in range(0, y_counter, 1):
            for c in range(0, z_counter, 1):
                vector = (x_min + x_step * a, y_min + y_step * b, z_min + z_step * c)
                vector_list.append(vector)

    coord_vectors = array(vector_list)

    return coord_vectors

# Run Script
if __name__ == '__main__':
    working_directory = get_working_directory()  # Ask for working directory containing ligand and receptor
    os.chdir(working_directory)  # change directory to input
    receptor = get_receptor()  # Receive file name of receptor
    flex_receptor = get_flex_receptor()
    ligand_dir = get_ligand_directory()  # Receive ligand directory name
    ligands = identify_ligands(ligand_dir) #make a list of the ligand files
    min_coordinates, max_coordinates = get_coordinates()  # Receive co-ordinates of box in format x y z
    box_size = get_box_size()  # Receive box size in format x y z
    seed_list = get_seeds()  # Receive list of seeds for docking experiments
    coord_matrix = get_covering_matrix(min_coordinates, max_coordinates, box_size)  # Create vector space for 3D iteration
    print('Starting docking.')

    for ligand_m in ligands:
        file_names = []  # Empty list for names of config files
        ligand_name = splitext(ligand_m)[0]
        # Generate config files for AutoDock Vina
        if ligand_name not in os.listdir():
            os.mkdir(ligand_name)

        for seed_m in seed_list:
            if f'conformations-{seed_m}' not in os.listdir(ligand_name):
                os.mkdir(os.path.join(ligand_name, f'conformations-{seed_m}'))

            if f'logs-{seed_m}' not in os.listdir(ligand_name):
                os.mkdir(os.path.join(ligand_name, f'logs-{seed_m}'))

        for iteration_counter, coordinates in enumerate(coord_matrix, start = 1):
            [config_writer(ligand_name, receptor, flex_receptor, coordinates, box_size, iteration_counter, ligand_dir, seeding=y) for y in seed_list]

        max_iterations = len(file_names)
        # Run AutoDock Vina using generated config files list
        for config_input in file_names:
            config_input_name = os.path.basename(config_input)
            file_n = config_input_name.split('-')  # extract iteration number and seed number from file name as a list
            iteration_num = file_n[0]  # store iteration number
            seed_num = file_n[-1].replace('.txt', '')  # store seed number
            log_path = os.path.join(ligand_name, f'logs-{seed_num}', f'{iteration_num}-log-{seed_num}.txt')
            with open(log_path, 'w') as log_file:
                result = subprocess.run(
                    ['vina_1.2.5_win.exe', '--config', config_input],
                    stdout = log_file,
                    stderr = log_file,
                    text = True
                )
            if result.returncode == 0:
                print(f'No. {iteration_num}/{max_iterations} space vector done.')
            else:
                print(f'No. {iteration_num}/{max_iterations} space vector resulted in error.')

        get_binding_data_csv(ligand_name)  # Collect binding affinity data and save as a CSV file

    print('Done with all processing.')
