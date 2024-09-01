
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import subprocess
import argparse

def run_script_for_each_file(folder_name):
    for filename in os.listdir(folder_name):
        print("Reached inside run script", filename)    
        if 'plt' in filename and 'proj' not in filename:
            file_path = os.path.join(folder_name, filename)
            subprocess.run(["python", "child_process.py", folder_name, filename])

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
    parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
    
    
    args = parser.parse_args()
    infolder = args.input_folder

    path_to_folder = os.path.join(data_home, infolder)
    print("paht to folder=", path_to_folder)
    run_script_for_each_file(path_to_folder)
