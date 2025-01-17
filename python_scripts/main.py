
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import subprocess
import argparse

def run_script_for_each_file(folder_name):
    # Check if the folder exists
    if not os.path.isdir(folder_name):
        print(f"Folder '{folder_name}' does not exist.")
        return
    print('here here')
    for filename in os.listdir(folder_name):
        # print(filename)
        # if(filename.split('_')[0]=='plt'):
        if('plt' in filename and 'proj' not in filename):
            file_path = os.path.join(folder_name, filename)
            print(file_path)
            subprocess.run(["python", "child_process.py", folder_name, filename])

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
    parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
    
    
    args = parser.parse_args()
    infolder = args.input_folder

    data_path = os.path.join('/g/data/jh2/av5889/quokka_myrepo/quokka/sims/', infolder)
    
   
    run_script_for_each_file(data_path)