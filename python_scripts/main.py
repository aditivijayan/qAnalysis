
import _init_
from constants import *
from set_path import *
from config import *
from functions import *
import subprocess
import argparse
import matplotlib.colors as mcolors
import matplotlib.cm


def run_script_for_each_file(folder_name, i , N):
     
    list_file = glob.glob(folder_name + "plt*/")
    list_of_lists_of_files = np.array_split(list_file, N)
    loop_over_list_file = list_of_lists_of_files[i]
    print("Reached here", folder_name,i)
    for filename in loop_over_list_file:
        subprocess.run(["python", "get_slices.py", folder_name, filename])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot slices for quokka plot files.')
    parser.add_argument('--input_folder', type=str, help='Path to input folder containing plt files')
    parser.add_argument('--N', type=int, help='Nodes to analyse files with.')
    parser.add_argument('--i', type=int, help='Go over every i-th file. Has to be < N.')


    args = parser.parse_args()
    infolder = args.input_folder
    i = args.i
    N = args.N
    if(i>=N):
        print("i should less than N! Choose i=0 to N-1!")
        print("Exiting now. Ciao!")
        exit()
    data_path = os.path.join(kronos, 'sims/',  infolder)
    data_path = os.path.join(orion_path,  infolder)
    run_script_for_each_file(data_path, i , N)
