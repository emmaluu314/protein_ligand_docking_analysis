import user_functions as uf
import score_functions as sf
import os
import pandas as pd

# User defined variables:
current_username = 'emmaluu'
remote_server_name = 'hive.hpc.ucdavis.edu'
main_server_path = "/quobyte/jbsiegelgrp/emmaluu/standard_tools/docking/cofactors/NAD_16/LDH/runs/" # main location of subfolders
main_local_path = '/Users/emluu/Documents/Siegel lab/standard/Rosetta Ligand/LDH/2025_10_30_16NAD_mutants/'
results_folder_name = "results/" #this is what your results folder is called. If you do not use a results folder, write "".
target_list_filename = 'jobs.txt' #name of file with a column of the subfolders you'd like to analyze
target_list_path = main_local_path + target_list_filename


def extract_scores_from_server(list_path, server_name, username, main_server_path, main_local_path):
    target_list = uf.compile_target_list(list_path)
    sf.pull_from_cluster(target_list, server_name, username, main_server_path, results_folder_name, main_local_path)
    # need to make it where I can have results flag

if __name__ == '__main__':
    # need to make a flag that indicates whether there is a results subfolder
    extract_scores_from_server(target_list_path, remote_server_name, current_username, main_server_path,
                               main_local_path)
