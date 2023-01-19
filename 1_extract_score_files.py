import user_functions as uf
import score_functions as sf

# User defined variables:
current_username = 'emmaluu'
remote_server_name = 'cacao.genomecenter.ucdavis.edu'
main_remove_server_path = '/share/siegellab/emmaluu/docking/template_for_automation/CAR/NMN/'
main_folder_path = '/Users/emluu/Documents/Siegel lab/Rosetta Ligand/ARPA-E/CAR/Cousins temp/'
target_list_path = main_folder_path + 'target_list.txt'


def extract_scores_from_server(list_path, server_name, username, main_server_path, main_local_path):
    target_list = uf.compile_target_list(list_path)
    sf.pull_from_cluster(target_list, server_name, username, main_server_path, main_local_path)
    # need to make it where I can have results flag

if __name__ == '__main__':
    # need to make a flag that indicates whether there is a results subfolder
    extract_scores_from_server(target_list_path, remote_server_name, current_username, main_remove_server_path,
                               main_folder_path)
