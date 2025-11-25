
def compile_target_list(target_list_path):
    target_file = open(target_list_path, "r")
    targets = target_file.read()
    target_file.close()
    target_list = targets.split("\n")
    return target_list

