import requests
import os

# Functions
def file_to_list(file_path):
    """Reads a text file and returns a list where each element is a line from the file."""
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        return [line.strip() for line in lines]
    except FileNotFoundError:
        print(f"Error: File not found at path: {file_path}")
        return None

prefix = 'AF-'
suffix = '-F1-model_v6.cif'
website = "https://alphafold.ebi.ac.uk/files/"
downloadPath = "" # edit this path
uniProtIDs_path = 'ALDH_af3_search.txt' # edit this with either local or full path to textfile

uniProtIDs = file_to_list(uniProtIDs_path)
for uniProtID in uniProtIDs:
    r = requests.get(website+prefix+uniProtID+suffix)
    with open(downloadPath+uniProtID+".cif", 'wb') as f:
        f.write(r.content)
