input_fasta_path = '/Users/emluu/Downloads/'
filename = 'aldh_cousins.fasta'


def clean_fasta_header(header_string):
    delimiter_string = '|'
    name_index = 1
    string_list = header_string.split(delimiter_string)
    new_header_string = '>'+string_list[name_index]+'\n'
    return new_header_string


def clean_fasta(fasta_path):
    with open(fasta_path, 'r') as original_fasta:
        fasta_text = original_fasta.readlines()

    for index in range(len(fasta_text)):
        line = fasta_text[index]
        if ">" in line:
            new_header_string = clean_fasta_header(line)
            fasta_text[index] = new_header_string

    with open(fasta_path, 'w') as new_fasta:
        new_fasta.writelines(fasta_text)


if __name__ == '__main__':
    clean_fasta(fasta_path=input_fasta_path+filename)
