# Small script to further filter contaminants or E value
# ***WARNING*** Assumes first column is query
# probably some bugs still
import argparse
import csv
import sys

input_path = ""
e_val = float
species_col = int
e_col = int
tax_path = ""
delim = int
header = int
contam_only = False
fasta_file = ""

version = True      # True for 3 false for 2

bacteria = {}
fungi = {}
insect = {}

contam_seq = []
filtered_seq = []
no_hit_seq = []

bacteria_db_path = "bacteria_db.txt"
fungi_db_path = "fungi_db.txt"
insects_db_path = "insects_db.txt"


def init_argparse():
    global input_path
    global e_val
    global species_col
    global e_col
    global tax_path
    global delim
    global header
    global fasta_file
    global contam_only
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest='input_path', help='Path to input file', type=str)
    parser.add_argument('-e', action='store', dest='e_val', help='Minimum E-Value to accept in format 3E-12',
                        type=str, default=10)
    parser.add_argument('-s', action='store', dest='species_col', help='Species column number, leftmost being 0', type=int)
    parser.add_argument('-c', action='store', dest='e_col', help='E-val column number, leftmost being 0', type=int)
    parser.add_argument('-d', action='store', dest='tax_path', help='Path to taxonomic databases', type=str)
    parser.add_argument('-p', action='store', dest='delim', help='Deliminator, either 1 for tab or 2 for comma',
                        type=int, default=1)
    parser.add_argument('-r', action='store_true', dest='header', help='Is there a header? Use this flag if there is',
                        default=False)
    parser.add_argument('-a', action='store', dest='fasta_path',  help="Path to fasta file you'd like to pull sequences"
                       "from.", type=str)
    parser.add_argument('-x', action='store_true', dest='contam_only', help="Use this flag if you just want the contaminants outputted",
                        default=False)
    args = parser.parse_args()
    input_path = args.input_path
    e_val = float(args.e_val)
    species_col = args.species_col
    e_col = args.e_col
    tax_path = args.tax_path
    delim = args.delim
    header = args.header
    contam_only = args.contam_only


    fasta_file = args.fasta_path

def parse_file():
    d = ""
    # col_num = int
    if (delim == 1):
        d = '\t'
    else:
        d = ','
    if (version):
        contams = open("entap_contaminants.tsv", 'w', newline='')
        filtered = open("entap_filtered.tsv", 'w', newline='')
        nohit = open("entap_removed.tsv", 'w', newline='')
    else:
        contams = open("entap_contaminants.tsv", 'wb')
        filtered = open("entap_filtered.tsv", 'wb')
        nohit = open("entap_removed.tsv", 'wb')

    writer_contams = csv.writer(contams, delimiter=d)
    writer_filtered = csv.writer(filtered, delimiter=d)
    writer_nohit = csv.writer(nohit, delimiter=d)
    with open(input_path, 'r') as file:
        read = csv.reader(file, delimiter=d)
        if header:
            next(file)
        query1 = ""
        query2 = ""
        line1 = []

        ### Use if you want column number ###
        # if (header):
        #     line = file.readline()
        #     col_num = line.count(d) + 1
        #
        # else:
        #     line = file.readline()
        #     col_num = line.count(d) + 1
        #     query1 = line[0:line1.find(d)]
        #     line1 = line.split(d)
        #     print(line1)

        for line2 in read:
            query2 = line2[0]
            try:
                query1 = line1[0]
            except:
                pass
            if is_contaminant(line2[species_col]):
                writer_contams.writerow(line2)
                contam_seq.append(line2[0])
                continue
            if query2 == query1:
                i = best_hit(line1[e_col], line2[e_col])
                if i == 0:
                    writer_nohit.writerow(line1)
                    no_hit_seq.append(line1[0])
                elif i == 1:
                    writer_nohit.writerow(line2)
                    no_hit_seq.append(line2[0])
                    continue
                elif i == 2:
                    writer_nohit.writerow(line1)
                    no_hit_seq.append(line1[0])
            else:
                writer_filtered.writerow(line1)
                try:
                    filtered_seq.append(line1[0])
                except: pass    # End

            line1 = line2
    contams.close()
    filtered.close()
    nohit.close()


def best_hit(e1, e2):
    # return 1 for hit 1, 2 for hit 2
    # 0 for none above minimum

    e1_val = float(e1)
    e2_val = float(e2)
    if e1_val > e_val and e2_val > e_val:
        return 0
    if e1_val < e2_val:
        return 1
    else:
        return 2


def is_contaminant(species):
    # Separate so can say where contaminant came from, if wanted...
    species_low = species.lower()
    if fungi.get(species_low):
        return True
    if bacteria.get(species_low):
        return True
    if insect.get(species_low):
        return True
    return False


def init_contaminant():
    global bacteria
    global insect
    global fungi

    # read into memory
    full_bacteria = tax_path + bacteria_db_path
    full_fungi = tax_path + fungi_db_path
    full_insect = tax_path + insects_db_path
    with open(full_bacteria, 'r') as bact_file:
        for line in bact_file:
            line = line.rstrip()
            bacteria[line.lower()] = 1
    with open(full_fungi, 'r') as fung_file:
        for line in fung_file:
            line = line.rstrip()
            fungi[line.lower()] = 1
    with open(full_insect, 'r') as insec_file:
        for line in insec_file:
            line = line.rstrip()
            insect[line.lower()] = 1


def write_files():
    filtered_file = open("entap_filtered.fasta", 'w')
    removed_file = open("entap_removed.fasta", 'w')
    contaminant_file = open("entap_contaminants.fasta", 'w')

    file_index = 0
    seq_found = False

    filtered_lines = []
    contam_line  = []
    no_hit_lines = []
    files = [filtered_lines, contam_line, no_hit_lines]

    if contam_only:
        with open(fasta_file, 'r') as fasta_input:
        # slow...fix todo
            for line in fasta_input:
                if not line:
                    continue
                if line[0] == '>':
                    seq_found = False
                    for c_seq in contam_seq:
                        if line.find(c_seq) == 1:
                            file_index = 1
                            files[file_index].append(line)
                            seq_found = True
                            contam_seq.remove(c_seq)
                            break
                else:
                    if not seq_found: continue
                    files[file_index].append(line)
        contaminant_file.writelines(contam_line)
    else:
        with open(fasta_file, 'r') as fasta_input:
            # VERY slow can change to class...
            # or reverse todo
            for line in fasta_input:
                if not line:
                    continue
                if line[0] == '>':
                    seq_found = False
                    for f_seq in filtered_seq:
                        if line.find(f_seq) == 1:
                            file_index = 0
                            # files[file_index].writelines(line)
                            files[file_index].append(line)
                            seq_found = True
                            filtered_seq.remove(f_seq)
                            break
                    if seq_found: continue
                    for c_seq in contam_seq:
                        if line.find(c_seq) == 1:
                            file_index = 1
                            # files[file_index].write(line)
                            files[file_index].append(line)
                            seq_found = True
                            contam_seq.remove(c_seq)
                            break
                    if seq_found: continue
                    for r_seq in no_hit_seq:
                        if line.find(r_seq) == 1:
                            file_index = 1
                            # files[file_index].write(line)
                            files[file_index].append(line)
                            no_hit_seq.remove(r_seq)
                            break
                else:
                    if not seq_found: continue
                    files[file_index].append(line)
        filtered_file.writelines(filtered_lines)
        removed_file.writelines(no_hit_lines)
        contaminant_file.writelines(contam_line)
    filtered_file.close()
    removed_file.close()
    contaminant_file.close()


def init_version():
    global version
    version = sys.version_info > (3,0)


def main():
    global bacteria
    global fungi
    global insect
    init_version()
    init_argparse()
    init_contaminant()
    parse_file()
    del bacteria, fungi, insect
    write_files()

if __name__ == "__main__":
    main()