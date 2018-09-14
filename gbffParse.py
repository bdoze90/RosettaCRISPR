"""THIS DOES NOT BELONG TO ROSETTA! TRANSFER TO CASPERgenome WHEN IT IS LINKED TO GITHUB."""

from SeqTranslate import SeqTranslate
from operator import itemgetter
import os,sys

input_gbff = sys.argv[1]
input_cspr = sys.argv[2]
output_csv = sys.argv[3]

grna_temp_storage = list()
# main storage vehicle
ScaffGeneDict = dict()


def parse_gbff(file_path):
    f = open(file_path)
    for line in f:
        if line.startswith("##"):
            print(line)
        else:
            istics = line[:-1].split("\t")
            # grab only the exon lines
            if istics[2] == 'exon':
                keywds = istics[8].split(";")
                geneid = keywds[1][keywds[1].find("=")+1:-5]
                # reports the gene identification, positions, and direction
                mistics = [geneid,int(istics[3]),int(istics[4]),istics[6]]
                # stores the data in a dict of lists for exon/gene and chromosome/scaffold lookup:
                if istics[0] not in ScaffGeneDict:
                    ScaffGeneDict[istics[0]] = [mistics]
                else:
                    ScaffGeneDict[istics[0]].append(mistics)
    f.close()
    # sort by the end location (for parsing index)
    for C_S in ScaffGeneDict:
        ScaffGeneDict[C_S] = sorted(ScaffGeneDict[C_S], key=itemgetter(1))


def sort_grnas_by_genes(cspr_file):
    S = SeqTranslate()
    f = open(cspr_file)
    for line in f:
        if line.startswith(">"):
            if grna_temp_storage:
                if cur_cs in ScaffGeneDict:
                    assign_to_genes(cur_cs)
                grna_temp_storage.clear()
            cur_cs = line[1:line.find("(")-1]
            print("running " + cur_cs)
        elif line.startswith("REPEAT"):
            if cur_cs in ScaffGeneDict:
                assign_to_genes(cur_cs)
            grna_temp_storage.clear()
            break
        else:
            grna_temp_storage.append(S.decompress_csf_tuple(line[:-1]))
    f.close()


# Function places the guideRNA sequences into a gene locus if its PAM location is within the exon start and end.
def assign_to_genes(cs_index):
    print("Number of genes in scaffold:")
    total = len(ScaffGeneDict[cs_index])
    print(total)
    index_start = -1
    progressindex = 0
    grna_index = 0
    for gene in ScaffGeneDict[cs_index]:
        while grna_temp_storage[grna_index][0] < (gene[2]-20):  # PAM site is before end of gene (+20 for intron PAMs)
            g_tup = grna_temp_storage[grna_index]
            if gene[1] < (g_tup[0]+20):  # PAM site is after start of gene (-20 for promoter/intron PAMs)
                index_start = grna_index
                item = list(g_tup) + [is_istop(g_tup,gene[1],gene[3],grna_temp_storage[grna_index][4])]
                gene.append(item)
            if grna_index < len(grna_temp_storage)-1:  # Still in the right scaffold
                grna_index += 1
            else:
                if index_start > 0:
                    grna_index = index_start
                break
        progressindex += 1
        progressBar(progressindex,total)
    generate_report(cs=cs_index)
    ScaffGeneDict[cs_index].clear()


def is_istop(grna,atg_pos, strand, gstrand):
    if strand == gstrand:  # checks to see if the PAM and the gRNA sequence are on the same strand
        codons = ['CAG','CAA','CGA']
    else:
        codons = ['CCA']
    locs = list()
    for codon in codons:
        qc = grna[1][2:11].find(codon)
        if qc != -1:  # Looking for stop codon somewhere between 3rd and 9th bp (0 indexing)
            locs.append(qc + 2)  # Fixes indexing of old system
    if locs:
        for loc in locs:
            if (atg_pos - loc-20+grna[0]) % 3 == 0 and grna[4]:
                return "YES"
            elif (atg_pos - loc - 20 + grna[0]) % 3 == 0 and not grna[4]:
                return "YES"
            else:
                return "off"
    else:
        return "NO"


def generate_report(cs):
    filename = output_csv
    if os.path.exists(filename):
        ap_stat = 'a'
    else:
        ap_stat = 'w'
    f = open(filename, ap_stat)
    f.write("SCAFFOLD: " + cs + "\n")
    for gene in ScaffGeneDict[cs]:
        f.write("GENE: " + gene[0] + "," + str(gene[1]) + "," + str(gene[2]) + "," + gene[3] + "\n")
        for grna in gene[4:]:
            for i in range(len(grna)):
                if i != 4:
                    f.write(str(grna[i]) + ",")
            if grna[4]:
                f.write("+")
            else:
                f.write("-")
            f.write("\n")
    print(" completed " + cs)
    f.close()


def progressBar(value, endvalue, bar_length=30):
    percent = float(value)/endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow+spaces,int(round(percent*100))))
    sys.stdout.flush()


parse_gbff(input_gbff)

sort_grnas_by_genes(input_cspr)






