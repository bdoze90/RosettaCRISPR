"""THIS DOES NOT BELONG TO ROSETTA! TRANSFER TO CASPERgenome WHEN IT IS LINKED TO GITHUB."""

from SeqTranslate import SeqTranslate
from operator import itemgetter

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
        ScaffGeneDict[C_S] = sorted(ScaffGeneDict[C_S], key=itemgetter(2))


def sort_grnas_by_genes(cspr_file):
    S = SeqTranslate()
    f = open(cspr_file)
    for line in f:
        if line.startswith(">"):
            if grna_temp_storage:
                assign_to_genes(cur_cs)
                grna_temp_storage.clear()
            cur_cs = line[1:-5]
            print("running " + cur_cs)
        else:
            grna_temp_storage.append(S.decompress_csf_tuple(line[:-1]))
    f.close()


def assign_to_genes(cs_index):
    for gene in ScaffGeneDict[cs_index]:
        grna_index = 0
        while grna_temp_storage[grna_index][0] < gene[2]:
            g_tup = grna_temp_storage[grna_index]
            if gene[1] < g_tup[0]:
                for item in g_tup:
                    gene.append(item)
                    gene.append(is_istop(g_tup, gene[1]))
            if grna_index < len(grna_temp_storage)-1:
                grna_index += 1
            else:
                break


def is_istop(grna,atg_pos):
    # noncoding change?
    codons = ['CAG','CAA','CGA']
    locs = list()
    for codon in codons:
        if 4 < grna[1].find(codon) < 9:
            locs.append(grna[1].find(codon))
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


def generate_report():
    f = open("kfd_concise_data.csv", "w")
    for cs in ScaffGeneDict:
        f.write("SCAFFOLD: " + cs + "\n")
        for gene in ScaffGeneDict[cs]:
            f.write("GENE " + gene[0] + "," + str(gene[1]) + "," + str(gene[2]) + "," + gene[3] + "\n")
            for grna in gene[4:]:
                for item in grna:
                    f.write(item + ",")
                f.write("\n")
    f.close()








parse_gbff("/Users/brianmendoza/Desktop/Kfedtschenkoi_382_v1.1.gene_exons.gff3")

sort_grnas_by_genes("/Users/brianmendoza/Desktop/kfdspCas9.cspr")

generate_report()





