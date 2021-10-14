from Bio.Seq import Seq
from Bio import SeqIO

file_list = ["bacterial1", "bacterial2", "bacterial3", "bacterial4", "mamalian1", "mamalian2", "mamalian3", "mamalian4"]
codon_list = ["TAA", "TAG", "TGA"]
start_list = ["ATG"]
pos_codons = "MFLIVSPTAYHQNKDECWRG"
pos_dicodons = "MFLIVSPTAYHQNKDECWRG"

# Functions
def find_orf_frames(seq):
    list = []
    for frame in range(0, 3):
        j = frame
        while j < len(seq):
            possible_codon = seq[j:j + 3]
            if possible_codon in start_list:
                for js in range(j + 3, len(seq), 3):
                    possible_codon = seq[js:js + 3]
                    if possible_codon in codon_list:
                        if len(seq[j:js + 3]) > 99:
                            list.append(Seq(seq[j:js]))
                        j = js + 3
                        break

            j = j + 3
    return list

def count(seq):
    newlist = [(x, 0) for x in pos_codons]
    size = 0
    for protein in seq:
        protein = protein.translate()
        for j in protein:
            size += 1
            for i, (letter, count) in enumerate(newlist):
                if (letter == j):
                    newlist[i] = (letter, count + 1)
                    break
    for i, (letter, count) in enumerate(newlist):
        newlist[i] = (letter, (count / size) * 100)
    return newlist

def di_count(seq):
    newlist = [(x + y, 0) for x in pos_dicodons for y in pos_codons]
    size = 0
    for protein in seq:
        protein = protein.translate()
        for j in range(0, len(protein), 2):
            size += 1
            dicodon = protein[j:j + 2]
            for i, (letter, count) in enumerate(newlist):
                if (letter == dicodon):
                    newlist[i] = (letter, count + 1)
                    break
    for i, (letter, count) in enumerate(newlist):
        newlist[i] = (letter, (count / size) * 100)
    return newlist


def compare(current, compare):
    i = 0
    for (x1, c1), (x2, c2) in zip(current, compare):
        i = i + abs(c1 - c2)
    return i / len(current)


# main calculations
codon_frequency = []
dicodon_frequency = []
for file in file_list:
    data = SeqIO.read("data/" + file + ".fasta", "fasta")
    sequence = data.seq
    orfs = find_orf_frames(sequence) + find_orf_frames(sequence.complement()[::-1])
    orf = open("orfs/" + file + "Orf.txt", "w")
    for protein in orfs:
        orf.write(str(protein) + '\n')
    codon_frequency.append(count(orfs))
    dicodon_frequency.append(di_count(orfs))

# writing results
f = open("codonFrequency.txt", "w")
for temp, file in zip(codon_frequency, file_list):
    f.write(file + "\n")
    for (codon, freq) in temp:
        f.write("%s %f \n" % (codon, freq))
    f.write("\n")
f.close()

f = open("dicodonFrequency.txt", "w")
for temp, file in zip(dicodon_frequency, file_list):
    f.write(file + "\n")
    for (dicodon, freq) in temp:
        f.write("%s %f \n" % (dicodon, freq))
    f.write("\n")
f.close()

r = open("codonDistanceMatrix.txt", "w")
for row, name in zip(codon_frequency, file_list):
    r.write("%s " % name)
    list = []
    for row2 in codon_frequency:
        list.append(compare(row, row2))
    for num in list:
        r.write("%1.4f " % num)
    r.write("\n")
r.close()

r = open("dicodonDistanceMatrix.txt", "w")
for row, name in zip(dicodon_frequency, file_list):
    r.write("%s " % name)
    list = []
    for row2 in dicodon_frequency:
        list.append(compare(row, row2))
    for num in list:
        r.write("%1.4f " % num)
    r.write("\n")
r.close()
