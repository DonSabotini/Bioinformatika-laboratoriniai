from Bio.Seq import Seq
from Bio import SeqIO

from pathlib import Path


def determineFormat(path):
    txt = Path(path).read_text().splitlines()
    quality = txt[3:][::4]
    low = ord(quality[0][0]);
    high = 34;
    for lis in quality:
        qua = list(lis)
        for letter in qua:
            if (ord(letter) < low):
                low = ord(letter)
            if (ord(letter) > high):
                high = ord(letter)
    print("ascii range [%d , %d]" % (low, high))
    print("formats:")
    Sanger_Phred33 = range(32, 74)
    Solexa_Solexa64 = range(58, 105)
    Illumina13_Phred64 = range(63, 105)
    Illumina15_Phred64 = range(66, 106)
    Illumina18_Phred33 = range(32, 75)
    if low in Sanger_Phred33 and high in Sanger_Phred33:
        print("Sanger Phred+33")
    if low in Solexa_Solexa64 and high in Solexa_Solexa64:
        print("Solexa Solexa+64")
    if low in Illumina13_Phred64 and high in Illumina13_Phred64:
        print("Illumina1.3 Phred+64")
    if low in Illumina15_Phred64 and high in Illumina15_Phred64:
        print("Illumina1.5 Phred+64")
    if low in Illumina18_Phred33 and high in Illumina18_Phred33:
        print("Illumina1.8 Phred+33")


def calculateCG(sequence):
    f = open("result.txt", "w")
    occurList = []
    idOccurList = []
    for record in sequence.records:
        cOccur = ((record.seq.count('C') + record.seq.count('G')) / len(record.seq)) * 100
        occurList.append((cOccur, 0))
        idOccurList.append((record.id, cOccur))
    occurList = list(set(occurList))
    occurList = [(x, y, []) for (x, y) in occurList]

    for (id, occur) in idOccurList:
        for i, (bOccur, count, idlist) in enumerate(occurList):
            if (bOccur == occur):
                idlist.append(id)
                occurList[i] = (bOccur, count + 1, idlist)
                break

    for (occur, count, idList) in occurList:
        f.write("%f;    %d;       %s\n" % (occur, count, str(idList)))
    f.close()


if __name__ == '__main__':
    sequence = SeqIO.parse('reads_for_analysis.fastq', 'fastq')
    determineFormat('reads_for_analysis.fastq')
    calculateCG(sequence)
