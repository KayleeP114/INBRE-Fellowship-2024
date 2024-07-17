from Bio import SeqIO

target = "'2ztb[1].cif' "
for record in SeqIO.parse(target, "pdb-atom"):
    print(record.seq)
