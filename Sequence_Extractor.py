from Bio import SeqIO

target = "K62.pdb"
for record in SeqIO.parse(target, "pdb-atom"):
    print(record.seq)
