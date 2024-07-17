from Bio import SeqIO

target = "2ztb.pdb"
for record in SeqIO.parse(target, "pdb-atom"):
    print(record.seq)
