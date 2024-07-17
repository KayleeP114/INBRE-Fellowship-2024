from Bio import SeqIO

target = "2ztb.cif"
for record in SeqIO.parse(target, "pdb-atom"):
    print(record.seq)
