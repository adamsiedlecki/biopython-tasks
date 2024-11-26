from Bio import SeqIO

# zad3
def count_atcg(seq_record):
    print("A: " + str(seq_record.seq.count("A")))
    print("T: " + str(seq_record.seq.count("T")))
    print("C: " + str(seq_record.seq.count("C")))
    print("G: " + str(seq_record.seq.count("G")))

# zad2
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    count_atcg(seq_record)



