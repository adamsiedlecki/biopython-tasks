from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
my_dna_seq = Seq("AAGAAATTCCAAGTCCAGGGATACACAAACAGGTGTACAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATATT")

# Policz wystąpienia każdego nukleotydu
print("A: " + str(my_dna_seq.count("A")))
print("T: " + str(my_dna_seq.count("T")))
print("C: " + str(my_dna_seq.count("C")))
print("G: " + str(my_dna_seq.count("G")))

# zawartość GC
gc_percent = "%.2f" % (gc_fraction(my_dna_seq) * 100)
print(f'gc_percent: {gc_percent}')

# Dokonaj transkrypcji sekwencji do sekwencji RNA.
messenger_rna = my_dna_seq.transcribe()
print(f'RNA: {messenger_rna}')

#Przetłumacz sekwencję DNA na sekwencję białkową.
protein_seq = my_dna_seq.translate()
print(f'protein_seq: {protein_seq}')

#Sekwencja antyrównoległa:
my_dna_reverse_compl = my_dna_seq.reverse_complement()
print(f'my_dna_reverse_compl: {my_dna_reverse_compl}')

#zapis do pliku
with open("sequence_analysis.txt", "w", encoding="utf-8") as f:
  print(f"Oryginalna sekwencja DNA: {my_dna_seq}", file=f)
  print("Liczba nukleotydów:")
  print(f"  A: {my_dna_seq.count('A')}", file=f)
  print(f"  T: {my_dna_seq.count('T')}", file=f)
  print(f"  C: {my_dna_seq.count('C')}", file=f)
  print(f"  G: {my_dna_seq.count('G')}", file=f)
  print(f"Zawartość GC: {gc_percent}%", file=f)
  print(f"Transkrybowany RNA: {messenger_rna}", file=f)
  print(f"Translowane białko: {protein_seq}", file=f)
  print(f"Odwrotne dopełnienie:: {my_dna_reverse_compl}", file=f)