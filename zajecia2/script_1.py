from Bio import Align

seq_1 = "AAAACCCCTTGGTTTT"
seq_2 = "AAACACCCTTTGTTTA"

aligner = Align.PairwiseAligner()
alignments = aligner.align(seq_1, seq_2)

for alignment in alignments:
    # Wyznacz wartość dopasowania pomiędzy sekwencjami.
    print(alignment.score)

for i in range(0, min(5, len(alignments))):
    # Wypisz kilka pierwszych przyrównań.
    print(alignments[i])

for i in range(0, 3):
    # Dla pierwszych trzech przyrównań wypisz ich macierze koordynat
    print(f'{i + 1}. Macierz koordynat')
    print(alignments[i].coordinates)