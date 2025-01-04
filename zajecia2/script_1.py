from Bio import Align

seq_1 = "AAAACCCCTTGGTTTT"
seq_2 = "AAACACCCTTTGTTTA"

aligner = Align.PairwiseAligner()
alignments = aligner.align(seq_1, seq_2)

# Wyznacz wartość dopasowania pomiędzy sekwencjami.
print(f'Score: {alignments[0].score}')

for i in range(0, min(3, len(alignments))):
    # Wypisz kilka pierwszych przyrównań.
    print(alignments[i])

for i in range(0, 3):
    # Dla pierwszych trzech przyrównań wypisz ich macierze koordynat
    print(f'{i + 1}. Macierz koordynat')
    print(alignments[i].coordinates)

# Wyjaśnij, jakie informacje można odczytać z macierzy koordynat
# (np. fragmenty pasujących regionów)
# (zapisz jako komentarz w kodzie)

# Z macierzy koordynat można odczytać dopasowania dwóch sekwencji.
# Macierz przyrównania dwóch sekwencji ma 2 wiersze o równej długości.
# W jednym wierszu są indeksy pierwszej sekwencji, a w drugim wierszu odpowiadające im indeksy drugiej sekwencji.
# W ten sposób można odczytać pasujące fragmenty oraz przerwy (luki).
# Luki można zidentyfikować po skokach indeksów.

# 5d Zaimplementuj funkcję liczącą przyrównania, które zawierają co najmniej jedno niedopasowanie

def policz_przyrownania_z_niedopasowaniami(alignments):
    niedopasowania = 0

    for alignment in alignments:
        if (alignment.counts().mismatches > 0):
            niedopasowania += 1

    return niedopasowania

print(f'Przyrownania z niedopasowaniami: {policz_przyrownania_z_niedopasowaniami(alignments)} / {len(alignments)}')