from Bio import Align

seq_1 = "AAGAAATTCCAAGTCCAGGGATACACAAACAGGTGTACAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATTTAATCTATACCCTAGGAAAATGCCAAAGTCACAATTGGGTGGATTGGGTGATTTTCCAGTAGAAAGAAAATTCCATCCCAT"
seq_2 = "AAGAAATTCACAGTCCAGGGATACACAAACAGGTGTACAAAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATTTAATCTATACCCTAGGAAAATGCCAAAGTCGGAATTGGGTGATTGGGTGTTTCCAGTAGAAAGAAAATTCCATCCCAT"

aligner = Align.PairwiseAligner()
alignments = aligner.align(seq_1, seq_2)

# 7. Znajdź najdłuższy fragment pasujących par nukleotydów
# spośród wszystkich możliwych przyrównania sekwencji

matches = []

for alignment in alignments:
    cords = alignment.coordinates
    for i in range(0, len(cords[0])-1):
        f1 = seq_1[cords[0][i]:cords[0][i+1]]
        f2 = seq_2[cords[1][i]:cords[1][i+1]]
        if f1 == f2:
            matches.append({'fragment': f1, 'alignment': alignment, 'cords': cords})
        else:
            print("para cords sie nie zgrała")

longest_match = max(matches, key=lambda match: len(match['fragment']))
print(longest_match['alignment'])
print("Najdłuższy fragment pasujących par nukleotydów:", longest_match['fragment'])
print("Długość fragmentu:", len(longest_match['fragment']))
print(longest_match['cords'])


