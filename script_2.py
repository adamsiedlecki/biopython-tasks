from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import nt_search


# zad3
def count_atcg(seq_record):
    return {"A": seq_record.seq.count("A"),
            "T": seq_record.seq.count("T"),
            "C": seq_record.seq.count("C"),
            "G": seq_record.seq.count("G")}


def get_gc_content(seq_record):
    return gc_fraction(seq_record) * 100


def rev_complement(seq_record):
    return seq_record.reverse_complement().seq


#zad 4
def find_motifs(seq_record):
    motifs = {"ATG": [], "TATA": [], "GAATTC": []}
    for motif in motifs:
        # https://biopython.org/docs/latest/api/Bio.SeqUtils.html#Bio.SeqUtils.nt_search
        result = nt_search(str(seq_record.seq), motif)  # tą metodę podpowiedział ChatGpt
        motifs[motif] = result[1::] # pierwsza pozycja na liście to string - motif

    return str(motifs)




# zad5
def seq_translation(record):
    min_pro_len = 1
    result = []
    counter = 1
    # kod pochodzi z Identifying open reading frames
    # https://biopython.org/docs/1.84/Tutorial/chapter_cookbook.html
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record) - frame) // 3)  # Multiple of three
            for pro in nuc[frame: frame + length].translate().split("*"):
                if len(pro) >= min_pro_len:
                    # print(
                    #     "%s...%s - length %i, strand %i, frame %i"
                    #
                    #     % (pro[:30], pro[-3:], len(pro), strand, frame)
                    # )
                    if strand == 1:
                        result.append(f'Frame{counter}: {len(pro)}')
                    else:
                        result.append(f'RevFrame{counter}: {len(pro)}')
                    counter += 1
    #print(str(result))
    return str(result)


def printValueToFile(val, file):
    print(val, end="", file=file)


#zapis do pliku CSV - zad 6
with open("sequence_analysis.csv", "w", encoding="utf-8") as f:
    print("SeqID,Nucleotide_Counts,GC_Content,Motif_Positions,Reverse_Complement,Translation_Lengths", file=f)

    # zad2
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        printValueToFile(seq_record.id, f)
        printValueToFile(", ", f)
        printValueToFile(count_atcg(seq_record), f)
        printValueToFile(", ", f)
        printValueToFile(get_gc_content(seq_record), f)
        printValueToFile(", ", f)
        printValueToFile(find_motifs(seq_record), f)
        printValueToFile(", ", f)
        printValueToFile(rev_complement(seq_record), f)
        printValueToFile(", ", f)
        printValueToFile(seq_translation(seq_record), f)
        print("", file=f) # nowa linia

