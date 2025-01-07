from primer3 import bindings

seq = "GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT\
AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA\
ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG\
CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG\
TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA\
TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA\
ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA\
TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA\
GATGTTTCCCTCTAGTAG"

# przykład z https://libnano.github.io/primer3-py/quickstart.html#workflow
primers = bindings.design_primers(
    seq_args={
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_INCLUDED_REGION': [70, 200]
    },
    global_args={
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 29,
        'PRIMER_OPT_TM': 57.0,
        'PRIMER_MIN_TM': 50.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [
            [70, 300]
        ],
    })

print(primers)

def printValueToFile(val, file):
    print(val, end="", file=file)


with open("PCR_primers.txt", "w", encoding="utf-8") as f:
    printValueToFile("Left primer Right primer Start End", file=f)
    print("", file=f)
    for i in range(0, len(primers["PRIMER_LEFT"])):
        start = primers["PRIMER_LEFT"][i]["COORDS"][0] # początek lewego primera
        end = primers["PRIMER_RIGHT"][i]["COORDS"][0] + primers["PRIMER_RIGHT"][i]["COORDS"][1] # koniec prawego primera
        if start < 100 and end > 250: # zgodnie z treścią zadania - "sekwencji od 100 do 250 nukleotydu"
            printValueToFile(f'{primers["PRIMER_LEFT"][i]["SEQUENCE"]} -- {primers["PRIMER_RIGHT"][i]["SEQUENCE"]} {start} {end}', f)
            print("", file=f)
