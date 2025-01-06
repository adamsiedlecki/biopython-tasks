import primer3
from primer3 import calc_tm

primery = []
with open('primers.txt', 'r', encoding='utf-8') as file:
    lines = file.readlines()
    primery = [line.strip() for line in lines]


def printValueToFile(val, file):
    print(val, end="", file=file)


with open("primer_analysis.csv", "w", encoding="utf-8") as f:
    print("Primer,Tm,Homodimer,Hairpin", file=f)

    # docs
    # https://libnano.github.io/primer3-py/api/bindings.html#primer3.bindings.calc_homodimer
    # https://libnano.github.io/primer3-py/api/bindings.html#primer3.bindings.calc_hairpin
    for primer in primery:
        printValueToFile(primer, f)
        printValueToFile(",", f)
        printValueToFile(f'{calc_tm(primer):.2f}', f)
        printValueToFile(",", f)
        printValueToFile(f'{"Tak" if primer3.calc_homodimer(primer).structure_found  else "Nie"}', f)
        printValueToFile(",", f)
        printValueToFile(f'{"Tak" if primer3.calc_hairpin(primer).structure_found else "Nie"}', f)
        print("", file=f)  # nowa linia

