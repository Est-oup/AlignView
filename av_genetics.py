"""NCBI genetic codes, vendored for offline use.

Source: https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
Retrieved: 2026-11-10; SHA256: 2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04
Codon order: T,C,A,G at each position. Internal translation uses ncbieaa;
context-dependent stop codons are not inferred without CDS context.
"""
from __future__ import annotations
from functools import lru_cache
from itertools import product

CODONS = tuple(a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG")
TABLES = {
    1: ('Standard', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '---M------**--*----M---------------M----------------------------'),
    2: ('Vertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG', '----------**--------------------MMMM----------**---M------------'),
    3: ('Yeast Mitochondrial', 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**----------------------MM---------------M------------'),
    4: ('Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate\n Mitochondrial; Mycoplasma; Spiroplasma', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--MM------**-------M------------MMMM---------------M------------'),
    5: ('Invertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG', '---M------**--------------------MMMM---------------M------------'),
    6: ('Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear', 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--------------*--------------------M----------------------------'),
    9: ('Echinoderm Mitochondrial; Flatworm Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG', '----------**-----------------------M---------------M------------'),
    10: ('Euplotid Nuclear', 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**-----------------------M----------------------------'),
    11: ('Bacterial, Archaeal and Plant Plastid', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '---M------**--*----M------------MMMM---------------M------------'),
    12: ('Alternative Yeast Nuclear', 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**--*----M---------------M----------------------------'),
    13: ('Ascidian Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG', '---M------**----------------------MM---------------M------------'),
    14: ('Alternative Flatworm Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG', '-----------*-----------------------M----------------------------'),
    15: ('Blepharisma Macronuclear', 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------*---*--------------------M----------------------------'),
    16: ('Chlorophycean Mitochondrial', 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------*---*--------------------M----------------------------'),
    21: ('Trematode Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG', '----------**-----------------------M---------------M------------'),
    22: ('Scenedesmus obliquus Mitochondrial', 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '------*---*---*--------------------M----------------------------'),
    23: ('Thraustochytrium Mitochondrial', 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--*-------**--*-----------------M--M---------------M------------'),
    24: ('Rhabdopleuridae Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG', '---M------**-------M---------------M---------------M------------'),
    25: ('Candidate Division SR1 and Gracilibacteria', 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '---M------**-----------------------M---------------M------------'),
    26: ('Pachysolen tannophilus Nuclear', 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**--*----M---------------M----------------------------'),
    27: ('Karyorelict Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--------------*--------------------M----------------------------'),
    28: ('Condylostoma Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**--*--------------------M----------------------------'),
    29: ('Mesodinium Nuclear', 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--------------*--------------------M----------------------------'),
    30: ('Peritrich Nuclear', 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '--------------*--------------------M----------------------------'),
    31: ('Blastocrithidia Nuclear', 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '----------**-----------------------M----------------------------'),
    32: ('Balanophoraceae Plastid', 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '---M------*---*----M------------MMMM---------------M------------'),
    33: ('Cephalodiscidae Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG', '---M-------*-------M---------------M---------------M------------'),
}

TABLES={identifier:(' '.join(value[0].split()),value[1],value[2]) for identifier,value in TABLES.items()}


def code_table(identifier=1):
    if identifier not in TABLES: raise ValueError(f"Table NCBI inconnue: {identifier}")
    return dict(zip(CODONS, TABLES[identifier][1]))

def start_codons(identifier=1):
    code_table(identifier)
    return {codon for codon, marker in zip(CODONS,TABLES[identifier][2]) if marker=="M"}

def stop_codons(identifier=1):
    code_table(identifier)
    return {codon for codon, marker in zip(CODONS,TABLES[identifier][2]) if marker=="*"}

def context_dependent(identifier=1):
    table=code_table(identifier)
    return any(table[codon]!="*" for codon in stop_codons(identifier))

@lru_cache(maxsize=8192)
def translate_codon(codon,identifier=1):
    table=code_table(identifier)
    if codon in table:return table[codon]
    alphabet={"R":"AG","Y":"CT","S":"CG","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
    possibilities={table.get("".join(parts),"X") for parts in product(*(alphabet.get(c,c) for c in codon))}
    return next(iter(possibilities)) if len(possibilities)==1 else "X"
