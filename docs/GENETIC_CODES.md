# Choose a translation table

AlignView vendors the 27 active tables in the NCBI gc.prt source retrieved on 10 November 2026. Table 1 (Standard) remains the default. Choose the code from the biological annotation, organism and compartment; it is not inferred automatically.

## Where to choose

Use **Advanced → Additional analyses → Sequences, motifs and translation → Translation and coding regions → Choose the translation table**, or the parameters screen before a translated calculation. The list shows each NCBI number and name.

- `--genetic-code ID` selects the query/internal translation code.
- `--reference-genetic-code ID` selects the reference code; unset means inherit the query code.
- `--orf-starts atg` keeps ATG-only initiation; `--orf-starts table` uses the NCBI-listed initiators.

```bash
alignview --list-genetic-codes
alignview query.fna refs.faa --mode blastx --genetic-code 4 --save --out code04_01 --no-menu
alignview query.faa refs.fna --mode tblastn --reference-genetic-code 11 --save --out code11_01 --no-menu
alignview --write-profile bacterial.json --genetic-code 11 --orf-starts table
```

Selections are recorded in profiles, projects, commands and successful-calculation provenance. BLAST and DIAMOND receive the applicable native query/database genetic-code flags. An engine version that does not support a requested code may reject it in its stderr log; AlignView does not silently substitute the standard table.

## What changes

TGA is a stop in table 1 and tryptophan W in tables 2 or 4. ATA is isoleucine in table 1 and methionine in table 2. Six-frame translation, translated MSA preparation, translated Needle/Water comparisons and translated dotplots use the selected code for each side.

An ambiguous IUPAC codon resolves to a residue only when all possible codons agree; otherwise it becomes X. Incomplete trailing codons are ignored. Six-frame exports retain stop symbols; MSA/dotplot preparation replaces stops with X as announced.

Fragment translation does not arbitrarily rewrite its first amino acid as methionine. For an ORF detected with table-listed starts, the alternative initiator is represented by M at the start of the exported protein. The detector remains a simple start-to-stop scan, without intron, RNA-editing or other contextual models.

Tables 27, 28 and 31 have context-dependent stops. Context-free translation uses their ncbieaa amino acid; the simple ORF detector explicitly refuses them because CDS context is unavailable. This does not prevent displaying or translating fragments with those tables.

Changing a code prepares subsequent calculations; it does not reinterpret existing results. Rerun the relevant action and rebuild HTML. Read historical parameters to identify the code actually used.

## Available tables

| NCBI ID | Name | Context-dependent stops |
| --- | --- | --- |
| 1 | Standard | no |
| 2 | Vertebrate Mitochondrial | no |
| 3 | Yeast Mitochondrial | no |
| 4 | Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma | no |
| 5 | Invertebrate Mitochondrial | no |
| 6 | Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear | no |
| 9 | Echinoderm Mitochondrial; Flatworm Mitochondrial | no |
| 10 | Euplotid Nuclear | no |
| 11 | Bacterial, Archaeal and Plant Plastid | no |
| 12 | Alternative Yeast Nuclear | no |
| 13 | Ascidian Mitochondrial | no |
| 14 | Alternative Flatworm Mitochondrial | no |
| 15 | Blepharisma Macronuclear | no |
| 16 | Chlorophycean Mitochondrial | no |
| 21 | Trematode Mitochondrial | no |
| 22 | Scenedesmus obliquus Mitochondrial | no |
| 23 | Thraustochytrium Mitochondrial | no |
| 24 | Rhabdopleuridae Mitochondrial | no |
| 25 | Candidate Division SR1 and Gracilibacteria | no |
| 26 | Pachysolen tannophilus Nuclear | no |
| 27 | Karyorelict Nuclear | yes |
| 28 | Condylostoma Nuclear | yes |
| 29 | Mesodinium Nuclear | no |
| 30 | Peritrich Nuclear | no |
| 31 | Blastocrithidia Nuclear | yes |
| 32 | Balanophoraceae Plastid | no |
| 33 | Cephalodiscidae Mitochondrial | no |

Sources: [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), [gc.prt data](https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt). av_genetics.py records the source, retrieval date and SHA256.

Internal translation and ORF detection remove alignment gaps first. Their coordinates refer to the prepared sequence rather than MSA columns; retain source inputs and preparation details to interpret those positions.
