# Glossary

| Term | Meaning in AlignView |
| --- | --- |
| FASTA | Headers beginning with >identifier followed by sequence lines; may contain multiple records |
| Query | Active sequence compared against references; also the target of interactive sequence analyses |
| Reference | A comparison sequence, not necessarily a validated taxonomic label |
| Hit | A retained reference; it may contain several local HSPs |
| HSP | One local aligned segment with coordinates, identity, scores and possibly a translated frame |
| Identity | Matching characters within the relevant alignment definition; local HSP identity is not whole-sequence identity |
| Coverage | Fraction of a sequence touched by the union of retained local segments; VSEARCH uses its own global-alignment coverage fields |
| E-value | Expected random matches of comparable score within the engine's search space, not a taxonomic probability |
| Bitscore | Engine-reported normalized local score; VSEARCH nucleotide searches do not calculate it |
| MSA | Multiple sequence alignment with comparable columns and gaps |
| Gap | Alignment placeholder, usually -, for an absent residue in a column |
| Consensus | A column-wise summary, not necessarily an observed biological sequence |
| Entropy | Diversity of observed characters; can reflect variation, ambiguity or heterogeneous alignment |
| p-distance | Pairwise difference fraction on comparable positions, without evolutionary correction |
| nt / aa | Nucleotide / amino acid units |
| Reading frame | Translation offset and orientation: +1,+2,+3,-1,-2,-3 |
| Genetic code | NCBI codon table selected for the query or references; not automatically inferred |
| ORF | Here, an exploratory start-to-stop segment under the selected code and start policy |
| K-mer | Exact word of length k |
| Newick | Parenthesized tree format with names, optional branch lengths and internal labels |
| SH-like support | Local FastTree support; not a tree bootstrap |
| Phylogram / cladogram | Branch-length-scaled tree / topology drawing without a distance scale |
| Dereplication | Merge identical input sequences while retaining counts |
| Cluster | Algorithmic similarity group; not automatically a biological species |
| Chimera | Candidate mosaic sequence screened by an algorithm such as UCHIME |
| Provenance | Inputs, hashes, settings, versions, commands and logs describing a result |
| SHA256 | File fingerprint used to detect changes, not an author's signature |
| Calculation profile | Reusable JSON parameters, without sequence data |
| Search preset / export profile | Threshold choices / selection of output formats |

A 100-nt HSP at 99% identity covers only about 10% of a 1000-nt query. Two HSPs covering 1–100 and 81–180 cover 180 distinct positions, not 200. The inclusive interval 3–62 becomes BED [2,62). BLASTX can align 60 query nucleotides to 20 reference amino acids, so the axes must not be treated as the same unit.
