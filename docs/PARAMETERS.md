# Understand and adjust scientific parameters

Before each calculation, use **Calculation parameters / Engine parameters**. F1 explains an option. **Inspect commands** uses the same argument builders as actual execution, without running the analysis. The complete [CLI reference](CLI.md) lists every flag and default.

## Calculation versus result filtering

Identity, coverage, bitscore and the retained-reference limit filter calculated local-search results. E-value, seeds, masks, scores, sensitivity and search-candidate limits change the calculation and require a rerun. A filter that changes the selected references can invalidate their derived MSA/tree.

The maximum retained references starts at 100. Standard/More permissive/Stricter presets use E-values 1e-5/1e-2/1e-10 and adjust DIAMOND sensitivity where applicable. These presets do not change the retained-reference maximum or MMseqs2 sensitivity. The automatic candidate limit is max(100, 10 × retained references), separately adjustable.

## Engine-specific settings

| Engine | Exposed controls | Interpretation |
| --- | --- | --- |
| BLAST | Task, seed size, masks, strand, scores, matrix, gaps, composition, HSP/candidate limits | Valid values and compatible score combinations depend on the BLAST task |
| DIAMOND | Sensitivity, matrix, gaps, masks, composition, HSP count, BLASTX minimum ORF | Separate from BLAST and MMseqs2 settings |
| MMseqs2 | Numeric sensitivity, masking, composition bias, candidates, explicit traceback | DIAMOND sensitivity does not configure this engine |
| MAFFT | auto/global/local/fast, gap opening, offset, refinements | Explicit refinements require global/local/fast rather than auto |
| MUSCLE 5 | align/super5, perturbation seed, guide-tree permutation, refinements | One selected run produces one MSA; no implicit ensemble |
| Clustal Omega | Iterations and full distance matrices | Full iterative matrices require iterations and may cost more |
| Needle/Water | Matrix and internal gap costs | These costs do not configure search or MSA engines |
| Needle | Free ends and separate terminal gap costs | AlignView penalizes ends by default; native EMBOSS defaults differ |
| FastTree | Alphabet-appropriate model, Gamma20, local support, seed | Gamma20 rescales after CAT optimization; local support is not a bootstrap |
| VSEARCH | Identity definition, strands, coverage, masks, heuristics, abundance and length limits, UCHIME thresholds | Global nucleotide metrics; no BLAST E-value or bitscore |

Unset/native means the option is not sent and the installed engine owns the default. AlignView choices are explicit: MAFFT auto; global/local strategies initially use 1000 refinements; DIAMOND sensitive; DIAMOND BLASTX minimum ORF 1; 4 threads where supported. Standard FastTree is single-threaded.

## Translation

`--genetic-code` selects the query/internal NCBI table; `--reference-genetic-code` selects the reference table and otherwise inherits the query table. `--frame` and `--sframe` choose internal complete-translation frames; they do not limit the frames searched by BLASTX. Tables and initiators are covered in [Genetic codes](GENETIC_CODES.md).

## Reusable profiles and resuming

```bash
alignview --write-profile tuned.json --mafft local --mafft-op 2 --mafft-maxiterate 100
alignview sequences.fasta --load-profile tuned.json --save --out run_01 --no-menu
alignview --resume run_01 --mafft-op 3 --msa align --save --out run_02 --no-menu
```

A JSON calculation profile contains engine settings and result filters, not sequence data or executable paths. Explicit flags override profile values. Invalid fields, types and ranges are rejected. Profile output never overwrites an existing file.

Resuming alone does not rerun engines. Explicit settings become future settings; explicit action flags request new calculations. Saved calculation blocks preserve the original parameters of still-existing results. The export profile chooses files and is distinct from a search threshold preset.

## A useful comparison procedure

Check inputs and alphabets first. Save a baseline and its provenance. Change one parameter or coherent group, inspect the command, rerun the relevant calculation, compare coverage/regions/gaps/MSA quality, and save a new project. More Matches alone do not establish better biological validity. Engine stderr logs remain the authority for unsupported matrix/cost combinations.
