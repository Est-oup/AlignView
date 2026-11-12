# Command-line reference

Generated from the actual parser by `python3 scripts/generate_cli_docs.py`. Do not edit the tables manually.

Use the installed `alignview` command, or replace it with `/path/to/AlignView/alignview.sh`.

## Common workflows

```bash
alignview --demo
alignview --docs quickstart
alignview aligned.fasta --msa view --html-out report.html --no-menu
alignview sequences.fasta --msa align --report --save --out project_01 --no-menu
alignview queries.fna refs.fna --vsearch search --query-id chosen_id --report --save --out search_01 --no-menu
alignview --resume project_01 --html-out updated.html --no-menu
```

## Storage and resuming

- Relative paths start from the launch directory; the launcher does not change it.
- `--save --out DIRECTORY` retains a project. `--out` alone does not save.
- `--html-out FILE.html` retains standalone HTML. `--export-dir DIRECTORY` retains exports of the selected profile.
- `--report` alone prepares temporary HTML. Add a destination or saving to retain it.
- Batch runs refuse an existing destination.
- Resume alone reruns no scientific engine. Only explicit actions (`--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif`, `--vsearch`) request calculations.
- Explicit parameters override profiles/saved future settings; historical settings remain associated with existing results.
- `--genetic-code` selects the query/internal table; `--reference-genetic-code` selects or inherits the reference table.
See [Usage](USAGE.md), [Parameters](PARAMETERS.md), [VSEARCH](VSEARCH.md), [Genetic codes](GENETIC_CODES.md) and [Files](FILES.md).

## General help

| Option | Default | Purpose and values |
| --- | --- | --- |
| `-h / --help` | — | show this help message and exit |

## Input and resuming

| Option | Default | Purpose and values |
| --- | --- | --- |
| `query` | — | Query FASTA, or a single FASTA to explore |
| `--query-id` | '' | FASTA identifier to use as the query when the first file contains several sequences |
| `references` | — | Reference FASTA |
| `--mode` | — | Search or alignment method Values: msa, vsearch, blastn, blastp, blastx, tblastn, tblastx, diamond-blastp, diamond-blastx, mmseqs-blastp, global-dna, global-aa. |
| `--resume` | — | Resume a saved project without rerunning its calculations |
| `--demo` | False | Explore a bundled example without external tools |

## Search and filtering

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--task` | '' | NCBI task compatible with the selected program |
| `--evalue` | 1e-05 | Maximum E-value (default: 1e-5) |
| `--identity` | 0 | Minimum local identity % (0) |
| `--coverage` | 0 | Query coverage from the union of HSPs % (0) |
| `--bitscore` | 0 | Minimum local bitscore (0) |
| `--hits` | 100 | Maximum retained references (default: 100) |
| `--preset` | 'standard' | Search preset: standard, sensible or strict; explicit options take priority Values: standard, sensitive, strict. |
| `--threads` | 4 | Number of threads (4) |
| `--sensitivity` | 'sensitive' |  Values: fast, sensitive, more-sensitive, very-sensitive, ultra-sensitive. |

## Alignments and analyses

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--msa` | 'none' |  Values: none, view, align, full, hits, all, global, fragments, hsp. |
| `--mafft` | 'auto' |  Values: auto, global, local, fast. |
| `--aligner` | 'mafft' | MSA engine: MAFFT, MUSCLE 5, Clustal Omega Values: mafft, muscle, clustalo. |
| `--tree` | False | Infer an exploratory FastTree tree from the MSA |
| `--pairwise` | 'none' |  Values: none, global, local. |
| `--ref` | 1 | Reference index (1) |
| `--frame` | — | Select a translation frame Values: -3, -2, -1, 1, 2, 3. |
| `--gapopen` | 10 |  |
| `--gapextend` | 0.5 |  |
| `--matrix` | '' | EMBOSS substitution matrix |
| `--free-ends` | False | Do not penalize terminal gaps in Needle |
| `--dotplot` | False | Generate a query / best-reference dotplot |
| `--kmer` | 9 | Dotplot word length (9) |
| `--motif` | '' | IUPAC motif in the query |
| `--vsearch` | 'none' | VSEARCH nucleotide analysis: search, dereplication, clustering or chimeras Values: none, search, dereplicate, cluster, chimeras. |

## Reports, exports and saving

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--report` | False | Prepare HTML and SVG figures; retain them using --save, --export-dir or --html-out |
| `--html-out` | '' | Save a standalone HTML here, even without --save; refuses to overwrite |
| `--no-open` | False | Do not automatically open HTML in the browser |
| `--out` | '' | Project directory with --save; relative to the launch directory, created at the end |
| `--save` | False | Explicitly save at the end, including batch runs |
| `--export-profile` | 'standard' | standard / publication / complete: files prepared when saving Values: standard, publication, complete. |
| `--export-dir` | '' | Explicitly export files to this directory, even without --save |

## Scientific engines and installation

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--executable` | [] | Exact external-tool path, repeatable, e.g. blastx=/opt/ncbi/bin/blastx |
| `--tools-dir` | [] | Add local executable directories to AlignView's PATH (repeatable) |
| `--module-policy` | 'auto' | auto: try HPC modules for missing commands; never: use existing paths and the managed environment Values: auto, never. |
| `--diagnose` | False | Report available executables and versions, then exit |
| `--setup` | — | Prepare engine installation using Conda/Mamba Values: alv, blast, clustalo, diamond, emboss, fasttree, full, mafft, minimal, mmseqs, muscle, vsearch. |
| `--install-tools` | False | Run the --setup plan in a dedicated environment |
| `--tools-prefix` | — | Dedicated Conda prefix (or ALIGNVIEW_TOOLS_PREFIX) |
| `--diagnose-json` | False | Machine-readable JSON diagnostics without loading modules |

## Interface and documentation

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--list-genetic-codes` | False | List available NCBI translation tables and exit |
| `--no-menu` | False | Without interactive menus |
| `--simple-ui` | False | Classic numbered menus |
| `--no-color` | False | Without ANSI colors |
| `--list-tools` | False | List features |
| `--docs` | — | Read an offline guide without opening data (default: index) Values: index, quickstart, usage, files, installation, parameters, html, exports, cli, troubleshooting, glossary, vsearch, translation, development, releasing, audit, getting-started, genetics. |
| `--version` | — | show program's version number and exit |

## Other scientific parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--genetic-code` | 1 | Query translation table (NCBI): Code used for internal query translation and translated query searches. Values: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| `--reference-genetic-code` | — | Reference translation table: Unset: use the query table; passed to TBLASTN/TBLASTX and internal reference translation. Values: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| `--orf-starts` | 'atg' | ORF start codons: atg: ATG only; table: all initiators from the selected NCBI table. Values: atg, table. |
| `--vsearch-id` | 0.97 | Minimum VSEARCH identity: Fraction from 0 to 1; distinct from the BLAST local identity filter. |
| `--vsearch-iddef` | 2 | VSEARCH identity definition: 0 to 4; 2 excludes terminal gaps. See the VSEARCH documentation. Values: 0, 1, 2, 3, 4. |
| `--vsearch-strand` | 'both' | VSEARCH strands: both: search and cluster both strands; de novo UCHIME uses the plus strand. Values: plus, both. |
| `--vsearch-query-cov` | 0.0 | Minimum VSEARCH query coverage: Fraction from 0 to 1, calculated by VSEARCH. |
| `--vsearch-target-cov` | 0.0 | Minimum VSEARCH reference coverage: Fraction from 0 to 1, calculated by VSEARCH. |
| `--vsearch-maxaccepts` | 0 | Accepted candidates before stopping: 0: no accepted-candidate limit; distinct from the reported reference limit. |
| `--vsearch-maxrejects` | 32 | Rejected candidates before stopping: 0: do not stop after rejections; also set maxaccepts and minwordmatches to 0 for exhaustive search. |
| `--vsearch-minwordmatches` | 10 | Minimum shared words for VSEARCH: 0 disables this prefilter; useful for short sequences, at greater computational cost. |
| `--vsearch-masking` | 'dust' | VSEARCH masking: dust: low complexity; none: no masking. Internal input sequences are normalized to uppercase. Values: none, dust. |
| `--vsearch-n-mismatch / --no-vsearch-n-mismatch` | False | Count ambiguous symbols as mismatches: Prevent N from counting as a match to every compatible nucleotide. |
| `--vsearch-sizein / --no-vsearch-sizein` | False | Read header abundances: Use existing ;size=N; annotations. Otherwise each input record counts as one sequence. |
| `--vsearch-min-length` | 1 | Minimum VSEARCH sequence length: Shorter sequences are excluded and counted in the summary. |
| `--vsearch-max-length` | 50000 | Maximum VSEARCH sequence length: Longer sequences are excluded and counted in the summary. |
| `--vsearch-abskew` | 2.0 | UCHIME abundance skew: Minimum parent abundance relative to the candidate chimera. |
| `--vsearch-minh` | 0.28 | Minimum UCHIME chimera score: Higher values make detection stricter, with fewer candidates detected. |
| `--search-candidates` | — | Engine candidate limit: Auto = max(100, 10 × retained references). This candidate limit can affect results and differs from the 100 displayed hits. |
| `--hsp-overlap` | 0.5 | Minimum HSP anchor overlap: Fraction of the anchor query region that a segment must cover in the exploratory MSA. |
| `--sframe` | — | Reference frame for Needle/Water: Unset: best HSP frame, otherwise +1. Applies to TBLASTN/TBLASTX. Values: -3, -2, -1, 1, 2, 3. |

## BLAST parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--blast-word-size` | — | BLAST seed size: Smaller seeds increase sensitivity and cost; allowed values depend on the task. |
| `--blast-reward` | — | DNA match score: Reward for an identical nucleotide; depends on mismatch and gap costs. |
| `--blast-penalty` | — | DNA mismatch score: Negative mismatch score; affects scoring and statistics. |
| `--blast-matrix` | — | BLAST protein matrix: Examples: BLOSUM62, BLOSUM45, PAM30. Gap costs must be compatible. |
| `--blast-gapopen` | — | BLAST gap opening: BLAST-specific cost, independent of Needle/Water. TBLASTX is ungapped. |
| `--blast-gapextend` | — | BLAST gap extension: BLAST-specific cost. Only supported matrix/cost combinations are accepted. |
| `--blast-masking` | 'native' | Masquage DUST / SEG: Mask low complexity. Disabling masking may reveal nonspecific matches. Values: native, yes, no. |
| `--blast-soft-masking` | 'native' | Masquage souple BLAST: true: mask seeds; false: also mask extension. The default depends on the program. Values: native, true, false. |
| `--blast-strand` | 'both' | DNA query strand: both: both strands; plus/minus: restrict the search to one strand. Values: both, plus, minus. |
| `--blast-max-hsps` | — | Maximum HSPs per reference: Limit reported segments and potentially the coverage calculated from them. |
| `--blast-comp-stats` | 'native' | BLAST composition statistics: 0: disabled; 1/2/3: NCBI corrections depending on the engine. Affects scores and E-values. Values: native, 0, 1, 2, 3. |

## DIAMOND parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--diamond-matrix` | — | DIAMOND matrix: Protein matrix; gap costs must be compatible. |
| `--diamond-gapopen` | — | DIAMOND gap opening: DIAMOND-specific cost, independent of other engines. |
| `--diamond-gapextend` | — | DIAMOND gap extension: Penalize gap extension. |
| `--diamond-masking` | 'native' | DIAMOND masking: 0: disabled; 1: tantan; seg: SEG (DIAMOND ≥ 2.0.12). Values: native, 0, 1, seg. |
| `--diamond-comp-stats` | 'native' | DIAMOND composition correction: 0: disabled; 1: composition correction; 3/4: matrix adjustments, depending on the version. Values: native, 0, 1, 3, 4. |
| `--diamond-max-hsps` | — | Maximum DIAMOND HSPs: 0: all; multiple HSPs may increase coverage and computational cost. |
| `--diamond-min-orf` | 1 | Minimum DIAMOND BLASTX ORF: AlignView uses 1 rather than DIAMOND's automatic choice. 0: native automatic behavior. |

## MMseqs2 parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--mmseqs-sensitivity` | — | MMseqs2 sensitivity (-s): From 1 (fast) to 7.5 (sensitive). Independent of DIAMOND sensitivity. |
| `--mmseqs-mask` | 'native' | MMseqs2 masking: 1: low complexity masked; 0: disabled. Values: native, 0, 1. |
| `--mmseqs-comp-bias` | 'native' | MMseqs2 composition correction: 0/1: disable/enable composition-bias correction. Values: native, 0, 1. |

## MAFFT parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--mafft-op` | — | MAFFT gap opening (--op): Group-to-group alignment gap penalty; independent of EMBOSS. |
| `--mafft-ep` | — | MAFFT offset (--ep): MAFFT offset parameter; different from EMBOSS gapextend. |
| `--mafft-maxiterate` | — | MAFFT refinement iterations: Global/local: AlignView defaults to 1000 iterations. Choose global, local or fast to set this number explicitly. |

## MUSCLE 5 parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--muscle-method` | 'align' | MUSCLE 5 algorithm: align: PPP; super5: intended for larger datasets. MUSCLE 3 is unsupported. Values: align, super5. |
| `--muscle-perturb` | 0 | MUSCLE perturbation seed: 0: standard HMM; another seed explores alignment variability. |
| `--muscle-perm` | 'none' | Guide-tree permutation: Change the guide-tree traversal. One permutation produces one MSA. Values: none, abc, acb, bca. |
| `--muscle-refineiters` | — | MUSCLE 5 refinement iterations: Refinement iteration count; cost and results may change. |

## Clustal Omega parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--clustalo-iterations` | 0 | Clustal Omega iterations: Combined guide-tree/HMM cycles. 0: no iterative refinement. |
| `--clustalo-full / --no-clustalo-full` | False | Full Clustal Omega distances: Full distance matrix instead of mBed; quadratic cost. |
| `--clustalo-full-iter / --no-clustalo-full-iter` | False | Full distances during iterations: Use a full matrix during refinement; requires iterations. |

## Needle parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--needle-endopen` | — | Needle terminal gap opening: Unset: same as AlignView gapopen; ignored when sequence ends are free. |
| `--needle-endextend` | — | Needle terminal gap extension: Unset: same as AlignView gapextend; ignored when sequence ends are free. |

## FastTree parameters

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--fasttree-model` | 'auto' | FastTree evolutionary model: auto: JC for DNA, JTT for proteins. GTR: DNA; JTT/WAG/LG: proteins. Values: auto, JC, GTR, JTT, WAG, LG. |
| `--fasttree-gamma / --no-fasttree-gamma` | False | Gamma20 reassessment: After CAT optimization, rescore and rescale branches under Gamma20; does not repeat topology inference. |
| `--fasttree-boot` | — | Local-support resampling: SH-like support, not a tree bootstrap. 0: disable local support. |
| `--fasttree-seed` | — | FastTree random seed: Set the random seed for stochastic calculations. |

## Reusable calculation profiles

| Option | Default | Purpose and values |
| --- | --- | --- |
| `--load-profile` | — | Load a calculation profile; explicit CLI options take priority |
| `--write-profile` | — | Write the calculation profile and exit without running an engine |
