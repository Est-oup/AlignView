# Calculation and settings review

Reviewed on 10 November 2026, with later workflow and English-only updates. This is a functional integration review, not a biological accuracy ranking of the supported tools.

## Calculation controls

| Earlier problem | Current behavior |
| --- | --- |
| Few engine-native options exposed | A typed registry shares settings across CLI, menus, profiles and provenance |
| Generic sensitivity appeared to apply to every engine | Separate DIAMOND sensitivity and numeric MMseqs2 sensitivity |
| Hidden candidate limit varied by engine | Visible, configurable candidate limit with engine-specific semantics |
| MMseqs2 alignments lacked an explicit traceback request | `-a` requests the gapped alignments consumed by AlignView |
| EMBOSS gap settings could be mistaken for BLAST/MSA settings | Separate engine-specific option names and labels |
| MAFFT refinement count was fixed | Strategy, gap controls and refinements are exposed; conflicting auto/refinement settings are rejected |
| HSP alignment ignored the requested strategy and fixed anchor overlap | Requested strategy is retained; anchor overlap is configurable and preparation is inspectable |
| Fragment addition displayed an ambiguous engine | MAFFT is explicitly required; reference alignment and fragment addition remain distinct steps |
| DIAMOND BLASTX fixed `--min-orf 1` invisibly | The setting is explicit and configurable; 0 restores native automatic behavior |
| Needle end penalties were unclear | Separate free-end selection and terminal costs; Water does not expose irrelevant end settings |
| FastTree model/support choices were hidden | Alphabet-compatible models, Gamma20, local support, resampling and seed controls |
| Future settings could overwrite an existing result's context | Successful calculations keep settings snapshots, versions and executed commands |
| Reruns reused error-log paths | Each execution receives a distinct log name |
| Failed MSA/tree reruns could damage existing results | Outputs are staged and validated before replacing a successful result |

## Interface references

The reviewed interfaces document useful patterns: [Jalview parameter dialogs](https://www.jalview.org/help/html/webServices/webServicesParams.html), [AliView configuration](https://www.ormbunkar.se/aliview/) and [UGENE MUSCLE controls](https://ugene.net/docs/extensions/umuscle/muscle-aligning/). AlignView applies parameter help, reusable profiles, task-specific controls and command previews. No local comparative user study was performed, and another application's MUSCLE settings are not assumed compatible with MUSCLE 5.

AlignView validates structured fields and generates previews from the same command builders used for execution. It does not accept arbitrary shell fragments as scientific parameters. Its supported controls are documented in [Parameters](PARAMETERS.md); not every native engine option is exposed.

## Recorded engine runs

Baseline and modified-setting runs are retained in [real-tools-results.json](real-tools-results.json) and [tuned-tools-results.json](tuned-tools-results.json).

| Engine | Recorded version | Modified settings exercised |
| --- | --- | --- |
| BLASTN | 2.17.0+ | Task, word size, reward/penalty, gaps, DUST, soft masking, strand, HSP and candidate limits |
| DIAMOND | 2.2.8 | Matrix, gaps, masking, composition, HSP and candidate limits |
| MMseqs2 | 18.8cc5c | Sensitivity 6.5, masking/correction disabled, 37 candidates and traceback |
| MAFFT | 7.526 | L-INS-i, gap opening 2, offset 0.1, 3 refinements |
| MUSCLE | 5.3 | Super5, abc permutation, perturbation 2, 3 refinements |
| Clustal Omega | 1.2.4 | One iteration and full distance matrices |
| EMBOSS Needle | 6.6.0.0 | EDNAFULL and separate terminal penalties |
| FastTree | 2.2.0 | LG, Gamma20, 50 local resamplings and seed 42 |

Fixtures check matches on identical sequences, valid alignments and resumable projects. Regression tests cover preview/execution agreement, engine-specific settings, HSP overlap effects, invalid profiles, CLI precedence, provenance after resuming, Slurm scripts and parameter access in a pseudo-terminal. Browser checks verify that summaries show executed search thresholds rather than newly prepared values.

```bash
python3 scripts/validate_real_tools.py \
  --require mafft blastn diamond mmseqs muscle clustalo needle FastTree
python3 scripts/validate_real_tools.py --tuned \
  --require mafft blastn diamond mmseqs muscle clustalo needle FastTree
python3 -m unittest discover -s tests -v
```

Those engine runs used an isolated prefix without replacing host software. They do not cover all matrix/gap combinations, genetic codes, engine releases, difficult datasets or the target cluster. Engine error logs remain authoritative for native incompatibilities.

Later updates added guided query selection, native VSEARCH workflows, 27 offline NCBI translation tables, explicit standalone HTML output and a single English interface/manual. Historical engine records describe their original runs and do not claim coverage of every newer option. See [VSEARCH](VSEARCH.md), [Genetic codes](GENETIC_CODES.md) and [Files](FILES.md).
