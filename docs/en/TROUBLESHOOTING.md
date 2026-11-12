# Troubleshooting

Start with `pwd`, `alignview --version`, `alignview --diagnose`, and `alignview --lang en --docs files`. Replace `alignview` with the portable launcher's absolute path if needed.

| Problem | Check and action |
| --- | --- |
| Launcher path missing | Use the launcher inside `AlignView/`, respecting filename case; update old aliases after a move |
| Permission denied | Run through bash or restore the executable bit |
| Missing av_* module | Copy the whole application folder or use the complete zipapp |
| No interactive terminal | Supply inputs and --no-menu, or open a real terminal |
| Several queries | Select one in the guided list, or supply --query-id in batch |
| No query selector | Ensure a reference file is loaded; a single FASTA has a different exploration workflow |
| Different lengths with --msa view | Create an MSA with --msa align rather than treating raw sequences as aligned |
| Protein/nucleotide mixture | Separate alphabets or prepare an appropriate biological translation |
| Missing engine | Check diagnostics, --tools-dir, --executable or the matching --setup profile |
| Engine fails despite being found | Read its stderr log and actual parameters; check supported score/matrix combinations |
| No matches | Check method, input alphabets, threshold and filters; absence of matches is not proof of no homology |
| Old hits after changing E-value | Rerun the search; this is a calculation parameter |
| Fewer than 100 references | 100 is a maximum; engine candidates, inputs and filters may yield fewer |
| No output after --report | Add --html-out, --save or a publication/complete export destination |
| No output after --out | Add --save; --out alone only specifies a destination |
| Existing destination refused | Choose a new file/directory; batch runs do not overwrite |
| Tree absent in HTML | Infer it from at least three aligned sequences, then rebuild HTML |
| Old project has only Newick | Resume with --html-out updated.html; no new FastTree run is needed |
| Browser buttons do nothing | Use a browser with JavaScript enabled rather than a restricted IDE preview |
| HTML is stale | Rebuild, open the new file and refresh; exported reports are snapshots |
| Unexpected translated residues | Check query/reference NCBI tables and frames; genetic code is not inferred from the organism |
| ORF rejection for table 27/28/31 | These have context-dependent stops; a simple detector lacks CDS context |
| VSEARCH finds no short-sequence matches | Check length limits, masking and the word prefilter; reducing identity alone does not disable it |
| VSEARCH on proteins | Unsupported: provide nucleotide sequences or select a protein-search engine |
| Chimera screening returns no chimera | Check abundance annotations and suitable input composition; this is not absolute validation |

MAFFT explicit refinement needs global/local/fast rather than auto. MUSCLE must be version 5. FastTree GTR is for nucleotides and WAG/LG for proteins. Unsupported external-engine combinations are explained by their logs.

On SSH/HPC, choose a suitable TMPDIR and final output location. Download standalone HTML using scp or a remote file explorer and open it locally. No server is required.

Resume a whole project, not an isolated HTML. Changed input hashes mean saved input copies were modified: use an intact copy or open the new data as a new analysis. A printed recovery path can be resumed before temporary-file cleanup, but forced termination may prevent recovery.

For a useful bug report, include AlignView/Python versions, the exact command or menu path, the error, the relevant engine log and a small shareable FASTA that reproduces it. Raw commands and logs may contain local paths; share them with the intended recipients.
