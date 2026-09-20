# Upstream issue triage — 2026-09-20

Scope: public open issues in `jiarong/VirSorter2`, reviewed against the current
fork and upstream source. An open issue is **not** proof that the maintainer has
not replied or fixed it. This fork does not close upstream issues or claim to
reproduce every reporter's dataset.

| Upstream issue | Evidence / decision | Fork action and acceptance boundary |
|---|---|---|
| [#233](https://github.com/jiarong/VirSorter2/issues/233), [#67](https://github.com/jiarong/VirSorter2/issues/67): interrupted resume | Maintainer says rerun same command. Our independent checkpoint/incomplete-output regression exposed unsafe assumptions in old/new nested DAG expansion. | Preserve journals and locks, completion-file checkpoints, explicit incomplete-file forcing; actual SIGTERM/SIGKILL tests. Fresh fork workdir required, not a claim that every historical reporter was affected by the same cause. |
| [#248](https://github.com/jiarong/VirSorter2/issues/248): merging large split sets | Public log shows merge failure; missing complete input/trace prevents establishing the reporter's root cause. | Same merge/checkpoint category is covered by recovery tests; **reporter-specific large-data issue remains unverified**, not marked fixed. |
| [#216](https://github.com/jiarong/VirSorter2/issues/216): missing previous classify config | Includes confusion between initial prediction and reclassification. | Reject cold `classify` before creating output/layout state; explain `all` first and reuse same `-w`. Does not invent a previous feature run. |
| [#249](https://github.com/jiarong/VirSorter2/issues/249): header-only score table | Maintainer explicitly explains that this means no viral sequences detected. | Not a scientific bug. Add `run-outcome.json` distinguishing completed zero predictions from nonzero predictions, only after native score/FASTA agreement. Missing/malformed outputs fail instead of reporting zero. |
| [#250](https://github.com/jiarong/VirSorter2/issues/250), [#59](https://github.com/jiarong/VirSorter2/issues/59): DRAM-v formatting / `potential_amg` | Cross-tool exception; posted examples alone do not identify an incorrect VS2 transform. | **Deferred** pending matched VS2/DRAM-v versions, complete inputs and minimal failing case. No category, coordinate or annotation rewriting to hide the exception. |
| [#213](https://github.com/jiarong/VirSorter2/issues/213): DRAM-v sequence lengths | Native DRAM-v route intentionally uses original sequence context; not equivalent to trimmed final FASTA. | Preserve upstream scientific behavior; portable suffix transforms only. Do not force lengths to match. |
| [#39](https://github.com/jiarong/VirSorter2/issues/39): integer sequence IDs | Maintainer says fixed. | Not claimed as a newly fixed fork issue. |
| [#43](https://github.com/jiarong/VirSorter2/issues/43): scratch fallback | Current source has both scratch and non-scratch hmmsearch branches. | Historical report alone does not justify a new patch. |
| [#198](https://github.com/jiarong/VirSorter2/issues/198): output not updating | Maintainer identifies active prophage extraction, not a hung process. | Do not impose arbitrary timeouts or kill active scientific jobs. |
| Dependency reports such as #236 / #227 / #230 | Missing packages and model/library compatibility vary by environment. | Controller modernization is separate from pinned scientific environment. No unvalidated sklearn/numpy replacement. Full scientific Linux/macOS install remains an acceptance gate. |

Additional source-level fix discovered while reviewing CLI failure paths:
`--provirus-off --max-orf-per-seq N --prep-for-dramv` previously constructed an
error message without rejecting the incompatible request. It now fails before
filesystem mutation. Non-finite `--min-score` and unsafe labels are rejected too.

## What is merged, and what is not claimed

The default fork branch contains tested controller, resume and portability work.
It is still experimental scientific software: Linux/macOS CI covers CLI/DAG,
synthetic interruption and byte-preserving transforms, not full biological
prediction parity or native Apple Silicon model availability. No release tag,
PyPI publication, upstream PR, database replacement or HPC rollout is implied.
