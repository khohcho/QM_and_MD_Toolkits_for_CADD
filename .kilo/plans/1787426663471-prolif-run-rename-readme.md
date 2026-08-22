# Plan: prolif.py → prolif_run.py rename + README updates (license removal, repo structure)

## Goal

1. Rename the ProLIF script without touching its content: `Post_MD_Analysis/prolif.py` → `Post_MD_Analysis/prolif_run.py`.
2. Update both README files to reflect (a) the removed project license and (b) the current repo structure.
3. Add a third-party tools & licenses table to the ProLIF README (ProLIF is Apache-2.0, not GPL-3.0 as previously claimed).

## Context / Decisions (confirmed with user)

- New filename: `prolif_run.py` (with `.py` extension).
- Root README: remove the now-misleading "open-source" wording.
- ProLIF README: remove the GPL-3.0 License section and replace it with a "Third-Party Tools & Licenses" table.
- License files: already gone from the repo (`LICENSE` deleted, root README license section removed in earlier commits). Do not re-add any license text.

## Task List

### 1. Rename the script (no content changes)

- `git mv Post_MD_Analysis/prolif.py Post_MD_Analysis/prolif_run.py`
- Do NOT edit the file contents. Verify with `git diff --stat` that only a rename is registered.

### 2. Update `Post_MD_Analysis/README.md`

- Change heading `## prolif.py` → `## prolif_run.py`.
- Change usage command `python prolif.py` → `python prolif_run.py`.
- Delete the `## License` section (lines 69-70) and its GPL-3.0 text.
- Add a new `## Third-Party Tools & Licenses` section with a table: tool, purpose in this pipeline, license.

  | Tool | Purpose in `prolif_run.py` | License |
  |---|---|---|
  | MDAnalysis | `.tpr`/`.xtc` handling, atom selections, topology repair | GPL-2.0-or-later |
  | ProLIF | IFP calculation, barcode, ligand-network plots | Apache-2.0 |
  | RDKit | Molecular representation used by ProLIF | BSD-3-Clause |
  | NumPy | Array math / axis ticks | BSD-3-Clause |
  | pandas | IFP DataFrames and CSV outputs | BSD-3-Clause |
  | matplotlib | Agg backend PNG rendering (barcode) | Matplotlib License (PSF-based) |
  | IPython | HTML display support for network/3D views | BSD-3-Clause |
  | py3Dmol | 3D visualization HTML | BSD-3-Clause |

- Keep the rest (features, installation, input files, known issues) unchanged.

### 3. Update root `README.md`

- Remove "open-source" from the first line: "A collection of quantum mechanics (QM) and molecular dynamics (MD) analysis toolkits for **Computer-Aided Drug Design (CADD)**."
- Replace the bare `### Post_MD_Analysis` section with a short structure section that matches the current repo:
  - List `Post_MD_Analysis/prolif_run.py` with a one-line description and link to `Post_MD_Analysis/README.md`.
  - One line noting QM toolkits will be added in the future (repo name covers QM & MD).
- Do NOT add any license section (project currently has no license).

## Validation

- `git status` shows: renamed `prolif_run.py` (detected as rename), modified READMEs; no other changes.
- `git diff --stat` / `git diff -M` confirms `prolif.py` content is byte-identical except the path.
- Grep the repo to confirm no stale references remain: `prolif.py`, `GPL`, `License`, `open-source`, `LICENSE`, `DOI`.
- Visually render-check both READMEs (headings, table formatting, links) — no broken relative links.

## Status

Steps 1-3 are already implemented and validated in the working tree (rename staged via `git mv`; README edits unstaged).

## 4. Commit (user requested)

The rename motivation (per user): running `python prolif.py` caused confusion with the `prolif` package imported inside the script, so the script was renamed.

- Stage only these files (never `.kilo/`):
  - `Post_MD_Analysis/README.md`
  - `README.md`
- The rename `Post_MD_Analysis/prolif.py -> Post_MD_Analysis/prolif_run.py` is already staged.
- Commit message (subject + body):

  `Rename prolif.py to prolif_run.py and update READMEs`

  `Renamed the ProLIF analysis script to prolif_run.py to avoid confusion with the prolif package imported inside the script when running it directly with python. Also removed the outdated GPL-3.0 license section, added a third-party tools & licenses table, and aligned the root README with the current repository structure.`

- After commit: verify `git log --oneline -1` and `git status` (only untracked `.kilo/` may remain).

## Status after commit (verified)

- Commit `4890f66` exists locally; `git status` reports "ahead of 'origin/main' by 1 commit".
- Remote: `origin` = `https://github.com/khohcho/QM_and_MD_Toolkits_for_CADD.git`, upstream = `origin/main`.
- The user's doubt is explained: changes are committed locally but not yet visible on GitHub.

## 5. Push (user requested)

- Run: `git push origin main`
- Never push or stage `.kilo/`.
- Verify: `git status` shows branch up to date with `origin/main`; `git log origin/main --oneline -1` equals `4890f66`.

## Out of Scope / Notes

- No changes to the script's code, outputs, or naming inside the script.
- No author/attribution header added to the script for now (keeps the rename a pure `git mv` and respects "no changes" constraint); can be added later as a separate commit or as an Author line in the README.
- No new LICENSE file; project remains license-free as the user removed it.
- Commit is now requested and covered by step 4 above.
