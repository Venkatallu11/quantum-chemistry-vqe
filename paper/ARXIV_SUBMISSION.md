# arXiv submission guide

**Read this caveat first.** `arxiv.org` and `info.arxiv.org` are blocked
from the environment these papers were prepared in (the HTTPS proxy
returns `403 CONNECT tunnel failed`). The requirements below are from
working knowledge of arXiv's conventions, **not** read off arXiv's live
help pages. Items marked **[VERIFY]** are ones where arXiv has changed
policy before, or where the exact limit matters — check them against the
submission form and the help pages as you go. Everything marked
**[CHECKED]** was verified against the files in this repository.

---

## The two papers

| | File | Topic |
|---|---|---|
| **I** | `h4_forte_benchmark.tex` | The pre-registered Forte Enterprise 1 hardware run: accurate raw prediction, failed correction, residual diagnostics, magnitude refit |
| **II** | `zne_interpolation_fallacy.tex` | Held-out interpolation does not validate zero-noise extrapolation; the 513x shot-noise amplification and its conditioning mechanism |

They cross-reference each other. Paper I is the empirical result; Paper II
is the methodological one and is likely the more broadly citable, because
it applies to anyone running ZNE rather than to one device.

**Submission order matters.** Post Paper I first, get its arXiv
identifier, then replace the `\bibitem{companion}` placeholder in Paper II
with the real `arXiv:XXXX.XXXXX` number before submitting Paper II.
Otherwise Paper II cites a manuscript nobody can find.

---

## Categories

- **Primary: `quant-ph`** for both.
- **Cross-list: `physics.chem-ph`** for Paper I (electronic structure).
- Paper II is arguably also `stat.AP` (it is a statistical conditioning
  argument), but `quant-ph` alone is the honest home for it. **[VERIFY]**
  cross-list rules on the submission form.

**Endorsement.** arXiv requires an endorsement for a first submission to a
category where you have no posting history. **[VERIFY]** whether your
account is already endorsed for `quant-ph` — if you have posted there
before, you are not asked again.

Note: Vadim Karpusenko stated in writing that he **cannot** endorse you —
*"I don't have publications on arXiv, I will not be able to endorse you
for quant-ph."* So if you need an endorser, he is not it. An endorser must
be an established arXiv author in that category. Options: a university
contact, a co-author on a prior paper, or an academic who knows the work.
Do not leave this to the last minute; it is the single most common thing
that blocks a first submission.

---

## Abstract field

arXiv's abstract metadata field has a length limit (**[VERIFY]** the exact
figure on the form; I have used 1920 characters as the working ceiling).

**[CHECKED]** Current plain-text lengths:
- Paper I: **1781** characters
- Paper II: **1865** characters

Both are under 1920. If the form rejects either, cut the final sentence of
the conclusions list — it is the least load-bearing.

Formatting rules for the abstract field, as distinct from the LaTeX
abstract:
- Plain text only. No LaTeX macros, no `\emph`, no `$...$`.
- Write `chi^2/dof`, `x` (not `\times`), `->`, `~` for approximately.
- Keep paragraph breaks out of it; use one block.
- A `python3` snippet to generate the plain-text version from the `.tex`
  is in the repo history, but it is easier to copy from the printed
  abstracts in the two `.txt` files alongside this guide.

---

## LaTeX submission

**[VERIFY]** all of this against arXiv's TeX help page, since the
toolchain is updated periodically.

- **Both papers use REVTeX 4.2** (`revtex4-2`), which arXiv supports.
- **Upload the `.tex` source, not a PDF.** arXiv compiles it themselves
  and a PDF-only submission is flagged. **[VERIFY]** — PDF-only is
  permitted in some circumstances but discouraged and it loses you the
  full-text indexing.
- **Include the `.bbl` file.** arXiv does not run BibTeX. Both papers
  avoid this entirely by using an inline `thebibliography` environment —
  **[CHECKED]**, so there is no `.bib` and no `.bbl` to worry about for
  either paper. This was deliberate.
- **Do not upload** `.aux`, `.log`, `.out`, `.synctex.gz`, `.bib`, or
  editor backup files.
- **Figures**: neither paper currently has any. If you add them, use
  PDF or EPS (vector), not PNG screenshots. **[VERIFY]** accepted formats.
- **Page size**: let REVTeX handle it; do not force `geometry` overrides.
- **No `\usepackage` of anything exotic.** Both papers use only
  `amsmath`, `amssymb`, `booktabs`, `hyperref`, `graphicx` — all standard.
  **[CHECKED]**.

### Compile before submitting

**I could not compile these — no `pdflatex` in the preparation
environment.** This is a real gap. Before uploading:

```bash
pdflatex h4_forte_benchmark.tex
pdflatex h4_forte_benchmark.tex      # twice, for references
pdflatex zne_interpolation_fallacy.tex
pdflatex zne_interpolation_fallacy.tex
```

Or paste into Overleaf, which is faster and shows errors inline. Fix
anything that breaks. **[CHECKED]** structurally: balanced environments,
balanced braces, even `$` count, no undefined `\ref`, no undefined or
uncited `\bibitem`, no leftover placeholders. That catches most but not
all compile errors.

---

## Licence

**[VERIFY]** the current list on the form. The usual choices:

- **CC BY 4.0** — most permissive, best for citation and reuse.
- **CC BY-SA 4.0**
- **CC BY-NC-SA 4.0**
- **arXiv non-exclusive licence to distribute** — the default; you keep
  copyright, arXiv may distribute. Choose this if you intend to submit to
  a journal that is fussy about prior CC licensing.

**Recommendation: CC BY 4.0** unless you have a specific journal in mind.
It maximises reach, which is the point of posting.

---

## Metadata fields

- **Title / Authors**: plain text. Author as `Venkata Rao Allu`. No
  affiliation markup in the metadata field.
- **Comments**: use this. Suggested —
  - Paper I: `14 pages, 2 tables. Raw hardware counts, frozen analysis script and full research log at https://github.com/Venkatallu11/quantum-chemistry-vqe`
  - Paper II: `9 pages, 3 tables. Code and result files at https://github.com/Venkatallu11/quantum-chemistry-vqe`
  - **[VERIFY]** page counts after you compile.
- **ACM-class / MSC-class**: leave blank; not conventional in `quant-ph`.
- **Report-no**: leave blank.
- **Journal-ref / DOI**: leave blank now. Fill them in later if either
  paper is published, by updating the arXiv entry.
- **Funding**: the IonQ Research Credits Program acknowledgement is in the
  body of both papers. **[CHECKED]**

---

## Before you click submit

1. **Compile both papers.** Not optional — see above.
2. **Send Paper I to Vadim Karpusenko, one week ahead.** This is a
   condition of the funding, not a courtesy. His words: *"yes, no approval
   needed. Send me the draft a week before posting so the IonQ-specific
   statements are accurate."* Send Paper II too; it names IonQ hardware
   and quotes his native-gate diagnosis.
3. **Confirm the `quant-ph` endorsement situation.** Vadim has already
   said he cannot endorse.
4. **Check the cost figure.** Both papers quote $3,882.21 as IonQ's
   *estimate*. `billed_cost_usd` was still `null` in the job metadata
   **[CHECKED]**. If the real billed figure has since appeared, update it
   and say "billed" instead of "estimated".
5. **File the qiskit-ionq issue** referenced in Paper I's practical notes,
   so the citation points at something real.
6. **Decide the author list.** Both are currently single-author. If anyone
   else contributed scientifically, add them now — it is awkward later.

---

## On AI disclosure

Both papers carry this in the acknowledgements, **[CHECKED]**:

> In accordance with arXiv's policy on generative AI tools, we disclose
> that a large language model (Anthropic Claude) was used as a coding and
> drafting assistant during this work: for software implementation,
> analysis scripting, and preparation of this manuscript. All scientific
> claims, experimental design decisions and interpretations are the
> author's own, and the author takes full responsibility for the content.

Keep it. Three reasons, in order of how much they should matter to you:

1. **arXiv policy requires it.** **[VERIFY]** the exact current wording,
   but the substance has been settled since 2023: authors must disclose
   generative-AI use and remain fully responsible for content. Undisclosed
   use is a moderation issue, and moderation problems attach to your name,
   not to the paper.
2. **It is independently discoverable.** Both papers give the repository
   as their data-availability statement, and that repository's commit
   history contains `Co-Authored-By: Claude Opus 5` on a large number of
   commits **[CHECKED]**. A moderator who follows the link finds it in
   under a minute. A disclosure that is contradicted by your own published
   artifacts is worse than no disclosure.
3. **It costs you nothing.** You designed the experiments, won the
   funding, set the pre-registration, ran the hardware and made every
   scientific judgement. Tool use is not an authorship question — nobody
   discloses that they used PySCF or NumPy, but generative AI has a
   specific policy, so it gets a specific sentence. The work is yours
   either way, and the papers stand on the data.

What actually protects the papers from a "written by AI" dismissal is not
a denial. It is that every number in them is reproducible from committed
code against committed raw hardware counts, that the pre-registration was
public before the data existed, and that both papers lead with results
that went against the author's own hypothesis. That is a stronger
provenance claim than most papers can make.

---

## A third paper, not written

There is probably a third paper in the practical material: single-qubit
gates were **58.9%** of this job's cost (12,582 1q vs 1,302 2q gates,
**[CHECKED]**), device debiasing multiplies cost by ~3.5x and is on by
default, and there are three qiskit-ionq SDK defects. That is genuinely
useful to practitioners and Vadim was receptive to the idea of practical
notes. It is too thin for `quant-ph` as a research paper but would work as
a short note or a blog-style writeup. Not drafted.

A fourth candidate — the multi-molecule generalisation results (LiH,
formaldehyde, H2O2, acetaldehyde) — is **not** currently publishable as a
chemical-accuracy claim. Those twelve numbers are simulator-only and come
from the shared-frame estimator, which the hardware run showed rests
mostly on oracle information. Posting them as hardware-validated chemical
accuracy would not survive the scrutiny Paper I invites.
