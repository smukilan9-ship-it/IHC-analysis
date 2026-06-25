# Spatial Biology for the Price of a Stain

### Recovering multiplex-grade immune maps from century-old brightfield pathology — and engineering the instrument that makes a cheap substitute trustworthy enough to trust with a clinical question

**A Moonshot Paper**

*Author: Mukilan Senthilkumar* · June 2026

> **Artifacts.** Working prototype (desktop + CLI application, three pipelines), full
> technical record (`ihc.md`), validation ledger, and reproducibility harnesses
> accompany this paper. Figures referenced as **[FIG-n]** are in the appendix.

---

## Abstract

The most important emerging tool in cancer medicine is *spatial* — knowing not just
which immune cells are present in a tumour, but where they sit relative to the tumour,
to each other, and to the molecular signals of immune exhaustion. This spatial
organisation is among the strongest predictors of whether a patient will respond to
immunotherapy. Today that information is locked behind multiplex imaging platforms
(CODEX, Xenium, MIBI, CyCIF) costing six figures in instrumentation and hundreds to
thousands of dollars per sample — beyond the reach of most hospitals on Earth.

Meanwhile, brightfield immunohistochemistry (IHC) — a ~100-year-old chromogenic stain —
runs in essentially every pathology laboratory in every country, at tens of dollars per
stain, and has already been performed on hundreds of millions of archived slides.

The gap between these two worlds is not fundamentally one of *resolution*.
It is one of **trust**. Every prior attempt to extract spatial, cell-resolution biology
from cheap brightfield IHC has overclaimed — inferring single-cell co-expression from
separate tissue slices, using uncalibrated statistics, accepting registrations that look
aligned but are tens of micrometres wrong — and has collapsed under scrutiny. The cheap
substitute failed not because the signal was absent, but because nobody could tell a real
signal from an artifact.

The contribution of this work is an instrument that closes that gap by being **trustworthy by
construction**: it recovers spatial immune information from ordinary brightfield IHC along
two engineered paths — same-section restaining for genuine single-cell co-expression, and
certified serial-section registration for population-level spatial association — and it is
**architecturally incapable of presenting a conclusion it cannot defend**. It fails closed.
It refuses uncertified registrations. It uses a null model proven correct by first
*breaking an earlier version*. It withdrew its own headline accuracy claim when that claim could not be backed.

This paper presents the problem the field has misunderstood, the first-principles
insight behind the instrument, its technical and statistical foundations, the validation
that exists today, the honest limits of what has been shown, and the future a trustworthy,
universally-affordable spatial-biology instrument makes possible.

---

## 1. The problem humanity has misunderstood

The conventional story of why spatial biology is expensive is a story about physics:
multiplex platforms resolve dozens of proteins on a single section, so of course they
cost more. The implied corollary is that cheap brightfield IHC is simply a lesser tool —
one or two markers, no spatial co-expression, fit for diagnosis but not for the frontier.

This framing is wrong, and the error is consequential.

The frontier question in cancer immunology is not "how many markers can you image at
once." It is "what is the **spatial relationship** between the immune system and the
tumour, at the resolution of individual cells." A cytotoxic T cell (CD8⁺) sitting inside
a tumour nest means something very different from one stranded in distant stroma. A T cell
expressing an exhaustion marker (such as TIM-3) *adjacent to* tumour means something
different again. These spatial relationships predict immunotherapy response better than
bulk marker counts. This is why multiplex exists.

But here is the misunderstanding: **most of that spatial information does not require
measuring all markers simultaneously.** It requires (a) knowing where each cell type is,
in a shared coordinate frame, and (b) a statistically honest way to ask whether two
populations are associated in space. Brightfield IHC already produces (a) — one marker at
a time, on tissue that is either the *same physical section* (restained) or an *adjacent
section* (serial). The missing piece has always been (b): a way to combine cheap,
single-marker images into trustworthy spatial conclusions.

The reason nobody trusts (b) is not that it is impossible. It is that it is **dangerously
easy to do wrong, and the wrong answer looks exactly like the right one.** Two cell
populations can appear "co-localised" simply because both prefer inflamed stroma — shared
real estate, not interaction. Two serial sections can appear perfectly registered to the
eye while being 30 µm misaligned — past the scale of the very biology you are measuring.
A membrane marker measured in the wrong cellular compartment can fabricate a strong,
publishable co-expression signal that does not exist. Each of these failure modes produces
a confident, beautiful, *false* result.

So the field bifurcated: expensive multiplex that you can trust, or cheap IHC that you
cannot. The trustworthiness gap — not the resolution gap — is what keeps spatial biology
out of most of the world's hospitals and out of the archive of slides already in storage.

**That gap is a software and methodology problem, not a hardware problem. And software
problems can be solved without building a new machine.**

---

## 2. Why existing solutions are insufficient

There are three families of existing approaches, and each is insufficient in a specific,
diagnosable way.

**Multiplex imaging (CODEX, Xenium, MIBI, CyCIF).** Scientifically excellent and the
right tool when available. But the cost and instrumentation requirements mean it will
never reach most of the world's pathology labs, and it cannot be applied retrospectively
to the enormous existing archive of brightfield slides. It answers the question for the
few. It is not a path to answering it for the many.

**Naïve "colocalisation" tools on brightfield IHC.** These overclaim. They take two
markers on two serial sections and report single-cell co-expression — a claim serial
sections physically cannot support, because adjacent sections contain *different cells*
separated by a 4–5 µm z-gap, and a marker like TIM-3 is not even restricted to the CD8
lineage. The output is confident and wrong. These tools are why the cheap substitute has
a bad reputation.

**Generic image-registration and spatial-statistics pipelines.** The components exist in
the literature (Ripley's K, mutual-information registration, KDE intensity nulls). But
assembled naïvely they fail silently on real serial tissue: registration metrics that
*alias* on repetitive tissue texture and report a 30 µm misalignment as near-zero; null
models that mistake shared tissue preference for biological association at rates
measured as high as 85–100%; multiple-testing across radii that manufactures significance.
These are known to fail because they were **built, measured failing, and
documented** (§4, §5).

The common thread: every cheap-substitute approach is **insufficiently adversarial toward
its own output**. It is built to produce a result, not to interrogate whether the result
is real. In a domain where the artifact is indistinguishable from the signal to the naked
eye, a tool that is not relentlessly self-skeptical is worse than no tool at all, because
it produces confident false discoveries that propagate into the literature and, eventually,
into clinical reasoning.

---

## 3. The first-principles insight

The central insight of this work is two-sided, and the two sides are inseparable:

**(i) The spatial immune signal that justifies multiplex is largely recoverable from
brightfield IHC — if, and only if, the recovery is made trustworthy by construction.**

This work decomposes the multiplex value proposition into what each cheap data type can *honestly*
deliver:

- **Same physical section, restained** (stain → image → strip → restain → image): the
  same cells appear in every capture. This supports the *strong* endpoint — genuine
  **single-cell co-expression** — on commodity chromogenic stains. This is a true
  multiplex substitute.
- **Adjacent serial sections, registered**: the cells differ between slices, so single-cell
  claims are forbidden. But the *populations* can be compared. With certified registration
  and a calibrated null, this supports **population- and architecture-level spatial
  association** — a large and clinically meaningful fraction of what multiplex is used for.

Critically, the design **refuses to blur these two**. The instrument enforces, in code, that a
serial-section dataset can never produce a single-cell co-expression claim. The honesty is
not a disclaimer in a footnote; it is a property of the data path.

**(ii) Trustworthiness is an instrument-design problem, and it can be engineered.**

The reason cheap substitutes failed is that they were optimised to produce conclusions.
This inverts the usual goal: the instrument is optimised to **refuse conclusions it cannot defend** —
a design discipline called *fail-closed scientific computing*. Concretely:

- Registration must be **certified** by expert anatomical landmarks with a measured,
  held-out error below an honest tolerance (≤5 µm). Uncertified pairs are *blocked* from
  producing statistics — not flagged, blocked.
- Co-expression requires **certified cell correspondence** between captures; uncertified
  tiles return `BLOCKED` with no statistics.
- The significance test uses a null model **calibrated to control false positives at the
  nominal rate**, validated by simulation — and when the first null failed that test, it
  was retired rather than shipped.
- Every result carries complete provenance, and claims the data cannot support are
  unavailable by construction.

Side (ii) is what makes side (i) believable. A cheap multiplex substitute that overclaims
is snake oil; the field has seen plenty. A cheap multiplex substitute that is *constitutionally
unable to overclaim* is a scientific instrument. The democratisation is the moon; the
fail-closed design is the engineering that earns the right to aim there.

---

## 4. Technical and scientific foundations

The prototype is a desktop and command-line application built on the open-source digital
pathology stack (QuPath with the InstanSeg `brightfield_nuclei` deep-learning model for
nucleus segmentation). It implements three pipelines.

### 4.1 Quantification (the validated core)
For each brightfield image, nuclei are segmented and each cell's DAB stain optical density
is measured and thresholded into marker-positive/negative, with per-stain thresholds and
reproducible export (CSV, GeoJSON, summary). This is the deterministic foundation both
spatial pipelines build on. The standard QuPath nuclear
measurement is deliberately preserved as the only user-facing default, having **withdrawn** an earlier, unvalidated
membrane-ring mode (§5).

### 4.2 Same-section restained co-expression (the true multiplex substitute)
Given multiple chromogenic captures of the *same physical section* (e.g. a hematoxylin
reference plus two marker images), the instrument segments nuclei **once** and reuses those
exact cell coordinates across all captures — so each cell can be classified for every
marker and assigned to one of four states (A-only, B-only, double-positive, double-negative).
Per-tile it reports the 2×2 contingency table, the double-positive count expected under
independence, the enrichment ratio, Fisher's exact odds ratio and p-value, and the phi
coefficient, with Benjamini–Hochberg FDR correction across tiles. Cells whose marker
compartment cannot be measured are labelled `UNMEASURED` and excluded from the denominator.

Two design choices are load-bearing and validated, not cosmetic:
- **Compartment correctness.** A membrane marker (CD8) must be measured in a Voronoi-clipped
  2 µm cytoplasmic ring, not the nucleus; a nuclear marker (FOXP3) in the nucleus. On a real
  tile, measuring CD8 in the *wrong* (nuclear) compartment fabricates a strong false
  co-expression (φ = 0.44, p ≈ 2×10⁻¹¹); the correct compartment shows none. The instrument
  makes the correct choice and demonstrates that the wrong one manufactures a discovery.
- **Fail-closed correspondence gate.** Co-expression statistics are produced only when the
  operator certifies that the captures share cell coordinates; otherwise the run is `BLOCKED`.

### 4.3 Serial-section spatial association (population-level)
Given two markers on adjacent sections, the instrument:
1. Segments both images.
2. **Certifies registration** via in-app expert landmark placement (§4.4).
3. Restricts analysis to the intersection of the two tissue masks (so a fold or tear in one
   section cannot masquerade as biological segregation).
4. Computes a **cross-type Ripley's K / L-function** describing association as a function of
   distance, under a **calibrated null model** (§4.5), with a single global significance
   call (the Diggle–Cressie–Loosmore–Ford rank envelope test) restricted to the
   biologically relevant 10–50 µm band.

The defensible output is *"projected cross-section spatial concordance between two
populations within anatomically corresponding regions"* — explicitly **not** single-cell
co-expression.

### 4.4 Landmark-certified registration (why it can be trusted)
Automated registration of serial sections is the central trap: individual nuclei do not
correspond across the z-gap, and every automated sub-5 µm quality
metric built here was *empirically falsified* — patch phase-correlation aliased so badly that a deliberate 30 µm shift read
as ~0.2 µm. So registration is made **operator-certified**: the expert places corresponding
anatomical landmarks (vessels, lumens, boundaries — never nuclei) on blinded grayscale images;
a similarity transform (rotation + uniform scale + translation, never a distance-distorting
warp) is fitted; and accuracy is measured by **leave-one-out** target registration error.
The verdict is one of four states — `CERTIFIED`, `LOCALLY_CERTIFIED` (only a coherent ROI
passes; analysis is clipped to it), `DEFORMED`, `NOT_CERTIFIABLE` — and statistics are
gated on it.

### 4.5 The calibrated null (the statistical heart)
Two cell populations can co-localise simply because both prefer the same tissue compartment.
A trustworthy test must subtract that shared preference. The production null is an
**intensity-reweighted inhomogeneous cross-K** (Baddeley–Møller–Waagepetersen) with
leave-one-out kernel intensity estimation and a parametric bootstrap that re-estimates
intensity per simulation. Reweighting each pair by the inverse of both populations'
local intensities makes the null expectation independent of the shared architecture — what
remains is interaction *beyond* shared preference. Homogeneous complete-spatial-randomness
is retained only as a *diagnostic baseline*; a result significant under it but not the
calibrated primary is flagged a **shared-preference artifact**, never a finding.

---

## 5. The honesty is demonstrated, not asserted

The claim that this instrument is trustworthy by construction is only credible if it can be
shown constraining its own author. It can, three times over:

1. **An earlier null model was broken — and reported.** A "three-null robustness"
   design that gated the production verdict was stress-tested over 500 simulated datasets
   carrying shared tissue preference but no real interaction. A correct test should fire ~5%
   of the time. It fired **85–100%** of the time — anti-conservative, exactly the failure
   it was meant to prevent. It was retired and replaced by the reweighted primary, which was then
   measured at **3.2%** false-positive under shared preference, **6.4%** under uniform
   randomness, and **100% / 99.2%** power to detect planted 7 µm / 25 µm attraction.

2. **A headline accuracy number was withdrawn.** A previously repeated "~90% agreement
   with manual counts" claim had no shipping artifact behind it. Rather than defend it, it was
   flagged UNVERIFIED and removed, and the instrument now ships without it.

3. **The instrument blocks its own flagship claim.** On the project's own CD8/TIM-3 serial-section
   cohort (seven usable pairs), **no pair passed registration certification** — held-out
   landmark errors ran from ~38 µm to ~92 µm against a 5 µm tolerance — so the application
   *refuses* to produce the spatial-association result the project most wanted. That outcome was kept
   rather than relaxing the threshold to manufacture a pass.

4. **And when it *does* run end-to-end on a real certified region, it still refuses the
   easy false positive.** On a real serial-section pair from a public benchmark
   (CIMA/ANHIR lung-lesion, Ki67 vs proSPC) that *locally* certified within a
   landmark-supported ROI, the naïve homogeneous-CSR diagnostic flagged a significant
   association (p = 0.001) — exactly the kind of result a less careful tool would headline.
   The calibrated primary null returned **p = 0.32: not significant.** The instrument
   correctly labelled the apparent association a *shared-tissue-preference artifact*, not a
   biological finding. This is the whole thesis in one run: the cheap substitute produces an
   answer, and the design is what stops the wrong answer from being believed.

The hard validation numbers behind these claims (full ledger accompanies the submission):

| What was tested | Result | What it establishes |
|---|---|---|
| Reweighted null vs. retired three-null design, 500 shared-preference simulations | retired nulls fired **85–100%** under a true null; production null **3.2%** | the calibration failure was real, measured, and fixed |
| Production null size / power, 500 runs | FP **3.2%** (shared) / **6.4%** (uniform); power **100%** @7 µm, **99.2%** @25 µm | controls false positives near nominal while staying sensitive |
| Cross-K estimator vs. R `spatstat` | match to **~10⁻¹⁴** (synthetic) and **1.4×10⁻¹⁰** (real CODEX spot); K-ratio exactly 1 | the core statistic equals the field-reference implementation |
| Global DCLF test under null / signal | FP **0.045**, mean p 0.515; clustered & separated patterns p ≈ **0.002** in correct direction | the single significance call is calibrated |
| Edge correction on/off, paired 500-run | FP **0.032 identical** either way | the uncorrected estimator is the right, justified choice |
| Nucleus segmentation vs. real expert masks (HNSCC, 268 tiles / 8 patients) | micro-F1 **0.776** (0.85 on a 3-tile pilot) | real, externally-validated segmentation — the *harder* number is reported |
| Registration certification on real expert landmarks (ANHIR/CIMA, 83 pairs) | **3** locally certified, **80** not certifiable, **0** globally certified; scrambled control rejected | the gate admits only measured spatial support |
| Real-data spatial statistic (Schürch CODEX, 258,385 cells) | CD8↔CD4 association in 20/40 spots; many CSR-only signals *disappear* under the calibrated null | known biology in direction, and *less* flattering after correction |
| Compartment correctness (same tile, real data) | wrong compartment fabricates φ = **0.44**, p ≈ 2×10⁻¹¹; correct compartment shows **none** | getting the biology of the measurement right is load-bearing, and the tool does |

This is the asymmetry that matters: an overclaiming project cannot afford to show its failures.
These failures are surfaced first, because they are the proof that the instrument does what it claims.

---

## 6. Honest limits of what has been shown

This work claims a **validated method and a working instrument that demonstrate feasibility** of the
multiplex substitute. It does **not** claim a completed biological result. Specifically:

- No CD8/TIM-3 pair has yet been certified end-to-end; the blocker is **data acquisition**
  (registrable, ideally restained same-section tissue), **not method**.
- Marker-positivity thresholds are not yet validated against expert cell-level ground truth;
  the same-section co-expression results to date are *machinery demonstrations*, not biology.
- The calibrated null assumes tissue architecture is coarser than its 75 µm bandwidth; where
  meaningful structure exists at cell scale, the test can become anti-conservative. This is
  disclosed and stamped in every result, and is on the roadmap to be *measured* rather than
  assumed.
- Public datasets pairing CD8 with TIM-3 on registrable tissue are effectively unavailable,
  which is why current validation rests on registration benchmarks and surrogate marker pairs.

These are the honest last mile. The instrument's design ensures none of them can be quietly
crossed: the conclusions they gate are unavailable until the evidence exists.

---

## 7. Long-term implications if this succeeds

If spatial, cell-resolution immune biology can be recovered trustworthily from brightfield IHC,
the consequences compound:

**The world's existing pathology infrastructure becomes a spatial-biology platform.** No new
instrument, no new reagent supply chain, no new training of the kind multiplex demands — the
stain and scanner are already installed in tens of thousands of labs, including in low- and
middle-income countries that will not buy a CODEX. The frontier tool of cancer immunology
becomes available where most of the world's cancer is actually treated.

**The archive becomes readable.** Hundreds of millions of brightfield slides already sit in
storage, many attached to patients with known treatment outcomes. A trustworthy substitute
turns that archive into a retrospective spatial-biology dataset of a scale no prospective
multiplex study could ever reach — potentially decoding *which spatial immune patterns predict
immunotherapy response* across enormous, already-existing cohorts.

**Trustworthy-by-construction becomes a template.** The fail-closed design — certified inputs,
calibrated nulls, claims gated on evidence, software that refuses to overclaim — is not specific
to IHC. It is a blueprint for scientific instruments in any domain where the artifact is
indistinguishable from the signal and the cost of a confident false positive is high. In an era
where AI systems increasingly produce confident outputs faster than humans can audit them, an
instrument engineered to *refuse* the conclusions it cannot defend is a primitive the rest of
biomedical AI will need.

---

## 8. The future this work is trying to build

A pathologist anywhere on Earth runs two ordinary, inexpensive stains on a tumour. The
instrument segments the cells, certifies that the tissue is well enough registered to be trusted,
asks — with a statistic calibrated to refuse shared-preference artifacts — whether cytotoxic
T cells and exhausted T cells occupy the same tumour neighbourhoods, and returns either a
defensible spatial answer or an honest "the data cannot support this." No six-figure machine.
No samples shipped to a distant core facility. And no confident lie.

That is the moonshot: not a cheaper imitation of multiplex, but a new category — **spatial
biology that is simultaneously universally affordable and constitutionally honest**. The first
of these makes the frontier reachable for everyone; the second is the only reason reaching it
that cheaply can be trusted with a human life.

The future is not discovered. It is built — and, if it is to be built in medicine, it must be
built to be unable to deceive.

---

## Acknowledgements and credited prior work

This work builds on, and credits, substantial open-source software and public data:
**QuPath** and the **InstanSeg** `brightfield_nuclei` model (segmentation); **R / spatstat**
(independent cross-validation of the spatial statistic); the **Schürch et al. 2020 colorectal
CODEX** dataset (Mendeley, CC BY 4.0) for real-data statistical validation; the **ANHIR/CIMA**
expert landmark benchmark (CC BY 4.0) for registration certification validation; the **HNSCC
mIF/mIHC comparison** dataset for real nuclear-mask segmentation validation; and the established
spatial-statistics literature (Ripley's K; Baddeley–Møller–Waagepetersen intensity-reweighted
inhomogeneous K; the Diggle–Cressie–Loosmore–Ford envelope test). All datasets and frameworks
are used under their respective licenses and cited as the independent references against which
the original contribution was tested.

The **original contribution of this work** is the synthesis: a fail-closed, certification-gated
instrument that recovers trustworthy spatial immune biology from commodity brightfield IHC along
two honest data paths, together with the calibrated reweighted null (and the documented failure
and redesign that produced it), the four-state landmark certification, and the enforced separation
between what serial-section and same-section data can each support.

---

*Appendix figures, the full technical record (`ihc.md`), the plain-language companion
(`learn.md`), the open-problems ledger (`problems.md`), and the validation harnesses accompany
this submission.*
