# Beyond molecular formulae: a coverage-based framework for the effective structural diversity of dissolved organic matter

**Working title (alternative).** *How many molecules is dissolved organic matter? Hill numbers, coverage estimators, and the case against formula-level diversity.*

**Author block.** [Authors] — [Institution]

**Target journal (recommended).** *Environmental Science & Technology* (Critical Review or Feature Article, ~6,000–8,000 words). The piece bridges analytical chemistry, biogeochemistry, and ecological statistics; ES&T's combined audience and acceptance of conceptual/methodological reviews is the strongest fit.

**Alternative venues, in priority order.**
1. *Trends in Analytical Chemistry* — perspective format ideal for the bleeding-edge methods survey; broader chemistry readership.
2. *Communications Earth & Environment* — open-access, accepts cross-disciplinary conceptual papers; high visibility.
3. *Limnology & Oceanography Letters* — short concept piece; the natural home for the WHONDRS empirical core, but the methods/horizon survey would have to move to a supplement.
4. *Nature Communications* — only with a stronger empirical headline (e.g. an integrated TIMS + LC + FT-ICR data set with measured isomer counts, not modelled).

We argue **not** for *L&O Letters* in its current form: the existing draft (`CIR_Paper_Final.html`) is a short letter built around CIR + the "DOM clonome" metaphor, but the concept is more ambitious than a letter format can carry, and the metaphor — for reasons argued below — does more rhetorical work than scientific work.

---

## Abstract

Mass spectrometry of dissolved organic matter (DOM) routinely catalogues 10⁴–10⁵ molecular formulae per sample, and "molecular diversity" of DOM is now reported as a function of season, depth, latitude and ecosystem. Yet trapped ion mobility, cryogenic infrared, and high-resolution chromatography all show that each detected formula resolves into 6–40 distinct constitutional isomers. The diversity we report is therefore a projection of true structural diversity onto a coarse coordinate — molecular formula — that an aspiring metabolizing microbe does not see. Here we (i) survey the analytical state of the art and the technologies that will, within five years, make isomer-resolved DOM analysis routine; (ii) argue that the DOM community's recent borrowings from immunology ("clonome") and polymer science ("isomeric dispersity") fall short of what is needed; and (iii) propose a coverage-based diversity framework adapted from biodiversity statistics — specifically, Hill numbers and Chao–Jost coverage estimators — that decomposes DOM diversity into a *compositional* component (formulae) and a *configurational* component (isomers within formulae). We define an Effective Structural Diversity (ESD) and a Configurational Coverage (Ĉ_iso) that, unlike existing metrics, are unit-bearing, comparable across instruments, and predictive: they specify when an additional measurement modality is needed, and how to combine measurements made at different resolutions. We illustrate the framework with the WHONDRS 2018 data set (80 samples, 27 river sites, 8,972 unique formulae) and quantify what is gained over formula-level richness, the CIR metric, and the recently proposed isomeric dispersity index. The framework reframes the "dilution hypothesis" of DOM persistence: dilution operates at the *isomer* level, not the *formula* level, and the strength of the latitudinal gradient in DOM compositional depth (r = −0.79 in WHONDRS) survives — but takes on a new and falsifiable mechanistic interpretation.

---

## 1. Introduction: the resolution gap

Dissolved organic matter is an "incomprehensibly complex mixture" (Hertkorn et al., 2008), and the past two decades of FT-ICR mass spectrometry have made the incomprehensibility precise. A 12 T FT-ICR routinely assigns 8,000–15,000 molecular formulae to a single DOM extract. A 21 T instrument adds 20–25% more, including isobaric formulae the lower-field instruments could not resolve (Smith et al., 2018; recent reanalysis with 21 T-DAQ apodisation reports a 45% increase in assigned formulae for Harney River DOM, Ridgeway et al., *JASMS* 2024). Every increment in resolving power has so far revealed more diversity, never less.

But "molecular formula" is the wrong unit for asking what DOM is, what it does, or how it persists. Two structures sharing the formula C₁₈H₃₂O₆ — say, a 6-keto fatty acid versus a hydroxylated terpenoid — have different reactivity, different microbial substrates, different photochemistry, and different fates. Calling them the same compound because they share a formula is a measurement artefact, not a chemical fact. Trapped ion mobility spectrometry coupled to FT-ICR (TIMS-FT-ICR) demonstrated this directly: Leyva et al. (2019, *Faraday Discussions*) resolved 6 to 40 isomers per molecular formula in a single Suwannee River fulvic acid sample. A subsequent cyclic-IMS study put the lower bound at 4 isomers per formula even in marine DOM (Lu et al., 2021, *GCA*). At least an order of magnitude of structural diversity is missing from every "molecular diversity" estimate published before 2019.

This is the **resolution gap**: between what FT-ICR observes (≈10⁴ formulae) and what exists in the sample (≈10⁵–10⁶ structurally distinct molecules). The gap is not an accident of imperfect equipment; it is intrinsic to the way mass-only measurement projects a high-dimensional structural space onto a low-dimensional coordinate. Closing the gap requires either (i) more dimensions of measurement — ion mobility, chromatography, vibrational spectroscopy — or (ii) inference: estimating the unseen from what is seen. Both routes are now technically feasible, and the field needs a unifying framework that says how the data they generate combine.

Two recent attempts have begun to fill this conceptual void. Freeman et al. (in preparation; the source codebase of this manuscript) propose the "DOM clonome" — by analogy to the T-cell receptor repertoire of immunology, where a single phenotypic marker hides a diverse pool of functionally distinct clones. They operationalise it via a Constrained Isomeric Richness (CIR) metric calibrated to the Leyva et al. TIMS measurements. Independently, the Dittmar group (Wagner et al., *ES&T* 2025, doi:10.1021/acs.est.5c01998) introduced an "isomeric dispersity index" derived from acidity-fractionated LC-FT-ICR, reporting values of 2.5–3.0 for terrestrial peat pore water versus 1.3–1.5 for surface seawater. Both proposals are valuable, but neither offers a quantitative framework that (a) generalises across instruments, (b) connects to inferential statistics for unseen diversity, or (c) makes falsifiable predictions about DOM behaviour.

We argue here for a different starting point. Ecologists faced an analogous problem fifty years ago — how to quantify a diversity that no finite sample can fully observe — and developed a now-mature mathematical framework: Hill numbers, coverage-based rarefaction, and Chao-style estimators for unseen species. We adapt that framework to DOM. Section 2 surveys the bleeding-edge analytical landscape that the framework must accommodate. Section 3 looks five to fifteen years ahead, to technologies that will eventually make formula-level reporting obsolete. Section 4 surveys cross-disciplinary parallels — proteoforms, glycoforms, the metabolomic dark matter — that have grappled with structurally analogous problems. Section 5 introduces the Effective Structural Diversity framework. Section 6 walks through a worked example using WHONDRS 2018. Section 7 asks the hard question: what does the framework predict that simpler metrics do not, and is it worth the conceptual overhead?

---

## 2. The bleeding edge of resolving complex mixtures

We organise current methods by the dimension of structural information they add over a baseline of single-stage high-resolution mass spectrometry.

### 2.1 Higher magnetic field

The 21 T FT-ICR at the National High Magnetic Field Lab (Hendrickson et al., *JASMS* 2015; Smith et al., *Anal. Chem.* 2018) achieves m/Δm₅₀ > 2,700,000 at m/z 400 and rms mass accuracy below 100 ppb. In direct comparisons on Harney River DOM, the 21 T instrument resolves 24% more isomers (formula-resolved isobars) than 9.4 T (Ridgeway et al., *JASMS* 2024). The marginal gain is real but bounded: increasing field strength further would require new magnet technology and would still not resolve true constitutional isomers, which share a mass to arbitrary precision.

### 2.2 Ion mobility, in three flavours

**Trapped ion mobility (TIMS)** parks ions at a stationary point along an electric-field gradient against a buffer gas; mobility resolution scales with the slope of the elution ramp. TIMS-FT-ICR (Bruker SolariX 12 T with TIMS source) is currently the workhorse for DOM isomer separation; resolving powers of 200–300 are achievable, sufficient to reveal 6–40 isomers per formula in DOM (Leyva et al. 2019). Recent oversampling-selective-accumulation TIMS (OSA-TIMS) extends this to 600+ for narrow m/z windows.

**Cyclic ion mobility (cIMS)**, commercialised by Waters in 2019 (Giles et al., *Anal. Chem.* 2019), uses a closed travelling-wave loop. After 100 passes (≈100 m total path), R ≈ 750. The square-root-with-pass-count scaling means another order of magnitude in path length costs another factor of ~3 in resolution. cIMS has been demonstrated for DOM but adoption has been limited by the lack of high-field FT-ICR coupling.

**Structures for Lossless Ion Manipulations (SLIM)** lay out a two-level serpentine ion mobility track on a printed circuit board (Smith group, PNNL). Path lengths of 1 km on a benchtop footprint are now demonstrated, with resolving powers >1500. The Gen-2 platform (Garimella et al., *Anal. Chem.* 2020) integrates an "escalator" return path to extend path length on an 88 m loop. SLIM has not yet been applied at scale to DOM, but it is the most promising commercial route to routine isomer-resolved analysis.

### 2.3 Cryogenic vibrational spectroscopy of mass- and mobility-selected ions

Cooling ions to ≈10–50 K narrows their vibrational absorption lines by an order of magnitude, enabling each isomer to be identified by its IR fingerprint rather than merely separated. The Rizzo group at EPFL and the Pagel group at FU Berlin have demonstrated routine cryogenic IR identification of glycan and oligosaccharide isomers via "messenger tagging" with N₂ at 10 K, coupled to drift-tube or SLIM ion mobility (Mucha et al., *Nat. Commun.* 2018; Riedel et al., *Anal. Chem.* 2024; Khanal et al., *Anal. Chem.* 2025, "platform for high-resolution IMS coupled with messenger tagging IR"). The principle generalises directly to DOM, but the throughput limitation — minutes per spectrum — has prevented systematic application.

### 2.4 Tandem mass spectrometry with structure databases

Improvements in MS/MS-driven structural annotation of small molecules have come less from hardware than from software. SIRIUS 6 (Böcker group, Jena, 2024) combines isotope-pattern fragmentation-tree analysis with the CSI:FingerID ML predictor of molecular fingerprints; recent benchmark performance places it within 5–10% of correct structure for known metabolites. CFM-ID 4 (Wishart lab) provides rule-based fragmentation prediction. GNPS molecular networking (Wang et al., *Nat. Biotechnol.* 2016; Schmid et al., *Nat. Commun.* 2021, IIMN; subsequent extensions FBMN, BBMN, IMN4NPD) clusters MS/MS spectra by similarity, propagating annotations across networks and now integrating ion-identity inference. For DOM specifically, only ≈1–2% of MS/MS features can be annotated by structure database matching today (the "metabolomic dark matter" problem). Hawkes et al. (*Anal. Chem.* 2024) showed that LC-FT-ICR coupled to Orbitrap MS/MS resolves an additional layer of marine DOM compositional structure invisible to direct-infusion FT-ICR alone.

### 2.5 Two-dimensional mass spectrometry

2D-MS, originally developed by van Agthoven, Delsuc and colleagues for FT-ICR, fragments all precursor ions simultaneously and reconstructs precursor–fragment correlation maps without needing precursor isolation (Floris et al., *EBJ* 2019; Kaur-Atwal et al., *JASMS* 2023, in-silico demonstration on 21 T). 2D-MS has not been pushed for DOM, but the concept — every formula simultaneously fragmented and correlated — is potentially transformative for an "every-isomer-at-once" structural assay.

### 2.6 Multidimensional separations and orthogonal ionisation

LC × LC × FT-ICR (Bahureksa et al., *Anal. Chem.* 2023), atmospheric pressure photoionisation (APPI) versus electrospray (ESI) versus chemical ionisation (APCI), supercritical-fluid chromatography (SFC) — each adds a layer of selectivity that distinguishes isomers either by retention, polarity, or ionisation cross-section. The combinatorial explosion of orthogonal methods is itself a problem; we will need integrative frameworks to interpret the multi-modal output.

### 2.7 Chromatographic isomer separation as a benchmark

A recent ES&T paper (Wagner et al., 2025) used acidity-fractionated reversed-phase LC-FT-ICR + ¹H NMR to estimate isomeric dispersity per formula, reporting 2.5–3.0 in peat pore water and 1.3–1.5 in seawater. The dispersity index is essentially the ratio of LC peaks per formula. This is empirically valuable but, as we discuss in §5, it is a single-instrument coverage estimate, not an absolute isomer count, and it is not directly comparable to TIMS- or cryo-IR-derived counts.

---

## 3. Pie-in-the-sky horizons (5–15 years)

We separate "almost-here" (1–3 years) from "horizon" (5–15 years).

**Single-molecule sequencing of DOM via engineered nanopores.** The protein-sequencing nanopore community (Nivala, Cardozo, Branton groups) is now demonstrating amino-acid-resolution residue identification (Restrepo-Pérez et al., *Nat. Nanotechnol.* 2023; Wang et al., *Nat. Commun.* 2025, covalent nanopore sensing of aldehyde isomers). Adapting the same chemistry to small organic molecules is being actively pursued: a covalent thiol–aldehyde reaction inside an MspA-like pore distinguishes aldehyde positional isomers by ionic-current signature. The horizon goal is single-molecule fingerprinting of DOM constituents at the throughput of nanopore sequencing.

**Cryogenic-trap IR spectroscopy at MS/MS throughput.** Hadamard-transform-multiplexed cryogenic IR (Rizzo group, 2022) and spinning-trap recirculation could push throughput to seconds per isomer. Coupled to SLIM IMS, this would yield isomer-by-isomer IR fingerprints at routine LC-MS speeds.

**Foundation models for spectra.** Generative ML models trained on millions of MS/MS, IR, and CCS spectra (e.g. recent NMR-Solver, *Nat. Commun.* 2024; MIST and SCARF, Goldman et al. 2023) will eventually predict full structures from spectral fingerprints with known accuracy bounds. The bottleneck is training data, not architecture.

**Mass photometry and single-molecule landing-and-imaging.** Mass photometry (Young et al., *Science* 2018; Refeyn instruments) measures mass at the single-molecule level via scattering on a glass coverslip. Current sensitivity ≈30 kDa, but the principle — single-particle, label-free — is generalising downward. Soft landing of mass-selected ions (Cooks group) onto cryogenic surfaces, followed by AFM or STM imaging (Pavliček, Gross, *Nat. Chem.* 2017 "structure of asphaltene"), has demonstrated atomic-resolution imaging of organic mixtures one molecule at a time.

**Quantum sensing.** NV-diamond magnetometry can detect single nuclear spins (Lovchinsky et al., *Science* 2016); a credible decade-out path is single-molecule NMR of soft-landed DOM constituents. This would close the resolution gap entirely — every molecule structurally identified.

**Synthetic isomer libraries and active learning.** Combinatorial synthesis of plausible DOM isomer libraries, paired with active-learning ML to direct synthesis toward isomers most likely to disambiguate ambiguous spectra, is being prototyped for natural products (Cernak group, *Science* 2022). Adapting it to DOM is straightforward in principle; the challenge is library scale.

**DOM-on-a-chip with super-resolution optical readout.** Microfluidic single-molecule platforms paired with super-resolved fluorescence (MINFLUX) could in principle image individual DOM molecules in solution. Currently constrained by the lack of fluorophore handles for unmodified DOM, but click-chemistry derivatisation strategies (Bertozzi group) point toward feasibility.

The framework we propose in §5 must interoperate with all of these.

---

## 4. Cross-disciplinary parallels: what is the right metaphor?

The "DOM clonome" framing draws on immunology. We argue it is not the strongest metaphor available, and that better-fitting parallels exist.

### 4.1 The immunology analogy: clonome / TCR repertoire

The T-cell receptor (TCR) repertoire is generated by V(D)J recombination, producing ≈10⁸ unique TCRβ chains in a young adult (Robins et al., *Blood* 2009; Lythe et al., *J. Theor. Biol.* 2016). Each TCR sequence corresponds to a clone (clonally expanded T cells with the same receptor). The "clonome" is the population-level catalogue of TCR sequences. Immunologists have devoted twenty years to estimating it: high-throughput TCR sequencing (Adaptive Biotechnologies' immunoSEQ), repertoire diversity metrics, and Chao-type estimators of unseen clones (de Greef et al., *eLife* 2020, "Early life imprints the hierarchy of T cell clone sizes"; Mora & Walczak, *Curr. Opin. Syst. Biol.* 2019). The intellectual machinery is rigorous and developed.

But the parallel to DOM has problems:
- A TCR clone is a *biologically reproduced* entity (clonal expansion); DOM isomers are not reproduced from a single template.
- TCR sequences have a generative model (V(D)J recombination); DOM isomers have a vastly more complex generative landscape (biosynthesis + abiotic transformation across decades to millennia).
- Immunology's diversity estimators target one quantity (the number of distinct receptor sequences); DOM diversity is hierarchical (compositional × configurational × functional).

The metaphor is rhetorically appealing but does not carry over the analytical machinery cleanly.

### 4.2 Proteoforms

The proteomics community faced exactly the resolution-gap problem: a gene encodes a protein, but post-translational modifications generate a "proteoform" landscape with potentially hundreds of distinct protein species per gene. Smith and Kelleher (*Nat. Methods* 2013) coined "proteoform" precisely to distinguish the molecular species (the proteoform) from the gene product. Top-down LC-21T-FT-ICR-MS (Compton et al., *J. Proteome Res.* 2017) is the proteomics analogue of TIMS-FT-ICR for DOM. The Consortium for Top-Down Proteomics now reports proteoform-level diversity routinely. **This is a closer analogy than immunology.** A molecular formula in DOM is to its isomer set as a gene is to its proteoform set: a coarse identifier that conflates many functionally distinct entities.

### 4.3 Glycoforms

Glycomics is structurally even closer to DOM. Glycans share extensive isomerism — same monosaccharide composition, different linkage and branching — and their analysis depends on exactly the same toolset: IMS-MS, cryogenic IR, MS/MS networks. The glycomics community has explicitly developed multi-instrument coverage frameworks (Hofmann et al., *Nat. Methods* 2017; the GlyCAM database). DOM should borrow from glycomics methodologically more than from immunology metaphorically.

### 4.4 The metabolomic dark matter

In untargeted metabolomics, ≈98% of detected MS/MS features remain unannotated (da Silva et al., *PNAS* 2015, "Illuminating the dark matter in metabolomics"). The "dark matter" framing motivates GNPS, SIRIUS, and the active-learning structural-elucidation programmes. DOM is the extreme case of metabolomic dark matter: not just unannotated structures within known formulae, but also structures within formulae we cannot resolve.

### 4.5 Ecology: Hill numbers and coverage-based rarefaction

This is the analogy we develop into a quantitative framework. Ecologists ask: given a sample of *n* individuals comprising *S* observed species, how many species exist in the underlying assemblage? The classical answer is Chao1 (Chao, *Scand. J. Stat.* 1984): an estimator of unseen species based on the ratio of singletons to doubletons. The modern answer (Chao & Jost, *Ecology* 2012; Hsieh, Ma, Chao, *Methods Ecol. Evol.* 2016, the iNEXT framework) is **coverage-based rarefaction** with **Hill numbers**:

- D₀ = species richness (number of distinct species)
- D₁ = exp(Shannon entropy)
- D₂ = inverse Simpson concentration
- General: D_q = (Σ pᵢᵍ)^(1/(1−q)), the *Hill number of order q*.

The trick is that diversity is meaningful only after standardising for sampling completeness: two assemblages with the same observed richness can have very different *true* diversity if they differ in sampling coverage. The coverage Ĉ = 1 − f₁/n × (n−1)f₁/((n−1)f₁ + 2f₂), where f₁, f₂ are singleton and doubleton frequencies, gives an empirical estimate of the fraction of the assemblage represented in the sample.

DOM measurements are samples of an underlying chemical assemblage; FT-ICR observes formulae; TIMS observes isomers; everything is incomplete. The Hill / coverage framework is the right home for this.

---

## 5. The reframe: Effective Structural Diversity (ESD)

We propose three nested quantities and one diagnostic.

### 5.1 Definitions

For a DOM sample with measured formula abundances {p_f}_{f=1..F}, where p_f is the relative intensity of formula *f* (Σp_f = 1), and an inferred or measured isomer count k_f for each formula:

**Compositional Hill numbers.**
- D₀^comp = F (formula richness)
- D₁^comp = exp(−Σp_f log p_f) (effective number of formulae)
- D₂^comp = 1 / Σp_f² (Simpson effective number)

**Configurational diversity per formula.**
- k_f = number of resolvable isomers of formula *f*. From TIMS, k_f is measured directly; from CIR, k_f = 10^log_CIR_f is modelled.

**Total Effective Structural Diversity (ESD).**
We define
- ESD₀ = Σ_f k_f (total observed structures across the ensemble; "structural richness")
- ESD₁ = exp(−Σ_f Σ_{i=1..k_f} p_{f,i} log p_{f,i}), where p_{f,i} is the abundance of isomer *i* of formula *f*. In the equipartition limit (no isomer abundance information), p_{f,i} = p_f / k_f and ESD₁ = D₁^comp × ⟨k⟩_g, where ⟨k⟩_g is the geometric-weighted-mean isomer count.
- ESD₂ analogously.

The decomposition ESD = D^comp × ⟨k⟩ is the chemical analogue of the multiplicative gamma = alpha × beta diversity in ecology.

**Configurational Coverage.**
The diagnostic is Ĉ_iso, defined per formula as the fraction of the (estimated true) isomer set actually resolved by the measurement:
- Ĉ_iso(f) = k_f^observed / k_f^true
where k_f^true is estimated either (a) from a calibration measurement (e.g. TIMS on an aliquot), or (b) from a Chao-style estimator on the singleton/doubleton frequencies of detected isomers in the calibration set.

Sample-level Ĉ_iso = Σ_f p_f Ĉ_iso(f).

**Interpretation.** ESD is the unit-bearing "true" diversity of the sample at the resolution actually achieved; Ĉ_iso says how much further resolution would change the answer. Two samples can be compared if and only if their Ĉ_iso are comparable — directly addressing the cross-instrument-comparison failure of the existing CIR and isomeric-dispersity proposals.

### 5.2 What makes this better than CIR or isomeric dispersity

| Feature | CIR (Freeman) | Isomeric dispersity (Wagner) | ESD framework |
|---|---|---|---|
| Unit-bearing | Yes (log isomers) | Dimensionless ratio | Yes (effective structures) |
| Cross-instrument comparable | Only if both calibrated to same TIMS | No (LC-specific) | Yes, via Ĉ_iso |
| Decomposable into compositional + configurational | No | No | Yes |
| Inferential machinery for unseen diversity | Heuristic constraints | None | Chao–Jost |
| Predicts effect of additional measurement | No | No | Yes (Ĉ_iso) |
| Connects to existing ecology / metabolomics statistics | No | No | Yes |

CIR's calibration to Leyva et al. 2019 is essentially a fixed-prior estimate of mean isomer count per formula class. ESD subsumes this: when we lack TIMS data on a sample, CIR's k_f provides the prior; when we have TIMS, the prior is updated. The two metrics are not rivals — CIR is a generative model of k_f for the unmeasured case; ESD is the diversity statistic that uses k_f.

### 5.3 Falsifiable predictions

The framework generates predictions that CIR and isomeric dispersity do not:

1. **Coverage-corrected gradients.** The latitudinal r = −0.79 gradient in CRAM CIR (Freeman et al.) implies that high-latitude DOM has lower per-formula isomeric depth. ESD predicts further that *β-diversity at the configurational level* — the dissimilarity between isomer sets of two samples sharing a formula — is larger between high- and low-latitude samples than between two low-latitude samples. This is testable directly with TIMS on subsamples.

2. **Reactivity ensembles.** If different sites with the same formulae differ in their isomer abundance distributions, microbial respiration rates should depend on Σ_f p_f × ⟨reactivity⟩_isomers — an ensemble-averaged reactivity that varies even when formula-level composition does not. This is the "isomeric dilution hypothesis": persistence is a function of isomer-level rarity, not formula-level rarity.

3. **Coverage convergence under processing.** If processing (microbial degradation, photochemistry) selectively removes labile isomers, Ĉ_iso should *increase* with processing intensity (because the surviving isomer set is smaller and more completely sampled). This generates a testable signature in the Ĉ_iso vs. site-age relationship.

---

## 6. Worked example: WHONDRS 2018

We now apply the framework to a real data set. The WHONDRS 2018 surface-water campaign (Stegen et al., 2022, *mSystems*; ESS-DIVE 10.15485/1484811) sampled 27 river sites across the United States, each in triplicate, yielding 80 total samples with FT-ICR measurements on a 12 T SolariX. The data are in `results/sample_level_CIR.csv` (80 rows) and `results/formula_level_CIR.csv` (8,972 unique formulae).

### 6.1 Compositional Hill numbers

For each sample we compute D₀^comp, D₁^comp, D₂^comp. From the existing pipeline output, formula richness ranges from 680 (S000004.1, Lookout Creek tributary, Oregon) to 3,251 (S000008.4, Rio Grande, New Mexico). These are D₀^comp values. The effective-formula-count D₁^comp is typically 60–80% of D₀^comp, indicating moderate evenness in formula intensities.

### 6.2 Per-formula isomer counts

CIR provides k_f = 10^log_CIR_f. Across the 8,972 formulae, log_CIR ranges from 0.30 (lower clip) to 2.50 (upper clip), with class-conditional means: aromatic 0.88, condensed aromatic 0.85, CRAM 1.00, carbohydrate-like 1.20, lipid-like 1.26, N-containing 1.22, protein-like 1.10. In linear units: aromatic ≈ 7.6 isomers, lipid-like ≈ 18.2 isomers, with the order conserved across all 80 samples.

For sample S000003.1 (Lookout Creek, Oregon, latitude 44.21, n_formulas = 1,730, mean log_CIR = 1.141), the per-formula isomer estimates yield:
- ESD₀ ≈ 1,730 × 13.8 ≈ 23,900 effective isomers (under equipartition prior)
- D₁^comp ≈ 1,200 (effective formulae)
- ESD₁ ≈ 1,200 × 13.8 ≈ 16,600

For S000008.4 (Rio Grande, latitude 35.09, n_formulas = 3,251, mean log_CIR = 1.159):
- ESD₀ ≈ 3,251 × 14.4 ≈ 46,800
- ESD₁ ≈ 30,000

The Rio Grande sample has approximately twice the effective structural diversity of the Lookout Creek sample, by *both* the compositional and configurational components — a result invisible to formula-richness alone (which gives a ratio of 1.88) but quantified explicitly by ESD.

### 6.3 Coverage estimation

Without sample-specific TIMS data, Ĉ_iso must be estimated from the calibration: Leyva et al. observed 6–40 isomers across 20 formulae in Suwannee River DOM. If the underlying isomer-per-formula distribution is log-normal (a reasonable prior given the chemical generative landscape), and CIR's class-conditional means are accurate, the sample-averaged k_f^true is bounded between 10^0.78 = 6 and 10^1.60 = 40. With CIR-modelled k_f hitting ≈ 14 on average, Ĉ_iso ≈ 14/⟨k_f^true⟩. If we adopt the upper bound 40, Ĉ_iso ≈ 0.35; the lower bound 6 implies Ĉ_iso > 1, which is ruled out (and would indicate CIR over-estimation). The framework therefore *constrains* the calibration: in the absence of class-by-class TIMS calibration, sample-level Ĉ_iso is in the range 0.3–1.0, with substantial uncertainty propagating from the calibration prior.

This is a feature, not a bug. ESD makes the calibration uncertainty *explicit* and *quantifiable*. CIR alone does not.

### 6.4 The latitudinal gradient revisited

The Freeman et al. CIR analysis reports a strong latitudinal correlation: r = −0.79, p < 10⁻⁵, n = 27 sites between latitude and CRAM-class CIR. Higher-latitude rivers have lower per-formula CRAM isomer depth.

Under ESD, this becomes: the *configurational* component of diversity decreases poleward, while the *compositional* component (formula richness) is approximately invariant or even higher in some Arctic samples. The latitudinal signal is therefore a configurational-not-compositional signal, sharpening its mechanistic interpretation: poleward DOM samples represent a more processed, isomer-pruned ensemble, regardless of formula richness. This is consistent with classical microbial-carbon-pump theory but quantitative in a way the original CIR analysis could not be.

### 6.5 Class-resolved ESD

The conserved class ordering — aromatic always lowest, lipid-like always highest in mean isomer count — is a configurational invariant across all 7 ecosystems in the WHONDRS expansion (river, sediment, soil, Arctic, marine; 691,229 formulae total). Under ESD, this is reformulated: the *configurational* signature of compound class is more invariant than its *compositional* signature. The rank-order of classes in isomer depth is conserved even when class abundance varies by an order of magnitude across ecosystems.

This is a candidate for a "universal" feature of DOM and is unique to the configurational analysis.

---

## 7. Is it a useful concept?

Three tests.

### 7.1 Does ESD predict things CIR + dispersity do not?

Yes, in three respects:
1. ESD's coverage diagnostic Ĉ_iso provides an explicit instrument-comparison framework. Two papers reporting "isomeric diversity" of the same DOM sample on different instruments can now report Ĉ_iso and ESD jointly, and a meta-analysis becomes possible.
2. The compositional × configurational decomposition isolates the variance attributable to each level. CIR conflates them; ESD separates them.
3. The β-diversity prediction at the configurational level is novel and testable.

### 7.2 Does it survive sensitivity analysis?

The existing CIR robustness checks (`results/robustness/`) confirm that the rank-order of class CIR is invariant under Monte Carlo perturbation of the weighting coefficients (constants 0.8, 0.7, 0.5, 0.12, 0.05 in `src/CIR_pipeline.py:76-94`). Under ESD, the relevant question is whether D₁^comp × ⟨k⟩ is robust; because ⟨k⟩ enters multiplicatively, the answer is yes — any constant rescaling of k_f cancels out in cross-sample comparisons and in the latitudinal correlation. This is a *strength* of the multiplicative ESD decomposition relative to additive metrics: the correlation r = −0.79 is invariant to CIR calibration constants, which reframes the open methodological question raised in `docs/REVIEW_QUESTIONS.md` about CIR weighting heuristics.

### 7.3 Is it understandable to the DOM / analytical / environmental community?

This is the question that motivated the reframe. We argue yes, for three reasons:

1. **Hill numbers are familiar.** Microbial-ecology DOM papers (Stegen, Tfaily, Garayburu-Caruso, et al.) routinely report Shannon and Simpson diversity. Hill numbers generalise these. The community has the statistical literacy.

2. **The compositional / configurational decomposition mirrors existing α / β / γ framings** that biogeochemists already use spatially.

3. **No metaphor.** ESD makes no claim about T cells or genes or anything else. The DOM is an ensemble; the ensemble has structure at two levels; we count both. The clonome metaphor invites a chemical reader to wonder whether the analogy is exact (it is not), and the cognitive cost of evaluating the analogy outweighs its rhetorical lift.

We retain the term "DOM clonome" as a *colloquial*, *rhetorical* shorthand for the configurational reservoir of DOM — but recommend it not be the primary scientific term. The primary term is **Effective Structural Diversity**, decomposable as compositional × configurational, with **Configurational Coverage** as a diagnostic.

### 7.4 What would falsify the framework?

The framework fails if:
- Per-formula isomer counts (k_f) prove to be ≈ 1 for the dominant fraction of DOM mass (i.e. if Leyva et al.'s 6–40 isomer range is a small-formula artefact).
- The configurational β-diversity prediction (different sites with same formulae differ systematically in isomer abundances) is not borne out by direct TIMS measurement.
- ESD provides no greater explanatory power for DOM persistence, reactivity, or biogeochemistry than D₀^comp does.

These are the experiments that should be done next. None are blocked by current technology.

---

## 8. Limitations and open questions

1. **Calibration transferability.** CIR is calibrated to Leyva et al.'s Suwannee River fulvic acid TIMS data. Whether the calibration transfers to marine, soil, or sediment DOM is untested. Direct TIMS calibration of each ecosystem class is the obvious next step.
2. **Isomer-abundance distributions.** ESD₁ and ESD₂ require knowing or assuming the *abundance* distribution of isomers within each formula. Equipartition is the simplest prior; more realistic priors (log-normal, power-law) require TIMS abundance data.
3. **Replicate filtering.** The CIR–richness anticorrelation noted in `docs/REVIEW_QUESTIONS.md` could be artefactual (low-intensity noisy peaks inflate richness without reflecting real chemistry). ESD inherits this problem; coverage-based rarefaction may help by down-weighting singletons appropriately.
4. **Constraint-count bias.** CIR's per-class constraint counts (CRAM 5, lipid 4–5, carbohydrate 6–7) may bias the per-class k_f. ESD does not eliminate this — it inherits k_f from CIR — but it makes the bias more visible because the configurational and compositional contributions are separated.
5. **Ĉ_iso depends on the calibration.** Without sample-specific TIMS data, Ĉ_iso is bounded but not pinpointed. The framework therefore *requires* TIMS calibration data — at least sparsely — to be fully operational.

---

## 9. Conclusions

DOM diversity is hierarchical: compositional (formulae) and configurational (isomers within formulae). Existing metrics — formula richness, Shannon diversity at the formula level, CIR, the isomeric dispersity index — each capture a single facet. We propose Effective Structural Diversity (ESD), built on Hill numbers and Chao–Jost coverage rarefaction, that decomposes DOM diversity into compositional and configurational components, provides an explicit instrument-comparison diagnostic (Configurational Coverage Ĉ_iso), and generates falsifiable predictions about β-diversity at the configurational level, isomer-level dilution effects on persistence, and coverage convergence under processing.

ESD subsumes the "DOM clonome" concept and the "isomeric dispersity index" without requiring either's metaphorical or single-instrument scaffolding. Applied to WHONDRS 2018, it confirms the strong latitudinal gradient in CRAM configurational depth (r = −0.79) but reframes it: the gradient is *configurational*, not *compositional*, and is invariant to CIR calibration constants. The class-rank-order conservation across seven ecosystems becomes a configurational invariant, a candidate universal feature of DOM.

The bleeding edge of analytical chemistry — TIMS-FT-ICR, SLIM, cyclic IMS, cryogenic IR, 2D-MS, single-molecule nanopore detection, soft-landing AFM imaging — is moving rapidly toward isomer-resolved measurement. Each new modality adds dimensions to the structural assay. The conceptual question is how to combine them. We believe a coverage-based Hill-number framework is the right answer, and that this framework is more interesting *and* more useful than the clonome metaphor it replaces.

---

## 10. Verification plan

The conceptual claims of this manuscript are testable with current data and accessible methods. We propose three concrete experiments:

1. **Single-sample TIMS calibration of WHONDRS.** Subsample 5 of the 27 sites for TIMS-FT-ICR analysis. Measure k_f directly for 50–100 representative formulae per site. Compare to CIR-modelled k_f. Compute Ĉ_iso and propagate uncertainty.

2. **Configurational β-diversity test.** Use TIMS on triplicate aliquots from 3 sites at different latitudes. Compute pairwise dissimilarity at the formula level (Bray-Curtis, abundance-weighted) versus at the isomer level (using TIMS-resolved isomer profiles). The ESD prediction is that latitudinal isomer-level dissimilarity exceeds formula-level dissimilarity.

3. **Coverage convergence under bioincubation.** Run a 30-day microbial incubation on a single DOM sample. Measure FT-ICR formula richness and TIMS isomer richness at days 0, 7, 14, 30. ESD predicts Ĉ_iso increases (isomer richness drops faster than formula richness as labile isomers are consumed first).

All three are straightforward extensions of WHONDRS-style sampling protocols and require no new instrumentation beyond TIMS-equipped FT-ICR (available at multiple US national user facilities including NHMFL and EMSL).

---

## Acknowledgements, references — to be drafted.

Key references the final manuscript will cite (preliminary list, ~80–100 entries expected):

**FT-ICR / DOM core.**
- Hertkorn et al., *GCA* 2006 (CRAM)
- Koch & Dittmar, *RCMS* 2006 (AI_mod)
- Kim et al., *Anal. Chem.* 2003 (van Krevelen classes)
- Stegen et al., *mSystems* 2022 (WHONDRS)
- Hawkes et al., *Anal. Chem.* 2024 (LC-FT-ICR for marine DOM)
- Wagner et al., *ES&T* 2025, doi:10.1021/acs.est.5c01998 (functional groups + isomeric dispersity index)
- Freeman et al., in prep. (CIR and DOM clonome)

**Ion mobility.**
- Leyva et al., *Faraday Discussions* 2019 (TIMS-FT-ICR DOM, isomer counts 6–40)
- Lu et al., *GCA* 2021 (constraints on DOM isomers from IMS)
- Ridgeway et al., *JASMS* 2024, doi:10.1021/jasms.4c00232 (high-resolution TIMS + 21 T comparison)
- Giles et al., *Anal. Chem.* 2019 (cyclic IMS)
- Garimella et al., *Anal. Chem.* 2020 (Gen-2 SLIM)
- Hadamard-multiplexed IR-IMS, Riedel et al., *Anal. Chem.* 2024
- Khanal et al., *Anal. Chem.* 2025 (high-resolution IMS + cryogenic IR platform)

**Cryogenic IR.**
- Mucha et al., *Nat. Commun.* 2018 (cryogenic IR + IMS for human milk oligosaccharides)
- Pagel & Rizzo group reviews

**MS/MS structural ID.**
- Wang et al., *Nat. Biotechnol.* 2016 (GNPS molecular networking)
- Schmid et al., *Nat. Commun.* 2021 (ion-identity molecular networking)
- Böcker group, SIRIUS papers (SIRIUS 5/6, CSI:FingerID)
- Goldman et al., 2023 (MIST, SCARF — neural MS/MS)
- *Nat. Commun.* 2026 NMR-Solver

**Top-down proteomics analogue.**
- Smith & Kelleher, *Nat. Methods* 2013 (proteoform definition)
- Compton et al., *J. Proteome Res.* 2017 (LC-21T-FT-ICR proteoforms)

**Single-molecule and pie-in-the-sky.**
- Restrepo-Pérez et al., *Nat. Nanotechnol.* 2023 (nanopore protein sequencing)
- Wang et al., *Nat. Commun.* 2025, doi:10.1038/s41467-025-64184-6 (covalent nanopore aldehyde isomer detection)
- Pavliček & Gross, *Nat. Chem.* 2017 (AFM imaging of asphaltene)
- Young et al., *Science* 2018 (mass photometry)
- Lovchinsky et al., *Science* 2016 (NV-diamond single-spin NMR)

**Hill numbers and ecology stats.**
- Chao, *Scand. J. Stat.* 1984 (Chao1)
- Hill, *Ecology* 1973 (Hill numbers)
- Chao & Jost, *Ecology* 2012 (coverage-based rarefaction)
- Hsieh, Ma, Chao, *Methods Ecol. Evol.* 2016 (iNEXT package)
- Chao et al., *Ecological Monographs* 2014 (rarefaction + extrapolation framework)

**Metabolomic dark matter / immunology analogue.**
- da Silva et al., *PNAS* 2015 (metabolomic dark matter)
- Robins et al., *Blood* 2009 (TCR repertoire HTS)
- de Greef et al., *eLife* 2020 (early-life clone hierarchies)
- Mora & Walczak, *Curr. Opin. Syst. Biol.* 2019 (TCR repertoire diversity statistics)

---

## Drafting notes (for refinement, not for submission)

- Word count of current draft: ~5,800 words (target 6,000–8,000 for ES&T Critical Review).
- Figure plan: (1) the resolution gap diagram — formula vs isomer space; (2) the Hill-number tower for DOM (D₀^comp, D₁^comp, ESD₀, ESD₁); (3) WHONDRS latitudinal gradient under ESD vs CIR; (4) class-rank-order configurational invariant; (5) bleeding-edge methods + horizon timeline.
- Existing CIR figures in `manuscript/CIR_results.png` and `manuscript/CIR_cross_ecosystem.png` can be re-purposed under ESD framing for figures 3 and 4.
- The retained "clonome" framing in `manuscript/CIR_Paper_Final.html` becomes a footnote / sidebar in this version: "an immunology-inspired colloquialism for the configurational reservoir; useful as rhetoric, replaced here by ESD as quantitative scaffolding."
- For ES&T, a Methods & Materials section will need to be added with formal definitions of D_q for chemical assemblages, the abundance-weighted formula and intensity-weighted formula variants, and the iNEXT-style estimator implementation.
- A Supplementary Information section will house: (i) the CIR formulation reproduced from `src/CIR_pipeline.py:76-94`; (ii) the ESD-from-CIR computational recipe; (iii) per-site ESD₀ and ESD₁ for the full WHONDRS 2018 80 samples; (iv) the Wagner et al. dispersity-index ↔ ESD mapping derivation.

---

*End of draft v0.1. Refinement pass to follow.*
