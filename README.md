# What does a weak burden result actually tell us?

A small, assistant-supported exploratory analysis prepared by Cathy Liu for the Spence Lab. It connects public gene-level loss-of-function (pLoF) burden estimates with the authors' GWAS intervals, then distinguishes tightly bounded small effects from unresolved estimates. It is an **analysis of processed summary statistics**, not a reanalysis of participant-level data, a new causal-gene method, or a validated clinical tool.

## Research use

At a GWAS interval, a weak burden association can mean either a precisely measured small effect or inadequate precision. This page lets the reader choose an explicit effect-size range, inspect every overlapping annotated gene, and see how the description changes with pointwise versus within-trait familywise intervals. It does not infer that the nearby gene is causal, that regulatory specificity explains the locus, or that a bounded pLoF association makes a therapeutic target safe.

The next practical use would be a region-review table for fine-mapping and variant-to-gene follow-up, or an explicit request for more pLoF carrier/replication evidence. No claim is made that this addresses an internally confirmed lab bottleneck.

## Sources

- Spence and colleagues, *Specificity, length and luck drive gene rankings in association studies*, Nature: https://doi.org/10.1038/s41586-025-09703-7
- Target-owned code/data: https://github.com/jeffspence/specificity_length_luck
- Pinned commit: `d85655034c18c2901913132aa90f1fd774c887ab`
- Upstream reuse: MIT, copyright 2025 Jeffrey Spence; the notice is retained in `UPSTREAM_LICENSE.txt`.

The canonical HGNC-indexed beta, standard-error and z-score matrices are used rather than the repository's smaller ENSG-indexed per-trait CSVs. The latter have mapping/version differences, so they are not mixed into this analysis. The genotype-frequency table is auxiliary metadata, not a carrier-count measurement or filtering rule. The source scripts were inspected and processed outputs were parsed; the authors' complete plotting/model-generation pipeline was **not** executed.

`analysis/build_data.py` contains exact pinned raw URLs and SHA-256 hashes. It reads:

| File | Fields used |
| --- | --- |
| `gamma_hats_by_trait_MAF1_hgncIDs.txt` | HGNC ID and three trait beta columns |
| `SE_hats_by_trait_MAF1_hgncIDs.txt` | HGNC ID and matching SE columns |
| `zscores_by_trait_MAF1_hgncIDs.txt` | Independent beta/SE consistency check |
| `plof_MAF1_hgncIDs.csv` | Effect-allele frequency, where unambiguous |
| `all_gwas_hits.txt` | Trait, SNP, chromosome, source interval, p value and locus ID |
| `genes.v46.gtf` | ENSG/HGNC/symbol, chromosome, start/end and biotype |
| `trait_lookup.tsv` | Source trait names |

The `.gtf` file is a processed tab-separated annotation table, **not** a conventional nine-column GTF. Its `GeneSymbol` field contains ENSG identifiers; its `gene` field contains display symbols.

## Coordinates and coverage

The repository-native annotation is GRCh37. GENCODE's official [v46lift37 header](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz) identifies the GRCh37 mapping. Independent UCSC `hg19/wgEncodeGencodeCompV46lift37` transcript-union checks match the source's one-based inclusive gene bounds for PCSK9, LDLR, NPR2 and HHIP. No liftover is performed. Membership requires matching chromosome and `gene.start <= locus.end && gene.end >= locus.start`; all source coordinates remain unchanged.

The page contains the 50 strongest supplied GWAS intervals per trait, ordered by source p value (source order breaks ties). The default pool is 20. A source p value recorded as zero is numerical underflow, not literal certainty. Neither burden significance nor the descriptive interval class filters gene membership.

The annotation contains 33,515 distinct ENSG records. All overlapping records are retained, including noncoding genes and missing burden estimates; this is not an assertion that all were assayed by a pLoF test. Fifteen assayed HGNC IDs each map to two ENSG records. Those ambiguous mappings are retained as unavailable rather than arbitrarily choosing one or duplicating an effect. Assayed genes absent from the annotation cannot be placed in intervals and are counted explicitly:

| Trait | Finite beta + positive SE | Missing from annotation | Ambiguous HGNC mapping | Displayed gene–locus rows |
| --- | ---: | ---: | ---: | ---: |
| LDL cholesterol (`30780_irnt`) | 18,505 | 752 | 15 | 2,963 |
| Standing height (`50_irnt`) | 18,524 | 758 | 15 | 5,247 |
| Red blood cell count (`30010_irnt`) | 18,514 | 756 | 15 | 3,948 |

Gene–locus rows can repeat a gene across overlapping intervals. Do not sum them and call the total unique genes. Every trait uses its full finite-estimate gene count, including genes outside displayed intervals and unmapped genes, for its Bonferroni family.

## Computation

For finite `beta` and positive finite `SE`, the normal-approximation interval is `beta ± z × SE`. Default `z = 1.959963984540054` gives pointwise 95% intervals. The optional within-trait familywise setting uses `z = Φ⁻¹(1 − 0.05 / (2 × nAssayed))`, approximately 4.692. This conservative Bonferroni setting addresses the family of assayed genes in one trait, **not multiplicity across all three traits** or misspecification of rare-carrier summary statistics.

For a chosen positive effect-size threshold `delta`:

- **Bounded-small:** `abs(beta) + z × SE <= delta`.
- **Above-threshold:** `abs(beta) − z × SE > delta`.
- **Unresolved:** neither rule is satisfied.
- **Unavailable:** missing/invalid estimate or ambiguous mapping.

The default `delta = 0.1` is an exploratory assumption in inverse-rank-normal transformed phenotype units. It is not a validated clinical threshold, a percentage, or native units such as mmol/L. These are descriptive compatibility-interval labels, not a calibrated equivalence test. Conventional two-sided burden p values are separate annotations and never determine the descriptive class.

The strongest supplied LDL interval contains 70 annotation records: 48 have usable estimates and 22 are unavailable. At `delta = 0.1`, four estimates are bounded-small with pointwise intervals; one remains so with the within-trait Bonferroni intervals. SORT1's beta is −0.0404364 with SE 0.106235, giving a pointwise interval of approximately −0.24865 to 0.16778. This interval remains compatible with effects outside ±0.1; it does not resolve the mechanism at that locus.

These computed summaries are specific to the shipped release and chosen threshold. They do not reproduce the paper's significance-filtered co-location counts: this analysis includes every overlapping annotation record and uses a different, explicitly descriptive question.

## Reproduce

Python 3.12+ and SciPy are sufficient; no credential or protected participant dataset is needed.

```sh
python -m venv .venv
.venv/bin/pip install -r analysis/requirements.txt
.venv/bin/python analysis/build_data.py --download --source-dir source-cache
.venv/bin/python -m unittest discover -s analysis -v
```

This downloads only the eight pinned source files (approximately 150 MB), not the full multi-gigabyte repository. SHA-256 verification fails closed if any bytes differ. `source-cache/` and the virtual environment are ignored by Git. Alternatively, use an existing exact upstream checkout:

```sh
python analysis/build_data.py --source-dir /path/to/specificity_length_luck
```

The command deterministically regenerates `data.json` and `data.js`; all metadata and sources are also embedded in the JSON. No timestamp, random seed or live query changes the result. `data.js` is a static browser wrapper around the same payload; there is no analytics service or server-side execution.

## Reuse and limits

A future adapter can accept canonical gene IDs, trait ID, beta, SE, optional EAF, genome build and intervals. It must verify uniqueness, positive SE, missingness, gene mapping and assembly compatibility before applying these functions. A new cohort needs its own schema/source adapter and must not be silently treated as identical to the present input. Independent replication and fine-mapping are still required before drawing biological or experimental conclusions.

Cathy's contribution is framing the question and integrating public evidence into a reviewable, interactive prototype with assistant-supported implementation. This artifact alone is not evidence of independent statistical-genetics or software-engineering competence. It contains only public processed research data, not personal participant data.
