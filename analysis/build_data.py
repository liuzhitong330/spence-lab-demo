#!/usr/bin/env python3
"""Rebuild the descriptive Spence-lab locus precision explorer.

Reads pinned public processed outputs, never individual-level UK Biobank data.
Use --source-dir for an existing checkout, or --download with a cache directory.
The output is deterministic; no current timestamps are embedded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import urllib.request

from scipy.stats import norm

COMMIT = "d85655034c18c2901913132aa90f1fd774c887ab"
REPOSITORY = "https://github.com/jeffspence/specificity_length_luck"
RAW = f"https://raw.githubusercontent.com/jeffspence/specificity_length_luck/{COMMIT}"
POINTWISE_Z = 1.959963984540054
TRAITS = [("30780_irnt", "LDL cholesterol"), ("50_irnt", "Standing height"),
          ("30010_irnt", "Red blood cell count")]
FILES = {
    "data/gamma_hats_by_trait_MAF1_hgncIDs.txt": "1a23756f86ee5d9fcd18573ef64e1b1aecdf2e643a085f7ed817c9146fbce36c",
    "data/SE_hats_by_trait_MAF1_hgncIDs.txt": "84fdfcf35b8cbdab9c69615e2480a18a6e79e140a3da3ad36e39d92090521d82",
    "data/zscores_by_trait_MAF1_hgncIDs.txt": "888738f57ab541535d55091e4da11be76ba616540ba49e1b7c3fe2fda20cf06f",
    "data/plof_MAF1_hgncIDs.csv": "08607ef5874df407ae289b6d45d0b367ef3481b6bf558c44f1a5a5add5c23919",
    "data/all_gwas_hits.txt": "4099521d287ee1b267ff36d64f4b17de54aa34e29f360f2d12564675240b5933",
    "data/genes.v46.gtf": "2341659f6db06e64141c9dcb808d6049547e9f696e61f6cba649a31080a04d7d",
    "data/trait_lookup.tsv": "e09cfddbaa67762209eef4d99e879f7a676c444d0c8f4208108c685fbf519717",
    "LICENSE": "54f92b728128a18b6184edd87ac63d1617c9457452fb0c58fd0b0667e9fb12fa",
}


def number(value):
    """Missing/non-finite input remains missing, never zero."""
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def valid_estimate(beta, se):
    return beta is not None and se is not None and math.isfinite(beta) and math.isfinite(se) and se > 0


def classify(beta, se, delta=0.1, z=POINTWISE_Z):
    """A descriptive interval classification, not a causal/equivalence test."""
    if not math.isfinite(delta) or delta <= 0 or not math.isfinite(z) or z <= 0:
        raise ValueError("delta and z must be positive finite numbers")
    if not valid_estimate(beta, se):
        return "unavailable"
    if abs(beta) + z * se <= delta:
        return "bounded-small"
    if abs(beta) - z * se > delta:
        return "above-threshold"
    return "unresolved"


def overlaps(gene_start, gene_end, locus_start, locus_end):
    """Closed-source-coordinate overlap; do not mix genome assemblies."""
    return gene_start <= locus_end and gene_end >= locus_start


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_inputs(root, download=False):
    for name, expected in FILES.items():
        path = root / name
        if not path.exists() and download:
            path.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(f"{RAW}/{name}", path)
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; supply --download or --source-dir")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"Source SHA-256 mismatch: {name}: {actual}")


def read_matrix(path):
    result = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert all(trait in reader.fieldnames for trait, _ in TRAITS)
        for row in reader:
            key = row["hgnc_id"]
            if key in result:
                raise ValueError(f"Duplicate matrix ID {key}: {path.name}")
            result[key] = {trait: number(row[trait]) for trait, _ in TRAITS}
    return result


def count_classes(genes, delta=0.1, z=POINTWISE_Z):
    counts = Counter(classify(g["beta"], g["se"], delta, z) for g in genes)
    return {key: counts[key] for key in ["bounded-small", "above-threshold", "unresolved", "unavailable"]}


def build(root):
    beta = read_matrix(root / "data/gamma_hats_by_trait_MAF1_hgncIDs.txt")
    se = read_matrix(root / "data/SE_hats_by_trait_MAF1_hgncIDs.txt")
    z_source = read_matrix(root / "data/zscores_by_trait_MAF1_hgncIDs.txt")
    if set(beta) != set(se) or set(beta) != set(z_source):
        raise ValueError("Canonical beta/SE/z gene sets differ")
    with (root / "data/plof_MAF1_hgncIDs.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len({r["hgnc_id"] for r in rows}):
        raise ValueError("Duplicate EAF ID")
    eaf = {r["hgnc_id"]: number(r["effect_allele_frequency"]) for r in rows}
    with (root / "data/genes.v46.gtf").open(newline="") as handle:
        annotation = list(csv.DictReader(handle, delimiter="\t"))
    if len(annotation) != len({r["GeneSymbol"] for r in annotation}):
        raise ValueError("Duplicate ENSG annotation")
    hgnc_count = Counter(r["hgnc_id"] for r in annotation)
    genes = []
    for row in annotation:
        chrom = row["chr"].removeprefix("chr")
        start, end = int(row["start"]), int(row["end"])
        if start < 1 or end < start:
            raise ValueError(f"Invalid annotation interval {row}")
        genes.append({"id": row["GeneSymbol"], "symbol": row["gene"] or row["GeneSymbol"],
                      "hgncId": row["hgnc_id"], "chr": chrom, "start": start, "end": end,
                      "biotype": row["conseq"]})
    with (root / "data/all_gwas_hits.txt").open(newline="") as handle:
        locus_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (root / "data/trait_lookup.tsv").open(newline="") as handle:
        trait_labels = dict(csv.reader(handle, delimiter="\t"))

    traits = []
    for trait_id, label in TRAITS:
        assayed = {key for key in beta if valid_estimate(beta[key][trait_id], se[key][trait_id])}
        z_family = float(norm.ppf(1 - 0.05 / (2 * len(assayed))))
        differences = []
        for key in assayed:
            released_z = z_source[key][trait_id]
            if released_z is not None:
                difference = abs(beta[key][trait_id] / se[key][trait_id] - released_z)
                differences.append(difference)
                if difference > 1e-4:
                    raise ValueError(f"Canonical beta/SE disagrees with released z: {trait_id} {key} {difference}")
        source_loci = [r for r in locus_rows if r["pheno"] == trait_id]
        if len(source_loci) < 50:
            raise ValueError(f"Fewer than 50 available loci for {trait_id}")
        # Decimal prevents tiny p values from collapsing to an artificial float tie.
        ranked = sorted(enumerate(source_loci), key=lambda item: (Decimal(item[1]["pval"]), item[0]))[:50]
        loci = []
        for rank, (_, row) in enumerate(ranked, start=1):
            chrom, start, end = str(row["chr"]), int(row["window_lower"]), int(row["window_upper"])
            if end < start:
                raise ValueError(f"Invalid locus window: {row}")
            member_genes = []
            for g in genes:
                if g["chr"] != chrom or not overlaps(g["start"], g["end"], start, end):
                    continue
                key = g["hgncId"]
                reason = None
                if key not in beta:
                    reason = "not_in_burden_matrix"
                elif hgnc_count[key] != 1:
                    reason = "ambiguous_hgnc_mapping"
                elif not valid_estimate(beta[key][trait_id], se[key][trait_id]):
                    reason = "missing_or_invalid_estimate"
                b = beta[key][trait_id] if reason is None else None
                s = se[key][trait_id] if reason is None else None
                freq = eaf.get(key) if reason is None else None
                z_value = b / s if reason is None else None
                member_genes.append({**g, "beta": b, "se": s, "eaf": freq,
                                     "z": z_value, "burdenP": float(2 * norm.sf(abs(z_value))) if z_value is not None else None,
                                     "availability": reason or "available"})
            member_genes.sort(key=lambda g: (g["start"], g["end"], g["id"]))
            loci.append({"id": f"{trait_id}:{row['name']}", "sourceId": row["name"], "rank": rank,
                         "chr": chrom, "start": start, "end": end, "snp": row["SNP"],
                         "pos": int(row["pos"]), "p": number(row["pval"]), "pText": row["pval"],
                         "genes": member_genes, "defaultCounts": count_classes(member_genes),
                         "familyCounts": count_classes(member_genes, z=z_family)})
        if len(loci) != len({x["id"] for x in loci}):
            raise ValueError("Duplicate locus IDs")
        traits.append({"id": trait_id, "name": label, "sourceName": trait_labels[trait_id],
                       "nAssayed": len(assayed), "zFamily": z_family,
                       "nSourceLoci": len(source_loci), "nMatrixRows": len(beta),
                       "nUnmappedAssayed": len(assayed - set(hgnc_count)),
                       "nAmbiguousAssayed": sum(hgnc_count[key] > 1 for key in assayed),
                       "maxZDiscrepancy": max(differences), "loci": loci})
    payload = {
        "meta": {"schemaVersion": 1, "repository": REPOSITORY, "commit": COMMIT,
                 "paper": "https://doi.org/10.1038/s41586-025-09703-7",
                 "analysisType": "Reanalysis of public processed gene-level summary statistics",
                 "assembly": "GRCh37", "coordinateConvention": "Unchanged source coordinates; inclusive interval overlap",
                 "units": "Inverse-rank-normal transformed phenotype units",
                 "pointwiseZ": POINTWISE_Z, "defaultDelta": 0.1, "defaultTrait": "30780_irnt",
                 "defaultLocus": traits[0]["loci"][0]["id"],
                 "defaultGene": "ENSG00000134243",
                 "defaultSelectionReason": "SORT1 illustrates imprecision inside the strongest supplied LDL GWAS interval; overlap is not a causal assignment.",
                 "nAnnotationGenes": len(annotation),
                 "coverageNote": "Every overlapping GRCh37 annotation record is retained. Finite assayed genes lacking annotation: LDL 752, height 758, RBC 756; each trait also has 15 ambiguous assayed HGNC mappings, retained as unavailable. Interval overlap does not nominate a causal gene.",
                 "mappingPolicy": "Retain every overlapping annotation row; ambiguous HGNC-to-ENSG mappings have unavailable estimates. Matrix genes without annotation cannot be located and are counted separately.",
                 "locusSelection": "Top 50 source GWAS intervals per trait, sorted by supplied p value with original source order breaking ties; no burden-based selection.",
                 "zeroPValueNote": "Some very small GWAS p values are recorded as 0 in the source. This is source numerical underflow, not a literal probability of zero.",
                 "assemblyEvidence": ["https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz", "https://api.genome.ucsc.edu/getData/track?genome=hg19;track=wgEncodeGencodeCompV46lift37;chrom=chr1;start=55505117;end=55530525"],
                 "sources": [{"path": name, "url": f"{RAW}/{name}", "sha256": digest} for name, digest in FILES.items()]},
        "traits": traits,
    }
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path.home() / ".cache/spence-demo-inputs")
    parser.add_argument("--download", action="store_true", help="Download missing pinned public inputs only; verify SHA-256")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    verify_inputs(args.source_dir, args.download)
    result = build(args.source_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(result, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    (args.output_dir / "data.json").write_text(encoded + "\n", encoding="utf-8")
    (args.output_dir / "data.js").write_text("window.SPENCE_DATA = " + encoded + ";\n", encoding="utf-8")
    print(json.dumps({t["id"]: {"nAssayed": t["nAssayed"], "loci": len(t["loci"]),
                     "zFamily": t["zFamily"], "geneLocusRows": sum(len(l["genes"]) for l in t["loci"]),
                     "unmappedAssayed": t["nUnmappedAssayed"], "ambiguousAssayed": t["nAmbiguousAssayed"]}
                      for t in result["traits"]}, indent=2))


if __name__ == "__main__":
    main()
