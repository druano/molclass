# Sharing variant classifications

Molclass shows on its *Previously classified* card how each participating center has classified a variant. To add a center, a file prepared as described below is needed.

## The file

One tab-separated text file (`.tsv`, UTF-8) with a header line and one row per classified variant:

| Column | Required | Example | Notes |
|---|---|---|---|
| `gene` | optional | `EGFR` | HGNC symbol |
| `chr` | yes | `7` | Without `chr`; `X`, `Y`, `MT` |
| `pos` | yes | `55259515` | 1-based |
| `ref` | yes | `T` | Reference bases, A/C/G/T |
| `alt` | yes | `G` | Alternative bases, A/C/G/T |
| `classification` | yes | `pathogenic` | One of the five values under [Classification](#classification) |
| `transcript` | recommended | `NM_005228.5` | Used to verify the variant position |
| `cdot` | recommended | `c.2573T>G` | Used to verify the variant position |

For example:

```tsv
gene	chr	pos	ref	alt	classification	transcript	cdot
EGFR	7	55259515	T	G	pathogenic	NM_005228.5	c.2573T>G
BRCA1	17	41209079	T	TG	pathogenic	NM_007294.4	c.5266dupC
```

## Genome build and file name

The coordinates are on GRCh37 (hg19) or GRCh38 (hg38). The build is part of the file name, and the two builds are not mixed in one file:

```
CENTER_BUILD_YYYY-MM-DD.tsv
```

with `BUILD` either `hg19` or `hg38`, for example `LUMC_hg19_2026-11-01.tsv` or `LUMC_hg38_2026-11-01.tsv`.

## Deletions and insertions

Deletions and insertions are written as in a VCF file, including the base before the event:

| | `pos` | `ref` | `alt` |
|---|---|---|---|
| Deletion | `100` | `ATG` | `A` |
| Insertion | `100` | `A` | `ATG` |

`-`, `.` or empty alleles are not accepted.

## Classification

The classification is exactly one of `benign`, `likely benign`, `VUS`, `likely pathogenic` or `pathogenic`.

- Variants that are unclassified, artifacts or still under discussion are left out.
- If other class names are used in the lab, a mapping of these names to the five values above is needed.
- One row per variant. If a variant was classified more than once, only the most recent classification is included.

## No patient data

> [!IMPORTANT]
> Only the columns above. No sample IDs, patient numbers, dates, VAF, read depth or any other information that relates to a patient or a sample.

## Sending and updates

- The file is sent through a secure file transfer (for example SURFfilesender), not as an email attachment.
- Each update contains the full list, not only the new variants: each file replaces the previous one.

## What happens next

Every position is checked against the reference genome and, where given, against the transcript and c. notation; rows that don't match are reported back. The classifications then appear on the *Previously classified* card under the name of the center, next to those of the Erasmus MC and the LUMC.

If an export can't easily be put in this format, a solution can be found together.
