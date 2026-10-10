import os
import re
import urllib.parse
from collections import Counter

import aiohttp
from jinja2 import BaseLoader, Environment
from markupsafe import escape
from pydantic_settings import BaseSettings

from inhouse_db import lookup, match, normalize_alleles, normalize_chr

from .source_result import Source

TEMPLATE = """
<table class='table table-sm mb-0 fs-6'>
    {% for center in centers %}
    {% for line in center.lines %}
    <tr>
        {% if loop.first %}<td rowspan='{{ center.lines|length }}'>{{ center.name }}</td>{% endif %}
        <td class='{{ line.style }}'>
            {% if line.check %}<small>Possible match, check {{ line.check }}:</small><br />{% endif %}
            {{ line.text }}
        </td>
    </tr>
    {% endfor %}
    {% endfor %}
</table>
"""

# green for benign, yellow for VUS and red for pathogenic
CLASSIFICATION_STYLES = {
    "benign": "table-success",
    "likely benign": "table-success",
    "vus": "table-warning",
    "likely pathogenic": "table-danger",
    "pathogenic": "table-danger",
}

GENOMIC_SNV_RE = re.compile(r"NC_0*(?P<chr>\d+)\.\d+:g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])")
REFERENCE_RE = re.compile("^[ACGTN]+$")

# how far around the variant to look, the same deletion/insertion can be written at several positions in a repeat
WINDOW = 100

# Mutalyzer and Ensembl are a lot slower than the database lookup
_HGVS_CACHE: dict[str, tuple] = {}
HGVS_CACHE_MAX = 1000


class Secrets(BaseSettings):
    inhouse_database: str = "databases/inhouse.db"


secrets = Secrets()


class InHouse(Source):
    def set_entries(self):
        self.entries = {
            ("transcript", "cdot"): self.transcript_cdot,
            ("chr", "pos", "ref", "alt"): self.chr_pos_ref_alt,
        }

    def is_complete(self) -> bool:
        if not os.path.exists(secrets.inhouse_database):
            print("Could not find in-house classification database:", secrets.inhouse_database)
            return False
        return True

    async def resolve_hg19(self, hgvs: str):
        """
        Returns (chr, pos, ref, alt) on the forward strand of GRCh37, or None.
        Mutalyzer is fast, but only gives the bases for substitutions and rejects intronic positions,
        for everything else Ensembl VEP is used, its vcf_string has the bases for deletions/insertions too.
        """
        enc_hgvs = urllib.parse.quote(hgvs, safe="")
        try:
            _, map_json = await self.async_get_json(
                f"https://mutalyzer.nl/api/map/?description={enc_hgvs}&reference_id=GRCH37&filter_out=true"
            )
            if isinstance(map_json, dict) and (m := GENOMIC_SNV_RE.fullmatch(map_json.get("genomic_description") or "")):
                chrom = {"23": "X", "24": "Y"}.get(m.group("chr"), m.group("chr"))
                return chrom, int(m.group("pos")), m.group("ref"), m.group("alt")
        except (TimeoutError, aiohttp.ClientError):
            self.timeout = False  # not fatal, Ensembl can still answer

        _, vep = await self.async_get_json(
            f"https://grch37.rest.ensembl.org/vep/human/hgvs/{enc_hgvs}?vcf_string=1&content-type=application/json",
            headers={"Accept": "application/json"},
            timeout=aiohttp.ClientTimeout(total=20),
        )
        if isinstance(vep, list) and vep and isinstance(vep[0].get("vcf_string"), str):
            parts = vep[0]["vcf_string"].split("-")
            if len(parts) == 4 and parts[1].isdigit():
                chrom, pos, ref, alt = parts
                return chrom, int(pos), ref, alt
        return None

    async def get_reference(self, chrom: str, start: int, end: int) -> str | None:
        """The forward strand of GRCh37 from start to end, needed to compare deletions/insertions"""
        try:
            _, sequence = await self.async_get_text(
                f"https://grch37.rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1?content-type=text/plain",
                headers={"Accept": "text/plain"},
                timeout=aiohttp.ClientTimeout(total=20),
            )
        except (TimeoutError, aiohttp.ClientError):
            self.timeout = False
            return None
        sequence = sequence.strip().upper()
        if len(sequence) != end - start + 1 or not REFERENCE_RE.match(sequence):
            return None
        return sequence

    async def transcript_cdot(self):
        """
        Resolves the transcript:cdot to hg19 genomic coordinates, the centers don't all report (the same) transcript
        """
        hgvs = f"{self.variant['transcript']}:{self.variant['cdot']}"
        resolved = _HGVS_CACHE.get(hgvs)
        if resolved is None:
            resolved = await self.resolve_hg19(hgvs)
            if resolved is None:
                self.log_warning(f"Could not map {hgvs} to hg19 (Mutalyzer and Ensembl VEP both failed)")
                self.html_text = "Could not determine the hg19 position of this variant."
                self.found = False
                self.complete = True
                return
            if len(_HGVS_CACHE) >= HGVS_CACHE_MAX:
                _HGVS_CACHE.clear()
            _HGVS_CACHE[hgvs] = resolved

        await self.report(*resolved)

    async def chr_pos_ref_alt(self):
        """For searches by position, these are expected to be hg19 and on the forward strand"""
        await self.report(self.variant["chr"], self.variant["pos"], self.variant["ref"], self.variant["alt"])

    async def report(self, chrom, pos, ref, alt):
        """
        Looks up the genomic variant, see inhouse_db.py for how variants are compared
        """
        self.complete = True
        chrom = normalize_chr(chrom)
        query = normalize_alleles(pos, ref, alt)
        if query is None:
            self.log_warning(f"Not a usable variant: {chrom}:{pos} {ref}>{alt}")
            self.html_text = "Could not determine the hg19 position of this variant."
            self.found = False
            return
        pos, ref, alt = query
        self.log_info(f"Looking up {chrom}:{pos} {ref or '-'}>{alt or '-'} (hg19)")

        center_names, rows, built = lookup(secrets.inhouse_database, chrom, pos - WINDOW, pos + WINDOW)
        if built:
            self.html_subtitle = f"<small>database of {escape(built)}</small>"

        # only a deletion/insertion can be written in more than one way, and only those need the reference
        sequence = None
        sequence_start = pos - WINDOW - 1
        if not (ref and alt) and any(row[1:4] != query and not (row[2] and row[3]) for row in rows):
            longest = max(len(ref), len(alt), *(len(row[2]) for row in rows))
            sequence = await self.get_reference(chrom, sequence_start, pos + WINDOW + longest)
            if sequence is None:
                self.log_warning("Could not get the reference sequence, deletions/insertions are compared less precise")

        centers = []
        for name in center_names:
            exact, possible = [], []
            for center, row_pos, row_ref, row_alt, classification in rows:
                if center != name:
                    continue
                result = match((row_pos, row_ref, row_alt), query, sequence, sequence_start)
                if result == "exact":
                    exact.append(classification)
                elif result == "possible":
                    possible.append((f"{chrom}:{row_pos} {row_ref or '-'}>{row_alt or '-'}", classification))

            # one colored line per classification, a center can have more than one for the same variant
            lines = []
            for classification, count in Counter(exact).most_common():
                text = f"{classification} ({count}x)" if count > 1 else classification
                lines.append({"text": text, "style": CLASSIFICATION_STYLES.get(classification.lower(), "")})
            for description, classification in sorted(set(possible)):
                style = CLASSIFICATION_STYLES.get(classification.lower(), "")
                lines.append({"text": classification, "style": style, "check": description})
            if not lines:
                lines.append({"text": "Not found", "style": "table-secondary"})
            centers.append({"name": name, "lines": lines, "found": bool(exact or possible)})

        self.found = any(center["found"] for center in centers)
        template = Environment(loader=BaseLoader(), autoescape=True).from_string(TEMPLATE)
        self.html_text = template.render(centers=centers)

    def get_name(self):
        return "Previously classified"
