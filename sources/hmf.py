import os

from jinja2 import BaseLoader, Environment
from pydantic_settings import BaseSettings
from util import reverse_complement

from .source_result import Source

SUMMARY_TABLE_TEMPLATE = """
<table class='table caption-top fs-6'>
    <caption>Sources</caption>
    <tr>
        <td>DoCM (<a href="http://www.docm.info/">http://www.docm.info/</a>)</td>
        <td class='table-{{ 'success' if "docm" in sources else 'danger' }}'>{{ 'Yes' if "docm" in sources else 'No' }}</td>
    </tr>
    <tr>
        <td>CIVIC (<a href="https://civicdb.org/">https://civicdb.org/</a>)</td>
        <td class='table-{{ 'success' if "vicc_civic" in sources else 'danger' }}'>{{ 'Yes' if "vicc_civic" in sources else 'No' }}</td>
    </tr>
    <tr>
        <td>CGI (Cancer Genome Interpreter)</td>
        <td class='table-{{ 'success' if "vicc_cgi" in sources else 'danger' }}'>{{ 'Yes' if "vicc_cgi" in sources else 'No' }}</td>
    </tr>
    <tr>
        <td>HMF WGS cohort</td>
        <td class='table-{{ 'success' if "hartwig_cohort" in sources else 'danger' }}'>{{ 'Yes' if "hartwig_cohort" in sources else 'No' }}</td>
    </tr>
    <tr>
        <td>HMF curated</td>
        <td class='table-{{ 'success' if "hartwig_curated" in sources else 'danger' }}'>{{ 'Yes' if "hartwig_curated" in sources else 'No' }}</td>
    </tr>
</table>
"""

class Secrets(BaseSettings):
    hmf_database: str = "databases/hmf_hotspots.tsv"

secrets = Secrets()

# (mtime, {(chr, pos, ref, alt): (gene, ensembl, pdot, sources)}), so the database is only read again when it changes
_hotspots: tuple[float, dict] = (0, {})


def get_hotspots() -> dict:
    """
    Returns the HMF database as {(chr, pos, ref, alt): (gene, ensembl, pdot, sources)}
    """
    global _hotspots
    mtime = os.path.getmtime(secrets.hmf_database)
    if _hotspots[0] != mtime:
        hotspots = {}
        with open(secrets.hmf_database) as f:
            for line in f:
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 8:
                    continue
                chrom, pos, ref, alt, gene, ensembl, pdot, sources = columns[:8]
                hotspots[(chrom, pos, ref, alt)] = (gene, ensembl, pdot, sources.strip())
        _hotspots = (mtime, hotspots)
    return _hotspots[1]


class HMF(Source):
    def set_entries(self):
        self.entries = {
            ("transcript", "cdot"): self.transcript_cdot,
            ("chr", "pos", "ref", "alt"): self.chr_pos_ref_alt,
        }

    def is_optional(self) -> bool:
        return True

    def is_complete(self) -> bool:
        if not os.path.exists(secrets.hmf_database):
            print("Could not find HMF hotspot database:", secrets.hmf_database)
            return False
        return True

    async def transcript_cdot(self):
        """
        The HMF database is on the forward strand of hg19, so use the genomic coordinates of the transcript:cdot
        """
        resolved = await self.resolve_coordinates()
        if resolved is None:
            self.hidden = True
            return  # chr_pos_ref_alt can still be tried when another source has the position

        self.report(self.lookup_variant(resolved.chr, resolved.pos, resolved.ref, resolved.alt))
        self.complete = True

    async def chr_pos_ref_alt(self):
        """
        For searches by position, and when the transcript:cdot could not be mapped.
        In the last case ref/alt are the bases of the cdot, so the reverse complement is tried too.
        """
        chrom = self.variant["chr"]
        pos = self.variant["pos"]
        ref = self.variant["ref"]
        alt = self.variant["alt"]

        info = self.lookup_variant(chrom, pos, ref, alt)
        if not info and ref and alt:
            info = self.lookup_variant(chrom, pos, reverse_complement(ref), reverse_complement(alt))
        self.report(info)

    def report(self, info):
        """
        Adds the meta data of the hotspot and shows which sources of the HMF database have it
        """
        if not info:
            # most variants are not a hotspot, don't use space on the page for a table of 'No'
            self.found = False
            self.hidden = True
            return

        sources = set()
        if info:
            gene, ensembl, pdot, sources = info
            sources = set(sources.split(","))

            if gene != "null":
                self.new_variant_data["gene"] = gene
            if ensembl != "null":
                self.new_variant_data["ensembl_transcript"] = ensembl
            if pdot not in ("null", "-"):
                self.new_variant_data["pdot"] = pdot

        template = Environment(loader=BaseLoader()).from_string(SUMMARY_TABLE_TEMPLATE)
        self.html_text = template.render(sources=sources)

    def lookup_variant(self, chrom, pos, ref, alt):
        """
        Searches the HMF database for the variant, returns (gene, ensembl, pdot, sources) or None
        """
        chrom = str(chrom)
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        return get_hotspots().get((chrom, str(pos), ref, alt))

    def get_name(self):
        return "HMF Hotspots Database"
