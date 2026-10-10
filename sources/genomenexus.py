import urllib.parse

from util import reverse_complement
from .source_result import Source, SourceURL

class GenomeNexus(Source):
    def set_entries(self):
        self.entries = {
            ("transcript", "cdot"): self.transcript_cdot,
            ("chr", "pos", "ref", "alt"): self.chr_pos_ref_alt
        }

    async def check_validity(self, query):
        enc_query = urllib.parse.quote(query)
        api_url = f"https://www.genomenexus.org/annotation/{enc_query}?fields=annotation_summary"
        resp, json = await self.async_get_json(api_url)
        return json.get("successfully_annotated", False)

    async def chr_pos_ref_alt(self):
        """
        Simply add a URL to the variant, won't work for more del, ins, dup, etc.
        Probably could add support for those, but would take other fields then pos, ref, alt.
        """
        chrom = self.variant["chr"]
        pos = self.variant["pos"]
        ref = self.variant["ref"]
        alt = self.variant["alt"]
        if not ref or not alt:
            self.found = False
            return

        # Check if annotation exists for ref/alt, if not try the reverse complement
        if not await self.check_validity(f"{chrom}:g.{pos}{ref}>{alt}"):
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
            if not await self.check_validity(f"{chrom}:g.{pos}{ref}>{alt}"):
                self.html_text = (f"No annotation found for either {chrom}:g.{pos}{self.variant['ref']}>{self.variant['alt']} "
                                  f"or {chrom}:g.{pos}{ref}>{alt} (reverse complement).")
                self.found = False
                return

        url = f"https://www.genomenexus.org/variant/{chrom}:g.{pos}{ref}>{alt}"
        self.html_links["main"] = SourceURL(f"{chrom}:g.{pos}{ref}>{alt} (hg19)", url)

    async def transcript_cdot(self):
        """
        Link to the variant by its genomic (hgvs) notation
        """
        resolved = await self.resolve_coordinates()
        if resolved is None:
            self.found = False
            return

        url = f"https://www.genomenexus.org/variant/{resolved.gdot}"
        self.html_links["main"] = SourceURL(f"{resolved.gdot} (hg19)", url)
        self.complete = True

    def get_name(self):
        return "Genome Nexus"