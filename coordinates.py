"""
Resolves a transcript:cdot to genomic coordinates, for every source that needs them.

Sources should not look these up themselves: here a variant is only looked up once,
also when several sources ask for it at the same time, and the result is kept for the next search.
Use Source.resolve_coordinates() from a source.
"""
import asyncio
import re
import urllib.parse
from dataclasses import dataclass

import aiohttp

ENSEMBL_URLS = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}

# Ensembl can take a while, don't let it depend on the timeout of the source that happens to ask first
MUTALYZER_TIMEOUT = aiohttp.ClientTimeout(total=5)
ENSEMBL_TIMEOUT = aiohttp.ClientTimeout(total=20)
CACHE_MAX = 1000

GENOMIC_DESCRIPTION_RE = re.compile(r"(?P<accession>NC_0*(?P<number>\d+)\.\d+):(?P<gdot>g\..+)")
GDOT_SNV_RE = re.compile(r"g\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])")
REFERENCE_RE = re.compile("^[ACGTN]+$")
ACCESSION_CHROMOSOMES = {"23": "X", "24": "Y", "12920": "MT"}


@dataclass(frozen=True)
class Coordinates:
    """
    A variant on the forward strand of the genome, pos/ref/alt as in a vcf
    (a deletion/insertion includes the base in front of it).
    gdot is the hgvs notation (7:g.55242465_55242479del).
    """
    build: str
    chr: str
    pos: int
    ref: str
    alt: str
    gdot: str


_coordinates: dict[tuple, Coordinates] = {}
_gdots: dict[tuple, str] = {}
_references: dict[tuple, str] = {}
_pending: dict[tuple, asyncio.Future] = {}


async def _once(cache: dict, key: tuple, fetch):
    """
    Returns the cached value for key, or awaits fetch() for it.
    If the same key is already being fetched, that result is awaited instead of fetching it again.
    None is not cached, so a failed lookup is tried again the next time.
    """
    if key in cache:
        return cache[key]
    if key in _pending:
        return await _pending[key]

    task = asyncio.ensure_future(fetch())
    _pending[key] = task
    try:
        result = await task
    finally:
        del _pending[key]

    if result is not None:
        if len(cache) >= CACHE_MAX:
            cache.clear()
        cache[key] = result
    return result


async def _get(session: aiohttp.ClientSession, url: str, timeout: aiohttp.ClientTimeout, as_json: bool = True):
    """Returns the parsed json (or text) of the url, None if the request or the parsing fails"""
    headers = {"Accept": "application/json" if as_json else "text/plain"}
    try:
        async with session.get(url, headers=headers, timeout=timeout) as response:
            if as_json:
                return await response.json(content_type=None)
            return await response.text()
    except (aiohttp.ClientError, asyncio.TimeoutError, ValueError):
        return None


def gdot_from_vcf(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    The hgvs genomic notation for a vcf style variant, for when Mutalyzer didn't give one.
    A duplication can't be recognised without the reference, it is written as an insertion.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0] and len(ref) + len(alt) > 2:
        ref, alt, pos = ref[1:], alt[1:], pos + 1

    if not ref:
        return f"{chrom}:g.{pos - 1}_{pos}ins{alt}"
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    span = str(pos) if len(ref) == 1 else f"{pos}_{pos + len(ref) - 1}"
    return f"{chrom}:g.{span}del" + (f"ins{alt}" if alt else "")


async def _mutalyzer_gdot(session: aiohttp.ClientSession, hgvs: str, build: str) -> str | None:
    """Returns the hgvs genomic notation (7:g.55242465_55242479del) according to Mutalyzer"""
    enc_hgvs = urllib.parse.quote(hgvs, safe="")
    url = f"https://mutalyzer.nl/api/map/?description={enc_hgvs}&reference_id={build.upper()}&filter_out=true"
    mapped = await _get(session, url, MUTALYZER_TIMEOUT)
    description = mapped.get("genomic_description") if isinstance(mapped, dict) else None
    if description and (m := GENOMIC_DESCRIPTION_RE.fullmatch(description)):
        chrom = ACCESSION_CHROMOSOMES.get(m.group("number"), m.group("number"))
        return f"{chrom}:{m.group('gdot')}"
    return None


async def _resolve(session: aiohttp.ClientSession, hgvs: str, build: str) -> Coordinates | None:
    """
    Mutalyzer is fast and gives the hgvs genomic notation, but only has the bases for substitutions
    and rejects intronic positions. For everything else Ensembl VEP is used,
    its vcf_string has the bases for deletions/insertions too.
    """
    enc_hgvs = urllib.parse.quote(hgvs, safe="")

    gdot = await _mutalyzer_gdot(session, hgvs, build)
    if gdot:
        chrom, description = gdot.split(":", 1)
        if snv := GDOT_SNV_RE.fullmatch(description):
            return Coordinates(build, chrom, int(snv.group("pos")), snv.group("ref"), snv.group("alt"), gdot)

    url = f"{ENSEMBL_URLS[build]}/vep/human/hgvs/{enc_hgvs}?vcf_string=1&content-type=application/json"
    vep = await _get(session, url, ENSEMBL_TIMEOUT)
    if isinstance(vep, list) and vep and isinstance(vep[0], dict) and isinstance(vep[0].get("vcf_string"), str):
        parts = vep[0]["vcf_string"].split("-")
        if len(parts) == 4 and parts[1].isdigit():
            chrom, pos, ref, alt = parts
            return Coordinates(build, chrom, int(pos), ref, alt, gdot or gdot_from_vcf(chrom, int(pos), ref, alt))
    return None


async def resolve(session: aiohttp.ClientSession, hgvs: str, build: str = "GRCh37") -> Coordinates | None:
    """Returns the genomic coordinates of a transcript:cdot (NM_005228.5:c.2573T>G), None if they can't be determined"""
    return await _once(_coordinates, (hgvs, build), lambda: _resolve(session, hgvs, build))


async def _resolve_gdot(session: aiohttp.ClientSession, hgvs: str, build: str) -> str | None:
    if (hgvs, build) in _coordinates:
        return _coordinates[(hgvs, build)].gdot
    gdot = await _mutalyzer_gdot(session, hgvs, build)
    if gdot is None:  # intronic for example, the full lookup knows what to do with those
        resolved = await resolve(session, hgvs, build)
        gdot = resolved.gdot if resolved else None
    return gdot


async def resolve_gdot(session: aiohttp.ClientSession, hgvs: str, build: str = "GRCh37") -> str | None:
    """
    Returns only the hgvs genomic notation of a transcript:cdot, for when the bases are not needed.
    For a deletion/insertion this saves the (slow) request to Ensembl.
    """
    return await _once(_gdots, ("gdot", hgvs, build), lambda: _resolve_gdot(session, hgvs, build))


async def _reference(session: aiohttp.ClientSession, chrom: str, start: int, end: int, build: str) -> str | None:
    url = f"{ENSEMBL_URLS[build]}/sequence/region/human/{chrom}:{start}..{end}:1?content-type=text/plain"
    sequence = await _get(session, url, ENSEMBL_TIMEOUT, as_json=False)
    sequence = (sequence or "").strip().upper()
    if len(sequence) != end - start + 1 or not REFERENCE_RE.match(sequence):
        return None
    return sequence


async def reference(session: aiohttp.ClientSession, chrom: str, start: int, end: int, build: str = "GRCh37") -> str | None:
    """Returns the forward strand of the reference genome from start to end, None if it can't be retrieved"""
    return await _once(_references, (build, str(chrom), start, end), lambda: _reference(session, str(chrom), start, end, build))
