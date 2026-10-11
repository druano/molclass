import os

import yaml
from jinja2 import BaseLoader, Environment
from markupsafe import escape
from pydantic_settings import BaseSettings

from .source_result import Source, SourceURL

TEMPLATE = """
{% if gene.remark %}<p class='small mb-1'>{{ gene.remark }}</p>{% endif %}
<table class='table table-sm mb-1 fs-6'>
    <tr>
        <td class='text-nowrap' title='{{ table.columns.report }}'>PA-verslag</td>
        <td class='text-break' style='white-space: pre-line'>{{ gene.report or '-' }}</td>
    </tr>
    <tr>
        <td class='text-nowrap' title='{{ table.columns.mtb }}'>MTB</td>
        <td class='text-break' style='white-space: pre-line'>{{ gene.mtb or '-' }}</td>
    </tr>
</table>
<small class='text-muted' title='{{ table.notes|join("\n") }}'>
    Alleen (waarschijnlijk) pathogene varianten, TMB &lt; 20 mut/Mb en niet evident subklonaal.
    {% if gene.changed %}Tekst gewijzigd ten opzichte van de vorige versie.{% endif %}
</small>
"""


class Secrets(BaseSettings):
    table3_database: str = "databases/table3.yml"


secrets = Secrets()

# (mtime, table), so the file is only read again when it changes
_table: tuple[float, dict] = (0, {})


def get_table() -> dict:
    """
    Returns 'Tabel 3' of the NVVP, VKGL and VKGN: for which genes a (likely) pathogenic variant
    found in a tumor is a reason to advise germline diagnostics. See the README for the format of the file.
    """
    global _table
    mtime = os.path.getmtime(secrets.table3_database)
    if _table[0] != mtime:
        with open(secrets.table3_database, encoding="utf-8") as f:
            table = yaml.safe_load(f)
        table["genes"] = {str(gene).upper(): advice or {} for gene, advice in (table.get("genes") or {}).items()}
        _table = (mtime, table)
    return _table[1]


class Table3(Source):
    def set_entries(self):
        self.entries = {
            ("transcript",): self.transcript,
            ("gene",): self.gene,
        }

    def is_complete(self) -> bool:
        if not os.path.exists(secrets.table3_database):
            print("Could not find Tabel 3:", secrets.table3_database)
            return False
        try:
            get_table()
        except (yaml.YAMLError, AttributeError, TypeError) as e:
            print("Tabel 3 is not in the expected format:", e)
            return False
        return True

    async def transcript(self):
        """
        Looks up the gene of the transcript, so the advice is shown without waiting for another source to supply the gene
        """
        gene = await self.resolve_gene()
        if gene:
            self.report(gene)
        else:
            # this is not complete, so gene() is tried when the meta data has a gene
            self.html_text = "<span class='text-muted'>Gen nog niet bekend...</span>"

    async def gene(self):
        self.report(self.variant["gene"])

    def report(self, symbol):
        """
        Shows the advice of the table for the gene, the text is shown as it is in the table (in Dutch)
        """
        table = get_table()
        gene = table["genes"].get(str(symbol).upper())

        # the name and version are on the link to the table, at the bottom of the card
        name = f"Tabel 3 ({table['version']})" if table.get("version") else "Tabel 3"
        if table.get("url"):
            self.html_title = ""
            self.html_links["main"] = SourceURL(escape(name), table["url"])
        else:
            self.html_subtitle = f"<small>{escape(name)}</small>"
        self.complete = True

        if gene is None:
            # most genes are not in the table, don't use space on the page to say so
            self.html_text = "Gen staat niet in Tabel 3."
            self.found = False
            self.hidden = True
            return

        template = Environment(loader=BaseLoader(), autoescape=True).from_string(TEMPLATE)
        self.html_text = template.render(table=table, gene=gene)

    def get_name(self):
        return "Tabel 3"
