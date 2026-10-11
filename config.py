from typing import List

from pydantic import Field
from pydantic_settings import BaseSettings

from sources import *
from util import get_release_tag, get_release_year


class Settings(BaseSettings):
    debug: bool = True
    port: int = 8080
    release_tag: str = Field(get_release_tag())
    release_year: int = Field(get_release_year())
    disabled_sources: str = Field(default="", env="DISABLED_SOURCES")

    sources: List[type] = [
        # is_banner()=True sources should be at the start, they are shown above the cards
        Table3,
        GenomeNexus,
        InHouse,
        GeneLists,
        Alamut,
        COSMIC,
        ClinVar,
        Franklin,
        CKB,
        TP53,
        LOVD,
        OncoKB,
        cBioPortal,
        Varsome,
        Mutalyzer,
        Mastermind,
        SpliceAI,
        BoostDM,
        HMF,
        # is_hidden()=True sources should be at the end, to not leave gaps in source layout
        Clingen,
        Ensembl,
        # GnomAD,
        dbSNP,
    ]

    # these are shown on top of each other in one column, at the place of the first of them in sources
    stacked_sources: List[type] = [
        Alamut,
        CKB,
        Franklin,
    ]

    def __init__(self):
        super().__init__()
        disabled_sources_list = self.disabled_sources.split(",")
        self.sources = [
            source
            for source in self.sources
            if source.is_complete(source)
            and source.__name__ not in disabled_sources_list
        ]


settings = Settings()
