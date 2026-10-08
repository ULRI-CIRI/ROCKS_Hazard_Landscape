"""PubChem API data retrieval and processing utilities."""

import re
from functools import reduce
from urllib.parse import parse_qs, urlparse

import polars as pl
import requests
from pydantic import BaseModel

PUBCHEM_NFPA_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/?heading=NFPA+Hazard+Classification"
PUBCHEM_GHS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/?heading=GHS+Classification"
PUBCHEM_VENDORS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{}/JSON/?heading=Chemical+Vendors"
PUBCHEM_NFPA_IMG_URL = "https://pubchem.ncbi.nlm.nih.gov/image/nfpa.cgi?code={}"

# Signal word severity ordering (higher index = higher severity)
_SIGNAL_PRIORITY = {"": 0, "Warning": 1, "Danger": 2}


# ----------------------------------------- Data models --------------------------------------------
class Reference(BaseModel):
    id: int | None
    source: str
    title: str
    url: str


class GhsPictogram(BaseModel):
    url: str
    label: str


class GhsData(BaseModel):
    compound_name: str = ""
    signal_word: str = ""
    pictograms: list[GhsPictogram] = list()
    hazard_statements: list[str] = list()
    references: list[Reference] = list()

    def has_data(self) -> bool:
        return (
            bool(self.signal_word)
            or len(self.pictograms) > 0
            or len(self.hazard_statements) > 0
            or len(self.references) > 0
        )


class NfpaDiamond(BaseModel):
    health: int
    fire: int
    instability: int
    special: set[str] = set()

    @classmethod
    def from_code(cls, code: str) -> "NfpaDiamond":
        """Create an NfpaDiamond from a PubChem NFPA code string (e.g. '322W' or '322')."""
        if len(code) < 3:
            return cls(health=0, fire=0, instability=0)

        special = code[3:].strip()
        return cls(
            health=_safe_char_to_int(code[0]),
            fire=_safe_char_to_int(code[1]),
            instability=_safe_char_to_int(code[2]),
            special={special} if special else set(),
        )

    @property
    def code(self) -> str:
        return f"{self.health}{self.fire}{self.instability}" + "".join(self.special)

    @property
    def img_url(self) -> str:
        return PUBCHEM_NFPA_IMG_URL.format(self.code)


class NfpaData(BaseModel):
    compound_name: str = ""
    diamonds: list[tuple[int, NfpaDiamond]] = list()
    references: list[Reference] = list()

    @property
    def consensus(self) -> NfpaDiamond | None:
        """Compute the consensus NFPA diamond from max ratings of all per-source diamonds."""
        if not self.has_data():
            return None

        # Get unique special values across all diamonds
        special = set()
        for _, d in self.diamonds:
            for s in d.special:
                special.add(s)

        # Compute the max value for each category (health, fire, instability)
        return NfpaDiamond(
            health=max(d.health for _, d in self.diamonds),
            fire=max(d.fire for _, d in self.diamonds),
            instability=max(d.instability for _, d in self.diamonds),
            special=special,
        )

    def has_data(self) -> bool:
        return len(self.diamonds) > 0 or len(self.references) > 0


class VendorRecord(BaseModel):
    source: str
    url: str
    sid: int
    source_record_id: str
    source_record_url: str


class VendorData(BaseModel):
    records: list[VendorRecord] = list()

    @property
    def _df(self) -> pl.DataFrame:
        """Return vendor records as a Polars DataFrame."""
        return pl.from_dicts([rec.model_dump() for rec in self.records])

    def to_summary(self) -> pl.DataFrame:
        """Collapse vendor records to one row per unique vendor with homepage URLs."""
        if self.is_empty():
            return pl.DataFrame()

        return (
            self._df.select(["source", "url"])
            .unique(subset=["source"], keep="first")
            .sort("source")
        )

    def to_detail(self) -> pl.DataFrame:
        """Return detailed vendor records with SID, name, registry ID, and source record URL."""
        if self.is_empty():
            return pl.DataFrame()

        return self._df.select(
            ["source", "source_record_id", "source_record_url"]
        ).sort("source")

    def is_empty(self) -> bool:
        return len(self.records) == 0


# -------------------------------- Parsing/processing helpers --------------------------------------
def _make_request(url: str) -> tuple[dict | list[dict], str]:
    """Make a GET request to the given URL."""
    response = requests.get(url)
    error = ""
    try:
        response.raise_for_status()
    except requests.HTTPError as e:
        error = f"HTTP error {response.status_code}: {e!s}"

    return response.json(), error


def _get_section_by_heading(
    sections: list[dict], heading: str, return_subsections: bool = True
) -> dict | list[dict]:
    """Return the first section matching the given heading, optionally returning its subsections."""
    section = next((s for s in sections if s.get("TOCHeading", "") == heading), {})
    return section if not return_subsections else section.get("Section", [])


def _deep_get_section_by_headings(sections: list[dict], headings: list[str]) -> dict:
    """Recursively navigate nested sections by a tree of headings."""
    return_subsections = [True] * (len(headings) - 1) + [False]
    return reduce(
        lambda sec, h: _get_section_by_heading(sec, h[0], return_subsections=h[1]),
        zip(headings, return_subsections),
        sections,
    )


def _safe_char_to_int(c: str, fallback: int = 0) -> int:
    """Safely convert a char string to integer or fallback to a provided value (default: 0)."""
    return int(c) if c.isdigit() else fallback


def _nfpa_code_from_url(url: str) -> str:
    """Extract the 'code' query parameter from a PubChem NFPA image URL."""
    codes = parse_qs(urlparse(url).query).get("code", [])
    return "".join(codes) if codes else ""


def _parse_nfpa_info(info: dict) -> tuple[int, NfpaDiamond] | None:
    """Parse a PubChem NFPA information dict into a (reference number, NfpaDiamond) tuple."""
    if info.get("Name", "") != "NFPA 704 Diamond":
        return None

    ref_num = info.get("ReferenceNumber", 0)
    markup = (
        info.get("Value", {}).get("StringWithMarkup", [{}])[0].get("Markup", [{}])[0]
    )
    img_url = markup.get("URL", "")
    code = _nfpa_code_from_url(img_url)
    if not code:
        return None

    return ref_num, NfpaDiamond.from_code(code)


def _parse_ghs_pictograms(blocks: list[dict]) -> dict[str, GhsPictogram]:
    """Parse a single PubChem GHS pictogram markup block into GhsPictogram objects."""
    pictograms = {}
    for block in blocks:
        for m in block.get("Markup", []):
            url = m.get("URL", "")
            label = m.get("Extra", "")
            if url:
                pictograms.update({url: GhsPictogram(url=url, label=label)})

    return pictograms


def _parse_ghs_hazard_statements(blocks: list[dict]) -> dict[str, str]:
    """Parse a list of PubChem GHS hazard statement markup blocks into cleaned strings."""
    statements = {}
    for block in blocks:
        raw = block.get("String", "")
        cleaned = _strip_pct(raw)
        if cleaned:
            statements.update({cleaned[:4]: cleaned})

    return statements


def _parse_ghs_signal_word(blocks: list[dict]) -> str:
    """Parse a list of PubChem GHS signal word markup blocks into a single string."""
    return blocks[0].get("String", "") if blocks else ""


def _parse_all_ghs_info(
    info_list: list[dict],
) -> tuple[str, list[GhsPictogram], list[str]]:
    """Parse all GHS information blocks into signal word, pictograms, and hazard statements."""
    signal_word = ""
    pictogram_set: dict[str, GhsPictogram] = {}
    hazard_statement_set: dict[str, str] = {}

    for info in info_list:
        name = info.get("Name", "")
        markup_blocks = info.get("Value", {}).get("StringWithMarkup", [])

        if name == "Pictogram(s)":
            pictograms = _parse_ghs_pictograms(markup_blocks)
            pictogram_set.update(pictograms)

        elif name == "Signal":
            signal_word_ = _parse_ghs_signal_word(markup_blocks)
            if _SIGNAL_PRIORITY.get(signal_word_, 0) > _SIGNAL_PRIORITY.get(
                signal_word, 0
            ):
                signal_word = signal_word_

        elif name == "GHS Hazard Statements":
            statements = _parse_ghs_hazard_statements(markup_blocks)
            hazard_statement_set.update(statements)

    return (
        signal_word,
        list(pictogram_set.values()),
        sorted(hazard_statement_set.values(), key=lambda s: s[:4]),
    )


def _parse_vendor_records(records: list[dict]) -> list[VendorRecord]:
    """Convert a list of PubChem vendor record dicts to a list of VendorRecord objects."""
    return [
        VendorRecord(
            source=r.get("SourceName", ""),
            url=r.get("SourceURL", ""),
            sid=r.get("SID", 0),
            source_record_id=r.get("RegistryID", ""),
            source_record_url=r.get("SourceRecordURL", ""),
        )
        for r in records
    ]


def _parse_references(references: list[dict]) -> list[Reference]:
    """Convert a list of PubChem reference dicts to a list of Reference objects."""
    return [
        Reference(
            id=r.get("ReferenceNumber", None),
            source=r.get("SourceName", ""),
            title=r.get("Name", ""),
            url=r.get("URL", ""),
        )
        for r in references
    ]


def _strip_pct(text: str) -> str:
    """Remove percentage annotations like '(48.1%)' from hazard statement strings."""
    return re.sub(r"\s*\(\d+\.?\d*%\)", "", text).strip()


# ---------------------------------- Data retrieval methods ----------------------------------------
def get_nfpa_data(cid: int) -> tuple[NfpaData, str]:
    """Retrieve NFPA 704 data for a compound by CID, returning per-source diamonds and consensus."""

    # Do API request
    data, error = _make_request(PUBCHEM_NFPA_URL.format(cid))
    if error:
        return (
            NfpaData(),
            f"Error retrieving NFPA data: {error}",
        )  # Return empty data on error

    # Navigate to NFPA data section
    record = data.get("Record", {})
    compound_name = record.get("RecordTitle", str(cid))
    nfpa_section = _deep_get_section_by_headings(
        record.get("Section", []),
        ["Safety and Hazards", "Hazards Identification", "NFPA Hazard Classification"],
    )

    # Parse NFPA diamonds from the "Information" list, filtering out any that fail to parse
    diamonds = [
        parsed
        for info in nfpa_section.get("Information", [])
        if (parsed := _parse_nfpa_info(info)) is not None
    ]

    # Parse references from the "Reference" list
    references = _parse_references(record.get("Reference", []))

    # Build NfpaData object and return it with an empty error string on success
    return NfpaData(
        compound_name=compound_name,
        diamonds=diamonds,
        references=references,
    ), ""


def get_ghs_data(cid: int) -> tuple[GhsData, str]:
    """Retrieve GHS classification data for a compound by CID."""

    # Do API request
    data, error = _make_request(PUBCHEM_GHS_URL.format(cid))
    if error:
        return (
            GhsData(),
            f"Error retrieving GHS data: {error}",
        )  # Return empty data on error

    # Navigate to GHS data section
    record = data.get("Record", {})
    compound_name = record.get("RecordTitle", str(cid))
    ghs_section = _deep_get_section_by_headings(
        record.get("Section", []),
        ["Safety and Hazards", "Hazards Identification", "GHS Classification"],
    )

    # Parse all GHS information blocks into signal word, pictograms, and hazard statements
    signal_word, pictograms, hazard_statements = _parse_all_ghs_info(
        ghs_section.get("Information", [])
    )

    # Parse references from the "Reference" list
    references = _parse_references(record.get("Reference", []))

    # Build GhsData object and return it with an empty error string on success
    return GhsData(
        compound_name=compound_name,
        signal_word=signal_word,
        pictograms=pictograms,
        hazard_statements=hazard_statements,
        references=references,
    ), ""


def get_vendor_data(cid: int) -> tuple[VendorData, str]:
    """Retrieve chemical vendor records for a compound by CID."""
    data, error = _make_request(PUBCHEM_VENDORS_URL.format(cid))
    if error:
        return (
            VendorData(),
            f"Error retrieving vendor data: {error}",
        )  # Return empty data on error

    categories = data.get("SourceCategories", {}).get("Categories", [])
    vendors_category = next(
        (c for c in categories if c.get("Category", "") == "Chemical Vendors"), {}
    )
    if not vendors_category:
        return VendorData(), ""  # Return empty data if no vendors category found

    vendors = vendors_category.get("Sources", [])
    return VendorData(
        records=_parse_vendor_records(vendors)
    ), ""  # Return empty error string on success


def merge_references(
    nfpa_refs: list[Reference], ghs_refs: list[Reference]
) -> pl.DataFrame:
    """Combine and deduplicate references from NFPA and GHS responses and build dataframe."""
    combined = pl.from_dicts([ref.model_dump() for ref in nfpa_refs + ghs_refs])
    if combined.is_empty():
        return combined

    return combined.unique(subset=["id", "title", "url"]).sort("id")
