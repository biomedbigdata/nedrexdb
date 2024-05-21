import gzip as _gzip
import re as _re
import sys as _sys
import itertools as _itertools
from csv import (
    DictReader as _DictReader,
    field_size_limit as _field_size_limit,
)

from Bio import SeqIO as _SeqIO, SeqRecord as _SeqRecord
from more_itertools import chunked as _chunked
from pymongo import UpdateOne
from tqdm import tqdm as _tqdm

from nedrexdb.db import MongoInstance
from nedrexdb.db.parsers import _get_file_location_factory
from nedrexdb.db.models.nodes.gene import Gene
from nedrexdb.db.models.nodes.protein import Protein
from nedrexdb.db.models.edges.protein_encoded_by_gene import (
    ProteinEncodedByGene,
)

get_file_location = _get_file_location_factory("uniprot")

_field_size_limit(_sys.maxsize)


class UniProtRecord:
    _CURLY_REGEX = _re.compile(r"{|}")
    _DESCRIPTION_CUTOFF_STRINGS = ["Contains:", "Includes:", "Flags:"]
    _CATEGORY_FIELDS = ["RecName:", "AltName:", "SubName:", ""]
    _SUBCATEGORY_FIELDS = [
        "Full=",
        "Short=",
        "EC=",
        "Allergen=",
        "Biotech=",
        "CD_antigen=",
        "INN=",
    ]
    _COMBINED_FIELDS = sorted(
        (" ".join(i) for i in _itertools.product(_CATEGORY_FIELDS, _SUBCATEGORY_FIELDS)),
        key=lambda i: len(i),
        reverse=True,
    )
    _COMBINATION_REGEX = _re.compile(r"|".join(_COMBINED_FIELDS))

    def __init__(self, record: _SeqRecord.SeqRecord):
        self._record = record

    def get_primary_id(self) -> str:
        return f"uniprot.{self._record.id}"

    def get_sequence(self) -> str:
        return str(self._record.seq)

    def get_display_name(self) -> str:
        return self._record.name

    def get_taxid(self) -> int:
        taxid = self._record.annotations["ncbi_taxid"][0]
        return int(taxid)

    def get_synonyms(self) -> list[str]:
        synonyms = self._record.description.split()
        cutoff = next(
            (val for val, item in enumerate(synonyms) if item in self._DESCRIPTION_CUTOFF_STRINGS),
            999_999,
        )
        synonyms = " ".join(synonyms[:cutoff])
        synonyms = self._COMBINATION_REGEX.split(synonyms)
        synonyms = [i.strip() for i in synonyms]
        synonyms = [i[:-1] if i.endswith(";") else i for i in synonyms]
        synonyms = [i for i in synonyms if i]
        return synonyms

    def get_gene_name(self) -> str:
        gene_name = self._record.annotations.get("gene_name", "")
        if not gene_name:
            pass
        else:
            if isinstance(gene_name, list):
                name = None
                for i in gene_name:
                    for k,v in i.items():
                        if not name:
                            name = v if isinstance(v, str) else v[0]
                        if k == "Name":
                            name = v
                gene_name = name
            if gene_name.startswith("Name="):
                gene_name = gene_name.replace("Name=", "").split(";", 1)[0]
                gene_name = self._CURLY_REGEX.split(gene_name)[0].strip()

        return gene_name


    def get_comments(self) -> str:
        return self._record.annotations.get("comment", "")

    def parse(self):
        p = Protein()

        p.primaryDomainId = self.get_primary_id()
        p.domainIds.append(p.primaryDomainId)

        p.displayName = self.get_display_name()
        p.synonyms = self.get_synonyms()
        p.comments = self.get_comments()
        p.geneName = self.get_gene_name()

        p.taxid = self.get_taxid()
        p.sequence = self.get_sequence()

        p.dataSources = ["uniprot"]

        return p


class IDMapRow:
    def __init__(self, row):
        self._row = row

    def get_source_domain_id(self) -> str:
        return f"uniprot.{self._row['UniProtKB-AC']}"

    def get_target_domain_ids(self) -> list[str]:
        genes = [f"entrez.{acc}" for acc in self._row["GeneID (EntrezGene)"].split(";") if acc.strip()]
        return genes

    def parse(self):
        pebg = ProteinEncodedByGene()
        pebg.sourceDomainId = self.get_source_domain_id()
        pebg.dataSources = ["uniprot"]

        genes = self.get_target_domain_ids()
        for gene in genes:
            pebg.targetDomainId = gene
            yield pebg


def _iter_gzipped_swiss(fname):
    with _gzip.open(fname, "rt") as f:
        for record in _SeqIO.parse(f, "swiss"):
            yield record


def parse_proteins():
    filenames = [get_file_location("trembl"), get_file_location("swissprot")]
    uniprot_records = _itertools.chain(*[_iter_gzipped_swiss(filename) for filename in filenames])
    updates = (UniProtRecord(record).parse().generate_update() for record in uniprot_records)

    for chunk in _tqdm(
        _chunked(updates, 1_000),
        desc="Parsing Swiss-Prot and TrEMBL",
        leave=False,
    ):
        MongoInstance.DB[Protein.collection_name].bulk_write(chunk)


def parse_idmap():
    fieldnames = (
        "UniProtKB-AC",
        "UniProtKB-ID",
        "GeneID (EntrezGene)",
        "RefSeq",
        "GI",
        "PDB",
        "GO",
        "UniRef100",
        "UniRef90",
        "UniRef50",
        "UniParc",
        "PIR",
        "NCBI-taxon",
        "MIM",
        "UniGene",
        "PubMed",
        "EMBL",
        "EMBL-CDS",
        "Ensembl",
        "Ensembl_TRS",
        "Ensembl_PRO",
        "Additional PubMed",
    )
    filename = get_file_location("idmapping")

    comment_char = "#"
    delimiter = "\t"

    gene_ids = {doc["primaryDomainId"] for doc in Gene.find(MongoInstance.DB)}
    protein_ids = {doc["primaryDomainId"] for doc in Protein.find(MongoInstance.DB)}

    # Pass one to get the relationships between NCBI genes and UniProt proteins
    with _gzip.open(filename, "rt") as f:
        filtered_f = filter(lambda row: row[0] != comment_char, f)
        reader = _DictReader(filtered_f, delimiter=delimiter, fieldnames=fieldnames)

        def record_gen():
            for row in reader:
                for pebg in IDMapRow(row).parse():
                    if pebg.sourceDomainId not in protein_ids:
                        continue
                    if pebg.targetDomainId not in gene_ids:
                        continue

                    yield pebg

        updates = (pebg.generate_update() for pebg in record_gen())
        for chunk in _tqdm(
            _chunked(updates, 1_000),
            desc="Parsing UniProt ID map (NCBI relationships)",
            leave=False,
        ):
            MongoInstance.DB[ProteinEncodedByGene.collection_name].bulk_write(chunk)

    # Pass two to add the ENSEMBL protein IDs
    with _gzip.open(filename, "rt") as f:
        filtered_f = filter(lambda row: row[0] != comment_char, f)
        reader = _DictReader(filtered_f, delimiter=delimiter, fieldnames=fieldnames)

        def update_gen():
            for row in reader:
                pdid = f"uniprot.{row['UniProtKB-AC']}"
                ensembl_ids = [f"ensembl.{i.strip()}" for i in row["Ensembl_PRO"].split(";") if i.strip()]
                for ensembl_id in ensembl_ids:
                    yield UpdateOne({"primaryDomainId": pdid}, {"$addToSet": {"domainIds": ensembl_id}}, upsert=False)

        for chunk in _tqdm(
            _chunked(update_gen(), 1_000), desc="Parsing UniProt ID map (ENSEMBL relationships)", leave=False
        ):
            MongoInstance.DB[Protein.collection_name].bulk_write(chunk)
