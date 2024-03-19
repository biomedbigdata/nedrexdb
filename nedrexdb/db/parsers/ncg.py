import gzip as _gzip
from collections import defaultdict as _defaultdict
from csv import DictReader as _DictReader
from itertools import chain as _chain
from pathlib import Path as _Path

from more_itertools import chunked as _chunked
from tqdm import tqdm as _tqdm

from nedrexdb.db import MongoInstance
from nedrexdb.db.models.edges.gene_associated_with_disorder import GeneAssociatedWithDisorder
from nedrexdb.db.models.nodes.disorder import Disorder
from nedrexdb.db.models.nodes.gene import Gene
from nedrexdb.db.parsers import _get_file_location_factory

get_file_location = _get_file_location_factory("ncg")


def _umls_to_nedrex_map() -> dict[str, list[str]]:
    d = _defaultdict(list)

    for dis in Disorder.find(MongoInstance.DB):
        umls_ids = [acc for acc in dis["domainIds"] if acc.startswith("umls.")]
        for umls_id in umls_ids:
            d[umls_id].append(dis["primaryDomainId"])

    return d


class NCGRow:
    def __init__(self, row):
        self._row = row

    def get_gene_id(self):
        return f"entrez.{self._row['entrez'].strip()}"

    def parse(self, ncg2mondo: dict[str, list[str]]) -> list[GeneAssociatedWithDisorder]:
        sourceDomainId = self.get_gene_id()
        asserted_by = ["ncg"]
        print(ncg2mondo)
        disorders = ncg2mondo[self._row["cancer_type"]]

        gawds = [
            GeneAssociatedWithDisorder(
                sourceDomainId=sourceDomainId, targetDomainId=disorder.replace(
                    "MONDO:", "mondo."),
                dataSources=asserted_by
            )
            for disorder in disorders
        ]

        return gawds


class NCGParser:
    def __init__(self, f: _Path, mapping: _Path):
        self.f = f

        if self.f.name.endswith(".gz") or self.f.name.endswith(".gzip"):
            self.gzipped = True
        else:
            self.gzipped = False

        import json
        self.ncg2mondo = json.load(open(mapping))["mondo_id"]

    def parse(self):
        if self.gzipped:
            f = _gzip.open(self.f, "rt")
        else:
            f = self.f.open()

        reader = _DictReader(f, delimiter="\t")

        genes = {gene["primaryDomainId"]
                 for gene in Gene.find(MongoInstance.DB)}

        updates = (NCGRow(row).parse(self.ncg2mondo) for row in reader)
        for chunk in _tqdm(_chunked(updates, 1_000), leave=False, desc="Parsing NCG"):
            chunk = list(_chain(*chunk))
            chunk = [gawd.generate_update()
                     for gawd in chunk if gawd.sourceDomainId in genes]

            if not chunk:
                continue

            MongoInstance.DB[GeneAssociatedWithDisorder.collection_name].bulk_write(
                chunk)

        f.close()


def parse_gene_disease_associations():
    fname = "/nfs/data3/nedrex_data/downloads/ncg/NCG_cancerdrivers_systemslevelproperties.tsv" #get_file_location("annotation")
    mapping_fname =  "/nfs/data3/nedrex_data/downloads/ncg/ncg2mondo.json" #get_file_location("mapping")
    NCGParser(fname, mapping_fname).parse()

parse_gene_disease_associations()