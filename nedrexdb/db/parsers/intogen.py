import gzip as _gzip
from csv import DictReader as _DictReader
from itertools import chain as _chain
from pathlib import Path as _Path

import requests
from more_itertools import chunked as _chunked
from tqdm import tqdm as _tqdm

from nedrexdb.db import MongoInstance
from nedrexdb.db.models.edges.gene_associated_with_disorder import GeneAssociatedWithDisorder
from nedrexdb.db.models.nodes.gene import Gene
from nedrexdb.db.parsers import _get_file_location_factory

get_file_location = _get_file_location_factory("intogen")


def biomart_symbol_transcript_to_entrez(symbol_list: list[str], filter_by: str = "hgnc_symbol", batch_size: int = 100):
    import xml.etree.ElementTree as ET
    query = ET.Element("Query", virtualSchemaName="default", formatter="CSV", header="0", uniqueRows="0", count="",
                       datasetConfigVersion="0.6")
    dataset = ET.SubElement(
        query, "Dataset", name="hsapiens_gene_ensembl", interface="default")
    ET.SubElement(dataset, "Filter", name=filter_by, value="{tr_ids}")
    ET.SubElement(dataset, "Attribute", name=filter_by)
    ET.SubElement(dataset, "Attribute", name="entrezgene_id")
    tree = ET.ElementTree(query)
    xml_string = f'<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>{
        ET.tostring(tree.getroot(), encoding="unicode")}'

    symbol2entrez: dict[str, str] = dict()
    for i in range(0, len(symbol_list), batch_size):
        response = requests.get(
            f'http://www.ensembl.org/biomart/martservice?query={format(xml_string.format(tr_ids=",".join(symbol_list[i:i + batch_size])))}')
        response.raise_for_status()
        symbol2entrez.update(dict(row.split(',', 1)
                             for row in response.content.decode('utf-8').splitlines()))
    return symbol2entrez


class IntOGenRow:
    def __init__(self, row):
        self._row = row

    def parse(self, intogen2mondo: dict[str, list[str]], symbol2entrez: dict[str, str]) -> list[
            GeneAssociatedWithDisorder]:
        sourceDomainId = symbol2entrez[self._row["SYMBOL"]]
        asserted_by = ["intogen"]
        disorders = intogen2mondo[self._row["CANCER_TYPE"]]

        gawds = [
            GeneAssociatedWithDisorder(
                sourceDomainId=sourceDomainId, targetDomainId=disorder.replace(
                    "MONDO:", "mondo."),
                dataSources=asserted_by
            )
            for disorder in disorders
        ]

        return gawds


class IntOGenParser:
    def __init__(self, f: _Path, mapping: _Path):
        self.f = f

        if self.f.name.endswith(".gz") or self.f.name.endswith(".gzip"):
            self.gzipped = True
        else:
            self.gzipped = False

        import json
        self.intogen2mondo = json.load(open(mapping))["mondo_id"]

    def parse(self):
        if self.gzipped:
            f = _gzip.open(self.f, "rt")
        else:
            f = self.f.open()

        reader = _DictReader(f, delimiter="\t")
        f_dict = [{"SYMBOL": row['SYMBOL'], "CANCER_TYPE": row['CANCER_TYPE']}
                  for row in reader]

        symbol2entrez = {gene["approvedSymbol"]: gene["primaryDomainId"]
                         for gene in Gene.find(MongoInstance.DB)}

        updates = (IntOGenRow(row).parse(self.intogen2mondo, symbol2entrez)
                   for row in f_dict)
        for chunk in _tqdm(_chunked(updates, 1_000), leave=False, desc="Parsing IntOGen"):
            chunk = list(_chain(*chunk))
            chunk = [gawd.generate_update() for gawd in chunk]

            if not chunk:
                continue

            MongoInstance.DB[GeneAssociatedWithDisorder.collection_name].bulk_write(
                chunk)

        f.close()


def parse_gene_disease_associations():
    fname = get_file_location("drivers")
    mapping_fname = get_file_location("mapping")
    IntOGenParser(fname, mapping_fname).parse()
