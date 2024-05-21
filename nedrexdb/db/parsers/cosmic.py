import gzip as _gzip
from csv import DictReader as _DictReader
from pathlib import Path as _Path

from more_itertools import chunked as _chunked
from tqdm import tqdm as _tqdm

from nedrexdb.db import MongoInstance
from nedrexdb.db.models.edges.variant_affects_gene import VariantAffectsGene
from nedrexdb.db.models.edges.variant_associated_with_disorder import VariantAssociatedWithDisorder
from nedrexdb.db.models.edges.gene_associated_with_disorder import GeneAssociatedWithDisorder
from nedrexdb.db.models.nodes.gene import Gene
from nedrexdb.db.models.nodes.genomic_variant import GenomicVariant
from nedrexdb.db.parsers import _get_file_location_factory
from nedrexdb.logger import logger

get_file_location = _get_file_location_factory("cosmic")
get_clinvar_file_location = _get_file_location_factory("clinvar")


# g_dot_re = re.compile("(..?):g\.(\d*)_?(\d*)(.*)")


def get_gdot2clinvar(fname: str) -> dict[str, str]:
    from nedrexdb.db.parsers.clinvar import ClinVarVCFParser
    vcf_parser = ClinVarVCFParser(fname)
    gdot2clinvar = {}
    for row in vcf_parser.iter_rows():
        full_gdot = row['INFO'].get('CLNHGVS')
        if full_gdot:
            _, gdot = row['INFO']['CLNHGVS'].split(':', maxsplit=1)
            gdot2clinvar[f"{row['CHROM']}:{gdot}"] = f"clinvar.{row['ID']}"
    return gdot2clinvar


def get_cancer2mondo(mapping_fname: _Path) -> dict[tuple: str]:
    mapping_columns = ['SITE_PRIMARY_COSMIC', 'SITE_SUBTYPE1_COSMIC',
                       'SITE_SUBTYPE2_COSMIC', 'SITE_SUBTYPE3_COSMIC', 'HISTOLOGY_COSMIC',
                       'HIST_SUBTYPE1_COSMIC', 'HIST_SUBTYPE2_COSMIC', 'HIST_SUBTYPE3_COSMIC']
    cancer2mondo = {}
    with open(mapping_fname, newline='') as mapping_file:
        reader = _DictReader(mapping_file, delimiter="\t")
        cancer2mondo = {tuple(
            row[column] for column in mapping_columns): row['mapped_curie'] for row in reader}
    return cancer2mondo


class COSMICRow:
    def __init__(self, row):
        self._row = row

    def get_HGVSG(self):
        return self._row["HGVSG"]

    def get_COSMIC(self):
        return f"cosmic.{self._row['GENOMIC_MUTATION_ID']}"

    def get_symbol(self):
        return self._row["Gene name"]

    def get_cancer_tuple(self) -> tuple:
        return tuple(
            self._row[column] for column in ['Primary site', 'Site subtype 1', 'Site subtype 2', 'Site subtype 3',
                                             'Primary histology', 'Histology subtype 1', 'Histology subtype 2',
                                             'Histology subtype 3'])

    def get_mutation_status(self):
        return self._row['Mutation somatic status']

    # def get_variant(self, gdot2clinvar) -> GenomicVariant:
    # match = g_dot_re.search(self._row["HGVSG"])

    # if variant_id:
    #     return GenomicVariant(primaryDomainId=variant_id, domainIds=[cosmic_id], dataSources=['COSMIC'])
    # if match:
    #     chrom, pos_start, pos_end, mut = match.group(1, 2, 3, 4)
    #     pos_start = int(pos_start)
    #     if pos_end: pos_end = int(pos_end)
    #     if mut == 'del':
    #         pos_start -= 1
    #         if not pos_end:
    #             pos_end = pos_start + 1
    #         genomic_variants = chr_pos_type2id[(chrom, pos_start, 'Deletion')]
    #         variant = [variant for variant in genomic_variants if
    #                    len(variant["referenceSequence"]) == pos_end + 1 - pos_start]
    #     elif mut == "dup":
    #         genomic_variants = chr_pos_type2id[(chrom, pos_start, 'Duplication')]
    #         variant = [variant for variant in genomic_variants if
    #                    len(variant["alternativeSequence"]) == pos_end + 2 - pos_start]
    #         if variant:
    #             breakpoint()
    #     elif mut.startswith('ins'):
    #         insertion = mut.strip('ins')
    #         genomic_variants = chr_pos_type2id[(chrom, pos_start, 'Insertion')]
    #         variant = [variant for variant in genomic_variants if
    #                    variant["referenceSequence"] + insertion == variant["alternativeSequence"]]
    #     elif mut.startswith('delins'):
    #         insertion = mut.strip('delins')
    #         if not pos_end:
    #             pos_end = pos_start + 1
    #         genomic_variants = chr_pos_type2id[(chrom, pos_start, 'Indel')]
    #         variant = [variant for variant in genomic_variants if
    #                    insertion == variant["alternativeSequence"] and
    #                    len(variant["referenceSequence"]) == pos_end + 1 - pos_start]
    #     else:
    #         if '>' not in mut:
    #             breakpoint()
    #         mut_from, mut_to = mut.split('>', 1)
    #         genomic_variants = chr_pos_type2id[(chrom, pos_start, 'Single Nucleotide Variant')]
    #         variant = [variant for variant in genomic_variants if
    #                    variant['referenceSequence'] == mut_from and
    #                    variant['alternativeSequence'] == mut_to]
    #     assert len(variant) <= 1, f"More than one matching variant found in Nedrex for {self._row['HGVSG']}"
    #     if variant:
    #         variant = variant[0]
    #         variant["domainIds"].append(cosmic_id)
    #         return GenomicVariant(**variant)
    # return None

    def parse(self, gdot2clinvar: dict[str, str], symbol2entrez: dict[str, str], cancer2mondo: dict[tuple, str]) -> \
            tuple[GenomicVariant, VariantAffectsGene, VariantAssociatedWithDisorder]:
        variant_id = gdot2clinvar.get(self.get_HGVSG())
        genomic_variant = None
        variant_gene = None
        variant_disorder = None
        if variant_id:
            cosmic_id = self.get_COSMIC()
            asserted_by = ["cosmic"]
            genomic_variant = GenomicVariant(primaryDomainId=variant_id, domainIds=[
                                             cosmic_id], dataSources=asserted_by)
            # data_update = genomic_variant.generate_dataSource_update()
            # MongoInstance.DB[GenomicVariant.collection_name].bulk_write([data_update])
            gene_id = symbol2entrez[self.get_symbol()]
            variant_gene = VariantAffectsGene(sourceDomainId=variant_id, targetDomainId=gene_id,
                                              dataSources=asserted_by)
            mondo_id = cancer2mondo.get(self.get_cancer_tuple())
            # if mondo_id:
            variant_disorder = VariantAssociatedWithDisorder(accession=cosmic_id, dataSources=asserted_by,
                                                             sourceDomainId=variant_id,
                                                             targetDomainId=mondo_id,
                                                             reviewStatus=self.get_mutation_status())
            # else:
            #     variant_disorder = None

        return genomic_variant, variant_gene, variant_disorder


class COSMICParser:
    COLUMN_NAMES = ['HGVSG', 'Gene name', 'GENOMIC_MUTATION_ID', 'Mutation somatic status', 'Primary site',
                    'Site subtype 1', 'Site subtype 2',
                    'Site subtype 3',
                    'Primary histology', 'Histology subtype 1', 'Histology subtype 2',
                    'Histology subtype 3']

    def __init__(self, f: _Path):
        self.f = f

        if self.f.name.endswith(".gz") or self.f.name.endswith(".gzip"):
            self.gzipped = True
        else:
            self.gzipped = False

    def parse(self, mapping_fname: _Path):
        if self.gzipped:
            f = _gzip.open(self.f, "rt")
        else:
            f = self.f.open()

        reader = _DictReader(f, delimiter="\t")
        f_dict = [{column: row[column]
                   for column in self.COLUMN_NAMES} for row in reader]
        f.close()

        all_symbols = {row['Gene name'] for row in f_dict}
        symbol2entrez = {gene["approvedSymbol"]: gene["primaryDomainId"] for gene in
                         Gene.find(MongoInstance.DB, {"approvedSymbol": {"$in": list(all_symbols)}})}
        non_approved_symbols = all_symbols - symbol2entrez.keys()
        for symbol in non_approved_symbols:
            genes = [gene["primaryDomainId"]
                     for gene in Gene.find(MongoInstance.DB, {"symbols": symbol})]
            assert len(genes) == 1, f"Multiple genes found for the symbol {symbol}"
            symbol2entrez.update({symbol: genes[0]})
        assert not (non_approved_symbols - symbol2entrez.keys()), \
            f"Not all symbols could be mapped: {non_approved_symbols - symbol2entrez.keys()}"

        # id2genomic_variant = {genomic_variant['primaryDomainId']: genomic_variant for genomic_variant in
        #                        GenomicVariant.find(MongoInstance.DB)}
        # chr_pos_type2id = defaultdict(list)
        # this_id =
        # genomic_variant['primaryDomainId'], genomic_variant
        #     chr_pos_type2id[
        #         (genomic_variant['chromosome'], genomic_variant['position'], genomic_variant['variantType'])].append(
        #         genomic_variant)
        gdot2clinvar = get_gdot2clinvar(
            get_clinvar_file_location("human_data"))
        cancer2mondo = get_cancer2mondo(mapping_fname)

        updates = (COSMICRow(row).parse(
            gdot2clinvar, symbol2entrez, cancer2mondo) for row in f_dict)
        for chunk in _tqdm(_chunked(updates, 10_000), leave=False, desc="Parsing COSMIC"):
            if not chunk:
                continue
            genomic_variant_updates, variant_gene_updates, variant_disorder_updates, gene_disorder_updates = [], [], [], []
            for genomic_variant, variant_gene, variant_disorder in chunk:
                if genomic_variant:
                    genomic_variant_updates.append(
                        genomic_variant.generate_update())
                    variant_gene_updates.append(variant_gene.generate_update())
                    if variant_disorder:
                        variant_disorder_updates.append(
                            variant_disorder.generate_update())
                        gene_disorder = variant_gene.copy()
                        gene_disorder.sourceDomainId = variant_gene.targetDomainId
                        gene_disorder.targetDomainId = variant_disorder.targetDomainId
                        gene_disorder_updates.append(GeneAssociatedWithDisorder(gene_disorder).generate_update())

            for this_collection_name, these_updates in zip([GenomicVariant.collection_name, VariantAffectsGene.collection_name, VariantAssociatedWithDisorder.collection_name, GeneAssociatedWithDisorder.collection_name],
                                                           [genomic_variant_updates, variant_gene_updates, variant_disorder_updates, gene_disorder_updates]):
                bulk_write_results = MongoInstance.DB[this_collection_name].bulk_write(
                    these_updates)
                if bulk_write_results.bulk_api_result['writeErrors'] or bulk_write_results.bulk_api_result['writeConcernErrors']:
                    print(bulk_write_results.bulk_api_result)
                    


def parse_gene_disease_associations():
    logger.info("Parsing COSMIC")
    fname = get_file_location("census")
    mapping_fname = get_file_location("mapping")
    COSMICParser(fname).parse(mapping_fname)
