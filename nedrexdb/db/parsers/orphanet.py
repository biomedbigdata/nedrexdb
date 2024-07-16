from collections import defaultdict as _defaultdict
from pathlib import Path as _Path
import xml.etree.cElementTree as _et
import openpyxl as _openpyxl
import zipfile as _zipfile


from nedrexdb.logger import logger
from nedrexdb.db.models.nodes.disorder import Disorder
from nedrexdb.db.models.nodes.gene import Gene
from nedrexdb.db.models.edges.gene_associated_with_disorder import (
    GeneAssociatedWithDisorder,
)

orphanet_path = "https://github.com/Orphanet/Orphadata_aggregated/blob/master/Genes%20associated%20with%20rare%20diseases/en_product*.xml"
nomenclature_pack = "https://github.com/Orphanet/Orphadata_aggregated/blob/fcb68b2c4b26ace67cc068a2569015c5b156e291/Rare%20diseases%20and%20classifications/Orphanet%20nomenclature%20files%20for%20coding/Orphanet_Nomenclature_Pack_EN.zip"

class OrphanetParser:
    def __init__(self, orphanet_path, nomenclature_pack):
        # file for disorder-genes associations
        self.associations_path = _Path(orphanet_path)
        if not self.associations_path.exists():
            raise Exception(f"{self.associations_path} does not exist")
        # file for Orphacode-icd10 mapping
        self.nomenclature_path = _Path(nomenclature_pack)
        if not self.nomenclature_path.exists():
            raise Exception(f"{self.nomenclature_path} does not exist")


    def get_dict_icd10_mondo(self):
        # get the mapping ICD10 to MONDO from the existing disorders
        icd10_mondo = _defaultdict(list)
        for item in Disorder.objects():
            for icd10_ids in [i for i in item.icd10]:
                for icd10_id in icd10_ids:
                    icd10_mondo[icd10_id].append(item.primaryDomainId)

        return icd10_mondo


    def get_dict_OrphaCode_icd10(self):
        target_file_name = "ORPHAnomenclature_MasterFile_*"
        try:
            with _zipfile.ZipFile(self.nomenclature_path, 'r') as zip_ref:
                # Check if the target file is present in the zip archive
                if target_file_name not in zip_ref.namelist():
                    print(f"The file '{target_file_name}' was not found in the zip archive.")
                # Read the content of the specified file
                with zip_ref.open(target_file_name, 'r') as file:
                    orpha_icd10 = _defaultdict(list)
                    workbook = _openpyxl.load_workbook(file)
                    sheet = workbook.active
                    # Get header names
                    headers = [cell.value for cell in sheet[1]]
                    for row in sheet.iter_rows(min_row=2, values_only=True):
                        # Assuming the column names are 'orpha' and 'ids'
                        row_data = dict(zip(headers, row))
                        # Splitting multiple Orpha codes
                        orpha = row_data['ORPHAcode']
                        icd10 = row_data['ICDcodes']
                        if icd10 != None:
                            orpha_icd10[orpha].append(icd10)

        except _zipfile.BadZipFile:
            print(f"The file '{self.nomenclature_path}' is not a valid zip file.")

        return orpha_icd10


    def get_OrphaCode(self):
        depth = 0
        OrphaCode_list = []

        for event, elem in _et.iterparse(self.associations_path, events=["start", "end"]):
            if not elem.tag == "OrphaCode":
                continue
            if event == "start":
                depth += 1
            if event == "end":
                depth -= 1
            if depth == 0 and event == "end":
                OrphaCode_list.append(elem.text)

        return OrphaCode_list


    def get_genes(self):
        list_genes = []
        depth = 0

        for event, elem in _et.iterparse(self.associations_path, events=["start", "end"]):
            if not elem.tag == "DisorderGeneAssociationList":
                continue
            if event == "start":
                depth += 1
            if event == "end":
                depth -= 1
            if depth == 0 and event == "end":
                genes = []
                if elem.attrib["count"] != 0:
                    for i in elem.iter('Symbol'):
                        genes.append(i.text)
                list_genes.append(genes)

        return list_genes


    @logger.catch
    
    def parse(self):

        logger.info("Parsing OrphaNet")
        logger.info("\tParsing disorder-gene associations from OrphaNet")

        orpha_icd10 = self.get_dict_OrphaCode_icd10()
        icd10_mondo = self.get_dict_icd10_mondo()

        ordered_OrphaCode = self.get_OrphaCode()
        ordered_associatedGenes = self.get_genes()

        dict_disorder_genes = {}
        for i in range(len(ordered_OrphaCode)):
            icd10 = orpha_icd10[ordered_OrphaCode[i]]
            mondo = icd10_mondo[icd10]
            dict_disorder_genes[mondo] = ordered_associatedGenes[i]

        for d, gs in dict_disorder_genes.items():
            for g in gs:
                # Query to see if a relationship is already recorded.
                gawd = GeneAssociatedWithDisorder.objects(
                    sourceDomainId = g,
                    targetDomainId = d
                )
                if not gawd:
                    # Create it.
                    gawd = GeneAssociatedWithDisorder(
                        sourceDomainId = g,
                        targetDomainId = d,
                        assertedBy = ["orphanet"]
                    )
                    gawd.save()
                else:
                    # Check that there is only one result.
                    assert len(gawd) == 1
                    gawd = gawd[0]
                    # Update
                    gawd.modify(
                        add_to_set__assertedBy = "orphanet"
                        )
