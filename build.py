#!/usr/bin/env python

import time

import click

import nedrexdb
from nedrexdb import config, downloaders
from nedrexdb.control.docker import NeDRexDevInstance, NeDRexLiveInstance
from nedrexdb.db import MongoInstance, mongo_to_neo, collection_stats, update_db_version
from nedrexdb.db.parsers import (
    biogrid,
    disgenet,
    ctd,
    drugbank,
    drug_central,
    hpo,
    hpa,
    iid,
    intact,
    uniprot,
    uniprot_signatures,
    ncbi,
    mondo,
    omim,
    reactome,
    go,
    clinvar,
    chembl,
    unichem,
    bioontology,
    sider,
    uberon,
    repotrial,
    cosmic,
    ncg,
    intogen,
    orphanet
)
from nedrexdb.post_integration import trim_uberon, drop_empty_collections


@click.group()
def cli():
    pass


@click.option("--conf", required=True, type=click.Path(exists=True))
@click.option("--download", is_flag=True, default=False)
@cli.command()
def update(conf, download):
    print(f"Config file: {conf}")
    print(f"Download updates: {download}")

    nedrexdb.parse_config(conf)

    version = config["db.version"]

    if version not in ["open", "licensed"]:
        raise Exception(f"invalid version {version!r}")

    dev_instance = NeDRexDevInstance()
    dev_instance.remove()
    dev_instance.set_up(use_existing_volume=False, neo4j_mode="import")
    MongoInstance.connect("dev")
    MongoInstance.set_indexes()

    if download:
        downloaders.download_all()

    # Parse sources contributing only nodes (and edges amongst those nodes)
    go.parse_go()
    mondo.parse_mondo_json()
    ncbi.parse_gene_info()
    uberon.parse()
    uniprot.parse_proteins()

    cosmic.parse_gene_disease_associations()
    ncg.parse_gene_disease_associations()

    # Sources that add node type but require existing nodes, too
    clinvar.parse()

    if version == "licensed":
        drugbank._parse_drugbank()  # requires proteins to be parsed first
    elif version == "open":
        drugbank.parse_drugbank()
        chembl.parse_chembl()

    uniprot_signatures.parse()  # requires proteins to be parsed first
    hpo.parse()  # requires disorders to be parsed first
    reactome.parse()  # requires protein to be parsed first
    bioontology.parse()  # requires phenotype to be parsed

    # Sources that add data to existing nodes
    drug_central.parse_drug_central()
    unichem.parse()
    repotrial.parse()

    # Sources adding edges.
    biogrid.parse_ppis()
    ctd.parse()
    disgenet.parse_gene_disease_associations()
    go.parse_goa()
    hpa.parse_hpa()
    iid.parse_ppis()
    intact.parse()

    if version == "licensed":
        omim.parse_gene_disease_associations()

    sider.parse()
    uniprot.parse_idmap()
    
    intogen.parse_gene_disease_associations()
    orphanet.parse_gene_disease_associations()    

    from nedrexdb.analyses import molecule_similarity

    molecule_similarity.run()

    # Post-processing
    trim_uberon.trim_uberon()
    drop_empty_collections.drop_empty_collections()

    # export to Neo4j
    mongo_to_neo.mongo_to_neo(dev_instance, MongoInstance.DB)

    # Profile the collections
    collection_stats.profile_collections(MongoInstance.DB)
    update_db_version.update_db_version(default_version="2.0.0")
    time.sleep(60)

    # remove dev instance and set up live instance
    dev_instance.remove()
    live_instance = NeDRexLiveInstance()
    live_instance.remove()
    live_instance.set_up(use_existing_volume=True, neo4j_mode="db")


@click.option("--conf", required=True, type=click.Path(exists=True))
@cli.command()
def restart_live(conf):
    print(f"Config file: {conf}")
    nedrexdb.parse_config(conf)

    live_instance = NeDRexLiveInstance()
    live_instance.remove()
    live_instance.set_up(use_existing_volume=True, neo4j_mode="db")


if __name__ == "__main__":
    cli()
