import subprocess as _sp
from pathlib import Path as _Path
import shutil as _shutil

from requests.exceptions import HTTPError as _HTTPError

from nedrexdb import config as _config
from nedrexdb.common import Downloader, change_directory
from nedrexdb.logger import logger

def download_orphanet():
    root = _Path(_config["db.root_directory"]) / _config["sources.directory"]
    target_dir = root / "orphanet"
    target_dir.mkdir(exist_ok=True, parents=True)

    # FIRST, download the mapping (zip file)
    orphanet = _config["sources.orphanet"]
    orphanet_mapping = orphanet["mapping"]

    url_mapping = orphanet_mapping["url"]
    zip_fname = target_dir / "Orphanet_Nomenclature_Pack_EN.zip"
    unzip_fname = target_dir / "Orphanet_Nomenclature_Pack_EN"
    target_mapping_fname = (unzip_fname / orphanet_mapping["filename"]).resolve()
 
    d = Downloader(
        url=url_mapping,
        target=zip_fname,
        username=None,
        password=None,
    )
    print("Downloader step passed", zip_fname)
    try:
        d.download()
    except _HTTPError as E:
        logger.warning(f"Unable to download DrugBank: {E}")
        return
    files = list(target_dir.iterdir())
    print(files)
    print("download step passed")
    # Unzip the zip
    with change_directory(target_dir):
        _sp.call(["unzip", f"{zip_fname.resolve()}"])
        zip_fname.unlink()

    # Move the target file from the unzipped directory to the desired target directory
    file = unzip_fname / orphanet_mapping["filename"]
    target_file_path = target_dir / orphanet_mapping["filename"]
    file.rename(target_file_path)

    # Delete the unzipped directory
    _shutil.rmtree(unzip_fname)
    

    # SECOND, download the data
    orphanet = _config["sources.orphanet"]
    orphanet_data = orphanet["data"]

    url_data = orphanet_data["url"]
    # zip_fname = target_dir / "all.zip"
    target_fname = (target_dir / orphanet_data["filename"]).resolve()

    d = Downloader(
        url=url_data,
        target=target_fname,
        username=None,
        password=None,
    )
    try:
        d.download()
    except _HTTPError as E:
        logger.warning(f"Unable to download DrugBank: {E}")
        return

    