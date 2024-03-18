import os as _os
import re as _re
import shutil as _shutil
import subprocess as _subprocess
from pathlib import Path as _Path
from urllib.request import urlretrieve as _urlretrieve

import requests  # type: ignore
from bs4 import BeautifulSoup

from nedrexdb import config as _config
from nedrexdb.common import change_directory as _cd
from nedrexdb.exceptions import (
    AssumptionError as _AssumptionError,
    ProcessError as _ProcessError,
)
from nedrexdb.logger import logger


def get_latest_biogrid_version() -> str:
    url = "https://wiki.thebiogrid.org/doku.php/statistics"
    response = requests.get(url)
    if response.status_code != 200:
        raise _ProcessError("got non-zero status code while scraping BioGRID to get latest version")

    soup = BeautifulSoup(response.text, features="lxml")
    result = soup.find_all(text=_re.compile("^Current Build Statistics"))
    text = str(result[0])
    version = text.replace("(", ")").split(")")[1]
    return version


def download_biogrid():
    version = get_latest_biogrid_version()

    url = (
        "https://downloads.thebiogrid.org/"
        "Download/BioGRID/Release-Archive/"
        f"BIOGRID-{version}/"
        f"BIOGRID-ORGANISM-{version}.tab3.zip"
    )

    zip_fname = url.rsplit("/", 1)[1]
    target_fname = _config.get("sources.biogrid.human_data.filename")
    biogrid_dir = _Path(_config.get("db.root_directory")) / _config.get("sources.directory") / "biogrid"

    biogrid_dir.mkdir(exist_ok=True, parents=True)

    with _cd(biogrid_dir):
        # NOTE: we have to remove the old file first, otherwise we get a bug where it gets deleted
        #       (because iterdir() will delete it).
        if _os.path.isfile(target_fname):
            _os.remove(target_fname)

        logger.debug("Downloading BioGRID v%s" % version)
        _urlretrieve(url, zip_fname)
        _subprocess.call(
            ["unzip", zip_fname],
            stdout=_subprocess.DEVNULL,
            stderr=_subprocess.DEVNULL,
        )
        _os.remove(zip_fname)

        counter = 0
        for f in _Path().iterdir():
            if "Homo_sapiens" in f.name:
                counter += 1
                _shutil.move(f.name, target_fname)
            else:
                _os.remove(f)

    if counter != 1:
        raise _AssumptionError("more than one BioGRID file containing 'Homo_sapiens' was found")
    
