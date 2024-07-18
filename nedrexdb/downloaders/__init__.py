import datetime as _datetime
import shutil as _shutil
from pathlib import Path as _Path

from nedrexdb import config as _config
from nedrexdb.common import Downloader
from nedrexdb.db import MongoInstance
from nedrexdb.downloaders.biogrid import download_biogrid as _download_biogrid
from nedrexdb.downloaders.orphanet import download_orphanet as _download_orphanet


class Version:
    def __init__(self, string):
        self.major, self.minor, self.patch = [int(i) for i in string.split(".")]

    def increment(self, level):
        if level == "major":
            self.major += 1
        elif level == "minor":
            self.minor += 1
        elif level == "patch":
            self.patch += 1

    def __repr__(self):
        return f"{self.major}.{self.minor}.{self.patch}"


def download_all(force=False):
    base = _Path(_config["db.root_directory"])
    download_dir = base / _config["sources.directory"]

    if force and (download_dir).exists():
        _shutil.rmtree(download_dir)
    download_dir.mkdir(exist_ok=True, parents=True)

    sources = _config["sources"]
    # Remove the source keys
    exclude_keys = {"directory", "username", "password"}

    metadata = {"source_databases": {}}

    metadata["source_databases"]["biogrid"] = {"date": f"{_datetime.datetime.now().date()}", "version": None}
    _download_biogrid()
    _download_orphanet()

    for source in filter(lambda i: i not in exclude_keys, sources):
        metadata["source_databases"][source] = {"date": f"{_datetime.datetime.now().date()}", "version": None}

        # Catch case to skip sources with bespoke downloaders.
        if source in {
            "biogrid",
            "drugbank",
            "orphanet",
        }:
            continue
        (download_dir / source).mkdir(exist_ok=True)

        data = sources[source]
        username = data.get("username")
        password = data.get("password")

        for _, download in filter(lambda i: i[0] not in exclude_keys, data.items()):
            url = download.get("url")
            filename = download.get("filename")
            if url is None:
                continue
            if filename is None:
                filename = url.rsplit("/", 1)[1]

            d = Downloader(
                url=url,
                target=download_dir / source / filename,
                username=username,
                password=password,
            )
            d.download()

    docs = list(MongoInstance.DB["metadata"].find())
    if len(docs) == 1:
        version = docs[0]["version"]
    elif len(docs) == 0:
        version = "0.0.0"
    else:
        raise Exception("should only be one document in the metadata collection")

    v = Version(version)
    v.increment("patch")

    metadata["version"] = f"{v}"

    MongoInstance.DB["metadata"].replace_one({}, metadata, upsert=True)
