from pathlib import Path as _Path

from nedrexdb import config as _config


def _get_file_location_factory(database):
    def inner(label):
        path = f"sources.{database}.{label}"
        details = _config[path]

        if "filename" in details:
            filename = details["filename"]
        else:
            filename = details["url"].rsplit("/", 1)[1]

        path = _Path(_config["db.root_directory"]) / _config["sources.directory"] / database / filename
        print(path)

        assert path.exists, f"{path} does not exist"
        return path

    return inner
