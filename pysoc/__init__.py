from datetime import datetime
from contextlib import ExitStack
import atexit
import importlib_resources


# The name of this package.
name = "PySOC"
# Brief description.
description = "PySOC spin-orbit coupling calculations"
# Whether this is a development version.
development = False
# Version information.
major_version = 2
minor_version = 3
revision = 0
version_number = "{}.{}.{}".format(major_version, minor_version, revision)
# The full version number of this package.
version = "{}{}".format(version_number, "-dev" if development else "")


# Program date (when we were last updated).
_last_updated_string = "30/05/2024"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")


# The name of the logger that will be used by pysoc.
logger_name = "PySOC"


def get_resource(name):
    """
    Get a pathlib path object to a package resource.
    """
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    ref = importlib_resources.files('pysoc') / name
    return file_manager.enter_context(importlib_resources.as_file(ref))