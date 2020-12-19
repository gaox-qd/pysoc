from datetime import datetime


# The name of this package.
name = "PySOC"
# Brief description.
description = "PySOC spin-orbit coupling calculations"
# Whether this is a development version.
development = False
# Version information.
major_version = 2
minor_version = 1
revision = 4
version_number = "{}.{}.{}".format(major_version, minor_version, revision)
# The full version number of this package.
version = "{}{}".format(version_number, "-dev" if development else "")


# Program date (when we were last updated).
_last_updated_string = "26/11/2020"
last_updated = datetime.strptime(_last_updated_string, "%d/%m/%Y")


# The name of the logger that will be used by pysoc.
logger_name = "PySOC"
