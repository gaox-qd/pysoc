# -*- mode: python ; coding: utf-8 -*-

import sys
from pathlib import Path
sys.path.insert(0,str(Path("./").resolve()))

#from base import binaries, datas
binaries = []
datas = [('../../pysoc/data', 'pysoc/data')]

script = "../../pysoc/program/main.py"
prog_name = "pysoc.exe"
package_name = "CentOS-6.10"

a = Analysis([script],
			 binaries=binaries,
			 datas=datas,
			 # 'pkg_resources.py2_warn' see https://github.com/pypa/setuptools/issues/1963
			 hiddenimports=['pkg_resources.py2_warn'],
			 hookspath=[],
			 runtime_hooks=[],
			 excludes=["openbabel"],
			 win_no_prefer_redirects=False,
			 win_private_assemblies=False,
			 cipher=None,
			 noarchive=False
)

pyz = PYZ(a.pure, a.zipped_data,
			 cipher=None
)

exe = EXE(pyz,
		a.scripts,
		[],
		exclude_binaries=True,
		name=prog_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True
)

import pysoc
coll = COLLECT(exe,
				a.binaries,
				a.zipfiles,
				a.datas,
				strip=False,
				upx=True,
				upx_exclude=[],
				name="{}.{}.{}".format("pysoc", pysoc.version, package_name)
)
