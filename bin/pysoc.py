#!/usr/bin/env python2

''' run soc calculation
'''
import re
import sys
import subprocess
import string
import fileinput

def shell_cmd(cmd):
   out = subprocess.call(cmd, shell=True)
   return out

# where are the soc scripts: soc.py and soc_td
scrip_soc = '/fsnfs/users/xinggao/work/gsh/thiothymine/gtsh/test_python/soc_tb/bin'

#please check and set the control parameters in init.py before run pysoc
shell_cmd(scrip_soc+'/soc.py')
shell_cmd(scrip_soc+'/soc_td')





