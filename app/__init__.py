#!/usr/bin/env python

import os
from flask import Flask

curdir = os.path.abspath(os.curdir)

app = Flask(__name__)
app.config['JOB_DIR'] = os.path.join(curdir, 'app', 'tmp')
app.config['CNS_EXE'] = os.path.join(curdir, 'app', 'bin', 'cns.exe')
app.config['CNS_INP'] = os.path.join(curdir, 'app', 'data', 'calc_pairwise_residue.inp')
app.config['TOPPAR'] = os.path.join(curdir, 'app', 'data', 'toppar')
app.config['ALLOWED_EXT'] = set(('pdb', 'ent'))

from app import views