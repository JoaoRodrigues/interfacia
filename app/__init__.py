#!/usr/bin/env python

import datetime
import json
import os
import shutil
import subprocess

from flask import Flask
from flask import abort, render_template, request, redirect, send_file, url_for
from werkzeug import secure_filename
from celery import Celery

from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

from py.MolProcesser import MolProcesser, MolError

curdir = os.path.abspath(os.curdir)

# Setup Flask App
app = Flask(__name__)
app.config['JOB_DIR'] = os.path.join(curdir, 'app', 'tmp')
app.config['CNS_EXE'] = os.path.join(curdir, 'app', 'bin', 'cns.exe')
app.config['ENERGY_INP'] = os.path.join(curdir, 'app', 'mm', 'calc_pairwise_residue.inp')
app.config['BUILDR_INP'] = os.path.join(curdir, 'app', 'mm', 'builder.inp')
app.config['TOPPAR'] = os.path.join(curdir, 'app', 'mm', 'toppar')
app.config['ALLOWED_EXT'] = set(('pdb', 'ent'))

# Setup Celery
app.config.update(
    CELERY_BROKER_URL='redis://localhost:6379',
    CELERY_RESULT_BACKEND='redis://localhost:6379'
)


celery_app = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery_app.conf.update(app.config)

@celery_app.task
def run_cns_analysis(dpath):
	"""
	Small function to run CNS in the background.
	Probably not scalable.
	"""

	os.chdir(dpath)
	builder = os.path.basename(app.config['BUILDR_INP'])
	builder = os.path.join(dpath, builder)
	energy = os.path.basename(app.config['ENERGY_INP'])
	energy = os.path.join(dpath, energy)

	if not os.path.exists('molecule.pdb') \
		or not os.path.exists(builder) or not os.path.exists(energy) \
		or not os.path.exists(os.path.dirname(app.config['TOPPAR'])):
			return False
	else:
		print("Launching {0}...".format(dpath))
		cmd = '{0} < {1} > builder.log\n{0} < {2} > energy.log'.format(app.config['CNS_EXE'], builder, energy)
		pid = subprocess.Popen(cmd, shell=True)
		status = pid.wait()
		return status

# Views
@app.route('/', methods=['GET', 'POST'])
def home():
	if request.method == 'POST':
		uploaded_pdb = request.files['structure']
		uploaded_ext = uploaded_pdb.filename.split('.')[-1]
		if uploaded_pdb:
			if uploaded_ext in app.config['ALLOWED_EXT']:
				# Validate structure
				try:
					mol = MolProcesser(uploaded_pdb)
				except MolError, e:
					# Flash error messages
					print "oh uh.."
					print e
					return '<h1>Error...</h1>'

				# Create job folder
				date_now = datetime.datetime.now()
				job_name = date_now.strftime('%Y-%m-%d-%H-%M-%S-%f')
				job_dir = os.path.join(app.config['JOB_DIR'], job_name)
				os.mkdir(job_dir)

				# Save PDB file with generic name
				local_fn = os.path.join(job_dir, 'molecule.pdb')
				with open(local_fn, 'w') as handle:
					handle.write(mol.tostring)

				# Copy CNS scripts & FF and send job to Celery
				shutil.copy(app.config['ENERGY_INP'], job_dir)
				shutil.copy(app.config['BUILDR_INP'], job_dir)
				shutil.copytree(app.config['TOPPAR'], os.path.join(job_dir, 'toppar'))
				task = run_cns_analysis.delay(job_dir)

				if not task:
					return "Error..."
				else:
					# print "Hello?"
					# return render_template("index.html", submit=True, job_name=job_name)
					return "Results: <a href={0}>{0}</a>".format(url_for('results', job_id=job_name))

	return render_template("index.html")

@app.route('/results/')
@app.route('/results/<job_id>')
def results(job_id=None):
	if not job_id:
		job_list = []
		for folder in os.listdir(app.config['JOB_DIR']):
			link = "<a href={0}>{0}</a><br />".format(url_for('results', job_id=folder))
			job_list.append(link)
		return ''.join(job_list)

	job_dir = os.path.join(app.config['JOB_DIR'], job_id)
	results_fn = os.path.join(job_dir, 'molecule_cmplt.pwr_ene')
	results_fe = os.path.exists(results_fn)
	if results_fe:
		results_fsize = os.path.getsize(results_fn)
	if not os.path.exists(job_dir):
		return 'Job not found...'
	elif os.path.exists(job_dir) and (not results_fe or not results_fsize):
		if os.path.isfile(os.path.join(job_dir, 'energy.log')):
			return '<pre>{0}</pre>'.format(open(os.path.join(job_dir, 'energy.log')).read())
		elif os.path.isfile(os.path.join(job_dir, 'builder.log')):
			return '<pre>{0}</pre>'.format(open(os.path.join(job_dir, 'builder.log')).read())
		else:
			return 'Job not started yet...'
	else:
		os.chdir(job_dir)

		# Build data structures for D3.js
		# Array with residue names and metadata (chain id, residue number, residue type)
		# 2D Array with interaction strenght per res-res interaction.
		# Order and sizes must match (i.e. 20x1 & 20x20x1)

		# First pass: extract data into graph
		nodes, edges = set(), {}
		with open(results_fn) as energetics:
			for line in energetics:
				if line.startswith('#') or not line.strip():
					continue

				(r_segm, r_resi, r_resn, _, l_segm, l_resi, l_resn,
				 _, _, _, ene_vdw, _, _, ene_elec ) = line.split()

				r_resi, l_resi = int(r_resi), int(l_resi)
				ene_vdw, ene_elec = float(ene_vdw), float(ene_elec)

				nodes.add((r_segm, r_resi, r_resn))
				nodes.add((l_segm, l_resi, l_resn))
				edges[((r_segm, r_resi, r_resn), (l_segm, l_resi, l_resn))] = (ene_vdw, ene_elec)
				edges[((l_segm, l_resi, l_resn), (r_segm, r_resi, r_resn))] = (ene_vdw, ene_elec)

		# Second pass: build JSON data structures
		nodes = sorted(nodes)

		n_nodes = len(nodes)
		# Could be made generic by a series of matrices
		vdw_graph = [[0 for _ in range(n_nodes)] for _ in range(n_nodes)]
		elec_graph = [[0 for _ in range(n_nodes)] for _ in range(n_nodes)]

		for e in edges:
			res_a, res_b = e
			i_a, i_b = nodes.index(res_a), nodes.index(res_b)
			vdw_graph[i_a][i_b] = edges[e][0]
			vdw_graph[i_b][i_a] = edges[e][0]
			elec_graph[i_a][i_b] = edges[e][1]
			elec_graph[i_b][i_a] = edges[e][1]

		nodes = [{'seg': i, 'resi': j, 'resn': k} for (i,j,k) in nodes]

		return render_template('results.html', vdw_graph=json.dumps(vdw_graph),
											   elec_graph=json.dumps(elec_graph),
											   nodes=json.dumps(nodes)
											  )


@app.route('/results/<job_id>/<path:filename>')
def serve_file(job_id=None, filename=None):

	job_dir = os.path.join(app.config['JOB_DIR'], job_id)
	filename = os.path.join(job_dir, filename)
	filename = filename if os.path.exists(filename) else False

	if not filename:
		abort(404)
	else:
		return send_file(filename)

