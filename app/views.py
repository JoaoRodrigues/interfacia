import datetime
import json
import os
import shutil
import subprocess

from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

from app import app
from app.py import MolProcesser
from flask import abort, render_template, request, redirect, send_file, url_for
from werkzeug import secure_filename

def run_cns_analysis(dpath):
	"""Small function to run CNS asynchronously"""

	os.chdir(dpath)
	cns_inp = os.path.basename(app.config['CNS_INP'])
	cns_inp = os.path.join(dpath, cns_inp)

	if not os.path.exists('molecule.pdb') \
		or not os.path.exists(cns_inp) \
		or not os.path.exists(os.path.dirname(app.config['TOPPAR'])):
			return False
	else:
		cmd = '{0} < {1} > cns.log'.format(app.config['CNS_EXE'], cns_inp)
		pid = subprocess.Popen(cmd, shell=True)
		return True

@app.route('/', methods=['GET', 'POST'])
def home():
	if request.method == 'POST':
		uploaded_pdb = request.files['pdbfile']
		uploaded_ext = uploaded_pdb.filename.split('.')[-1]
		if uploaded_pdb:
			if uploaded_ext in app.config['ALLOWED_EXT']:
				# Create job folder
				date_now = datetime.datetime.now()
				job_name = date_now.strftime('%Y-%m-%d-%H-%M-%S-%f')
				job_dir = os.path.join(app.config['JOB_DIR'], job_name)
				os.mkdir(job_dir)

				# Save PDB file with generic name
				local_fn = os.path.join(job_dir, 'molecule.pdb')
				uploaded_pdb.save(local_fn)

				# Copy CNS scripts & FF and send job to Celery
				shutil.copy(app.config['CNS_INP'], job_dir)
				shutil.copytree(app.config['TOPPAR'], os.path.join(job_dir, 'toppar'))
				status = run_cns_analysis(job_dir)
				if not status:
					return "Error..."
				else:
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
	results_fn = os.path.join(job_dir, 'molecule.pwr_ene')
	results_fe = os.path.exists(results_fn)
	if results_fe:
		results_fsize = os.path.getsize(results_fn)
	if not os.path.exists(job_dir):
		return 'Job not found...'
	elif os.path.exists(job_dir) and (not results_fe or not results_fsize):
		return 'Job Running...'
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
				if not line.startswith('# '):
					continue

				(_, r_segm, r_resi, r_resn, _, l_segm, l_resi, l_resn,
				 _, _, _, ene_vdw, _, _, ene_elec, _, _, ene_total ) = line.split()
				
				r_resi, l_resi = int(r_resi), int(l_resi)
				ene_vdw, ene_elec, ene_total = float(ene_vdw), float(ene_elec), float(ene_total)

				nodes.add((r_segm, r_resi, r_resn))
				nodes.add((l_segm, l_resi, l_resn))
				edges[((r_segm, r_resi, r_resn), (l_segm, l_resi, l_resn))] = (ene_vdw, ene_elec, ene_total)
				edges[((l_segm, l_resi, l_resn), (r_segm, r_resi, r_resn))] = (ene_vdw, ene_elec, ene_total)

		# Second pass: build JSON data structures
		nodes = sorted(nodes)

		n_nodes = len(nodes)
		# Could be made generic by a series of matrices
		vdw_graph = [[0 for _ in range(n_nodes)] for _ in range(n_nodes)]
		elec_graph = [[0 for _ in range(n_nodes)] for _ in range(n_nodes)]
		total_graph = [[0 for _ in range(n_nodes)] for _ in range(n_nodes)]

		for e in edges:
			res_a, res_b = e
			i_a, i_b = nodes.index(res_a), nodes.index(res_b)
			vdw_graph[i_a][i_b] = edges[e][0]
			vdw_graph[i_b][i_a] = edges[e][0]
			elec_graph[i_a][i_b] = edges[e][1]
			elec_graph[i_b][i_a] = edges[e][1]
			total_graph[i_a][i_b] = edges[e][2]
			total_graph[i_b][i_a] = edges[e][2]
		
		nodes = [{'seg': i, 'resi': j, 'resn': k} for (i,j,k) in nodes]

		return render_template('results.html', vdw_graph=json.dumps(vdw_graph),
											   elec_graph=json.dumps(elec_graph),
											   total_graph=json.dumps(total_graph),
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

