#!/usr/bin/env python

# To avoid the creation of a temporary PDB file
try:
	from cStringIO import StringIO
except ImportError:
	from StringIO import StringIO

# import PDB / mmCIF IO modules
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionException

# Define custom exceptions to facilitate error handling
class MolError(Exception):
	"""Custom exception to report validation problems"""
	pass


# Main processer class
class MolProcesser(object):
	"""Class to read and validate a molecular structure"""

	def __init__(self, pdbdata):
		# in toppar/: egrep -i '^RESI' *top | awk '{ print $2 }' | sort -u
		self.topology = set(("A", "ACE", "ADE", "AG", "AG1", "AL", "AL3", "ALA", "ALY",
							"AR", "ARG", "AS", "ASH", "ASN", "ASP", "AU", "AU1", "AU3",
							"BR", "BR1", "C", "CA", "CA2", "CD", "CD2", "CHX", "CL", "CL1",
							"CO", "CO2", "CO3", "CR", "CR2", "CR3", "CS", "CS1", "CSP",
							"CU", "CU1", "CU2", "CYF", "CYM", "CYS", "CYT", "DA", "DC",
							"DG", "DT", "DUM", "F", "F1", "FE", "FE2", "FE3", "G", "GLH",
							"GLN", "GLU", "GLY", "GUA", "HG", "HG1", "HG2", "HIS", "HO",
							"HO3", "HYP", "I", "I1", "ILE", "IR", "IR3", "K", "K1", "KR",
							"LEU", "LI1", "LYS", "M3L", "MET", "MG", "MG2", "MLY", "MLZ",
							"MN", "MN2", "MN3", "MO", "MO3", "MSE", "NA", "NA1", "NEP",
							"NI", "NI2", "OS", "OS4", "PB", "PB2", "PHE", "PO4", "PRO",
							"PT", "PT2", "PTR", "PUR", "SEP", "SER", "SO4", "SR", "SR2",
							"THR", "THY", "TOP", "TRP", "TYP", "TYR", "TYS", "U", "U3",
							"U4", "URI", "V", "V2", "V3", "VAL", "WO4", "XE", "YB", "YB2",
							"YB3", "ZN", "ZN2"))

		self.parser = PDBParser(QUIET=1)
		self.structure = self.parse(pdbdata)
		self.validate()

	@property
	def tostring(self):
		"""Returns the validated structure as a string"""

		stream = StringIO()

		io = PDBIO()
		io.set_structure(self.structure)
		io.save(file=stream)

		stream.seek(0)
		contents = stream.read()
		stream.close()

		return contents

	def parse(self, handle):
		"""Parses a PDB using Biopython"""

		try:
			s = self.parser.get_structure('xyz', handle)
		except PDBConstructionException, e:
			raise MolError(e)
		else:
			return s

	def validate(self):
		"""Top-level validation routine"""

		s = self.structure

		# Check for multiple models
		# Keep first model only
		n_models = sum(1 for _ in s.get_models())
		if n_models > 1:
			model_lst = s.child_list[1:]
			for model in model_lst:
				s.detach_child(model.id)

		# Check for multiple chains
		# Convert chain identifier to segid element and delete chain identifier (CNS)
		n_chains = sum(1 for _ in s.get_chains())
		if n_chains < 2:
			raise MolError('Structure contains only one chain')
		else:
			for chain in s[0].child_list:
				chain_id = chain.id
				chain.id = ' '
				res_lst = chain.child_list[:]
				for res in res_lst:
					# Make residue-level checks here to avoid looping again
					# Remove HETATM
					if res.id[0] != ' ':
						chain.detach_child(res.id)
						continue

					# Check for insertion codes
					if res.id[2] != ' ':
						msg = 'Residue insertions are unsupported: {0}.{1}.{2}.{3}'
						raise MolError(msg.format(chain_id, res.id[1], res.resname, res.id[2]))

					# Check for unknown residues based on name
					if res.resname.strip() not in self.topology:
						msg = 'Unsupported or unknown residue name: {0!r} ({1}.{2})'
						raise MolError(msg.format(res.resname, chain_id, res.id[1]))

					# Add segment id
					res.segid = chain_id

		# Check for double occuppancies
		# Keep highest occupancy atom
		atom_lst = list(s.get_atoms())
		for atom in atom_lst:
			if atom.is_disordered():
				residue = atom.parent
				sel_atom = atom.selected_child
				sel_atom.altloc = ' '
				sel_atom.disordered_flag = 0
				residue.detach_child(atom.id)
				residue.add(sel_atom)

if __name__ == '__main__':

	import sys
	import unittest

	if not sys.argv[1:]:
		# Test the Class
		# Dummy 1: Double Occ, Multiple Models, HETATMs (should pass)
		dummy_1 = """
MODEL        1
ATOM      1  N   VAL A   1      16.783  48.812  26.447  1.00 30.15      A    N
ATOM      2  CA  VAL A   1      17.591  48.101  25.416  1.00 27.93      A    C
ATOM      3  C   VAL A   1      16.643  47.160  24.676  1.00 25.52      A    C
ATOM      4  O   VAL A   1      16.213  46.129  25.220  1.00 26.86      A    O
ATOM      5  CB  VAL A   1      18.813  47.329  25.940  1.00 28.10      A    C
ATOM      6  CG1 VAL A   1      19.512  46.623  24.799  1.00 29.27      A    C
ATOM      7  CG2 VAL A   1      19.838  48.188  26.670  1.00 28.76      A    C
ATOM      8  N   ILE B  10      16.327  47.509  23.466  1.00 23.19      B    N
ATOM      9  CA AILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     10  CA BILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     11  C   ILE B  10      16.364  45.902  21.716  1.00 19.99      B    C
ATOM     12  O   ILE B  10      16.847  46.590  20.823  1.00 18.74      B    O
ATOM     13  CB  ILE B  10      14.415  47.516  21.981  1.00 22.04      B    C
ATOM     14  CG1 ILE B  10      13.710  48.294  23.129  1.00 22.39      B    C
ATOM     15  CG2 ILE B  10      13.485  46.762  21.052  1.00 19.05      B    C
ATOM     16  CD1 ILE B  10      12.870  49.416  22.486  1.00 23.20      B    C
HETATM   17  O   HOH A   2      51.263  27.344   1.862  1.00 39.29      A    O
HETATM   18  O   HOH A   2      45.191  34.800 -10.616  1.00 45.76      A    O
ENDMDL
MODEL        2
ATOM      1  N   VAL A   1      16.783  48.812  26.447  1.00 30.15      A    N
ATOM      2  CA  VAL A   1      17.591  48.101  25.416  1.00 27.93      A    C
ATOM      3  C   VAL A   1      16.643  47.160  24.676  1.00 25.52      A    C
ATOM      4  O   VAL A   1      16.213  46.129  25.220  1.00 26.86      A    O
ATOM      5  CB  VAL A   1      18.813  47.329  25.940  1.00 28.10      A    C
ATOM      6  CG1 VAL A   1      19.512  46.623  24.799  1.00 29.27      A    C
ATOM      7  CG2 VAL A   1      19.838  48.188  26.670  1.00 28.76      A    C
ATOM      8  N   ILE B  10      16.327  47.509  23.466  1.00 23.19      B    N
ATOM      9  CA AILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     10  CA BILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     11  C   ILE B  10      16.364  45.902  21.716  1.00 19.99      B    C
ATOM     12  O   ILE B  10      16.847  46.590  20.823  1.00 18.74      B    O
ATOM     13  CB  ILE B  10      14.415  47.516  21.981  1.00 22.04      B    C
ATOM     14  CG1 ILE B  10      13.710  48.294  23.129  1.00 22.39      B    C
ATOM     15  CG2 ILE B  10      13.485  46.762  21.052  1.00 19.05      B    C
ATOM     16  CD1 ILE B  10      12.870  49.416  22.486  1.00 23.20      B    C
HETATM   17  O   HOH A   2      51.263  27.344   1.862  1.00 39.29      A    O
HETATM   18  O   HOH A   2      45.191  34.800 -10.616  1.00 45.76      A    O
ENDMDL
END
	"""
		# Dummy 2: One Chain
		dummy_2 = """
MODEL        1
ATOM      1  N   VAL A   1      16.783  48.812  26.447  1.00 30.15      A    N
ATOM      2  CA  VAL A   1      17.591  48.101  25.416  1.00 27.93      A    C
ATOM      3  C   VAL A   1      16.643  47.160  24.676  1.00 25.52      A    C
ATOM      4  O   VAL A   1      16.213  46.129  25.220  1.00 26.86      A    O
ATOM      5  CB  VAL A   1      18.813  47.329  25.940  1.00 28.10      A    C
ATOM      6  CG1 VAL A   1      19.512  46.623  24.799  1.00 29.27      A    C
ATOM      7  CG2 VAL A   1      19.838  48.188  26.670  1.00 28.76      A    C
ENDMDL
END
	"""
		# Dummy 3: Insertion Codes
		dummy_3 = """
MODEL        1
ATOM      1  N   VAL A   1      16.783  48.812  26.447  1.00 30.15      A    N
ATOM      2  CA  VAL A   1      17.591  48.101  25.416  1.00 27.93      A    C
ATOM      3  C   VAL A   1      16.643  47.160  24.676  1.00 25.52      A    C
ATOM      4  O   VAL A   1      16.213  46.129  25.220  1.00 26.86      A    O
ATOM      5  CB  VAL A   1      18.813  47.329  25.940  1.00 28.10      A    C
ATOM      6  CG1 VAL A   1      19.512  46.623  24.799  1.00 29.27      A    C
ATOM      7  CG2 VAL A   1      19.838  48.188  26.670  1.00 28.76      A    C
ATOM      8  N   ILE B  10A     16.327  47.509  23.466  1.00 23.19      B    N
ATOM      9  CA  ILE B  10A     15.435  46.625  22.666  1.00 21.26      B    C
ATOM     10  C   ILE B  10A     16.364  45.902  21.716  1.00 19.99      B    C
ATOM     11  O   ILE B  10A     16.847  46.590  20.823  1.00 18.74      B    O
ATOM     12  CB  ILE B  10A     14.415  47.516  21.981  1.00 22.04      B    C
ATOM     13  CG1 ILE B  10A     13.710  48.294  23.129  1.00 22.39      B    C
ATOM     14  CG2 ILE B  10A     13.485  46.762  21.052  1.00 19.05      B    C
ATOM     15  CD1 ILE B  10A     12.870  49.416  22.486  1.00 23.20      B    C
ATOM      8  N   ALA B  10B     16.327  47.509  23.466  1.00 23.19      B    N
ATOM      9  CA  ALA B  10B     15.435  46.625  22.666  1.00 21.26      B    C
ATOM     10  C   ALA B  10B     16.364  45.902  21.716  1.00 19.99      B    C
ATOM     11  O   ALA B  10B     16.847  46.590  20.823  1.00 18.74      B    O
ATOM     12  CB  ALA B  10B     14.415  47.516  21.981  1.00 22.04      B    C
ENDMDL
END
	"""

		# Dummy 4: Unrecognized Residue (HOH)
		dummy_4 = """
MODEL        1
ATOM      1  N   VAL A   1      16.783  48.812  26.447  1.00 30.15      A    N
ATOM      2  CA  VAL A   1      17.591  48.101  25.416  1.00 27.93      A    C
ATOM      3  C   VAL A   1      16.643  47.160  24.676  1.00 25.52      A    C
ATOM      4  O   VAL A   1      16.213  46.129  25.220  1.00 26.86      A    O
ATOM      5  CB  VAL A   1      18.813  47.329  25.940  1.00 28.10      A    C
ATOM      6  CG1 VAL A   1      19.512  46.623  24.799  1.00 29.27      A    C
ATOM      7  CG2 VAL A   1      19.838  48.188  26.670  1.00 28.76      A    C
ATOM     17  O   HOH A   2      51.263  27.344   1.862  1.00 39.29      A    O
ATOM      8  N   ILE B  10      16.327  47.509  23.466  1.00 23.19      B    N
ATOM      9  CA AILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     10  CA BILE B  10      15.435  46.625  22.666  0.50 21.26      B    C
ATOM     11  C   ILE B  10      16.364  45.902  21.716  1.00 19.99      B    C
ATOM     12  O   ILE B  10      16.847  46.590  20.823  1.00 18.74      B    O
ATOM     13  CB  ILE B  10      14.415  47.516  21.981  1.00 22.04      B    C
ATOM     14  CG1 ILE B  10      13.710  48.294  23.129  1.00 22.39      B    C
ATOM     15  CG2 ILE B  10      13.485  46.762  21.052  1.00 19.05      B    C
ATOM     16  CD1 ILE B  10      12.870  49.416  22.486  1.00 23.20      B    C
ENDMDL
END
"""

		class TestMolecularProcesser(unittest.TestCase):

			def test_validinput(self):
				"""Process a valid structure with double occupancies, multiple models, and HETATM"""

				stream = StringIO()
				stream.write(dummy_1)
				stream.seek(0)
				mol = MolProcesser(stream)

				n_models = sum(1 for _ in mol.structure.get_models()) #1
				n_chains = sum(1 for _ in mol.structure.get_chains()) #2
				n_resids = sum(1 for _ in mol.structure.get_residues()) #2
				n_atoms = sum(1 for _ in mol.structure.get_atoms()) #15
				has_docc = sum(1 for a in mol.structure.get_atoms() if a.is_disordered())
				has_hatm = sum(1 for r in mol.structure.get_residues() if r.id[0] != ' ')

				self.assertEqual(n_models, 1)
				self.assertEqual(n_chains, 2)
				self.assertEqual(n_resids, 2)
				self.assertEqual(n_atoms, 15)
				self.assertEqual(has_docc, 0)
				self.assertEqual(has_hatm, 0)

			def test_single_chain(self):
				"""Process a single-chain structure"""

				stream = StringIO()
				stream.write(dummy_2)
				stream.seek(0)

				self.assertRaisesRegexp(MolError, 'contains only one chain', MolProcesser, stream)

			def test_insertion_codes(self):
				"""Process structure with insertion codes"""

				stream = StringIO()
				stream.write(dummy_3)
				stream.seek(0)

				self.assertRaisesRegexp(MolError, 'Residue insertions are unsupported', MolProcesser, stream)

			def test_unknown_residue(self):
				"""Process structure with unknown residues"""

				stream = StringIO()
				stream.write(dummy_4)
				stream.seek(0)

				self.assertRaisesRegexp(MolError, 'Unsupported or unknown residue name', MolProcesser, stream)

			def test_to_string(self):
				"""Write structure as string"""

				stream = StringIO()
				stream.write(dummy_1)
				stream.seek(0)

				mol = MolProcesser(stream)
				n_models = sum(1 for _ in mol.structure.get_models()) #1
				n_chains = sum(1 for _ in mol.structure.get_chains()) #2
				n_resids = sum(1 for _ in mol.structure.get_residues()) #2
				n_atoms = sum(1 for _ in mol.structure.get_atoms()) #15
				has_docc = sum(1 for a in mol.structure.get_atoms() if a.is_disordered())
				has_hatm = sum(1 for r in mol.structure.get_residues() if r.id[0] != ' ')

				stream_2 = StringIO()
				stream_2.write(mol.tostring)
				stream_2.seek(0)

				p = PDBParser(QUIET=1)
				mol_2 = p.get_structure('xyz', stream_2)

				n_models_2 = sum(1 for _ in mol_2.get_models()) #1
				n_resids_2 = sum(1 for _ in mol_2.get_residues()) #2
				n_atoms_2 = sum(1 for _ in mol_2.get_atoms()) #15
				has_docc_2 = sum(1 for a in mol_2.get_atoms() if a.is_disordered())
				has_hatm_2 = sum(1 for r in mol_2.get_residues() if r.id[0] != ' ')

				self.assertEqual(n_models, n_models_2)
				self.assertEqual(n_resids, n_resids_2)
				self.assertEqual(n_atoms, n_atoms_2)
				self.assertEqual(has_docc, has_docc_2)
				self.assertEqual(has_hatm, has_hatm_2)

		unittest.main(verbosity=0)

	else:
		for structure in sys.argv[1:]:
			with open(structure) as handle:
				mol = MolProcesser(handle)

				io = PDBIO()
				io.set_structure(mol.structure)
				io.save('{0}_val.pdb'.format(structure.split('.')[0]))
