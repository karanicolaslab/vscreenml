#! /usr/bin/env python

#This code is a python implementation of the atom counts features used in 
#Ballester PJ, Mitchell JB. A machine learning approach to predicting protein-ligand binding affinity with applications to molecular docking. Bioinformatics. 2010; 26:1169-75.
#Inspiration for the CartesianPoint,Atom and Loader classes gotten from Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.


import numpy as np
#import h5py
import fnmatch
#from __future__ import print_function
import math
import os
import sys
import textwrap

class CartesianPoint:

	x=88888.0
	y=88888.0
	z=88888.0
    
	def __init__ (self, x, y ,z):
		self.x = x
		self.y = y
		self.z = z
        
	def distance_to(self,a_point):
		return math.sqrt(math.pow(self.x - a_point.x,2) + math.pow(self.y - a_point.y,2) + math.pow(self.z- a_point.z,2))

class Atom:
        
	def __init__ (self):
		self.record_name = ""
		self.atom_name = ""
		self.residue_name = ""
		self.coordinates = CartesianPoint(88888.0,88888.0,88888.0)
		self.elem_sym = ""
		self.pdb_index = ""
		#self.line=""
		#self.atom_type=""
		self.atom_no = 0
		self.res_id = 0
		self.chain_id = ""

	def read_pdb_line(self, line):
		self.line = line
		self.atom_name = line[12:16].strip()      
		self.coordinates = CartesianPoint(float(line[30:38]), float(line[38:46]), float(line[46:54]))

                if self.elem_sym == "": 
			# guessing elem from name
			first_two_letters = self.atom_name[0:2].strip().upper()
			if first_two_letters=='BR':
				self.elem_sym='BR'
			elif first_two_letters=='AL':
				self.elem_sym='AL'
			elif first_two_letters=='CL':
				self.elem_sym='CL'
			elif first_two_letters=='BI':
				self.elem_sym='BI'
			elif first_two_letters=='AS':
				self.elem_sym='AS'
			elif first_two_letters=='AG':
				self.elem_sym='AG'
			elif first_two_letters=='LI':
				self.elem_sym='LI'
			elif first_two_letters=='MG':
				self.elem_sym='MG'
			elif first_two_letters=='MN':
				self.elem_sym='MN'
			elif first_two_letters=='RH':
				self.elem_sym='RH'
			elif first_two_letters=='ZN':
				self.elem_sym='ZN'
			elif first_two_letters=='FE':
				self.elem_sym='FE'

			else: #So,we use just the first letter.
			# Remove any number from elem_sym
				self.elem_sym = self.atom_name
				self.elem_sym = self.elem_sym.replace('0','')
				self.elem_sym = self.elem_sym.replace('1','')
				self.elem_sym = self.elem_sym.replace('2','')
				self.elem_sym = self.elem_sym.replace('3','')
				self.elem_sym = self.elem_sym.replace('4','')
				self.elem_sym = self.elem_sym.replace('5','')
				self.elem_sym = self.elem_sym.replace('6','')
				self.elem_sym = self.elem_sym.replace('7','')
				self.elem_sym = self.elem_sym.replace('8','')
				self.elem_sym = self.elem_sym.replace('9','')
				self.elem_sym = self.elem_sym.replace('@','')

				self.elem_sym = self.elem_sym[0:1].strip().upper()
		elem_type = self.elem_sym[0:2].strip().upper()
		if elem_type  == 'C':
			self.atom_no = 6
		elif elem_type  == 'N':
                        self.atom_no = 7
		elif elem_type  == 'O':
                        self.atom_no = 8
		elif elem_type  == 'F':
                        self.atom_no = 9
		elif elem_type  == 'P':
                        self.atom_no = 15
		elif elem_type  == 'S':
                        self.atom_no = 16
		elif elem_type  == 'SD':
                        self.atom_no = 16
		elif elem_type  == 'SG':
                        self.atom_no = 16
		elif elem_type  == 'CL':
                        self.atom_no = 17
		elif elem_type  == 'BR':
                        self.atom_no = 35
		elif elem_type  == 'I':
                        self.atom_no = 53


		self.pdb_index = line[6:11].strip()
		self.residue_name = line[17:20]
		self.residue_name = " " + self.residue_name[-3:] # this only uses the rightmost three characters, essentially removing unique rotamer identification
        
		try: self.res_id = int(line[22:26]) # because it's possible the pdbqt might not have any resid entries.
		except: pass
        
		self.chain_id = line[21]
		if self.residue_name.strip() == "": self.residue_name = " MOL"
		self.record_name = line[0:6]

class Loader:

	def __init__ (self):
		self.all_atoms = {}
		self.non_protein_atoms = {}
		self.ligand_com = CartesianPoint(88888.0,88888.0,88888.0)
		self.max_x = -8888.88
		self.min_x = 8888.88
		self.max_y = -8888.88
		self.min_y = 8888.88
		self.max_z = -8888.88
		self.min_z = 8888.88

		self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
     
	def PDBLoad(self, file_name, min_x=-8888.88, max_x=8888.88, min_y=-8888.88, max_y=8888.88, min_z=-8888.88, max_z=8888.88):

		auto_index = 1

		self.__init__()
        
		# Now load the file into a list
		file = open(file_name,"r")
		lines = file.readlines()
		file.close()
        
		atom_already_loaded = [] # going to keep track of atomname_resid_chain pairs, to make sure redundants aren't loaded. This basically
                                 # gets rid of rotomers, I think.
		print file_name
		for t in range(0,len(lines)):
			#print "OK"
			line=lines[t]
            
			if line[:3] == "END":
				#print "OK"
				t = textwrap.wrap("WARNING: END or ENDMDL encountered in " + file_name + ". Everyline after this will be ignored. If your PDB file has multiple pose or is an ensemble structure, split them up into individual pose or model", 80)
				print "\n".join(t) + "\n"
				print line
				break
                    
			if len(line) >= 7:
				#print "OK"
				if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
					temp_atom = Atom()
					temp_atom.read_pdb_line(line)
					#print "OK"
					if temp_atom.coordinates.x > min_x and temp_atom.coordinates.x < max_x and temp_atom.coordinates.y > min_y and temp_atom.coordinates.y < max_y and temp_atom.coordinates.z > min_z and temp_atom.coordinates.z < max_z:
						#print "OK"
						if self.max_x < temp_atom.coordinates.x: self.max_x = temp_atom.coordinates.x
						if self.max_y < temp_atom.coordinates.y: self.max_y = temp_atom.coordinates.y
						if self.max_z < temp_atom.coordinates.z: self.max_z = temp_atom.coordinates.z
                        
						if self.min_x > temp_atom.coordinates.x: self.min_x = temp_atom.coordinates.x
						if self.min_y > temp_atom.coordinates.y: self.min_y = temp_atom.coordinates.y
						if self.min_z > temp_atom.coordinates.z: self.min_z = temp_atom.coordinates.z

						key = temp_atom.atom_name.strip() + "_" + str(temp_atom.res_id) + "_" + temp_atom.residue_name.strip() + "_" + temp_atom.chain_id.strip() # this string unique identifies each atom
                        

                        
						if not key in atom_already_loaded or not temp_atom.residue_name.strip() in self.protein_resnames: # so either the atom hasn't been loaded, or else it's a non-protein atom
                                                                                                            # so note that non-protein atoms can have redundant names, but protein atoms cannot.
                                                                                                            # This is because protein residues often contain rotamers
							atom_already_loaded.append(key) # so each atom can only be loaded once. No rotamers.
							self.all_atoms[auto_index] = temp_atom # So you're actually reindexing everything here.
							if not temp_atom.residue_name[-3:] in self.protein_resnames: self.non_protein_atoms[auto_index] = temp_atom;#print "OK"
							
							auto_index = auto_index + 1


class AtomCountFeaturizer:
	def __init__ (self,complex_pdb):
		self.complex_pdb = complex_pdb
		self.Complex_PDB = Loader()
		self.Complex_PDB.PDBLoad(complex_pdb)



	def calc_feature(self):

		feature_dict = {
				'6.6' : 0,
				'7.6'  :0,
				'8.6' :0,
				'16.6' :0,
				'6.7' :0,
				'7.7' :0,
				'8.7' :0,
				'16.7' :0,
				'6.8' :0,
				'7.8' :0,
				'8.8' :0,
				'16.8' :0,
				'6.9' :0,
				'7.9' :0,
				'8.9' :0,
				'16.9' : 0,
				'6.15' :0,
				'7.15' :0,
				'8.15' :0,
				'16.15' :0,
				'6.16' :0,
				'7.16' :0,
				'8.16' :0,
				'16.16' :0,
				'6.17' :0,
				'7.17' :0,
				'8.17' :0,
				'16.17' :0,
				'6.35' :0,
				'7.35' :0,
				'8.35' :0,
				'16.35' :0,
				'6.53' :0,
				'7.53' :0,
				'8.53' :0,
				'16.53' :0
				}

		feature_list = ['6.6','7.6','8.6','16.6','6.7','7.7','8.7','16.7','6.8','7.8','8.8','16.8',
				'6.9','7.9','8.9','16.9','6.15','7.15','8.15','16.15','6.16','7.16','8.16',
				'16.16','6.17','7.17','8.17','16.17','6.35','7.35','8.35','16.35','6.53','7.53','8.53','16.53']

		elem_interest = [6,7,8,9,15,16,17,35,53]
		for jatom in self.Complex_PDB.non_protein_atoms:
			for iatom in self.Complex_PDB.all_atoms:
				if self.Complex_PDB.all_atoms[iatom].res_id != self.Complex_PDB.non_protein_atoms[jatom].res_id:
					dist = 0.0
					dist = self.Complex_PDB.non_protein_atoms[jatom].coordinates.distance_to(self.Complex_PDB.all_atoms[iatom].coordinates)
					#12 Arngstrom distance
					prot_atom_no = 0
					lig_atom_no = 0
					lig_atom_no = self.Complex_PDB.non_protein_atoms[jatom].atom_no
					prot_atom_no = self.Complex_PDB.all_atoms[iatom].atom_no
					if dist < 12:
						if (lig_atom_no in elem_interest) and (prot_atom_no in elem_interest):
							#key_list = sorted([lig_atom_no,prot_atom_no])
							key_list = [lig_atom_no,prot_atom_no]
							key = str(key_list[1])+'.'+str(key_list[0])
							feature_dict[key] = feature_dict.get(key, 0) + 1

		for feat in feature_list:
			print feat,
		print 'PDB'
		for feat in feature_list:
			print str(feature_dict[feat]),
		print self.complex_pdb
		return feature_dict


					

def exec_func(File):
	rf = AtomCountFeaturizer(File).calc_feature()
	return rf

DIRIN='/fccc/users/karanicolaslab/adeshiy'

INFILE=sys.argv[1]

exec_func(INFILE)
