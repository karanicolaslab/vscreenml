#! /usr/bin/env python
import numpy as np
#import h5py
import fnmatch
#from __future__ import print_function
import math
import os
import sys
import textwrap

class point:

	X=88888.0
	Y=88888.0
	Z=88888.0
    
	def __init__ (self, X, Y ,Z):
		self.X = X
		self.Y = Y
		self.Z = Z

	def copy_of(self):
		return point(self.X, self.Y, self.Z)

	def print_coord(self):
		print str(self.X)+"\t"+str(self.Y)+"\t"+str(self.Z)
        
        
	def distance_to(self,apoint):
		return math.sqrt(math.pow(self.X - apoint.X,2) + math.pow(self.Y - apoint.Y,2) + math.pow(self.Z- apoint.Z,2))

	def description(self):
		return str(self.X) + " " + str(self.Y) + " " + str(self.Z)

	def Magnitude(self):
		return self.distance_to(point(0,0,0))


class atom:
        
	def __init__ (self):
		self.RecordName = ""
		self.AtomName = ""
		self.ResidueName = ""
		self.Coordinates = point(88888.0,88888.0,88888.0)
		self.ElemSym = ""
		self.PDBIndex = ""
		#self.line=""
		#self.atomtype=""
		self.atomno = 0
		self.IndeciesOfAtomsConnecting=[]
		self.Charge = 0
		self.ResId = 0
		self.ChainId = ""
		self.Structure = ""
		self.Comment = ""
        
	def copy_of(self):
		atomcopy = atom()
		atomcopy.RecordName = self.RecordName
		atomcopy.AtomName = self.AtomName 
		atomcopy.ResidueName = self.ResidueName
		atomcopy.Coordinates = self.Coordinates.copy_of()
		atomcopy.ElemSym = self.ElemSym 
		atomcopy.PDBIndex = self.PDBIndex 
                atomcopy.atomno = self.atomno 
		#atomcopy.line= self.line
		#atomcopy.AtomType= self.AtomType
		atomcopy.IndeciesOfAtomsConnecting = self.IndeciesOfAtomsConnecting[:]
		atomcopy.Charge = self.Charge 
		atomcopy.ResId = self.ResId 
		atomcopy.ChainId = self.ChainId 
		atomcopy.Structure = self.Structure 
		atomcopy.Comment = self.Comment
        
		return atomcopy


	def ReadPDBLine(self, Line):
		self.line = Line
		self.AtomName = Line[12:16].strip()      
		self.Coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))

                if self.ElemSym == "": # try to guess at element from name
			first_two_letters = self.AtomName[0:2].strip().upper()
			if first_two_letters=='BR':
				self.ElemSym='BR'
			elif first_two_letters=='AL':
				self.ElemSym='AL'
			elif first_two_letters=='CL':
				self.ElemSym='CL'
			elif first_two_letters=='BI':
				self.ElemSym='BI'
			elif first_two_letters=='AS':
				self.ElemSym='AS'
			elif first_two_letters=='AG':
				self.ElemSym='AG'
			elif first_two_letters=='LI':
				self.ElemSym='LI'
			elif first_two_letters=='MG':
				self.ElemSym='MG'
			elif first_two_letters=='MN':
				self.ElemSym='MN'
			elif first_two_letters=='RH':
				self.ElemSym='RH'
			elif first_two_letters=='ZN':
				self.ElemSym='ZN'
			elif first_two_letters=='FE':
				self.ElemSym='FE'

			else: #So,we use just the first letter.
			# Remove any number from ElemSym
				self.ElemSym = self.AtomName
				self.ElemSym = self.ElemSym.replace('0','')
				self.ElemSym = self.ElemSym.replace('1','')
				self.ElemSym = self.ElemSym.replace('2','')
				self.ElemSym = self.ElemSym.replace('3','')
				self.ElemSym = self.ElemSym.replace('4','')
				self.ElemSym = self.ElemSym.replace('5','')
				self.ElemSym = self.ElemSym.replace('6','')
				self.ElemSym = self.ElemSym.replace('7','')
				self.ElemSym = self.ElemSym.replace('8','')
				self.ElemSym = self.ElemSym.replace('9','')
				self.ElemSym = self.ElemSym.replace('@','')

				self.ElemSym = self.ElemSym[0:1].strip().upper()
		elemtype = self.ElemSym[0:2].strip().upper()
		if elemtype  == 'C':
			self.atomno = 6
		elif elemtype  == 'N':
                        self.atomno = 7
		elif elemtype  == 'O':
                        self.atomno = 8
		elif elemtype  == 'F':
                        self.atomno = 9
		elif elemtype  == 'P':
                        self.atomno = 15
		elif elemtype  == 'S':
                        self.atomno = 16
		elif elemtype  == 'SD':
                        self.atomno = 16
		elif elemtype  == 'SG':
                        self.atomno = 16
		elif elemtype  == 'CL':
                        self.atomno = 17
		elif elemtype  == 'BR':
                        self.atomno = 35
		elif elemtype  == 'I':
                        self.atomno = 53


		self.PDBIndex = Line[6:11].strip()
		self.ResidueName = Line[17:20]
		self.ResidueName = " " + self.ResidueName[-3:] # this only uses the rightmost three characters, essentially removing unique rotamer identification
        
		try: self.ResId = int(Line[22:26]) # because it's possible the pdbqt might not have any resid entries.
		except: pass
        
		self.ChainId = Line[21]
		if self.ResidueName.strip() == "": self.ResidueName = " MOL"
		self.RecordName = Line[0:6]

class Loader:

	def __init__ (self):
		self.AllAtoms = {}
		self.NonProteinAtoms = {}
		self.ligand_CoM = point(88888.0,88888.0,88888.0)
		self.max_x = -8888.88
		self.min_x = 8888.88
		self.max_y = -8888.88
		self.min_y = 8888.88
		self.max_z = -8888.88
		self.min_z = 8888.88

		self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
     
	def PDBLoad(self, FileName, min_x=-8888.88, max_x=8888.88, min_y=-8888.88, max_y=8888.88, min_z=-8888.88, max_z=8888.88):

		autoindex = 1

		self.__init__()
        
		# Now load the file into a list
		file = open(FileName,"r")
		lines = file.readlines()
		file.close()
        
		atom_already_loaded = [] # going to keep track of atomname_resid_chain pairs, to make sure redundants aren't loaded. This basically
                                 # gets rid of rotomers, I think.
		print FileName
		for t in range(0,len(lines)):
			#print "OK"
			line=lines[t]
            
			if line[:3] == "END":
				#print "OK"
				t = textwrap.wrap("WARNING: END or ENDMDL encountered in " + FileName + ". Everyline after this will be ignored. If your PDB file has multiple pose or is an ensemble structure, split them up into individual pose or model", 80)
				print "\n".join(t) + "\n"
				print line
				break
                    
			if len(line) >= 7:
				#print "OK"
				if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
					TempAtom = atom()
					TempAtom.ReadPDBLine(line)
					#print "OK"
					if TempAtom.Coordinates.X > min_x and TempAtom.Coordinates.X < max_x and TempAtom.Coordinates.Y > min_y and TempAtom.Coordinates.Y < max_y and TempAtom.Coordinates.Z > min_z and TempAtom.Coordinates.Z < max_z:
						#print "OK"
						if self.max_x < TempAtom.Coordinates.X: self.max_x = TempAtom.Coordinates.X
						if self.max_y < TempAtom.Coordinates.Y: self.max_y = TempAtom.Coordinates.Y
						if self.max_z < TempAtom.Coordinates.Z: self.max_z = TempAtom.Coordinates.Z
                        
						if self.min_x > TempAtom.Coordinates.X: self.min_x = TempAtom.Coordinates.X
						if self.min_y > TempAtom.Coordinates.Y: self.min_y = TempAtom.Coordinates.Y
						if self.min_z > TempAtom.Coordinates.Z: self.min_z = TempAtom.Coordinates.Z

						key = TempAtom.AtomName.strip() + "_" + str(TempAtom.ResId) + "_" + TempAtom.ResidueName.strip() + "_" + TempAtom.ChainId.strip() # this string unique identifies each atom
                        
						#if key in atom_already_loaded and TempAtom.ResidueName.strip() in self.protein_resnames: # so this is a protein atom that has already been loaded once
							#self.printout("Warning: Duplicate protein atom detected: \"" + TempAtom.line.strip() + "\". Not loading this duplicate.")
							#print ""
                        
						if not key in atom_already_loaded or not TempAtom.ResidueName.strip() in self.protein_resnames: # so either the atom hasn't been loaded, or else it's a non-protein atom
                                                                                                            # so note that non-protein atoms can have redundant names, but protein atoms cannot.
                                                                                                            # This is because protein residues often contain rotamers
							atom_already_loaded.append(key) # so each atom can only be loaded once. No rotamers.
							self.AllAtoms[autoindex] = TempAtom # So you're actually reindexing everything here.
							if not TempAtom.ResidueName[-3:] in self.protein_resnames: self.NonProteinAtoms[autoindex] = TempAtom;#print "OK"
							
							autoindex = autoindex + 1
	def printout(self, thestring):
                lines = textwrap.wrap(thestring, 80)
                for line in lines:
                        print line

        def CoM(self,pointlist):
                x_ = 0
                y_ = 0
                z_ = 0
                for j in pointlist:
                        x_ += j.X
                        y_ += j.Y
                        z_ += j.Z
                x_ = x_/float(len(pointlist))
                y_ = y_/float(len(pointlist))
                z_ = z_/float(len(pointlist))

                return(x_,y_,z_)


class RFfeature:
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
		for jatom in self.Complex_PDB.NonProteinAtoms:
			for iatom in self.Complex_PDB.AllAtoms:
				if self.Complex_PDB.AllAtoms[iatom].ResId != self.Complex_PDB.NonProteinAtoms[jatom].ResId:
					dist = 0.0
					dist = self.Complex_PDB.NonProteinAtoms[jatom].Coordinates.distance_to(self.Complex_PDB.AllAtoms[iatom].Coordinates)
					#12 Arngstrom distance
					prot_atomno = 0
					lig_atomno = 0
					lig_atomno = self.Complex_PDB.NonProteinAtoms[jatom].atomno
					prot_atomno = self.Complex_PDB.AllAtoms[iatom].atomno
					if dist < 12:
						if (lig_atomno in elem_interest) and (prot_atomno in elem_interest):
							#key_list = sorted([lig_atomno,prot_atomno])
							key_list = [lig_atomno,prot_atomno]
							key = str(key_list[1])+'.'+str(key_list[0])
							feature_dict[key] = feature_dict.get(key, 0) + 1

		for feat in feature_list:
			print feat,
		print 'PDB'
		for feat in feature_list:
			print str(feature_dict[feat]),
		print self.complex_pdb
		return feature_dict


class ComplexChainSelector:

	def __init__ (self,complex_pdb):
		self.complex_pdb = complex_pdb
		self.Complex_PDB = Loader()
		self.Complex_PDB.PDBLoad(complex_pdb)

	def ChainSelect(self):	
		#Center of mass of ligand
		x_ = 0.0
                y_ = 0.0
                z_ = 0.0
		ligandkey = ''
		for jatom in self.Complex_PDB.NonProteinAtoms:
			ligandkey = self.Complex_PDB.NonProteinAtoms[jatom].ResidueName
                        x_ += self.Complex_PDB.NonProteinAtoms[jatom].Coordinates.X
                        y_ += self.Complex_PDB.NonProteinAtoms[jatom].Coordinates.Y
                        z_ += self.Complex_PDB.NonProteinAtoms[jatom].Coordinates.Z
                x_ = x_/float(len(self.Complex_PDB.NonProteinAtoms))
                y_ = y_/float(len(self.Complex_PDB.NonProteinAtoms))
                z_ = z_/float(len(self.Complex_PDB.NonProteinAtoms))

                self.Complex_PDB.ligand_CoM = point(x_,y_,z_)

		#Translation of protein such center of mass is at (0,0,0)
		#iTempPoint = point(88888.0,88888.0,88888.0)
		chain = []
		for iatom in self.Complex_PDB.AllAtoms:
			dist = 0.0
			dist = self.Complex_PDB.AllAtoms[iatom].Coordinates.distance_to(self.Complex_PDB.ligand_CoM)
			if dist < 12:
				chain.append(self.Complex_PDB.AllAtoms[iatom].ChainId)
		chainUniq = list(set(chain))
		print self.complex_pdb[-30:-26],chainUniq,ligandkey
					

def exec_func(File):
	rf = RFfeature(File).calc_feature()
	return rf

DIRIN='/fccc/users/karanicolaslab/adeshiy'

INFILE=sys.argv[1]

exec_func(INFILE)
