#!/usr/bin/env python
#Author: Andreas Stroehlein
#Email: astroehlein@unimelb.edu.au
#Microbiological Diagnostic Unit Public Health Laboratory
#https://biomedicalsciences.unimelb.edu.au/departments/microbiology-Immunology/research/services/microbiological-diagnostic-unit-public-health-laboratory
#For info on Meningococcal Deduced Vaccine Antigen Reactivity (MenDeVAR) Index see https://doi.org/10.1128/JCM.02161-20

import os
import sys
from pkg_resources import resource_filename

###############################################################

#MenDeVar config file
MDVRCONFIG = "mendevar_config.txt"

#Output naming for Bexsero vaccine MenDeVAR column (change nomenclature here if reports require different wording)
N_BEXSERO = {
	"green" : "Exact match", 
	"amber" : "Cross-reactive", 
	"red" : "No match",
	"grey" : "Insufficent data", 
	}

#Output naming for Trumenba vaccine MenDeVAR column (change nomenclature here if reports require different wording)
N_TRUMENBA = {
	"green" : "Exact match", 
	"amber" : "Cross-reactive", 
	"red" : "No match",
	"grey" : "Insufficent data", 
	}

###############################################################

class MenDeVAR:

	def __init__(self, vacc, mdvr, atg_type, antigen):
		self.vacc = vacc
		self.mdvr = mdvr
		self.atg_type = atg_type
		self.antigen = antigen

def readMendevarConf():
	DBpath = resource_filename(__name__, 'db')
	mdvrconf = os.path.join(DBpath, MDVRCONFIG)

	with open(mdvrconf, 'r') as mdvrf:
		lines = mdvrf.readlines()[1:]
		lines = [line.rstrip() for line in lines]

	bexsero = {
		"green" : {
			"fhbp" :    set(),
			"nhba" :    set(),
			"nada" :    set(),
			"poravr2" : set()
		},

		"amber" : {
			"fhbp" : set(),
			"nhba" : set(),
			"nada" : set()
		},
		 
		"red" : {
			"fhbp" : set(),
			"nhba" : set(),
			"nada" : set()
		},
	}

	trumenba = {

		"green" : {
			"fhbp" : set()
		},

		"amber" : {
			"fhbp" : set()
		}
	}


	for cfg_line in lines:
		line = cfg_line.split()
		entry = MenDeVAR(line[0], line[1], line[3], line[2])
	
		if entry.vacc == "bexsero":
			if entry.mdvr == "green":
				if entry.atg_type == "fhbp":
					bexsero["green"]["fhbp"].add(entry.antigen)

				elif entry.atg_type == "nhba":
					bexsero["green"]["nhba"].add(entry.antigen)

				elif entry.atg_type == "nada":
					bexsero["green"]["nada"].add(entry.antigen)

				elif entry.atg_type == "poravr2":
					bexsero["green"]["poravr2"].add(entry.antigen)

				else:
					raise Exception("Corrupt config file: Entry in column 3 should be one of \"fhbp\", \"nhba\", \"nada\" and \"poravr2\", but is " + entry.atg_type)
					sys.exit(1)	

			elif entry.mdvr == "amber":
				if entry.atg_type == "fhbp":
					bexsero["amber"]["fhbp"].add(entry.antigen)

				elif entry.atg_type == "nhba":
					bexsero["amber"]["nhba"].add(entry.antigen)

				elif entry.atg_type == "nada":
					bexsero["amber"]["nada"].add(entry.antigen)

				else:
					raise Exception("Corrupt config file: Entry in column 3 should be one of \"fhbp\", \"nhba\" and \"nada\", but is " + entry.atg_type)
					sys.exit(1)	

			elif entry.mdvr == "red":
				if entry.atg_type == "fhbp":
					bexsero["red"]["fhbp"].add(entry.antigen)

				elif entry.atg_type == "nhba":
					bexsero["red"]["nhba"].add(entry.antigen)

				elif entry.atg_type == "nada":
					bexsero["red"]["nada"].add(entry.antigen)

				else:
					raise Exception("Corrupt config file: Entry in column 3 should be one of \"fhbp\", \"nhba\" and \"nada\", but is " + entry.atg_type)
					sys.exit(1)	

			else:
				raise Exception("Corrupt config file: Entry in column 2 should be one of \"green\", \"amber\" and \"red\", but is " + entry.mdvr)
				sys.exit(1)				


		elif entry.vacc == "trumenba":
			if entry.mdvr == "green":
				if entry.atg_type == "fhbp":
					trumenba["green"]["fhbp"].add(entry.antigen)

				else:
					raise Exception("Corrupt config file: Entry in column 3 should be \"fhbp\", but is " + entry.atg_type)
					sys.exit(1)	

			elif entry.mdvr == "amber":
				if entry.atg_type == "fhbp":
					trumenba["amber"]["fhbp"].add(entry.antigen)

				else:
					raise Exception("Corrupt config file: Entry in column 3 should be \"fhbp\", but is " + entry.atg_type)
					sys.exit(1)

			else:
				raise Exception("Corrupt config file: Entry in column 2 should be one of \"green\" and \"amber\", but is " + entry.mdvr)
				sys.exit(1)


		else:
			raise Exception("Corrupt config file: Entry in column 1 should be one of \"bexsero\" and \"trumenba\" but is " + entry.vacc)
			sys.exit(1)

	return bexsero, trumenba

def createMendevar(config, poravr, fhbp, nhba, nada):
	# NOTE can there be an "/" in any of poravr, fhbp, nhba, nada?
	# See "'/'.join" lines in meningotype.py
	
	bexsero = config[0]
	trumenba = config[1]

	#separate VR1 from VR2
	if poravr != '-':
		poravr2 = poravr.split(',')[1]
		# This removes the subvariants indicated by "-" which is not in concordance with online MenDeVAR implementation, e.g. porA 4 is an "Exact match", but 4-13 isn't.
		# To be consistent with online tool, subvariants are matched against the config sets (which at this point only contain main variants, i.e. no variants containing "-"
		# poravr2 = poravr.split(',')[1].split('-')[0] 
	else:
		poravr2 = '-'
	
	###Bexsero logic###
	# green: isolate contains ≥1 exact sequence match to antigenic variants found in the vaccine.
	if (
		fhbp in bexsero["green"]["fhbp"]
		or nhba in bexsero["green"]["nhba"]
		or nada in bexsero["green"]["nada"]
		or poravr2 in bexsero["green"]["poravr2"]
	):
		bexsero = N_BEXSERO["green"]

	# amber: isolate contains ≥1 antigenic variant deemed cross-reactive to vaccine variants through experimental studies.
	elif (
		fhbp in bexsero["amber"]["fhbp"]
		or nhba in bexsero["amber"]["nhba"]
		or nada in bexsero["amber"]["nada"]
	):
		bexsero = N_BEXSERO["amber"]

	# red: all the isolate's antigenic variants have been deemed not cross-reactive to vaccine variants through experimental studies.
	elif (
		fhbp in bexsero["red"]["fhbp"]
		and nhba in bexsero["red"]["nhba"]
		and nada in bexsero["red"]["nada"]
		and poravr2 != 4
	):
		bexsero = N_BEXSERO["red"]

	# grey: isolate contains antigens for which there is insufficient data from or are yet to be tested in experimental studies.
	else:
		bexsero = N_BEXSERO["grey"]		


	###Trumenba logic###
	# green: isolate contains ≥1 exact sequence match to antigenic variants found in the vaccine.
	if fhbp in trumenba["green"]["fhbp"]:
		trumenba = N_TRUMENBA["green"]

	# amber: isolate contains ≥1 antigenic variant deemed cross-reactive to vaccine variants through experimental studies.
	elif fhbp in trumenba["amber"]["fhbp"]:
		trumenba = N_TRUMENBA["amber"]

	# red: all the isolate's antigenic variants have been deemed not cross-reactive to vaccine variants through experimental studies.
	# NOTE this seems to contradict "missing" from online, but online status "missing" also does not create a red tag
	elif fhbp == "-":
		trumenba = N_TRUMENBA["red"]

	# grey: isolate contains antigens for which there is insufficient data from or are yet to be tested in experimental studies.
	else:
		trumenba = N_TRUMENBA["grey"]			

	return bexsero, trumenba
