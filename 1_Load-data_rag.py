################################################
######### Create tumors data structure #########
################################################

import os
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.stats import wasserstein_distance as wass_dist
import scipy.stats
import pandas as pd
import seaborn as sns
from lifelines import KaplanMeierFitter
import math
from statannot import add_stat_annotation
import matplotlib.colors as colors
from datetime import timedelta

tumors = []

#General tumor class and functions for all data

class Tumor():
	def __init__(self, line):
		self.TumorID = line.TumorID
		self.Labels = []
		self.Days = []
		self.DaysToResolve = []
		self.DaysPostTreatment = []
		self.Areas = []
		self.ROIs = []
		self.ImageData = []
		self.Resolve = line.Resolve
		self.Respond = False
		self.Merge = line.Merge
		self.Treatment = line.Treatment
		self.Cage = line.Cage
		self.Mouse = line.Mouse
		self.Ear = line.Ear
		self.Date = line.Date
	def add_point(self, line):
		self.Labels.append(line.Label)
		self.Days.append(line.Day)
		self.Areas.append(line.Area)
		self.ROIs.append(line.ROI)
	def add_image_data(self, ROI):
		self.ImageData.append(ROI)
	def add_days_to_resolve(self, input):
		self.DaysToResolve = input
	def add_respond(self, input):
		self.Respond = input
	def add_days_post_treatment(self, input):
		self.DaysPostTreatment = input

def get_days_to_resolve(arg1):
	DaysToResolve=[]
	for i in range(len(arg1.Days)):
		tmp = arg1.Days[i]-max(arg1.Days)
		DaysToResolve.append(tmp)
	return DaysToResolve

def get_days_post_treatment(arg1,arg2):
	DaysPostTreatment=[]
	for i in range(len(arg1.Days)):
		tmp = arg1.Days[i]-arg2
		DaysPostTreatment.append(tmp)
	return DaysPostTreatment

def find_T_pheno(arg1):
	pheno_str = ["excluded", "inflamed", "desert", "uncertain"]
	if np.isnan(arg1.Pos[0])==False: 
		ex_pheno = True
		inf_pheno = True
		des_pheno = False
		j=0
		if np.median(arg1.Tcell) <= 1540.124: # for 25% use 1540.124, for 50% use 1600.1584
			des_pheno=True
		while arg1.Scaled_Pos[j]<0.25:
			normed = (arg1.Tcell[j]-min(arg1.Tcell))/(max(arg1.Tcell)-min(arg1.Tcell))
			#print normed
			if normed>0.6:
				ex_pheno=False
			if normed<0.4:
				inf_pheno=False
			j=j+1
		core_index=int(len(arg1.Tcell)*0.5)
		core_median=np.median(arg1.Tcell[0:core_index])
		if des_pheno == True and ex_pheno == False:
			return pheno_str[2]
		if des_pheno == True and ex_pheno == True:
			return pheno_str[0]
		elif ex_pheno == True and core_median <= 6259.61875:#6259.61875: # from median of inflamed core medians
			return pheno_str[0];
		elif inf_pheno == True or core_median > 6259.61875:#6259.61875:
			return pheno_str[1];
		elif ex_pheno == False and inf_pheno == False:
			return pheno_str[3];

def find_core(arg1):
	if np.isnan(arg1.Pos[0])==False: 
		j=0
		core_index=int(len(arg1.Tcell)*0.5)
		core_median=np.median(arg1.Tcell[0:core_index])
		return core_median;

def find_ex_ratio(arg1):
	if np.isnan(arg1.Pos[0])==False: 
		j=0
		while arg1.Scaled_Pos[j]<0.2:
			normed = (arg1.Tcell[j]-min(arg1.Tcell))/(max(arg1.Tcell)-min(arg1.Tcell))
			j=j+1
		return normed;

def normalize_list(arg1):
	norm_vals = []
	minval = min(arg1)
	norm_max = max(arg1)-minval
	for i in range(len(arg1)):
		normed = (arg1[i]-minval)/(norm_max)
		norm_vals.append(normed)
	return norm_vals

def get_respond(arg1):
	if arg1.Areas[-1] < max(arg1.Areas)*0.8:
		return True
	elif arg1.Areas[-1] >= max(arg1.Areas)*0.8:
		return False

#####################################################################################
################################# 20200114 Data #####################################
#####################################################################################

#important information specific to the experiment!

tmpts = [6,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,26,27,30,32,34,37,39,40]
num_cages = 1
num_mice = 24
treatment_group = ['Control', 'aPD-L1', 'aTGFb', 'Combo']
implantation_day=0
treatment_initiated=2

class Line():
	#This class is useful for reading in data from CSV files but is not used after that
	#Probably worth eliminating this class in the future
	def __init__(self, input):
		self.Label = input[1]
		self.Area = float(input[2])*(4.679E-6) #conversion to mm2 for tiled 30x images
		self.ROI = int(input[1].split(":")[1].split("-")[0]);
		self.Day = int(input[1].split("D")[1].split("p")[0])
		self.Mean = float(input[3])
		self.Stdev = float(input[4])
		self.Min = float(input[5])
		self.Max = float(input[6])
		self.TumorID = int(input[7])
		self.Resolve = input[8]
		self.Respond = False
		self.Merge = input[9]
		self.Treatment = input[10]
		self.Cage = int(input[11])
		self.Mouse = int(input[12])
		self.Ear = input[13]
		self.Date = input[14]
	def get_Day(self, input):
		try: 
			day = int(input[1].split("D")[1].split("p")[0])
			self.Day = day
		except: 
			day = int(input[1].split("_d")[1].split("_")[0])
			self.Day = day

class ROI():
	def __init__(self, input):
		self.Day = int(input[0].split("D")[1].split("_")[0].split("p")[0])
		self.Cage = int(1);
		self.Mouse = int(input[0].split("_m")[1].split("e")[0].split("_")[0]);
		self.Ear = input[0].split("_e")[1].split("/")[0];
		self.ROI = int(input[1][0]);
		self.Area = float(input[1][1]);
		self.X = float(input[1][2]);
		self.Y = float(input[1][3]);
		self.Size = float(input[1][4]);
		self.Pos = input[2];
		self.Scaled_Pos = [];
		self.Tumor = input[3];
		self.Scaled_Tumor = [];
		self.Tcell = input[5];
		self.Scaled_Tcell = [];
		self.T_pheno = [];
		self.Ex_ratio = [];
		self.Core = [];
	def add_Scaled_Pos(self, input):
		self.Scaled_Pos = input;
	def add_T_pheno(self, input):
		self.T_pheno = input;
	def add_Ex_ratio(self, input):
		self.Ex_ratio = input;
	def add_Core(self, input):
		self.Core = input;
	def add_Scaled_Tcell(self, input):
		self.Scaled_Tcell = input;
	def add_Scaled_Tumor(self, input):
		self.Scaled_Tumor = input;
	def __eq__(self, other): 
		if not isinstance(other, ROI):
			return NotImplemented
		return self.Day == other.Day and self.Cage == other.Cage and self.Mouse == other.Mouse and self.Ear==other.Ear and self.ROI == other.ROI


path = './'

all_csvs=glob.glob(path+"20200114_Combo_withTcells_Lupe/Ear Folders/*.csv")

all_data=[]
for filename in all_csvs:
	with open(filename, 'r') as csvfile:
		file = list(csv.reader(csvfile))
		#print filename
		for i in range(len(file)):
			#add placeholder cage number
			file[i].insert(11,1)
			#file[i].insert(11,filename[:2].strip("C"))
			file[i].insert(12,filename.split("/")[-1].split("_")[0].split("m")[1])
			file[i].insert(13,filename[-5])
			file[i].insert(14,"20200114_Combo_withTcells_Lupe")
		all_data.append(file[1:])


all_image_data=[]
# Set path for ROI T cell and tumor intensity data CSV files
for c in all_csvs:
	path_2 = c[:-4]+"/"
	print(path_2)
	#Read in T cell and tumor intensity data to ROI objects
	for filename in glob.glob(path_2+"/D*.csv"):
		with open(filename, 'r') as csvfile:
			file = list(csv.reader(csvfile))
			#print filename
			tmp_lst_3=[]
			i=0
			while i<len(file):
				#print("pass! "+ filename)
				tmp_lst=[]
				#add day
				tmp_lst.append(filename)
				#add ROI num, areas, T cell and Tumor intensities
				tmp_lst.append(file[i])
				#add scaled position
				for j in [1, 2, 3, 4]:
					tmp_lst_2=[]
					if i+j < len(file):
						for k in range(len(file[i+j])): 
							if file[i+j][k] == "":
								continue
							else:
								tmp_lst_2.append(float(file[i+j][k]))
						tmp_lst.append(tmp_lst_2)
				tmp_lst_3.append(tmp_lst)
				#skip to next ROI
				if file[i+1][0].isnumeric():
					i=i+5
				else:
					i=i+4
		all_image_data.append(tmp_lst_3)


rois = []
#populate rois with all data
for j in range(len(all_image_data)):
	for i in range(len(all_image_data[j])):
		if not math.isnan(all_image_data[j][i][2][0]): 
			if all_image_data[j][i][3][0] != 0.0 and all_image_data[j][i][4][0] != 0.0:
				x = ROI(all_image_data[j][i])
				scaledpos=normalize_list(all_image_data[j][i][2])
				x.add_Scaled_Pos(scaledpos)
				scaledtumor=normalize_list(all_image_data[j][i][3])
				x.add_Scaled_Tumor(scaledtumor)
				scaledtcell=normalize_list(all_image_data[j][i][4])
				x.add_Scaled_Tcell(scaledtcell)
				T_pheno = find_T_pheno(x)
				x.add_T_pheno(T_pheno)
				Ex_ratio = find_ex_ratio(x)
				x.add_Ex_ratio(Ex_ratio)
				Core = find_core(x)
				x.add_Core(Core)
				rois.append(x)


###Generate tumors list with consolidated growth and imaging data
lines = []
for j in range(len(all_data)):
	tmp=[]
	for i in range(len(all_data[j])):
		if float(all_data[j][i][2]) > 1:
			x = Line(all_data[j][i])
			tmp.append(x)
	lines.append(tmp)


catchall = []
for h in range(len(lines)):
	tmp = []
	x=1
	z=1
	for i in range(len(lines[h])):
		if lines[h][i].TumorID == x: #true
			y=Tumor(lines[h][i])
			y.add_point(lines[h][i])
			for j in range(len(rois)):
				if rois[j].Cage==lines[h][i].Cage and rois[j].Mouse==lines[h][i].Mouse and rois[j].Ear==lines[h][i].Ear:
					if rois[j].Day==lines[h][i].Day and rois[j].ROI==lines[h][i].ROI:
						if len(y.ImageData)>0 and rois[j] == y.ImageData[-1]:
							continue
						else:
							y.add_image_data(rois[j])
			tmp.append(y)
			x+=1
			z=x-1
		elif lines[h][i].TumorID == z:
			tmp[z-1].add_point(lines[h][i])
			for j in range(len(rois)):
				if rois[j].Cage==lines[h][i].Cage and rois[j].Mouse==lines[h][i].Mouse and rois[j].Ear==lines[h][i].Ear:
					if rois[j].Day==lines[h][i].Day and rois[j].ROI==lines[h][i].ROI:
						#if rois[j] == tmp[z-1].ImageData[-1]:
						#	continue
						#else:
						tmp[z-1].add_image_data(rois[j])
			z=x-1
		else:
			continue
	catchall.append(tmp)



for i in range(len(catchall)):
	for j in range(len(catchall[i])):
		if catchall[i][j].Resolve == "TRUE":
			tmp_days=catchall[i][j].Days;
			index = tmpts.index(tmp_days[-1:][0])+1;
			tmp_days.append(tmpts[index])
			tmp_areas=catchall[i][j].Areas;
			tmp_areas.append(0.001)
			tmp_days.insert(0,implantation_day);
			tmp_areas.insert(0,0.001);
			get_dtr=get_days_to_resolve(catchall[i][j])
			catchall[i][j].add_days_to_resolve(get_dtr)
			get_res=get_respond(catchall[i][j])
			catchall[i][j].add_respond(get_res)
			get_dpt=get_days_post_treatment(catchall[i][j],treatment_initiated)
			catchall[i][j].add_days_post_treatment(get_dpt)
			tumors.append(catchall[i][j])
		elif catchall[i][j].Resolve == "FALSE":
			tmp_days=catchall[i][j].Days;
			tmp_areas=catchall[i][j].Areas;
			get_dtr=get_days_to_resolve(catchall[i][j])
			catchall[i][j].add_days_to_resolve(get_dtr)
			get_res=get_respond(catchall[i][j])
			catchall[i][j].add_respond(get_res)
			get_dpt=get_days_post_treatment(catchall[i][j],treatment_initiated)
			catchall[i][j].add_days_post_treatment(get_dpt)
			tmp_days.insert(0,implantation_day);
			tmp_areas.insert(0,0.001);
			tumors.append(catchall[i][j])


#####################################################################################
################################# 20191104 Data #####################################
#####################################################################################

#important information specific to the experiment!

tmpts = [7,8,9,10,11,13,14,15,16,17,19,20,21,28,30,33,35,37,40,42]
num_cages = 1
num_mice = 15
treatment_group = ['Control', 'Combo']
implantation_day=0
treatment_initiated=2

class Line():
	#This class is useful for reading in data from CSV files but is not used after that
	#Probably worth eliminating this class in the future
	def __init__(self, input):
		self.Label = input[1]
		self.Area = float(input[2])*(10**-6) #scaling factor for µm2 scaled 10x images to mm2
		self.ROI = int(input[1].split(":")[1].split("-")[0]);
		self.Day = int(input[1].split("D")[1].split("p")[0])
		self.Mean = float(input[3])
		self.Stdev = float(input[4])
		self.Min = float(input[5])
		self.Max = float(input[6])
		self.TumorID = int(input[7])
		self.Resolve = input[8]
		self.Respond = False
		self.Merge = input[9]
		self.Treatment = input[10]
		self.Cage = int(input[11])
		self.Mouse = int(input[12])
		self.Ear = input[13]
		self.Date = input[14]
	def get_Day(self, input):
		try: 
			day = int(input[1].split("D")[1].split("p")[0].split(".")[0])
			self.Day = day
		except: 
			day = int(input[1].split("D")[1].split("p")[0])
			self.Day = day



class ROI():
	def __init__(self, input):
		self.Day = int(input[0].split("D")[1].split("p")[0])
		self.Cage = int(1);
		self.Mouse = int(input[0].split("_m")[1].split("_")[0]);
		self.Ear = input[0].split("_e")[1].split("/")[0];
		self.ROI = int(input[1][0]);
		self.Area = float(input[1][1]);
		self.X = float(input[1][2]);
		self.Y = float(input[1][3]);
		self.Size = float(input[1][4]);
		self.Pos = input[2];
		self.Scaled_Pos = [];
		self.Tumor = input[3];
		self.Scaled_Tumor = [];
		self.Tcell = input[4];
		self.Scaled_Tcell = [];
		self.T_pheno = [];
		self.Ex_ratio = [];
		self.Core = [];
	def add_Scaled_Pos(self, input):
		self.Scaled_Pos = input;
	def add_T_pheno(self, input):
		self.T_pheno = input;
	def add_Ex_ratio(self, input):
		self.T_pheno = input;
	def add_Ex_ratio(self, input):
		self.Ex_ratio = input;
	def add_Core(self, input):
		self.Core = input;
	def add_Scaled_Tcell(self, input):
		self.Scaled_Tcell = input;
	def add_Scaled_Tumor(self, input):
		self.Scaled_Tumor = input;
	def __eq__(self, other): 
		if not isinstance(other, ROI):
			return NotImplemented
		return self.Day == other.Day and self.Cage == other.Cage and self.Mouse == other.Mouse and self.Ear==other.Ear and self.ROI == other.ROI

##Parse data into a list of Tumor objects called "tumors". Use tumors[0].__dict__ to see attributes

path = './'

all_csvs=glob.glob(path+"20191104_Lupe_selectedIMPLANTATION/Ear Folders/*.csv")

all_data=[]
for filename in all_csvs:
	with open(filename, 'r') as csvfile:
		file = list(csv.reader(csvfile))
		#print filename
		for i in range(len(file)):
			#Cage
			#file[i].insert(11,filename[:2].strip("C"))
			file[i].insert(11,"1")
			#Mouse
			file[i].insert(12,filename.split("/")[-1].split("_")[0].strip("m"))
			#Ear
			file[i].insert(13,filename[-5])
			file[i].insert(14,"20191104_Lupe_selectedIMPLANTATION")
		all_data.append(file[1:])


all_image_data=[]
# Set path for ROI T cell and tumor intensity data CSV files
for c in all_csvs:
	path_2 = c[:-4]+"/"
	print(path_2)
	#Read in T cell and tumor intensity data to ROI objects
	for filename in glob.glob(path_2+"/D*.csv"):
		with open(filename, 'r') as csvfile:
			file = list(csv.reader(csvfile))
			#print filename
			tmp_lst_3=[]
			i=0
			while i<len(file):
				#print("pass! "+ filename)
				tmp_lst=[]
				#add day
				tmp_lst.append(filename)
				#add ROI num, areas, T cell and Tumor intensities
				tmp_lst.append(file[i])
				#add scaled position
				for j in [1, 2, 3, 4]:
					tmp_lst_2=[]
					if i+j < len(file):
						for k in range(len(file[i+j])): 
							if file[i+j][k] == "":
								continue
							else:
								tmp_lst_2.append(float(file[i+j][k]))
						tmp_lst.append(tmp_lst_2)
				tmp_lst_3.append(tmp_lst)
				#skip to next ROI
				if file[i+1][0].isnumeric():
					i=i+5
				else:
					i=i+4
		all_image_data.append(tmp_lst_3)



rois = []
#populate rois with all data
for j in range(len(all_image_data)):
	for i in range(len(all_image_data[j])):
		if not math.isnan(all_image_data[j][i][2][0]): 
			if all_image_data[j][i][3][0] != 0.0 and all_image_data[j][i][4][0] != 0.0:
				x = ROI(all_image_data[j][i])
				scaledpos=normalize_list(all_image_data[j][i][2])
				x.add_Scaled_Pos(scaledpos)
				scaledtumor=normalize_list(all_image_data[j][i][3])
				x.add_Scaled_Tumor(scaledtumor)
				scaledtcell=normalize_list(all_image_data[j][i][4])
				x.add_Scaled_Tcell(scaledtcell)
				T_pheno = find_T_pheno(x)
				x.add_T_pheno(T_pheno)
				Ex_ratio = find_ex_ratio(x)
				x.add_Ex_ratio(Ex_ratio)
				Core = find_core(x)
				x.add_Core(Core)
				rois.append(x)


lines = []
for j in range(len(all_data)):
	tmp=[]
	for i in range(len(all_data[j])):
		if float(all_data[j][i][2]) > 1:
			x = Line(all_data[j][i])
			tmp.append(x)
	lines.append(tmp)



catchall = []
for h in range(len(lines)):
	tmp = []
	x=1
	z=1
	for i in range(len(lines[h])):
		if lines[h][i].TumorID == x: #true
			y=Tumor(lines[h][i])
			y.add_point(lines[h][i])
			for j in range(len(rois)):
				if rois[j].Cage==lines[h][i].Cage and rois[j].Mouse==lines[h][i].Mouse and rois[j].Ear==lines[h][i].Ear:
					if rois[j].Day==lines[h][i].Day and rois[j].ROI==lines[h][i].ROI:
						if len(y.ImageData)>0 and rois[j] == y.ImageData[-1]:
							continue
						else:
							y.add_image_data(rois[j])
			tmp.append(y)
			x+=1
			z=x-1
		elif lines[h][i].TumorID == z:
			tmp[z-1].add_point(lines[h][i])
			for j in range(len(rois)):
				if rois[j].Cage==lines[h][i].Cage and rois[j].Mouse==lines[h][i].Mouse and rois[j].Ear==lines[h][i].Ear:
					if rois[j].Day==lines[h][i].Day and rois[j].ROI==lines[h][i].ROI:
						#if rois[j] == tmp[z-1].ImageData[-1]:
						#	continue
						#else:
						tmp[z-1].add_image_data(rois[j])
			z=x-1
		else:
			continue
	catchall.append(tmp)


for i in range(len(catchall)):
	for j in range(len(catchall[i])):
		if catchall[i][j].Resolve == "TRUE":
			tmp_days=catchall[i][j].Days;
			index = tmpts.index(tmp_days[-1:][0])+1;
			tmp_days.append(tmpts[index])
			tmp_areas=catchall[i][j].Areas;
			tmp_areas.append(0.001)
			tmp_days.insert(0,implantation_day);
			tmp_areas.insert(0,0.001);
			get_dtr=get_days_to_resolve(catchall[i][j])
			catchall[i][j].add_days_to_resolve(get_dtr)
			get_dpt=get_days_post_treatment(catchall[i][j],treatment_initiated)
			catchall[i][j].add_days_post_treatment(get_dpt)
			get_res=get_respond(catchall[i][j])
			catchall[i][j].add_respond(get_res)
			tumors.append(catchall[i][j])
		elif catchall[i][j].Resolve == "FALSE":
			tmp_days=catchall[i][j].Days;
			tmp_areas=catchall[i][j].Areas;
			get_dtr=get_days_to_resolve(catchall[i][j])
			catchall[i][j].add_days_to_resolve(get_dtr)
			get_dpt=get_days_post_treatment(catchall[i][j],treatment_initiated)
			catchall[i][j].add_days_post_treatment(get_dpt)
			get_res=get_respond(catchall[i][j])
			catchall[i][j].add_respond(get_res)
			tmp_days.insert(0,implantation_day);
			tmp_areas.insert(0,0.001);
			tumors.append(catchall[i][j])


medians=[]
for i in range(len(tumors)):
	for j in range(len(tumors[i].ImageData)):
		if tumors[i].ImageData[j].Day in [6,7]:
			medians.append(np.median(tumors[i].ImageData[j].Tcell))

np.percentile(medians,25)



inflamed_medians=[]
for i in range(len(tumors)):
	for j in range(len(tumors[i].ImageData)):
		if tumors[i].ImageData[j].T_pheno in ["inflamed"]:
			core_index=int(len(tumors[i].ImageData[j].Tcell)*0.5)
			core_median=np.median(tumors[i].ImageData[j].Tcell[0:core_index])
			inflamed_medians.append(np.median(tumors[i].ImageData[j].Tcell))


np.percentile(inflamed_medians,25)
