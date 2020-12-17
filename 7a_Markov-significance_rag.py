
import quantecon as qe

phenos=["desert", "excluded", "inflamed","resolved","euthanized"]
tmpts = [7,10,11,12,13,15,17,18,19,20,21,24,25,26,28,29]
tmpts_index=17

transitions=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Combo","Control"]: #and tumors[i].Respond in [True, False]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			#print(tumors[i].ImageData[j].T_pheno)
			if j>3 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-3].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-4].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-4].T_pheno)
			elif j>2 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-3].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-3].T_pheno)
			elif j>1 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-2].T_pheno)
			elif j>0 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-1].T_pheno)
			elif j<2 and tumors[i].ImageData[j].T_pheno in ["uncertain",None]:
				score=phenos.index("desert")
			elif tumors[i].ImageData[j].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j].T_pheno)
			if math.isnan(float(score)):
				print(score)
			tmp.append(score)
		if tmp != [] and len(tmp)<tmpts_index and tumors[i].Resolve=="FALSE":
			#to_add=[tmp[-1]] * (tmpts_index-len(tmp))
			to_add=[4]
			tmp=tmp+to_add
			transitions.append(tmp)
			tag.append(tumors[i].Treatment+" "+str(tumors[i].Resolve))
		elif tmp != [] and len(tmp)<tmpts_index and tumors[i].Resolve=="TRUE":
			#to_add=[tmp[-1]] * (tmpts_index-len(tmp))
			to_add=[3]
			tmp=tmp+to_add
			transitions.append(tmp)
			tag.append(tumors[i].Treatment+" "+str(tumors[i].Resolve))


############ Generate Transition Matrix ################

all_vals=[]
for i in range(len(transitions)):
	for j in range(len(transitions[i])):
		all_vals.append(transitions[i][j])


import random

random.shuffle(all_vals)

T=all_vals

#create matrix of zeros

M = [[0]*5 for _ in range(5)]

for (i,j) in zip(T,T[1:]):
    M[i][j] += 1

#now convert to probabilities:
for row in M:
    n = sum(row)
    if n > 0:
        row[:] = [f/sum(row) for f in row]

#print M:

for row in M:
    print(row)

val=0
for i in range(len(M)):
	for j in range(len(M[i])):
		val+=M[i][j]


#########################################################
################### Show Scatter Plot ###################
#########################################################

trans=["D-D","D-Ex","D-I","D-R","D-Eu","Ex-D","Ex-Ex","Ex-I","Ex-R","Ex-Eu","I-D","I-Ex","I-I","I-R","I-Eu"]
cols=["CR_PR","SD_PD","Combo","Control","R1","R2","R3","R4","R5","R6","R6","R8","R9","R10"]
CR_PR=[0.375745527,0.280318091,0.302186879,0.039761431,0.001988072,0.01854067,0.666866029,0.219497608,0.088516746,0.006578947,0.013678554,0.123595506,0.684904739,0.153395213,0.024425989]
SD_PD=[0.517160686,0.197347894,0.198907956,0,0.086583463,0.076376554,0.708703375,0.147868561,0,0.06705151,0.064397906,0.14921466,0.642931937,0,0.143455497]
Combo=[0.360,0.298,0.285,0.013,0.044,0.036,0.617,0.285,0.031,0.031,0.016,0.120,0.702,0.115,0.472]
Control=[0.324,0.317,0.278,0.007,0.073,0.039,0.699,0.193,0.016,0.053,0.026,0.148,0.684,0.037,0.105]
R1=[0.173669468,0.356862745,0.375910364,0.043697479,0.049859944,0.171552383,0.371144532,0.351771603,0.04690288,0.058628601,0.163507708,0.364417488,0.373768006,0.04447814,0.053828658]
R2=[0.161904762,0.366386555,0.369187675,0.051540616,0.050980392,0.170532756,0.36375223,0.364007137,0.043079276,0.058628601,0.160222391,0.369471822,0.367955522,0.045489007,0.056861259]
R3=[0.178711485,0.346218487,0.377030812,0.040896359,0.057142857,0.154434251,0.366462793,0.377166157,0.04841998,0.05351682,0.170626896,0.372851365,0.355915066,0.042214358,0.058392315]
R4=[0.166946779,0.349579832,0.385994398,0.040896359,0.056582633,0.170744139,0.364169215,0.367227319,0.04587156,0.051987768,0.158240647,0.375631951,0.362992922,0.044489383,0.058645096] 
R5=[0.169747899,0.366946779,0.370308123,0.040896359,0.05210084,0.162079511,0.361620795,0.372833843,0.049949032,0.05351682,0.162032356,0.362740142,0.37057634,0.043983822,0.060667341]
R6=[0.166946779,0.367507003,0.363585434,0.044257703,0.057703081,0.164373089,0.362640163,0.370540265,0.044852192,0.057594292,0.166329626,0.363245703,0.37512639,0.041456016,0.053842265]
R7=[0.157983193,0.354621849,0.37254902,0.053781513,0.061064426,0.162120826,0.371144532,0.367320928,0.040785113,0.058628601,0.171341926,0.365428355,0.365681072,0.047763457,0.049785191]
R8=[0.165826331,0.350140056,0.379831933,0.047619048,0.056582633,0.158766565,0.374362895,0.363149847,0.049184506,0.054536188,0.169362993,0.367795753,0.366784631,0.040192113,0.05586451]
R9=[0.17535014,0.37535014,0.359103641,0.034733894,0.055462185,0.162589195,0.368501529,0.374362895,0.044342508,0.050203874,0.166835187,0.363751264,0.360465116,0.048281092,0.060667341]
R10=[0.173109244,0.382633053,0.347338936,0.044257703,0.052661064,0.163608563,0.363404689,0.369775739,0.046381244,0.056829766,0.165065723,0.358442872,0.375631951,0.043983822,0.056875632]

l=[CR_PR,SD_PD,Combo,Control,R1,R2,R3,R4,R5,R6,R6,R8,R9,R10]

vals=[]
series=[]
labels=[]

for i in range(len(l)):
	for j in range(len(l[i])):
		vals.append(l[i][j])
		series.append(trans[j])
		labels.append(cols[i])

for_df=np.array([vals,series,labels])

DF = pd.DataFrame(data=for_df.transpose(), columns=["Value","Transition","Sample"])
DF['Value'] = pd.to_numeric(DF['Value'])

colors=["cyan","red","magenta","blue","lightgray","lightgray","lightgray","lightgray","lightgray","lightgray","lightgray","lightgray","lightgray"]

ax=sns.stripplot(data=DF, x="Transition", y="Value", hue="Sample", jitter=0.3, alpha=0.8, linewidth=0.5, palette=colors)
plt.show()


