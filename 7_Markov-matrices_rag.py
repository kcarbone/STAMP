##########################################################################
################### Create markov transition matrices ####################
##########################################################################

phenos=["desert", "excluded", "inflamed","resolved","euthanized"]
tmpts = [6,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,26,27,30,32,34,37,39,40]
tmpts_index=24

transitions=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Control","Combo"] and tumors[i].Respond in [False]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			#print(tumors[i].ImageData[j].T_pheno)
			if j>0 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-3].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-4].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-4].T_pheno)
			elif j>0 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-3].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-3].T_pheno)
			elif j>0 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-2].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-2].T_pheno)
			elif j>0 and tumors[i].ImageData[j].T_pheno in ["uncertain",None] and tumors[i].ImageData[j-1].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j-1].T_pheno)
			elif j<2 and tumors[i].ImageData[j].T_pheno in ["uncertain",None]:
				score=phenos.index("desert")
			elif tumors[i].ImageData[j].T_pheno not in ["uncertain", None]:
				score=phenos.index(tumors[i].ImageData[j].T_pheno)
			value=score
			tmp.append(value)
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

#Check that markov probabilities sum to == 1

val=0
for i in range(len(M)):
	for j in range(len(M[i])):
		val+=M[i][j]

#Edited M
m_R=[[0.3757455268389662, 0.2803180914512923, 0.30218687872763417, 0.039761431411530816, 0.0019880715705765406],
[0.01854066985645933, 0.6668660287081339, 0.21949760765550239, 0.08851674641148326, 0.006578947368421052],
[0.013678553981436248, 0.12359550561797752, 0.684904738641915, 0.1533952125061065, 0.024425989252564728],
[0.0, 0.0, 0.0, 1.0, 0.0], 
[0.0, 0.0, 0.0, 0.0, 1.0]]

m_NR=[[0.517160686427457, 0.19734789391575663, 0.19890795631825273, 0.0, 0.08658346333853355],
[0.0763765541740675, 0.7087033747779752, 0.14786856127886322, 0.0, 0.06705150976909414],
[0.06439790575916231, 0.14921465968586387, 0.6429319371727749, 0.0, 0.14345549738219895],
[0.0, 0.0, 0.0, 1.0, 0.0], 
[0.0, 0.0, 0.0, 0.0, 1.0]]


a_R=np.array(m_R)
a_NR=np.array(m_NR)

rat=a_R/a_NR



