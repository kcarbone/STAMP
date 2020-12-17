####################################################################################
######################### Generate data lists for fitting ##########################
####################################################################################

phenos=["desert", "excluded", "inflamed"]
tmpts = [6,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,26,27,30,32,34,37,39,40]
tmpts_index=24


tmp_list_NR=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Combo"] and tumors[i].Respond in [False]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			score=tumors[i].ImageData[j].Ex_ratio
			tmp.append(score)
		tmp_list_NR.append(tmp)

x_vals_NR=[]
y_vals_NR=[]
d_vals_NR=[]
i_vals_NR=[]
for i in range(len(tmp_list_NR)):
	if len(tmp_list_NR[i]) > 1:
		t_step=100/(len(tmp_list_NR[i])-1)
		tmp_d=[]
		tmp_i=[]
		for j in range(len(tmp_list_NR[i])):
			x_vals_NR.append(j*t_step)
			y_vals_NR.append(tmp_list_NR[i][j])
			tmp_d.append(j*t_step)
			tmp_i.append(tmp_list_NR[i][j])
		d_vals_NR.append(tmp_d)
		i_vals_NR.append(tmp_i)



tmp_list_R=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Combo"] and tumors[i].Respond in [True]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			score=tumors[i].ImageData[j].Ex_ratio
			tmp.append(score)
		tmp_list_R.append(tmp)

x_vals_R=[]
y_vals_R=[]
d_vals_R=[]
i_vals_R=[]
for i in range(len(tmp_list_R)):
	if len(tmp_list_R[i]) > 1:
		t_step=100/(len(tmp_list_R[i])-1)
		tmp_d=[]
		tmp_i=[]
		for j in range(len(tmp_list_R[i])):
			x_vals_R.append(j*t_step)
			y_vals_R.append(tmp_list_R[i][j])
			tmp_d.append(j*t_step)
			tmp_i.append(tmp_list_R[i][j])
		d_vals_R.append(tmp_d)
		i_vals_R.append(tmp_i)


tmp_list_cNR=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Control"] and tumors[i].Respond in [False]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			score=tumors[i].ImageData[j].Ex_ratio
			tmp.append(score)
		tmp_list_cNR.append(tmp)


x_vals_cNR=[]
y_vals_cNR=[]
d_vals_cNR=[]
i_vals_cNR=[]
for i in range(len(tmp_list_cNR)):
	if len(tmp_list_cNR[i]) > 1:
		t_step=100/(len(tmp_list_cNR[i])-1)
		tmp_d=[]
		tmp_i=[]
		for j in range(len(tmp_list_cNR[i])):
			x_vals_cNR.append(j*t_step)
			y_vals_cNR.append(tmp_list_cNR[i][j])
			tmp_d.append(j*t_step)
			tmp_i.append(tmp_list_cNR[i][j])
		d_vals_cNR.append(tmp_d)
		i_vals_cNR.append(tmp_i)




tmp_list_cR=[]
for i in range(len(tumors)):
	if tumors[i].Treatment in ["Control"] and tumors[i].Respond in [True]:
		tmp=[]
		for j in range(len(tumors[i].ImageData)):
			score=tumors[i].ImageData[j].Ex_ratio
			tmp.append(score)
		tmp_list_cR.append(tmp)


x_vals_cR=[]
y_vals_cR=[]
d_vals_cR=[]
i_vals_cR=[]
for i in range(len(tmp_list_cR)):
	if len(tmp_list_cR[i]) > 1:
		t_step=100/(len(tmp_list_cR[i])-1)
		tmp_d=[]
		tmp_i=[]
		for j in range(len(tmp_list_cR[i])):
			x_vals_cR.append(j*t_step)
			y_vals_cR.append(tmp_list_cR[i][j])
			tmp_d.append(j*t_step)
			tmp_i.append(tmp_list_cR[i][j])
		d_vals_cR.append(tmp_d)
		i_vals_cR.append(tmp_i)

################################################################################
#################### Plot spline firt for tumor trajectories ###################
################################################################################

from scipy.interpolate import interp1d

cR_splines=[]
fig, ax=plt.subplots()
for i in range(len(d_vals_cR)):
	offset=0
	#sns.lineplot(d_vals_cR[i+offset],i_vals_cR[i+offset], linewidth=1,alpha=0.2, color="blue")
	if len(d_vals_cR[i+offset])>3:
		fit=interp1d(d_vals_cR[i+offset], i_vals_cR[i+offset], kind='cubic')
	else:
		fit=interp1d(d_vals_cR[i+offset], i_vals_cR[i+offset], kind='linear')
	xnew = np.linspace(0, 100, num=500, endpoint=True)
	ynew=fit(xnew)
	plt.plot(xnew,ynew,linewidth=1,alpha=0.1,color="blue")
	cR_splines.append(ynew.tolist())

#fig.show()

cR_avg_spline=[]
for i in range(500):
	val=np.median([item[i] for item in cR_splines])
	cR_avg_spline.append(val)

#fig, ax=plt.subplots()
plt.plot(xnew,cR_avg_spline,linewidth=2,alpha=1,color="navy")
fig.show()




from scipy.interpolate import interp1d

R_splines=[]
fig, ax=plt.subplots()
for i in range(len(d_vals_R)):
	offset=0
	#sns.lineplot(d_vals_R[i+offset],i_vals_R[i+offset], linewidth=1,alpha=0.2)
	if len(d_vals_R[i+offset])>3:
		fit=interp1d(d_vals_R[i+offset], i_vals_R[i+offset], kind='cubic')
	else:
		fit=interp1d(d_vals_R[i+offset], i_vals_R[i+offset], kind='linear')
	xnew = np.linspace(0, 100, num=500, endpoint=True)
	ynew=fit(xnew)
	plt.plot(xnew,ynew,linewidth=1,alpha=0.1,color="cyan")
	R_splines.append(ynew.tolist())

#fig.show()

R_avg_spline=[]
for i in range(500):
	val=np.median([item[i] for item in R_splines])
	R_avg_spline.append(val)

#fig, ax=plt.subplots()
plt.plot(xnew,R_avg_spline,linewidth=2,alpha=1, color="darkturquoise")
fig.show()



from scipy.interpolate import interp1d

NR_splines=[]
fig, ax=plt.subplots()
for i in range(len(d_vals_NR)):
	offset=0
	#sns.lineplot(d_vals_NR[i+offset],i_vals_NR[i+offset], linewidth=1,alpha=0.2)
	if len(d_vals_NR[i+offset])>3:
		fit=interp1d(d_vals_NR[i+offset], i_vals_NR[i+offset], kind='cubic')
	else:
		fit=interp1d(d_vals_NR[i+offset], i_vals_NR[i+offset], kind='linear')
	xnew = np.linspace(0, 100, num=500, endpoint=True)
	ynew=fit(xnew)
	plt.plot(xnew,ynew,linewidth=1,alpha=0.1,color="magenta")
	NR_splines.append(ynew.tolist())

#fig.show()

NR_avg_spline=[]
for i in range(500):
	val=np.median([item[i] for item in NR_splines])
	NR_avg_spline.append(val)

#fig, ax=plt.subplots()
plt.plot(xnew,NR_avg_spline,linewidth=2,alpha=1, color="mediumvioletred")
fig.show()




from scipy.interpolate import interp1d

cNR_splines=[]
fig, ax=plt.subplots()
for i in range(len(d_vals_cNR)):
	offset=0
	#sns.lineplot(d_vals_cNR[i+offset],i_vals_cNR[i+offset], linewidth=1,alpha=0.2)
	if len(d_vals_cNR[i+offset])>3:
		fit=interp1d(d_vals_cNR[i+offset], i_vals_cNR[i+offset], kind='cubic')
	else:
		fit=interp1d(d_vals_cNR[i+offset], i_vals_cNR[i+offset], kind='linear')
	xnew = np.linspace(0, 100, num=500, endpoint=True)
	ynew=fit(xnew)
	plt.plot(xnew,ynew,linewidth=1,alpha=0.1, color="red")
	cNR_splines.append(ynew.tolist())

#fig.show()

cNR_avg_spline=[]
for i in range(500):
	val=np.median([item[i] for item in cNR_splines])
	cNR_avg_spline.append(val)

#fig, ax=plt.subplots()
plt.plot(xnew,cNR_avg_spline,linewidth=2,alpha=1,color="firebrick")
fig.show()




fig, ax=plt.subplots()
ax.set_ylim([0,1])
ax.plot(xnew, cR_avg_spline,linewidth=2,alpha=0.5,color="b")
ax.plot(xnew, cNR_avg_spline,linewidth=2,alpha=0.5,color="r")
ax.plot(xnew, R_avg_spline,linewidth=2,alpha=0.5,color="c")
ax.plot(xnew, NR_avg_spline,linewidth=2,alpha=0.5,color="m")
fig.show()