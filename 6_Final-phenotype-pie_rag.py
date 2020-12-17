##########################################################################
#################### Plot final phenotype pie charts #####################
##########################################################################

phenos=["desert", "excluded", "inflamed"]
respond=[True, False]
treatment=["Control","aPD-L1","aTGFb","Combo"]

labels=[]
counts=[]

for c in treatment:
	for r in respond:
		for p in phenos:
			count=0
			for t in tumors:
				if len(t.ImageData)>=3:
				#ALSO INTERESTNG TO COMAPER TO t.ImageData[1]
					if t.Treatment==c and t.Respond==r and t.ImageData[-1].T_pheno==p:
						count+=1
			counts.append(count)
			labels.append(c+" "+str(r)+" "+p)
			print(c+"\t Respond: "+str(r)+"\t Phenotype: "+p+"\t Count = "+str(count))

control_total_r = counts[0:3]
control_total_nr = counts[3:6]
pdl1_total_r = counts[6:9]
pdl1_total_nr = counts[9:12]
tgfb_total_r = counts[12:15]
tgfb_total_nr = counts[15:18]
combo_total_r = counts[18:21]
combo_total_nr = counts[21:24]

pie_size=max(sum(control_total_r),sum(control_total_nr),
	sum(pdl1_total_r),sum(pdl1_total_nr),
	sum(tgfb_total_r),sum(tgfb_total_r),
	sum(combo_total_r),sum(combo_total_nr))

control_sizes_r = [x/sum(control_total_r)*360.0 for x in control_total_r]
control_sizes_nr = [x/sum(control_total_nr)*360.0 for x in control_total_nr]
tgfb_sizes_r = [x/sum(tgfb_total_r)*360.0 for x in tgfb_total_r]
tgfb_sizes_nr = [x/sum(tgfb_total_nr)*360.0 for x in tgfb_total_nr]
pdl1_sizes_r = [x/sum(pdl1_total_r)*360.0 for x in pdl1_total_r]
pdl1_sizes_nr = [x/sum(pdl1_total_nr)*360.0 for x in pdl1_total_nr]
combo_sizes_r = [x/sum(combo_total_r)*360.0 for x in combo_total_r]
combo_sizes_nr = [x/sum(combo_total_nr)*360.0 for x in combo_total_nr]

colors = ['gold','red', 'green']
explode = (0, 0, 0)  # explode 1st slice

sizes=[control_sizes_r, control_sizes_nr,
tgfb_sizes_r, tgfb_sizes_nr,
pdl1_sizes_r, pdl1_sizes_nr,
combo_sizes_r, combo_sizes_nr]

totals=[control_total_r, control_total_nr,
tgfb_total_r, tgfb_total_nr,
pdl1_total_r, pdl1_total_nr,
combo_total_r, combo_total_nr]

titles=["control_r", "control_nr",
"tgfb_r", "tgfb_nr",
"pdl1_r", "pdl1_nr",
"combo_r", "combo_nr"]

rc = {"axes.spines.left" : False,
      "axes.spines.right" : False,
      "axes.spines.bottom" : False,
      "axes.spines.top" : False,
      "xtick.bottom" : False,
      "xtick.labelbottom" : False,
      "ytick.labelleft" : False,
      "ytick.left" : False}
plt.rcParams.update(rc)

fig, ax = plt.subplots(4,2);
for i in range(len(sizes)):
	ax=fig.add_subplot(4, 2, i+1)
	ax.set_title(titles[i])
	_=ax.pie(sizes[i], radius=(sum(totals[i])/pie_size*1.5), explode=explode, colors=colors, shadow=False, startangle=140);

fig.show();