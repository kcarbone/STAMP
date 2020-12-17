##########################################################################
######## Plot individual tumor growth curves for treatment groups ########
##########################################################################

##Generate lists based on treatment group with optional cutoff value

Control = []
aPD_L1 = []
aTGFb = []
Combo = []

for i in range(len(tumors)):
	if tumors[i].Treatment == 'Control': #and tumors[i].Date == '20190207_Lupe_kppgfp_WT':
		Control.append(tumors[i])
	elif tumors[i].Treatment == 'aPD-L1': #and tumors[i].Date == '20190207_Lupe_kppgfp_WT':
		aPD_L1.append(tumors[i])
	elif tumors[i].Treatment == 'aTGFb': #and tumors[i].Date == '20190207_Lupe_kppgfp_WT':
		aTGFb.append(tumors[i])
	elif tumors[i].Treatment == 'Combo': #and tumors[i].Date == '20190207_Lupe_kppgfp_WT':
		Combo.append(tumors[i])



to_plot = Combo
n=len(to_plot)

title = "Individual Tumor Growth: "+to_plot[1].Treatment+" (n="+str(n)+")"

fig=plt.figure()
ax=fig.add_subplot(1, 1, 1)
for c in range(len(to_plot)):
	tmp_tumor=to_plot[c]
	ax=fig.add_subplot(1, 1, 1)
	if tmp_tumor.Resolve == "FALSE":
		if tmp_tumor.Merge == "FALSE":
			_=ax.plot(tmp_tumor.Days, tmp_tumor.Areas, color='gray', alpha=0.5);
		elif tmp_tumor.Merge == "TRUE":
			_=ax.plot(tmp_tumor.Days, tmp_tumor.Areas, color='gray', alpha=0.5);
	elif tmp_tumor.Resolve == "TRUE":
		if tmp_tumor.Merge == "FALSE":
			_=ax.plot(tmp_tumor.Days, tmp_tumor.Areas, color='firebrick', alpha=0.5);
		elif tmp_tumor.Merge == "TRUE":
			_=ax.plot(tmp_tumor.Days, tmp_tumor.Areas, color='firebrick', alpha=0.5);

_=ax.set_yscale('log', basey=10);
_=ax.set_ylim(10**-3, 10**3);
#_=ax.set_xlim([10.5,24.5]);
_=plt.title(title);
_=plt.tight_layout();
_=plt.subplots_adjust(top=0.87)	
#plt.suptitle(title);
fig.show()
