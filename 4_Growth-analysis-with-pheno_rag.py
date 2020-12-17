##################################################################################
#### Plot tumors that meet a condition with scatter points colored by T_pheno ####
##################################################################################

if show_plot_5 == True:
	title = "Growth Curve (Combo)"
	to_plot = []
	for i in range(len(tumors)):
		if len(tumors[i].ImageData) >= 1:
			if tumors[i].Treatment == "Combo" and tumors[i].Resolve=='TRUE':
				to_plot.append(tumors[i])
	
	#Uncomment this to plot all tumors at once
	#to_plot = tumors
	n=len(to_plot)
	
	fig=plt.figure()
	for c in range(len(to_plot)):
		tmp_tumor=to_plot[c]
		ax=fig.add_subplot(1, 1, 1)
		if tmp_tumor.Resolve == "TRUE":
			if tmp_tumor.Merge == "FALSE":
				_=ax.plot(tmp_tumor.DaysToResolve, tmp_tumor.Areas, color='black', alpha=0.5, zorder=1, lw=0.5);
			elif tmp_tumor.Merge == "TRUE":
				_=ax.plot(tmp_tumor.DaysToResolve, tmp_tumor.Areas, color='black', alpha=0.5, zorder=1, lw=0.5);
		elif tmp_tumor.Resolve == "FALSE":
			if tmp_tumor.Merge == "FALSE":
				_=ax.plot(tmp_tumor.DaysToResolve, tmp_tumor.Areas, color='black', alpha=0.5, zorder=1, lw=0.5);
			elif tmp_tumor.Merge == "TRUE":
				_=ax.plot(tmp_tumor.DaysToResolve, tmp_tumor.Areas, color='black', alpha=0.5, zorder=1, lw=0.5);
	
	ax=fig.add_subplot(1, 1, 1)
	for c in range(len(to_plot)):
		tmp_tumor=to_plot[c]
		color="black"
		#_=ax.scatter(tmp_tumor.Days[-1], tmp_tumor.Areas[-1], color=color, alpha=0.3, zorder=2)
		for i in range(len(tmp_tumor.ImageData)+1):
			color="lightgray"
			if i > 0:
				if tmp_tumor.ImageData[i-1].T_pheno == "excluded":
					color="red"
				if tmp_tumor.ImageData[i-1].T_pheno == "inflamed":
					color="green"
				if tmp_tumor.ImageData[i-1].T_pheno == "desert":
					color="gold"
				if tmp_tumor.ImageData[i-1].T_pheno == "uncertain":
					color="lightgray"	
			_=ax.scatter(tmp_tumor.DaysToResolve[i], tmp_tumor.Areas[i], color=color, alpha=0.3, s=30, zorder=2)
	
	_=ax.set_yscale('log', basey=10);
	_=ax.set_ylim(10**-3, 10**3);
	#_=ax.set_xlim([min(tmpts)-1,max(tmpts)+1]);
	_=ax.set_xlim([-40,0]);
	#_=ax.set_xlim([0,40]);
	_=plt.title(title+" n=" +str(n));
	_=plt.tight_layout();
	_=plt.subplots_adjust(top=0.87)	
	#plt.suptitle(title);
	fig.show()


