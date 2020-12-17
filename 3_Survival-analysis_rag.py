from lifelines import KaplanMeierFitter

######################################################################################
################  Output is a kaplan meyer curve of individual mice ##################
######################################################################################

#array of durations
T = []
#boolean array whether "death" was observed
E = []
#group descriptor, i.e. treatment
G = []

#colors: Combo
#C=['chartreuse',
#"palegreen",
#'darkgreen',
#"seagreen",
#"mediumspringgreen",
#'lawngreen',
#'lightgreen',
#'g',
#'mediumseagreen',
#"mediumaquamarine",
#"darkolivegreen",
#'forestgreen',
#"green",
#"springgreen",
#"aquamarine",
#"darkseagreen",
#'limegreen',
#'lime',
#"palegreen",
#'darkgreen',
#"seagreen",
#"mediumspringgreen",
#'lawngreen',
#'lightgreen',
#'g',
#'mediumseagreen',
#"mediumaquamarine",
#"darkolivegreen",
#'forestgreen',
#"green",
#"springgreen",
#"aquamarine",
#"darkseagreen",
#'limegreen',
#'lime',
#"red"]

#colors: Control
#C = ['darkturquoise',
#'deepskyblue',
#'royalblue',
#'navy',
#'blue',
#'skyblue',
#'dodgerblue',
#'darkblue',
#'slateblue',
#'mediumblue',
#'darkslateblue',
#'teal',
#'cyan',
#'lightblue',
#'steelblue',
#'cornflowerblue',
#'midnightblue',
#'b',
#'mediumslateblue',
#"red"]

#colors: aPD-L1
#C = ['mediumpurple',
#'darkorchid',
#'plum',
#'m',
#'mediumvioletred',
#'rebeccapurple',
#'darkviolet',
#'violet',
#'fuchsia',
#'deeppink',
#'blueviolet',
#'mediumorchid',
#'purple',
#'magenta',
#'hotpink',
#'pink',
#'indigo',
#'thistle',
#'darkmagenta',
#'orchid',
#"red"]

#colors: aTGFb
C= ['rosybrown',
"firebrick",
"red", 
'darksalmon',
'sienna',
'sandybrown',
'gold',
'lightcoral',
'maroon',
'coral',
'darkorange',
"orange", 
'crimson',
'indianred',
'darkred',
'salmon',
"orangered", 
'goldenrod',
'brown',
'r',
'tomato',
"lightsalmon",
'saddlebrown']

colors=[]

for t in range(len(tumors)):
	if tumors[t].Treatment == "Combo":
		T.append(max(tumors[t].DaysPostTreatment))
		G.append(str(tumors[t].Cage)+" "+str(tumors[t].Mouse)+" "+tumors[t].Date) #+tumors[t].Ear
		if tumors[t].Resolve == 'TRUE':
			E.append(True)
		else:
			E.append(False)

df = pd.DataFrame(list(zip(T, E, G)), columns =['T', 'E', 'G']) 

ax = plt.subplot(111)

kmf = KaplanMeierFitter()

for name, grouped_df in df.groupby('G'):
    kmf.fit(grouped_df["T"], grouped_df["E"], label=name)
    kmf.plot(ax=ax, ci_show=False, ci_legend=False)
    for i in range(len(plt.gca().get_lines())):
    	plt.gca().get_lines()[i].set_color(C[i])

#plt.legend(title="Mouse",frameon=False)
plt.legend("")
plt.xlim(0,42)
plt.ylim(0,1.1)
plt.show()

##############################################################################
#############  Kaplan meyer curve to compare treatment groups ################
##############################################################################


#array of durations
T = []
#boolean array whether "death" was observed
E = []
#group descriptor, i.e. treatment
G = []

#colors
C=["firebrick","navy","gold","orangered"]
#C=["navy","blue",'steelblue','royalblue']
#C= ["firebrick","red",'crimson','maroon']
#C=['goldenrod','gold']
#C=['orangered','darkorange']

#for t in range(len(tumors)):
#	T.append(max(tumors[t].DaysPostTreatment))
#	G.append(tumors[t].Treatment)
#	if tumors[t].Resolve == "TRUE":
#		E.append(True)
#	else:
#		E.append(False)


# Euthanized mice are coded as persisiting tumors 

for t in range(len(tumors)):
	if tumors[t].Resolve == "TRUE":
		T.append(max(tumors[t].Days))
		E.append(True)
		G.append(tumors[t].Treatment)
	if tumors[t].Resolve == "FALSE":
		T.append(42)
		E.append(False)
		G.append(tumors[t].Treatment)

df = pd.DataFrame(list(zip(T, E, G)), columns =['T', 'E', 'G']) 

ax = plt.subplot(111)

kmf = KaplanMeierFitter()

i=0
for name, grouped_df in df.groupby('G'):
    kmf.fit(grouped_df["T"], grouped_df["E"], label=name)
    kmf.plot(ax=ax, ci_show=True, ci_legend=False, color=C[i])
    i+=1

plt.legend(title="Treatment",frameon=False)
plt.xlim(0,42)
plt.ylim(0,1.1)
plt.show()



######Stats:

from lifelines.statistics import logrank_test


ex = df['G'] == 'Combo'
con = df['G'] == 'Control'
T_exp, E_exp = df.loc[ex, 'T'], df.loc[ex, 'E']
#T_con, E_con = df.loc[~ex, 'T'], df.loc[~ex, 'E']
T_con, E_con = df.loc[con, 'T'], df.loc[con, 'E']

results = logrank_test(T_exp, T_con, event_observed_A=E_exp, event_observed_B=E_con)
results.print_summary()
print(results.p_value)        
print(results.test_statistic) 

