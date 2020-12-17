############################################################
####### Compare maximum tumor size for C.R. vs. P.R. #######
############################################################

from statannot import add_stat_annotation

control_R_max_area=[]
control_all_max_area=[]
combo_R_max_area=[]
combo_all_max_area=[]
for i in range(len(tumors)):
	if tumors[i].Respond in [True] and tumors[i].Treatment=="Control":
		control_R_max_area.append(max(tumors[i].Areas))
	if tumors[i].Respond in [False] and tumors[i].Treatment=="Control":
		control_all_max_area.append(max(tumors[i].Areas))
	if tumors[i].Respond in [True] and tumors[i].Treatment=="Combo":
		combo_R_max_area.append(max(tumors[i].Areas))
	if tumors[i].Respond in [False] and tumors[i].Treatment=="Combo":
		combo_all_max_area.append(max(tumors[i].Areas))

lst1=control_all_max_area+control_R_max_area+combo_all_max_area+combo_R_max_area
lst2=["Control S.D. + P.D."]*len(control_all_max_area)+["Control C.R. + P.R."]*len(control_R_max_area)+["Combo S.D. + P.D."]*len(combo_all_max_area)+["Combo C.R. + P.R."]*len(combo_R_max_area)

df = pd.DataFrame(list(zip(lst2,lst1)))

box_pairs=[
(("Control C.R. + P.R."),("Control S.D. + P.D.")),
(("Combo C.R. + P.R."),("Combo S.D. + P.D.")),
(("Combo C.R. + P.R."),("Control C.R. + P.R.")),
(("Combo S.D. + P.D."),("Control S.D. + P.D.")),
]

fig, ax = plt.subplots(1, 1)
#ax.set_ylim([0,100])
ax.set_yscale('log', basey=10)
sns.boxplot(x=df[0],y=df[1],dodge=True, fliersize=0,palette=["r","b","magenta","cyan"])
#statistical annotation (text_format: 'star', 'full', 'simple', test='t-test_ind', 'Mann-Whitney'
add_stat_annotation(ax, x=df[0],y=df[1],
                    box_pairs=box_pairs,
                    test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
sns.stripplot(x=df[0],y=df[1], dodge=True, jitter=0.3, size=2, alpha=1,palette=["r","b","magenta","cyan"])
#ax.set_ylim([0,100])
ax.set_yscale('log', basey=10)
fig.show()
