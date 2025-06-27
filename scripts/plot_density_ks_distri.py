import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if len(sys.argv)!=4: #if the input arguments not 3, showing the usage.
	print("Please try this! Usage:python3 plot_density_ks_distri.py <input XX.synteny.blocks.ks.info> <name of the plot> <name of the title> ")
	sys.exit()


df = pd.read_csv(sys.argv[1],delimiter ='\t',header = 0).fillna(0)

sns.histplot(df['Average Ks'], kde=True, binwidth=0.05, color = 'grey', stat="density", label = "Ks distribution")
# To show the rug plot by rug =True
#sns.distplot(df.loc[df['Average Ks']>0]['Average Ks'], rug =True,bins = 40,hist= True, kde= True, color = 'grey',hist_kws= {'edgecolor':'black'},label = "Ks distribution")

plt.xlabel('Ks')
plt.ylabel('Density')
plt.title(sys.argv[3])
plt.savefig(sys.argv[2], format="pdf")


