#%%
from tkinter import TRUE
from scipy.stats import spearmanr
import os
import sys
import pandas as pd
sys.path.append('../')
sys.path.append('./')
from utils.load_data import dataset
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import Formula
import argparse
import pingouin as pg
import matplotlib.pyplot as plt
from statsmodels.stats import multitest
import numpy as np
import seaborn as sns

# pre setting
PWD = '/data/CC_chort_analyses_basic'
ds = dataset()
def pd2ri(df):
	with localconverter(ro.default_converter + pandas2ri.converter):
		r_df = ro.conversion.py2rpy(df)
	return r_df
def ri2pd(df):
	with localconverter(ro.default_converter + pandas2ri.converter):
		r_df = ro.conversion.rpy2py(df)
	return r_df
maaslin2 = importr('Maaslin2')
base = importr('base')
TMP = 'tmp_paired'
# end

# case-control time point 
def main(level):
	output_root_par = os.path.join(PWD,'Section_A')
	os.system(f'mkdir {output_root_par}')
	output_root = os.path.join(output_root_par,'Q6_Masslin2_case_control')
	os.system(f'mkdir {output_root}')
	raw_output = os.path.join(output_root,f'{level}_raw') 
	raw_output2 = os.path.join(output_root,f'{level}_raw2') 
	plots_output = os.path.join(output_root,f'{level}_plots')
	lm_plots_output = os.path.join(output_root,f'{level}_lm_plots')
	
	os.system(f'mkdir {plots_output}')
	os.system(f'mkdir {lm_plots_output}')
	
	if level == 'genus':
		df = ds.genus
	elif level == 'asv':
		df = ds.asv
	elif level == 'phylum':
		df = ds.phylum
	meta_df = ds.meta[['case','age_days_npp_collection','antibiotics_use_at_collection_analysis']]
	meta_df['age_days_npp_collection'] = meta_df['age_days_npp_collection'].astype(int)
	df_r = pd2ri(df)
	meta_df_r = pd2ri(meta_df)
	#run model 1
	all_vars = ['case','age_days_npp_collection','antibiotics_use_at_collection_analysis']
	maaslin2.Maaslin2(input_data = df_r,
					input_metadata = meta_df_r,
					output = raw_output,
					normalization='TSS',
					transform='LOG',
					plot_scatter='FALSE',
					fixed_effects = ro.StrVector(all_vars),
					reference = ro.StrVector(['antibiotics_use_at_collection_analysis','no']))
	res = pd.read_csv(os.path.join(raw_output,'significant_results.tsv'),sep='\t',index_col=0)
	#res = res[res.metadata=='case']
	if res.shape[0]==0:
		print('No significant results detected')
		return 
	res['FDR'] = multitest.multipletests(res['pval'],method='fdr_bh')[1]
	res.to_csv(os.path.join(output_root,f'{level}_sig.csv'))
	#sel significant taxa	
	sig = res.index.tolist()
	tax_lst = sig
	## plot
	# change to relative abundance
	for c in df.columns:
		ser = df[c]
		ser_sum = ser.sum()
		ser = ser.apply(lambda x:x/ser_sum)
		df[c] = ser
	
	## draw box plots
	for t in tax_lst:
		df_t_nan = df.loc[t].copy().to_frame()
		#df_t_nan = df_t_nan.replace(0., np.nan)
		plot_df = pd.concat([meta_df, df_t_nan], axis=1)
		sns.violinplot(data=plot_df, x='case',y=t,palette={'Positive':"#E97B69",'Negative':"#8FA3D1"})
		sns.stripplot(data=plot_df,x='case',y=t,color='black',alpha=0.25)	
		#sns.stripplot(x="case", y=t, data=plot_df, color=".25")
		plt.tight_layout()
		plt.savefig(os.path.join(plots_output, f'{t}.svg'))
		plt.close()
	## draw scatter plots
	lab = ds.lab
	lab['16S_copies_log10'] = [np.log10(x) for x in lab.templateDNA_16SqPCR_copies]
	lab = lab[['16S_copies_log10']]
	for t in tax_lst:
		# first plot
		df_t = df.loc[t].to_frame()
		plot_df = pd.concat([meta_df, df_t, lab], axis=1)
		sns.lmplot(x='16S_copies_log10',y=t,hue='case', data=plot_df,order=2,\
			scatter_kws={"s": 10}, \
				palette={'Positive':"#E97B69",'Negative':"#8FA3D1"})
		corr = spearmanr(plot_df['16S_copies_log10'],plot_df[t])[0]
		plt.title('Spearman = {}'.format(round(corr,3)))
		#plt.legend()
		plt.tight_layout()
		plt.savefig(os.path.join(lm_plots_output, f'{t}_16S_copy.svg'))
		plt.close()
		### second plot	
		sns.lmplot(x='age_days_npp_collection',y=t,hue='case', data=plot_df,order=2,
				scatter_kws={"s": 10}, \
				palette={'Positive':"#E97B69",'Negative':"#8FA3D1"})
		corr = spearmanr(plot_df['age_days_npp_collection'],plot_df[t])[0]
		plt.title('Spearman = {}'.format(round(corr,3)))
		#plt.legend()
		plt.tight_layout()
		plt.savefig(os.path.join(lm_plots_output, f'{t}_age.svg'))
		plt.close()
	return

if __name__ == "__main__":
	#main('phylum')	
	#main('genus')	
	main('asv')	
