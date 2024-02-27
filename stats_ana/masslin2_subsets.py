#%%
import os
import sys
sys.path.append('../')
sys.path.append('./')
from utils.load_data import dataset
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import Formula
import argparse
import pandas as pd
import pingouin as pg
import seaborn as sns
import math

import matplotlib.pyplot as plt
import numpy as np

CF = 0.1
# pre setting
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
# end

# case-control time point 
case_toi = ['Adenovirus',	'Influenza_A_B_C',	'Parainfluenza_1_2_3_4',\
				'RSV_A_B',	'Rhinovirus',	'Cytomegalovirus',	'S_pneumoniae',
				'H_influenzae_spp',	'S_aureus',	'M_catarrhalis',	'H_influenzae_B',	'k_pneumoniae']
 
#%%
def main(level,y,case_control):
	
	ds = dataset()
	if level == 'genus':
		df = ds.genus
	elif level == 'phylum':
		df = ds.phylum
	elif level == 'asv':
		df = ds.asv
	else:
		print('Wrong level')
	
	meta = ds.meta
	meta = meta[meta[y]==case_control]
	df = df[meta.index]
	# define y
	if y in ['H_influenzae_']:
		ser = meta.loc[:,meta.columns.str.startswith(y)]
		ser.replace('Positive',1,inplace=True)
		ser.replace('Negative',0,inplace=True)
		ser_sum = ser.apply(np.sum, axis=1)
		ser_sum = pd.Series(['Positive' if x > 0 else 'Negative' for x in ser_sum], 
		index=ser_sum.index,name=ser_sum.name)
		assert meta.index.tolist() == ser_sum.index.tolist()
		
		meta[f'{y}combined'] = ser_sum
		y = f'{y}combined'
		print(meta[y])
	if y in ['S_pneumoniae_abs']:
		meta['S_pneumoniae_abs'] = meta['S_pneumoniae_abs'].astype(float)

	meta_cols = ['age_days_npp_collection',y,'antibiotics_use_at_collection_analysis','case']
	meta_df = meta[meta_cols]
	meta_df['age_days_npp_collection'] = meta_df['age_days_npp_collection'].astype(int)
	
	# define output
	output_root_par = 'Section_B'
	os.system(f'mkdir {output_root_par}')
	output_root = os.path.join(output_root_par,f'New_ana_Masslin2_FTD_{level}_{CF}')
	os.system(f'mkdir {output_root}')
	
	output_root_case = os.path.join(output_root, y) 
	os.system(f'mkdir {output_root_case}')
	
	output_y = os.path.join(output_root_case, case_control)
	os.system(f'mkdir {output_y}')
	
	raw_output = os.path.join(output_y, 'masslin2_raw')
	
	plots_output = os.path.join(output_y, f'{level}_plots')
	os.system(f'mkdir {plots_output}')
	
	# model
	df_r = pd2ri(df)
	meta_df_r = pd2ri(meta_df)

	
	if len(meta_df['antibiotics_use_at_collection_analysis'].unique()) != 1:
		all_vars = ['case','age_days_npp_collection','antibiotics_use_at_collection_analysis']
		maaslin2.Maaslin2(input_data = df_r,
					input_metadata = meta_df_r,
					output = raw_output,
					normalization='TSS',
					transform='LOG',
					plot_scatter = 'FALSE',
					min_prevalence=CF,
					fixed_effects = ro.StrVector(all_vars),
					reference = ro.StrVector(['antibiotics_use_at_collection_analysis','no']))
	else:
		all_vars = ['case','age_days_npp_collection']
		maaslin2.Maaslin2(input_data = df_r,
					input_metadata = meta_df_r,
					output = raw_output,
					normalization='TSS',
					transform='LOG',
					plot_scatter = 'FALSE',
					min_prevalence=CF,
					fixed_effects = ro.StrVector(all_vars))
	
	res = pd.read_csv(os.path.join(raw_output,'significant_results.tsv'),sep='\t',index_col=0)
	if res.shape[0]==0:
		print('No significant results detected')
	res = res[res.metadata=='case']
	res.to_csv(os.path.join(output_y,f'{level}_sig.csv'))
	#os.system(f'rm -rf {raw_output}')	
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
	if y == 'S_pneumoniae_abs':
		return
	'''
	for t in tax_lst:
		df_t = df.loc[t].to_frame()
		
		plot_df = pd.concat([meta_df, df_t], axis=1)
		#plot_df = plot_df.replace(0,np.nan)
	
		sns.violinplot(data=plot_df, x=y,y=t,palette={'Positive':"#E97B69",'Negative':"#8FA3D1"})
		sns.stripplot(data=plot_df,x=y,y=t,color='black',alpha=0.25)	
		plt.tight_layout()
		plt.savefig(os.path.join(plots_output, f'{t}.svg'))
		plt.close()
		print(plot_df)
		plot_df[t] = plot_df[t]+0.001
		print(plot_df[t].min())
		plot_df[t] = [math.log10(x) for x in plot_df[t]]
		print(plot_df)
		sns.violinplot(data=plot_df, x=y,y=t,palette={'Positive':"#E97B69",'Negative':"#8FA3D1"})
		sns.stripplot(data=plot_df,x=y,y=t,color='black',alpha=0.25)	
		plt.tight_layout()
		plt.savefig(os.path.join(plots_output, f'{t}_log10.svg'))
		plt.close()
	'''	

	return 

for y in case_toi:
	for level in ['asv']:
		for case_control in ['Positive','Negative']:
		#for case_control in ['Ambulatory']:
			main(level,y,case_control)	