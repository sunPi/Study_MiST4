import pandas as pd
import os
import sys
import re

def to_bin(response):
	response_bin = []

	for r in response :
		if r == 'clinical-benefit':
			response_bin.append(1)
		else:
			response_bin.append(0)

	return response_bin

def prep_response(df):
	df['Clinical outcome'] = to_bin(df['Clinical outcome'])
	df.sort_values(by="trial_number", inplace=True)
	df.reset_index(drop=True, inplace=True)
	df.drop('response', axis=1, inplace=True)
	id21 = df[(df['trial_number'] == 'M4-001-021')].index
	df.drop(id21, inplace=True)
	df.reset_index(drop=True, inplace=True)
	df.rename(columns={'Clinical outcome':'response'}, inplace=True)
	df.rename(columns={'trial_number':'SampleID'}, inplace=True)
	df['index'] = [1,3,6,7,9,11,12,13,14,15,17,18,19,22,23,25,26,4,5,8,10,24]
	df.sort_values(by="index", inplace=True)
	df.reset_index(drop=True, inplace=True)
	df.drop("index", axis=1,inplace=True)
	return df

def prep_counts(df):
	df = df.transpose()
	df.columns = df.iloc[0,:]
	df.drop('Taxonomy', inplace=True)
	df.drop('Tax_detail', inplace=True)
	df.index.name = 'SampleID'
	df.reset_index(inplace=True)
	df.sort_values(by="SampleID", inplace=True)
	df.reset_index(drop=True, inplace=True)
	id2 = df[(df['SampleID'] == 'MIST4_002')].index
	df.drop(id2, inplace=True)
	df.reset_index(drop=True, inplace=True)
	id20 = df[(df['SampleID'] == 'MIST4_020')].index
	df.drop(id20, inplace=True)
	df.reset_index(drop=True, inplace=True)
	rid = df['SampleID'].tail(1).index
	df.iloc[rid,0] = 'MIST4_023'
	df.reset_index(drop=True, inplace=True)
	df.sort_values(by="SampleID", inplace=True)
	df.reset_index(drop=True, inplace=True)
	return df, df['SampleID']



# os.chdir("/home/jr429/Documents/UoL/Cancer_Studies_PhD/Study_MiST4/datasets/raw")
abun = "./asv_table.species.absolute.xls"
abun = sys.argv[1]
resp = "../MIST4 best response for correlative analysis.xlsx"

m = re.search("(?<=table\.).*", abun)

df_name = os.path.splitext("abun_" +m.group(0))[0]

print("Processing " + df_name + "...")

abun = pd.read_csv(abun, sep = "\t")
resp = pd.read_excel(resp)

left_df = prep_response(resp)
right_df, sid = prep_counts(abun)
left_df['SampleID'] = sid
final_df = pd.merge(left_df, right_df,on='SampleID')

save_path = "./processed/" + df_name + ".xlsx"
final_df.to_excel(str(save_path), index=False, sheet_name=str(df_name))

#
# # id_tax = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,22,24,25,26,23]
# # id_response = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,21,22,23,24,25,26]
# # fin_id = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,22,23,24,25,26]
# # col_id = set(id_1).intersection(id_2)
# # len(col_id)
# # dif = set(id_1) ^ set(id_2)
# # dif
