import warnings
warnings.simplefilter('ignore')
import argparse
import logging
import sys
import os
import pandas as pd
import scanpy as sc

# Configure logging
logger = logging.getLogger(__name__)

def get_args():
	'''
	Get arguments from command line
	'''

	parser = argparse.ArgumentParser(
		description = "This script takes STARsolo SJ files and outputs junction read counts",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument('-i', '--experiment', type = str, help = 'Experiment table', required = True)
	parser.add_argument('-o', '--out', type = str, help = 'Output junction file', required = True)
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

	args = parser.parse_args()
	return(args)

def load_experiment_table(experiment_table):
	'''
	Load experiment table and returns a DataFrame object
	'''

	experiment_table_df = pd.read_csv(experiment_table, sep = "\t")
	return(experiment_table_df)

def make_sjpath_list(experiment_table_df):
	'''
	Generate a list of SJ file paths
	'''

	sjpath_list = experiment_table_df["SJ"].tolist()
	return(sjpath_list)

def load_sj(sj_file):
	'''
	Load SJ file and returns a DataFrame object
	'''

	adata = sc.read_mtx(sj_file + "/matrix.mtx")
	var = pd.read_csv(sj_file + "/barcodes.tsv", header = None)
	obs = pd.read_csv(sj_file + "/features.tsv", sep = "\t", header = None, usecols = [0, 1, 2])
	obs["SJ"] = obs[0].astype(str) + ":" + obs[1].astype(str) + "-" + obs[2].astype(str)
	adata.obs.index = obs["SJ"]
	adata.var.index = var[0]

	return(adata)

def make_grouppath_list(experiment_table_df):
	'''
	Generate a list of group file paths
	'''

	grouppath_list = experiment_table_df["barcode"].tolist()
	return(grouppath_list)

def load_group(grouppath):
	'''
	Load group file and returns a DataFrame object
	'''

	group_df = pd.read_csv(grouppath, sep = "\t")
	group_df = group_df.drop_duplicates()
	return(group_df)

def make_sample_group_dict(group_df):
	'''
	Generate a dictionary of sample and group
	'''

	sample_group_dict = {}
	group_list = group_df["group"].unique().tolist()
	for i in range(len(group_list)):
		sample_group_dict[group_list[i]] = group_df[group_df["group"] == group_list[i]]["barcode"].tolist()

	return(sample_group_dict)

def grouping_read_count_each(adata, sample_group_dict, group, j):
	'''
	Group junction read counts
	'''

	# Sample list
	if group not in sample_group_dict:
		count_df = pd.DataFrame()
	else:
		sample_list = list(sample_group_dict[group])
		adata = adata[:, sample_list]
		count_sum = adata.X.sum(axis = 1).flatten().tolist()[0]
		# Summarize junction read counts
		count_df = pd.DataFrame({"SJ": list(adata.obs.index), j: count_sum})
		count_df = count_df[count_df[j] != 0]

	return(count_df)

def formatting_output(count_df):
	'''
	Formatting output junction file
	'''

	output_df = count_df.astype(int)
	output_df = output_df.reset_index()
	output_df.columns = ["ID"] + list(output_df.columns[1:])
	output_df["chr"] = output_df["ID"].str.split(":").str[0]
	output_df["start"] = output_df["ID"].str.split(":").str[1].str.split("-").str[0].astype(int) - 1
	output_df["start"] = output_df["start"].astype(str)
	output_df["end"] = output_df["ID"].str.split(":").str[1].str.split("-").str[1].astype(int) + 1
	output_df["end"] = output_df["end"].astype(str)
	output_df["ID"] = output_df["chr"] + ":" + output_df["start"] + "-" + output_df["end"]
	output_df = output_df[["chr", "start", "end", "ID"] + list(output_df.columns[1:-3])]

	return(output_df)

def main():
	'''
	Main function
	'''
	args = get_args()

	# Set up logging
	logging.basicConfig(
		format="[%(asctime)s] %(levelname)7s %(message)s",
		level=logging.DEBUG if args.verbose else logging.INFO,
	)
	logger.info("Starting junction read count calculation")
	logger.debug(args)

	# Parse arguments
	experiment_table = args.experiment
	output_path = args.out

	# Make directory
	logger.info("Making output directory ...")
	os.makedirs(os.path.dirname(output_path), exist_ok = True)

	# Load experiment table
	logger.info("Loading experiment table ...")
	experiment_table_df = load_experiment_table(experiment_table)

	# Make a list of SJ file paths
	logger.info("Making SJ file paths ...")
	sjpath_list = make_sjpath_list(experiment_table_df)

	# Load SJ files
	logger.info("Loading SJ files ...")
	adata_list = []
	for x in sjpath_list:
		logger.debug(f"Loading {x}...")
		adata = load_sj(x)
		adata_list.append(adata)

	# Make a list of group file paths
	logger.info("Making group file paths ...")
	grouppath_list = make_grouppath_list(experiment_table_df)

	# Load group files
	logger.info("Loading group files ...")
	group_df_list = []
	group_list = []
	sample_group_dict_list = []
	for x in grouppath_list:
		logger.debug(f"Loading {x}...")
		group_df = load_group(x)
		group_df_list.append(group_df)
		sample_group_dict = make_sample_group_dict(group_df)
		sample_group_dict_list.append(sample_group_dict)
		group_list += list(sample_group_dict.keys())
	group_list = sorted(list(set(group_list)))

	logger.info("Grouping junction read counts ...")
	sj_grouped_df = pd.DataFrame()
	for group in group_list:
		logger.debug(f"Group: {group}")
		sj_tmp_df = pd.DataFrame()
		for j in range(len(adata_list)):
			logger.debug(f"Sample: {j}")
			tmp_df = grouping_read_count_each(adata_list[j], sample_group_dict_list[j], group, j)
			if tmp_df.empty:
				logger.warning(f"No junction read counts found for group {group} in sample {j}.")
				continue
			sj_tmp_df = pd.merge(sj_tmp_df, tmp_df, on="SJ", how="outer").fillna(0) if not sj_tmp_df.empty else tmp_df

		sj_tmp_df = sj_tmp_df.set_index("SJ")
		# Summarize junction read counts
		sj_tmp_df = pd.DataFrame({group: sj_tmp_df.sum(axis = 1)})
		if sj_grouped_df.empty:
			sj_grouped_df = sj_tmp_df
		else:
			sj_grouped_df = pd.merge(sj_grouped_df, sj_tmp_df, on = "SJ", how = "outer")
			sj_grouped_df = sj_grouped_df.fillna(0)

	# Formatting output junction file
	logger.info("Formatting output junction file ...")
	output_df = formatting_output(sj_grouped_df)

	# Write output junction file
	logger.info("Writing output junction file ...")
	output_df.to_csv(output_path, sep = "\t", index = False)

	logger.info("All processes completed.")

if __name__ == '__main__':

	main()
