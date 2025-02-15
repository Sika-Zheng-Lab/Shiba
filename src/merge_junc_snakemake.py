import argparse
import sys
import pandas as pd
import logging

# Configure logging
logger = logging.getLogger(__name__)

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Merge junction read counts"
	)

	parser.add_argument("--exonexon", type = str, help = "Exon-exon Junction files", nargs = "+")
	parser.add_argument("--exonintron", type = str, help = "Exon-intron Junction files", nargs = "+")
	parser.add_argument("--output", type = str, help = "Output name")
	parser.add_argument("-v", "--verbose", action = "store_true", help = "Verbose output")

	args = parser.parse_args()
	return(args)

def merge_exonexon(filelist):

	result_df = pd.DataFrame()
	for i in range(len(filelist)):
		logger.debug(filelist[i].split("/")[-1] + "...")
		junc_df = pd.read_csv(
			filelist[i],
			sep = "\t",
			header = None,
			dtype = str
		)

		sample = filelist[i].split("/")[-1].rstrip("_exon-exon.junc")
		junc_df = junc_df.iloc[:, [0, 1, 2, 4, 10]]
		junc_df.columns = ["chr", "start", "end", "count", "block"]
		junc_df = junc_df.reset_index()
		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]
		junc_df["count"] = junc_df["count"].astype("int32")
		junc_df["blockSize1"] = junc_df["block"].str.split(",", expand = True)[0].astype("int32")
		junc_df["start"] = junc_df["start"].astype("int32") + junc_df["blockSize1"]
		junc_df["blockSize2"] = junc_df["block"].str.split(",", expand = True)[1].astype("int32")
		junc_df["end"] = junc_df["end"].astype("int32") - junc_df["blockSize2"] + 1
		junc_df["ID"] = junc_df["chr"].astype(str) + ":" + junc_df["start"].astype(str) + "-" + junc_df["end"].astype(str)
		junc_df["sample"] = sample
		junc_df = junc_df[["ID", "sample", "count"]]
		# Group by ID and sample
		junc_df = junc_df.groupby(["ID", "sample"], as_index = False).sum()
		result_df = pd.concat([result_df, junc_df], axis = 0) if not result_df.empty else junc_df

	result_df["count"] = result_df["count"].astype("int32")
	# Check if there are duplicated junctions
	if result_df.duplicated(subset=["ID", "sample"]).any():
		duplicates = result_df[result_df.duplicated(subset=["ID", "sample"], keep=False)]
		logger.debug(f"Duplicated junctions found: {duplicates}")
		logger.debug("Duplicated junctions occur when the same junction is detected in both strands.")
		logger.debug("The duplicated junctions will be summed and merged.")
		# Group by ID and sample and sum counts
		result_df = result_df.groupby(["ID", "sample"], as_index=False).sum()

	result_df = result_df.pivot(
		index = "ID",
		columns = "sample",
		values = "count"
	).fillna(0).reset_index()
	result_df = result_df.rename(columns = {"index": "ID"})

	result_df["chr"] = result_df["ID"].str.split(":", expand = True)[0]
	result_df["start"] = result_df["ID"].str.split(":", expand = True)[1].str.split("-", expand = True)[0].astype("int32")
	result_df["end"] = result_df["ID"].str.split(":", expand = True)[1].str.split("-", expand = True)[1].astype("int32")
	result_df["chr-start"] = result_df["chr"] + "-" + result_df["start"].astype(str)
	result_df["chr-end"] = result_df["chr"] + "-" + result_df["end"].astype(str)
	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID", "chr-start", "chr-end", "mean"]]
	for j in col:
		result_df = result_df.astype({j: "int32"})
	col = ["chr", "start", "end", "ID"] + col
	result_df = result_df[col]

	return(result_df)

def merge_exonintron(filelist):

	result_df = pd.DataFrame()
	for i in range(len(filelist)):
		logger.debug(filelist[i].split("/")[-1] + "...")
		junc_df = pd.read_csv(
			filelist[i],
			sep = "\t",
			comment = "#",
			dtype = str
		)

		sample = filelist[i].split("/")[-1].rstrip("_exon-intron.junc")
		junc_df = junc_df.iloc[:, [1, 2, 3, 0, 6]]
		junc_df.columns = ["chr", "start", "end", "ID", "count"]
		junc_df["sample"] = sample
		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]
		result_df = pd.concat([result_df, junc_df], axis = 0) if result_df is not None else junc_df

	result_df["count"] = result_df["count"].astype("int32")
	result_df = result_df.pivot(
		index = ["chr", "start", "end", "ID"],
		columns = "sample",
		values = "count"
	).fillna(0).reset_index()
	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID"]]
	for j in col:
		result_df = result_df.astype({j: "int32"})

	return(result_df)

def main():

	args = get_args()
	# Set up logging
	logging.basicConfig(
		format="[%(asctime)s] %(levelname)7s %(message)s",
		level=logging.DEBUG if args.verbose else logging.INFO,
	)
	logger.info("Starting merge junctions")
	logger.debug(args)

	# Parse args
	exonexon = args.exonexon
	exonintron = args.exonintron
	output = args.output

	# exon-exon junctions
	logger.info("Merge exon-exon junction count...")
	exon_exon_junc_df = merge_exonexon(exonexon)

	# exon-intron junctions
	logger.info("Merge exon-intron junction count...")
	exon_intron_junc_df = merge_exonintron(exonintron)

	# Combine and save results
	logger.info("Combine and save results...")
	result_df = pd.concat(
		[exon_exon_junc_df, exon_intron_junc_df]
	).sort_values(["chr", "start"])
	result_df.to_csv(
		args.output,
		sep = "\t",
		index = False
	)

	junc_num = str(result_df.count()[0])
	logger.debug(f"Total number of junctions: {junc_num}")
	logger.info("Merge junctions completed")

if __name__ == '__main__':

	main()
