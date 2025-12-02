import argparse
import sys
import os
from lib import expression
import pandas as pd
import logging

# Configure logging
logger = logging.getLogger(__name__)

def str2bool(v):
    if v == "True":
        return True
    elif v == "False":
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Calculate TPM and CPM"
	)

	parser.add_argument("--countfiles", type = str, help = "Text file of read counts to be processed", nargs = "+")
	parser.add_argument("--onlypsi", type = str2bool, help = "Onlypsi", nargs = "?", const = True, default = False)
	parser.add_argument("--onlypsi-group", type = str2bool, help = "Onlypsi-group", nargs = "?", const = True, default = False)
	parser.add_argument("--reference-gtf", type = str, help = "Reference GTF file")
	parser.add_argument("--output", type = str, help = "Output directory")
	parser.add_argument("--excel", help = "Make result files in excel format", type = str2bool, nargs = "?", const = True, default = False)
	parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
	args = parser.parse_args()
	return(args)

def merge_table(countfiles, reference_gtf):

	count_all_df = pd.DataFrame()
	for i in range(len(countfiles)):
		sample = countfiles[i].split('/')[-1].rstrip("_counts.txt")
		count_df = pd.read_csv(
			countfiles[i],
			sep = "\t",
			skiprows = 1,
			usecols = [0, 5, 6],
			dtype = {
				0: "str",
				5: "int32",
				6: "int32"
			}
		)

		count_df.columns = count_df.columns[0:2].tolist() + [sample]
		if count_all_df.empty:
			count_all_df = count_df
		else:
			count_all_df = pd.merge(
				count_all_df,
				count_df[["Geneid", sample]],
				on = "Geneid"
			)

	count_all_df = count_all_df.rename(columns = {"Geneid": "gene_id"})
	count_all_df = count_all_df.sort_values("gene_id")
	return(count_all_df)

def main():
	## Main
	args = get_args()

	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Starting TPM and CPM calculation")
	logger.debug(args)

	# Make gene dictionary from GTF
	logger.info("Making gene ID to gene name mapping from GTF...")
	gene_dict = expression.gene_id_to_name(args.reference_gtf)
	if not gene_dict:
		logger.warning("No gene ID to gene name mapping found in GTF.")
	else:
		logger.info(f"Mapped {len(gene_dict)} gene IDs to gene names.")

	# Merge count tables
	logger.info("Merge count tables...")
	merged_table_df = merge_table(args.countfiles, args.reference_gtf)
	count_df = merged_table_df.drop(columns = ["Length"])
	try:
		count_df['gene_name'] = count_df['gene_id'].map(gene_dict)
	except Exception as e:
		logger.warning(f"Failed to map gene names: {e}")
		count_df['gene_name'] = count_df['gene_id']
	# Make sure gene_id is the first column and gene_name is the second column
	cols = count_df.columns.tolist()
	cols = [cols[0], cols[-1]] + cols[1:-1]
	count_df = count_df[cols]
	# Save count files
	if args.onlypsi == False and args.onlypsi_group == False:
		count_df.to_csv(
			os.path.join(args.output, "counts.txt"),
			sep = "\t",
			index = False
		)

	# Calculate TPM and CPM
	logger.info("Calculate TPM...")
	tpm_df = expression.ExpressionProcessor(merged_table_df.copy()).TPM()
	try:
		tpm_df['gene_name'] = tpm_df['gene_id'].map(gene_dict)
	except Exception as e:
		logger.warning(f"Failed to map gene names: {e}")
		tpm_df['gene_name'] = tpm_df['gene_id']
	# Make sure gene_id is the first column and gene_name is the second column
	cols = tpm_df.columns.tolist()
	cols = [cols[0], cols[-1]] + cols[1:-1]
	tpm_df = tpm_df[cols]
	# Save TPM
	tpm_df.to_csv(os.path.join(args.output, "TPM.txt"), sep = "\t", index = False)

	logger.info("Calculate CPM...")
	cpm_df = expression.ExpressionProcessor(merged_table_df.copy()).CPM()
	try:
		cpm_df['gene_name'] = cpm_df['gene_id'].map(gene_dict)
	except Exception as e:
		logger.warning(f"Failed to map gene names: {e}")
		cpm_df['gene_name'] = cpm_df['gene_id']
	# Make sure gene_id is the first column and gene_name is the second column
	cols = cpm_df.columns.tolist()
	cols = [cols[0], cols[-1]] + cols[1:-1]
	cpm_df = cpm_df[cols]
	# Save CPM
	cpm_df.to_csv(os.path.join(args.output, "CPM.txt"), sep = "\t", index = False)

	# Excel file
	if args.excel:
		logger.info("Exporting results to Excel...")
		# Import StyleFrame
		from styleframe import StyleFrame, Styler, utils
		# StyleFrame
		style = Styler(
			horizontal_alignment = utils.horizontal_alignments.left,
			border_type = utils.borders.default_grid,
			wrap_text = False
		)
		with StyleFrame.ExcelWriter(os.path.join(args.output, "TPM_CPM.xlsx")) as writer:
			tpm_sf = StyleFrame(tpm_df)
			tpm_sf.set_column_width(columns = tpm_df.columns, width = 20)
			tpm_sf.apply_column_style(cols_to_style = tpm_df.columns, styler_obj = style, style_header = True)
			tpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "C2", sheet_name = "TPM")
			cpm_sf = StyleFrame(cpm_df)
			cpm_sf.set_column_width(columns = cpm_df.columns, width = 20)
			cpm_sf.apply_column_style(cols_to_style = cpm_df.columns, styler_obj = style, style_header = True)
			cpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "C2", sheet_name = "CPM")

	logger.info("TPM and CPM calculation completed")

if __name__ == '__main__':
	main()
