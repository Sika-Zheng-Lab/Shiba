import argparse
import os
import sys
import subprocess
import logging
import pandas as pd
from lib import expression, general
from lib.general import str2bool

# Configure logging
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Shared helper functions
# ---------------------------------------------------------------------------

def prepare_output_dir(output_dir):
    os.makedirs(f"{output_dir}/logs", exist_ok=True)

def run_featurecounts_single(bam, annotation, output, threads, fmt="GTF", long_read=False, extra_args=None, log_file=None):
    """Run featureCounts for a single BAM file.

    Used by both the all-in-one pipeline and the 'featurecounts' subcommand.
    """
    # Check if BAM is paired-end
    paired_flag = expression.is_paired_end(bam)
    paired_option = ["-p", "-B"] if paired_flag else []
    # Check if long read mode is enabled
    longread_option = ["-L"] if long_read else []
    # Build the command
    featurecounts_command = [
        "featureCounts", "-a", annotation, "-o", output, "-T", str(threads),
    ]
    if fmt == "SAF":
        featurecounts_command += ["-F", "SAF", "--fracOverlapFeature", "1.0", "-O"]
    else:
        featurecounts_command += ["-t", "exon", "-g", "gene_id"]
    featurecounts_command += paired_option + longread_option
    if extra_args:
        featurecounts_command += extra_args
    featurecounts_command += [bam]
    # Remove empty strings
    featurecounts_command = list(filter(None, featurecounts_command))
    return_code = general.execute_command(featurecounts_command, log_file)
    if return_code != 0:
        logger.error("Error executing featureCounts. Exiting...")
        sys.exit(1)

def process_samples(experiment_file, reference_gtf, output_dir, processors):
    count_all_df = pd.DataFrame()
    with open(experiment_file, "r") as experiment:
        for line in experiment:
            line = line.strip()
            if not line or line.startswith("sample"):
                continue
            # Parse the line
            try:
                sample, bam_file, _group, technology = line.split(maxsplit=3)
            except ValueError:
                sample, bam_file, _group = line.split(maxsplit=2)
                technology = "short"
            bam_index = f"{bam_file}.bai"

            logger.info(f"Processing sample: {sample}")
            logger.debug(f"BAM file: {bam_file}")
            if technology.lower() == "long":
                logger.debug(f"{sample} will be processed as a long read sequencing experiment.")

            # Ensure BAM index exists
            if not os.path.isfile(bam_index):
                logger.error(f"BAM index file not found for sample : {sample}")
                logger.error("Please create an index file using 'samtools index' for all BAM files.")
                sys.exit(1)
            else:
                logger.debug(f"Found BAM index for {bam_file}")

            # Run featureCounts
            counts_file = f"{output_dir}/{sample}_counts.txt"
            long_read = technology.lower() == "long"
            featurecounts_log = f"{output_dir}/logs/featureCounts.log"
            run_featurecounts_single(bam_file, reference_gtf, counts_file, processors, fmt="GTF", long_read=long_read, log_file=featurecounts_log)

            # Simplify counts file
            logger.info(f"Simplifying counts file for sample {sample}")
            counts_file_df = pd.read_csv(counts_file, sep="\t", comment="#")
            counts_file_df = counts_file_df.iloc[:, [0, 5, 6]]
            counts_file_df.columns = ["GeneID", "Length", sample]
            count_all_df = pd.merge(count_all_df, counts_file_df[["GeneID", sample]], on = "GeneID") if not count_all_df.empty else counts_file_df

    count_all_df = count_all_df.rename(columns = {"GeneID": "gene_id"})
    count_all_df = count_all_df.sort_values("gene_id")
    return count_all_df

def run_deseq2(src_path, experiment_file, counts_file, refgroup, altgroup, output_dir):
    if refgroup != "NA" and altgroup != "NA":
        logger.info("Running differential expression analysis using DESeq2...")
        deseq2_command = [
            "Rscript", f"{src_path}/deseq2.R", experiment_file, counts_file, refgroup, altgroup, f"{output_dir}/DEG.txt"
        ]
        return_code = general.execute_command(deseq2_command, f"{output_dir}/logs/DESeq2.log")
        if return_code != 0:
            logger.error("DESeq2 failed")
            sys.exit(1)

def compute_tpm_cpm(merged_table_df, gene_dict, output_dir, excel=False, save_counts=True):
    """Compute TPM and CPM from a merged count table. Shared by all-in-one and tpm subcommand."""
    count_df = merged_table_df.drop(columns=["Length"])
    try:
        count_df['gene_name'] = count_df['gene_id'].map(gene_dict)
    except Exception as e:
        logger.warning(f"Failed to map gene names: {e}")
        count_df['gene_name'] = count_df['gene_id']
    # Make sure gene_id is the first column and gene_name is the second column
    cols = count_df.columns.tolist()
    cols = [cols[0], cols[-1]] + cols[1:-1]
    count_df = count_df[cols]
    if save_counts:
        count_df.to_csv(os.path.join(output_dir, "counts.txt"), sep="\t", index=False)

    # Calculate TPM
    logger.info("Calculating TPM...")
    tpm_df = expression.ExpressionProcessor(merged_table_df.copy()).TPM()
    try:
        tpm_df['gene_name'] = tpm_df['gene_id'].map(gene_dict)
    except Exception as e:
        logger.warning(f"Failed to map gene names: {e}")
        tpm_df['gene_name'] = tpm_df['gene_id']
    cols = tpm_df.columns.tolist()
    cols = [cols[0], cols[-1]] + cols[1:-1]
    tpm_df = tpm_df[cols]
    tpm_df.to_csv(os.path.join(output_dir, "TPM.txt"), sep="\t", index=False)

    # Calculate CPM
    logger.info("Calculating CPM...")
    cpm_df = expression.ExpressionProcessor(merged_table_df.copy()).CPM()
    try:
        cpm_df['gene_name'] = cpm_df['gene_id'].map(gene_dict)
    except Exception as e:
        logger.warning(f"Failed to map gene names: {e}")
        cpm_df['gene_name'] = cpm_df['gene_id']
    cols = cpm_df.columns.tolist()
    cols = [cols[0], cols[-1]] + cols[1:-1]
    cpm_df = cpm_df[cols]
    cpm_df.to_csv(os.path.join(output_dir, "CPM.txt"), sep="\t", index=False)

    # Excel file
    if excel:
        logger.info("Exporting results to Excel...")
        from styleframe import StyleFrame, Styler, utils
        style = Styler(
            horizontal_alignment=utils.horizontal_alignments.left,
            wrap_text=False
        )
        with StyleFrame.ExcelWriter(os.path.join(output_dir, "TPM_CPM.xlsx")) as writer:
            tpm_sf = StyleFrame(tpm_df)
            tpm_sf.set_column_width(columns=tpm_df.columns, width=20)
            tpm_sf.apply_column_style(cols_to_style=tpm_df.columns, styler_obj=style, style_header=True)
            tpm_sf.to_excel(writer, index=False, columns_and_rows_to_freeze="C2", sheet_name="TPM")
            cpm_sf = StyleFrame(cpm_df)
            cpm_sf.set_column_width(columns=cpm_df.columns, width=20)
            cpm_sf.apply_column_style(cols_to_style=cpm_df.columns, styler_obj=style, style_header=True)
            cpm_sf.to_excel(writer, index=False, columns_and_rows_to_freeze="C2", sheet_name="CPM")

# ---------------------------------------------------------------------------
# Subcommand: featurecounts — run featureCounts for a single BAM
# (replaces expression_featureCounts_snakemake.py)
# ---------------------------------------------------------------------------

def cmd_featurecounts(args):
    run_featurecounts_single(args.bam, args.gtf, args.output, args.threads, fmt="GTF", long_read=args.long_read)
    logger.info("Done.")

# ---------------------------------------------------------------------------
# Subcommand: tpm — merge count files and compute TPM/CPM
# (replaces tpm_snakemake.py)
# ---------------------------------------------------------------------------

def merge_count_files(countfiles):
    """Merge per-sample featureCounts output files into a single DataFrame."""
    count_all_df = pd.DataFrame()
    for i in range(len(countfiles)):
        sample = countfiles[i].split('/')[-1].rstrip("_counts.txt")
        count_df = pd.read_csv(
            countfiles[i],
            sep="\t",
            skiprows=1,
            usecols=[0, 5, 6],
            dtype={0: "str", 5: "int32", 6: "int32"}
        )
        count_df.columns = count_df.columns[0:2].tolist() + [sample]
        if count_all_df.empty:
            count_all_df = count_df
        else:
            count_all_df = pd.merge(
                count_all_df,
                count_df[["Geneid", sample]],
                on="Geneid"
            )
    count_all_df = count_all_df.rename(columns={"Geneid": "gene_id"})
    count_all_df = count_all_df.sort_values("gene_id")
    return count_all_df

def cmd_tpm(args):
    # Make gene dictionary from GTF
    logger.info("Making gene ID to gene name mapping from GTF...")
    gene_dict = expression.gene_id_to_name(args.reference_gtf)
    if not gene_dict:
        logger.warning("No gene ID to gene name mapping found in GTF.")
    else:
        logger.info(f"Mapped {len(gene_dict)} gene IDs to gene names.")

    # Merge count tables
    logger.info("Merging count tables...")
    merged_table_df = merge_count_files(args.countfiles)

    # Save counts only when not in onlypsi mode
    save_counts = not args.onlypsi and not args.onlypsi_group

    compute_tpm_cpm(merged_table_df, gene_dict, args.output, excel=args.excel, save_counts=save_counts)
    logger.info("TPM and CPM calculation completed")

# ---------------------------------------------------------------------------
# Subcommand: deseq2 — run DESeq2 for differential expression
# (replaces deseq2_snakemake.py)
# ---------------------------------------------------------------------------

def cmd_deseq2(args):
    # Load experiment table
    df = pd.read_csv(args.experiment_table, sep="\t")
    count_df = df.groupby("group").count()
    count_reference = count_df["sample"][args.reference]
    count_alternative = count_df["sample"][args.alternative]

    if count_reference >= 2 and count_alternative >= 2:
        src_dir = os.path.dirname(os.path.abspath(__file__))
        run_command = [
            "Rscript", os.path.join(src_dir, "deseq2.R"),
            args.experiment_table, args.count, args.reference, args.alternative, args.output
        ]
        return_code = general.execute_command(run_command)
        if return_code != 0:
            logger.error("Error executing the command. Exiting...")
            sys.exit(1)
    else:
        logger.error("Error: The number of samples is less than 2.")
        with open(args.output, "w") as f:
            f.write("Error: The number of samples is less than 2.\n")
    logger.info("Done.")

# ---------------------------------------------------------------------------
# Default: all-in-one pipeline (original expression.py behavior)
# ---------------------------------------------------------------------------

def cmd_all(args):
    # Prepare output directory
    prepare_output_dir(args.output)

    # Make gene dictionary from GTF
    logger.info("Making gene ID to gene name mapping from GTF...")
    gene_dict = expression.gene_id_to_name(args.reference)
    if not gene_dict:
        logger.warning("No gene ID to gene name mapping found in GTF.")
    else:
        logger.info(f"Mapped {len(gene_dict)} gene IDs to gene names.")

    # Process samples and generate count files
    count_all_df = process_samples(args.input, args.reference, args.output, args.processors)

    compute_tpm_cpm(count_all_df, gene_dict, args.output, excel=args.excel, save_counts=True)

    # Run DESeq2 for differential expression analysis
    run_deseq2(os.path.dirname(__file__), args.input, f"{args.output}/counts.txt", args.refgroup, args.altgroup, args.output)

    # Cleanup
    for count_file in os.listdir(args.output):
        if count_file.endswith("_counts.txt") or count_file.endswith("_counts.txt.summary"):
            os.remove(os.path.join(args.output, count_file))

    logger.info("RNA expression analysis completed!")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="RNA expression analysis using featureCounts and DESeq2."
    )
    subparsers = parser.add_subparsers(dest="subcommand")

    # Default (all-in-one) arguments — used when no subcommand is given
    parser.add_argument("-i", "--input", help="Experiment table")
    parser.add_argument("-g", "--reference", help="Reference GTF file")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-r", "--refgroup", default="NA", help="Reference group for differential expression analysis")
    parser.add_argument("-a", "--altgroup", default="NA", help="Alternative group for differential expression analysis")
    parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
    parser.add_argument("--excel", help="Make result files in excel format", type=str2bool, nargs="?", const=True, default=False)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")

    # Subcommand: featurecounts
    fc_parser = subparsers.add_parser("featurecounts", help="Run featureCounts for a single BAM file")
    fc_parser.add_argument("-b", "--bam", type=str, required=True, help="Input BAM file")
    fc_parser.add_argument("-g", "--gtf", type=str, required=True, help="Input GTF file")
    fc_parser.add_argument("-o", "--output", type=str, required=True, help="Output count file")
    fc_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    fc_parser.add_argument("-l", "--long-read", action="store_true", help="Long read mode")
    fc_parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")

    # Subcommand: tpm
    tpm_parser = subparsers.add_parser("tpm", help="Merge count files and calculate TPM/CPM")
    tpm_parser.add_argument("--countfiles", type=str, nargs="+", required=True, help="Count files to merge")
    tpm_parser.add_argument("--onlypsi", type=str2bool, nargs="?", const=True, default=False, help="Only PSI mode")
    tpm_parser.add_argument("--onlypsi-group", type=str2bool, nargs="?", const=True, default=False, help="Only PSI group mode")
    tpm_parser.add_argument("--reference-gtf", type=str, required=True, help="Reference GTF file")
    tpm_parser.add_argument("--output", type=str, required=True, help="Output directory")
    tpm_parser.add_argument("--excel", help="Make result files in excel format", type=str2bool, nargs="?", const=True, default=False)
    tpm_parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")

    # Subcommand: deseq2
    deseq2_parser = subparsers.add_parser("deseq2", help="Run DESeq2 for differential expression")
    deseq2_parser.add_argument("--count", type=str, required=True, help="Read count file")
    deseq2_parser.add_argument("--experiment-table", type=str, required=True, help="Experiment table")
    deseq2_parser.add_argument("--reference", type=str, required=True, help="Reference group")
    deseq2_parser.add_argument("--alternative", type=str, required=True, help="Alternative group")
    deseq2_parser.add_argument("--output", type=str, required=True, help="Output file")
    deseq2_parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")

    return parser.parse_args()

def main():
    args = parse_args()

    # Set up logging
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)7s %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    logger.debug(args)

    if args.subcommand == "featurecounts":
        logger.info("Running featureCounts...")
        cmd_featurecounts(args)
    elif args.subcommand == "tpm":
        logger.info("Starting TPM and CPM calculation")
        cmd_tpm(args)
    elif args.subcommand == "deseq2":
        logger.info("Running DESeq2...")
        cmd_deseq2(args)
    else:
        # Default: all-in-one pipeline (original behavior)
        logger.info("Starting RNA expression analysis...")
        cmd_all(args)

if __name__ == "__main__":
    main()
