import pysam
import logging
logger = logging.getLogger(__name__)

def is_paired_end(bam_file, max_reads=100000):
	"""
	Determine if a BAM file is paired-end using pysam only.
	Reads up to `max_reads` records to check for any paired read.
	"""
	paired_count = 0
	total_checked = 0
	with pysam.AlignmentFile(bam_file, "rb") as bam:
		for read in bam:
			total_checked += 1
			if read.is_paired:
				paired_count += 1
				if paired_count >= 10:
					return True
			if total_checked >= max_reads:
				break
	logger.debug(f"Checked {total_checked} reads in {bam_file}. Paired count: {paired_count}")
	return False

class ExpressionProcessor:
	def __init__(self, df):
		self.df = df

	def TPM(self):
		samples = self.df.columns[2:]
		for s in samples:
			self.df[s] = (self.df[s] / self.df["Length"]) * 1000
			self.df[s] = (self.df[s] / sum(self.df[s])) * 10**6
		self.df = self.df.drop(columns=["Length"])
		self.df = self.df.round(2)
		return self.df

	def CPM(self):
		samples = self.df.columns[2:]
		for s in samples:
			self.df[s] = (self.df[s] / sum(self.df[s])) * 10**6
		self.df = self.df.drop(columns=["Length"])
		self.df = self.df.round(2)
		return self.df

# Usage:
# processor = ExpressionProcessor(df)
# df_tpm = processor.TPM()
# df_cpm = processor.CPM()

def gene_id_to_name(gtf):
	gene_dict = {}
	try:
		with open(gtf, "r") as gtf_file:
			for line in gtf_file:
				if line.startswith("#"):
					continue
				fields = line.strip().split("\t")
				if fields[2] != "gene":
					continue
				attributes = fields[8].split("; ")
				gene_id = ""
				gene_name = ""
				for attribute in attributes:
					if attribute.startswith("gene_id"):
						gene_id = attribute.split(" ")[1].replace('"', '')
					elif attribute.startswith("gene_name"):
						gene_name = attribute.split(" ")[1].replace('"', '')
				if gene_id and gene_name:
					gene_dict[gene_id] = gene_name
		return gene_dict
	except FileNotFoundError:
		logger.error(f"GTF file not found: {gtf}")
		sys.exit(1)
