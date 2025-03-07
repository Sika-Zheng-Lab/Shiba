import pysam

def is_paired_end(bam_file):
	"""
	Determine if a BAM file is paired-end.
	Returns True if paired-end, False otherwise.
	"""
	with pysam.AlignmentFile(bam_file, "rb") as bam:
		for read in bam:
			if read.is_paired:
				return True
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
