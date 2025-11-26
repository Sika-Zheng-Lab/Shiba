import argparse
import os
import sys
import pandas as pd
import logging
import scipy.stats as stats
import numpy as np
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

# Configure logging
logger = logging.getLogger(__name__)

def get_args():
    ## Get arguments from command line

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "Principal Component Analysis for matrix of gene expression and splicing"
    )
    parser.add_argument('--input-tpm', type=str, help='Input TPM file')
    parser.add_argument('--input-psi', type=str, help='Input PSI file')
    parser.add_argument('-g', '--genes', type=int, help='Number of highly-variable genes or splicing events to calculate PCs', default=3000)
    parser.add_argument('-o', '--output', type=str, help='Output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args

def load_tpm_table(tpm_file: str) -> pd.DataFrame:
    '''
    Load TPM table from input file

    Args:
    - tpm_file (str): input file containing TPM values

    Returns:
    - tpm_df (pd.DataFrame): dataframe containing TPM values
    '''

    tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)
    # Drop 'gene_name' column if exists
    if 'gene_name' in tpm_df.columns:
        tpm_df = tpm_df.drop(columns = ["gene_name"])
    # Drop rows with all zeros
    logger.debug("Dropping rows with all zeros...")
    tpm_df = tpm_df.loc[(tpm_df != 0).any(axis=1)]
    # Drop rows with NaN
    logger.debug("Dropping rows with NaN...")
    tpm_df = tpm_df.dropna()
    return tpm_df

def logit_conversion(psi: float, epsilon: float = 1e-6) -> float:
    '''
    Convert PSI value to logit scale

    Args:
    - psi (float): PSI value

    Returns:
    - logit_psi (float): logit-transformed PSI value
    '''
    
    # Avoid division by zero or log of zero
    if psi <= 0:
        psi = epsilon
    elif psi >= 1:
        psi = 1 - epsilon
    logit_psi = np.log(psi / (1 - psi))
    return logit_psi

def load_psi_table(psi_file: str) -> pd.DataFrame:
    '''
    Load PSI table from input file

    Args:
    - psi_file (str): input file containing PSI values

    Returns:
    - psi_df (pd.DataFrame): dataframe containing PSI values
    '''

    psi_df = pd.read_csv(psi_file, sep="\t", index_col=0)
    psi_df = psi_df.drop(columns = ["pos_id"])
    # KNN imputation when psi_df has less than 6000 rows without NaN values
    if psi_df.dropna().shape[0] < 6000:
        logger.warning("PSI table has less than 6000 rows without NaN values. Performing KNN imputation...")
        psi_df = psi_df.dropna(axis=0, thresh=psi_df.shape[1]*0.5) # Drop rows with more than 50% NaN
        logger.info("Number of rows after dropping rows with more than 50% NaN: {}".format(psi_df.shape[0]))
        # Check if there are any columns with all NaN
        column_all_nan = psi_df.columns[psi_df.isnull().all()]
        logger.debug("Columns with all NaN: {}".format(column_all_nan))
        logger.info("Dropping columns with all NaN...")
        psi_df = psi_df.drop(columns=column_all_nan)
        imputer = KNNImputer(n_neighbors=5)
        psi_df = pd.DataFrame(imputer.fit_transform(psi_df), index=psi_df.index, columns=psi_df.columns)
        logger.debug("Number of NaN values after KNN imputation: {}".format(psi_df.isnull().sum().sum()))
        logger.info("KNN imputation completed!")
    else:
        # Drop rows with NaN
        logger.debug("Number of NaN values: {}".format(psi_df.isnull().sum().sum()))
        logger.debug("Dropping rows with NaN...")
        psi_df = psi_df.dropna()
    # Logit conversion
    logger.info("Performing logit conversion...")
    psi_df = psi_df.applymap(logit_conversion)
    return psi_df

def mtx2pca(df, genes) -> pd.DataFrame:
    '''
    Perform PCA on the input dataframe

    Args:
    - df (pd.DataFrame): input dataframe
    - genes (int): number of highly-variable genes to calculate PCs

    Returns:
    - feature_df (pd.DataFrame): dataframe containing principal components
    - contribution_df (pd.DataFrame): dataframe containing the contribution of each principal component
    '''

    # Keep rows of top n highly-variable genes
    if df.shape[0] > genes:
        df = df.loc[df.var(axis=1).sort_values(ascending=False).index[:genes]]
    # Z-score normalization across samples
    normalized_df = df.T.apply(stats.zscore, ddof = 1)
	# PCA
    pca = PCA()
    pca.fit(normalized_df)
    # Feature
    feature = pca.transform(normalized_df)
    feature_df = pd.DataFrame(feature, columns=["PC{}".format(x + 1) for x in range(len(feature))])
    feature_df.index = df.columns
    # Contribution
    contribution_df = pd.DataFrame(pca.explained_variance_ratio_, index=["PC{}".format(x + 1) for x in range(len(feature))])
    return(feature_df, contribution_df)

def main():

    # Get arguments
    args = get_args()
    # Set up logging
    logging.basicConfig(
        format = "[%(asctime)s] %(levelname)7s %(message)s",
        level = logging.DEBUG if args.verbose else logging.INFO
    )
    logger.info("Starting PCA analysis...")
    logger.debug(args)

    # Load input files
    logger.info("Loading input files...")
    logger.debug(f"TPM file: {args.input_tpm}")
    tpm_df = load_tpm_table(args.input_tpm)
    logger.debug(f"PSI file: {args.input_psi}")
    psi_df = load_psi_table(args.input_psi)

    # Check if sample size is less than 2
    if tpm_df.shape[1] < 2:
        logger.info("PCA analysis cannot be performed because number of samples is less than 2. Exiting...")
        sys.exit(0)
    if psi_df.shape[1] < 2:
        logger.info("PCA analysis cannot be performed because number of samples is less than 2. Exiting...")
        sys.exit(0)

    # Perform PCA
    logger.info("Performing PCA...")
    logger.debug("Calculating PCA for TPM...")
    tpm_feature_df, tpm_contribution_df = mtx2pca(tpm_df, args.genes)
    logger.debug("Calculating PCA for PSI...")
    psi_feature_df, psi_contribution_df = mtx2pca(psi_df, args.genes)
    # Save output
    logger.info("Saving output...")
    logger.debug(f"Output directory: {args.output}")
    os.makedirs(args.output, exist_ok=True)
    logger.debug("Saving TPM PCA results...")
    tpm_feature_df.to_csv(os.path.join(args.output, "tpm_pca.tsv"), sep="\t")
    tpm_contribution_df.to_csv(os.path.join(args.output, "tpm_contribution.tsv"), sep="\t", header=False)
    logger.debug("Saving PSI PCA results...")
    psi_feature_df.to_csv(os.path.join(args.output, "psi_pca.tsv"), sep="\t")
    psi_contribution_df.to_csv(os.path.join(args.output, "psi_contribution.tsv"), sep="\t", header=False)

    logger.info("PCA analysis completed!")

if __name__ == '__main__':

    main()
