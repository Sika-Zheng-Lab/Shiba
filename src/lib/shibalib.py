# Modules used in psi.py and scpsi.py

import warnings
# warnings.simplefilter('ignore')
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import concurrent.futures
import logging
logger = logging.getLogger(__name__)

def read_events(event_path) -> dict:
    """
    Reads alternative splicing events from text files and returns a dictionary of dataframes.

    Args:
    - event_path (str): Path to the directory that contains text files of alternative splicing events.

    Returns:
    - event_df_dict (dict): A dictionary of dataframes containing alternative splicing events.
    """

    event_types = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]
    event_df_dict = {}
    for event in event_types:
        event_df_dict[event] = pd.read_csv(
            f"{event_path}/EVENT_{event}.txt",
            sep="\t",
            dtype="str"
        )
    return event_df_dict

def read_events_sc(event_path) -> dict:
    """
    Reads alternative splicing events from text files and returns a dictionary of dataframes for single cell data.

    Args:
    - event_path (str): Path to the directory that contains text files of alternative splicing events.

    Returns:
    - event_df_dict (dict): A dictionary of dataframes containing alternative splicing events.
    """

    event_types = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]
    event_df_dict = {}
    for event in event_types:
        event_df_dict[event] = pd.read_csv(
            f"{event_path}/EVENT_{event}.txt",
            sep="\t",
            dtype="str"
        )
    return event_df_dict

def read_junctions(junction_path) -> pd.DataFrame:
    """
    Reads junctions from a file specified in the command line arguments.

    Args:
    - junction_path (str): Path to the junction file.

    Returns:
    - pd.DataFrame: A pandas DataFrame containing the junction information.
    """

    junc_df = pd.read_csv(
        junction_path,
        sep = "\t",
        dtype = "str"
    )
    # Change dtype of junction read counts
    junc_df.iloc[:, 4:] = junc_df.iloc[:, 4:].astype(int)
    return(junc_df)

def read_group(group_path) -> pd.DataFrame:
    """
    Reads group information from a file specified in the command line arguments.

    Args:
    - group_path (str): Path to the group file.

    Returns:
    - pd.DataFrame: A pandas DataFrame containing the group information.
    """

    group_df = pd.read_csv(
        group_path,
        sep = "\t",
        dtype = "str",
        usecols = ["sample", "group"]
    )
    return(group_df)

def set_group(group_df, onlypsi_group, reference, alternative) -> list:
    """
    Sets the group information based on the command line arguments.

    Args:
    - group_df (pd.DataFrame): A pandas DataFrame containing the group information.
    - onlypsi_group (bool): Whether to use only the group specified in the command line arguments.
    - reference (str): The name of the reference group.
    - alternative (str): The name of the alternative group.

    Returns:
    - list: A list containing the group information in the order [reference, alternative].
    """

    if onlypsi_group or len(group_df['group'].unique().tolist()) == 1:
        group_list = sorted(group_df['group'].unique().tolist())
    else:
        group_list = [reference, alternative]
    return(group_list)

def make_sample_list(junc_df) -> list:
    """
    Returns a list of sample names from the junction DataFrame.

    Args:
    - junc_df (pd.DataFrame): A pandas DataFrame containing junction information.

    Returns:
    - list: A list of sample names.
    """

    sample_list = list(junc_df.columns[4:])
    return(sample_list)

def sample_in_group_list(group_df, group_list) -> list:
    """
    Returns a list of sample names that belong to the groups specified in group_list.

    Args:
    - group_df (pd.DataFrame): A pandas DataFrame containing the group information.
    - group_list (list): A list containing the group information in the order [reference, alternative].

    Returns:
    - list: A list of sample names that belong to the groups specified in group_list.
    """

    group1 = group_list[0]
    group2 = group_list[1]
    sample_group1 = list(group_df[group_df['group'] == group1]["sample"])
    sample_group2 = list(group_df[group_df['group'] == group2]["sample"])
    sample_in_group_list = sample_group1 + sample_group2
    return(sample_in_group_list)

def sum_reads(onlypsi_group, junc_df, group_df, group_list) -> pd.DataFrame:
    """
    Sum reads for each junction and return a pandas DataFrame.

    Args:
    - onlypsi_group (bool): Whether to use only the group specified in the command line arguments.
    - junc_df (pd.DataFrame): A pandas DataFrame containing junction information.
    - group_df (pd.DataFrame): A pandas DataFrame containing the group information.
    - group_list (list): A list containing the group information in the order [reference, alternative].

    Returns:
    - pd.DataFrame: A pandas DataFrame containing the sum of reads for each junction.
    """

    # Melt dataframe
    junc_df = pd.melt(
        junc_df,
        id_vars = ["ID"],
        value_vars = group_df["sample"].tolist(),
        var_name = "sample",
        value_name = "reads"
    )

    # Add group information
    junc_df = pd.merge(
        junc_df,
        group_df,
        on = "sample"
    )

    # Sum reads for each group using groupby
    junc_df = junc_df.groupby(["ID", "group"], as_index = False)["reads"].sum()

    # Pivot table
    junc_df = junc_df.pivot(
        index = "ID",
        columns = "group",
        values = "reads"
    ).reset_index()

    junc_df['chr'] = "chr" + junc_df['ID'].str.split(":", expand = True)[0]
    junc_df['start'] = junc_df['ID'].str.split(":", expand = True)[1].str.split("-", expand = True)[0]
    junc_df['end'] = junc_df['ID'].str.split(":", expand = True)[1].str.split("-", expand = True)[1]
    col = ["chr", "start", "end", "ID"] + group_list
    junc_df = junc_df[col]
    return(junc_df)

def junc_dict(junc_df) -> dict:
    """
    Make dictionary for junction read counts for each sample.

    Args:
    - junc_df (pd.DataFrame): DataFrame containing junction read counts for each sample.

    Returns:
    - dict: A dictionary containing junction read counts for each sample.
    """

    junc_dict_all = {}
    sample_list = [i for i in list(junc_df.columns) if i not in ["chr", "start", "end", "ID"]]
    for i in sample_list:
        junc_dict = dict(zip(junc_df["ID"], junc_df[i]))
        junc_dict_all[i] = junc_dict
    return(junc_dict_all)

def make_junc_set(junc_df) -> set:
    """
    Make a set of all junctions.

    Args:
    - junc_df (pd.DataFrame): DataFrame containing junction read counts for each sample.

    Returns:
    - set: A set containing all junctions.
    """

    junc_set = set(junc_df["ID"])
    return(junc_set)

def event_for_analysis_se(event_df, junc_set) -> pd.DataFrame:
    """
    Select SE events for analysis based on whether they contain junctions in the junction set.

    Args:
    - event_df (pd.DataFrame): DataFrame containing SE event information.
    - junc_set (set): Set containing all junctions.

    Returns:
    - pd.DataFrame: DataFrame containing selected SE events for analysis.
    """

    event_list = []
    intron_a_values = event_df.intron_a.values
    intron_b_values = event_df.intron_b.values
    intron_c_values = event_df.intron_c.values
    for index in range(event_df.shape[0]):
        intron_a = intron_a_values[index]
        intron_b = intron_b_values[index]
        intron_c = intron_c_values[index]
        ## Select events with junctions in junction list
        if (intron_a in junc_set) or (intron_b in junc_set) or (intron_c in junc_set):
            event_list.append(index)
    event_df = event_df[event_df.index.isin(event_list)]
    return(event_df)

def event_for_analysis_mse(event_df, junc_set) -> pd.DataFrame:
    """
    Select MSE events for analysis based on whether they contain junctions in the junction set.

    Args:
    - event_df (pd.DataFrame): DataFrame containing MSE event information.
    - junc_set (set): Set containing all junctions.

    Returns:
    - pd.DataFrame: DataFrame containing selected MSE events for analysis.
    """

    event_list = []
    intron_values = event_df.intron.values
    for index in range(event_df.shape[0]):
        intron = intron_values[index]
        intron_set = set(intron.split(";"))
        ## Select events whose at least one junction is in junction list
        if len(intron_set & junc_set) > 0:
            event_list.append(index)
    event_df = event_df[event_df.index.isin(event_list)]
    return(event_df)

def event_for_analysis_five_three_afe_ale(event_df, junc_set) -> pd.DataFrame:
    """
    Select FIVE, THREE, AFE, and ALE events for analysis based on whether they contain junctions in the junction set.

    Args:
    - event_df (pd.DataFrame): DataFrame containing FIVE, THREE, AFE, and ALE event information.
    - junc_set (set): Set containing all junctions.

    Returns:
    - pd.DataFrame: DataFrame containing selected FIVE, THREE, AFE, and ALE events for analysis.
    """

    event_list = []
    intron_a_values = event_df.intron_a.values
    intron_b_values = event_df.intron_b.values
    for index in range(event_df.shape[0]):
        intron_a = intron_a_values[index]
        intron_b = intron_b_values[index]
        ## Select events with junctions in junction list
        if (intron_a in junc_set) or (intron_b in junc_set):
            event_list.append(index)
    event_df = event_df[event_df.index.isin(event_list)]
    return(event_df)

def event_for_analysis_mxe(event_df, junc_set) -> pd.DataFrame:
    """
    Select MXE events for analysis based on whether they contain junctions in the junction set.

    Args:
    - event_df (pd.DataFrame): DataFrame containing MXE event information.
    - junc_set (set): Set containing all junctions.

    Returns:
    - pd.DataFrame: DataFrame containing selected MXE events for analysis.
    """

    event_list = []
    intron_a1_values = event_df.intron_a1.values
    intron_a2_values = event_df.intron_a2.values
    intron_b1_values = event_df.intron_b1.values
    intron_b2_values = event_df.intron_b2.values
    for index in range(event_df.shape[0]):
        intron_a1 = intron_a1_values[index]
        intron_a2 = intron_a2_values[index]
        intron_b1 = intron_b1_values[index]
        intron_b2 = intron_b2_values[index]
        ## Select events with junctions in junction list
        if (intron_a1 in junc_set) or (intron_a2 in junc_set) or (intron_b1 in junc_set) or (intron_b2 in junc_set):
            event_list.append(index)
    event_df = event_df[event_df.index.isin(event_list)]
    return(event_df)

def event_for_analysis_ri(event_df, junc_set) -> pd.DataFrame:
    """
    Select RI events for analysis based on whether they contain junctions in the junction set.

    Args:
    - event_df (pd.DataFrame): DataFrame containing RI event information.
    - junc_set (set): Set containing all junctions.

    Returns:
    - pd.DataFrame: DataFrame containing selected RI events for analysis.
    """

    event_list = []
    intron_a_values = event_df.intron_a.values
    for index in range(event_df.shape[0]):
        intron_a = intron_a_values[index]
        chr = str(intron_a_values[index].split(":")[0])
        intron_a_start = int(intron_a_values[index].split(":")[1].split("-")[0])
        intron_a_end = int(intron_a_values[index].split(":")[1].split("-")[1])
        intron_a_start_junc = chr + ":" + str(intron_a_start) + "-" + str(intron_a_start + 1)
        intron_a_end_junc = chr + ":" + str(intron_a_end - 1) + "-" + str(intron_a_end)
        ## Select events with junctions in junction list
        if (intron_a in junc_set) or (intron_a_start_junc in junc_set) or (intron_a_end_junc in junc_set):
            event_list.append(index)
    event_df = event_df[event_df.index.isin(event_list)]
    return(event_df)

def col_se(sample_id, group_or_not) -> list:
    """
    Returns a list of column names for output files for SE events.

    Args:
    - sample_id (list): List of sample IDs.
    - group_or_not (bool): Whether it is a group analysis or not.

    Returns:
    - list: List of column names for output files.
    """

    col = []
    for i in sample_id:
        col += [i + "_junction_a", i + "_junction_b", i + "_junction_c", i + "_PSI"]
    col = ["event_id", "pos_id", "exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name", "label"] + col
    return(col)

def se(junc_dict_all, sample_id, event_df, num_process, minimum_reads, k) -> list:
    """
    Calculate PSI for each sample.

    Args:
    - junc_dict_all (dict): A dictionary containing junction counts for each sample.
    - sample_id (list): A list of sample IDs.
    - event_df (pd.DataFrame): A pandas DataFrame containing information about alternative splicing events.
    - num_process (int): The number of processes.
    - minimum_reads (int): The minimum number of reads for each junction.
    - k (int): The index of the current process.

    Returns:
    - list: A list of lists containing PSI values for each sample and information about alternative splicing events.
    """

    event_l = []
    # Sample ID
    sample_size = len(sample_id)
    # event list
    AS_event_l = list(event_df["event_id"])
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    pos_values = event_split_df.pos_id.values
    exon_values = event_split_df.exon.values
    intron_a_values = event_split_df.intron_a.values
    intron_b_values = event_split_df.intron_b.values
    intron_c_values = event_split_df.intron_c.values
    strand_values = event_split_df.strand.values
    gene_id_values = event_split_df.gene_id.values
    gene_name_values = event_split_df.gene_name.values
    label_values = event_split_df.label.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index], pos_values[index], exon_values[index], intron_a_values[index], intron_b_values[index], intron_c_values[index], strand_values[index], gene_id_values[index], gene_name_values[index], label_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_b_count = junc_list[intron_b_values[index]]
            except:
                intron_b_count = 0
            try:
                intron_c_count = junc_list[intron_c_values[index]]
            except:
                intron_c_count = 0
            # Check minimum junction read count
            if (intron_a_count + intron_b_count >= minimum_reads*2) or (intron_c_count >= minimum_reads):
                # PSI
                if ((intron_a_count + intron_b_count) / 2 + intron_c_count) != 0:
                    psi = ((intron_a_count + intron_b_count) / 2) / ((intron_a_count + intron_b_count) / 2 + intron_c_count)
                else:
                    psi = np.nan
            else:
                psi = np.nan
            psi_list += [intron_a_count, intron_b_count, intron_c_count, psi]
        event_l += [psi_list]
    return(event_l)

def col_ind(sample_list) -> list:
    """
    Returns a list of column names for output files.

    Args:
    - sample_list: a list of sample IDs

    Returns:
    - col: a list of column names for output files, including "event_id" and each sample ID followed by "_PSI"
    """

    col = ["event_id"] + [i + "_PSI" for i in sample_list]
    return(col)

def se_ind(junc_dict_all, event_df, sample_id, num_process, k) -> list:
    """
    Calculates PSI for each sample.

    Args:
    - junc_dict_all: a dictionary containing junction counts for each sample
    - event_df: a pandas DataFrame containing information about each event
    - sample_id: a list of sample IDs
    - num_process (int): The number of processes.
    - k: an integer representing the index of the current process

    Returns:
    - event_l: a list of lists containing the event ID and PSI values for each sample
    """

    event_l = []
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    intron_a_values = event_split_df.intron_a.values
    intron_b_values = event_split_df.intron_b.values
    intron_c_values = event_split_df.intron_c.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_b_count = junc_list[intron_b_values[index]]
            except:
                intron_b_count = 0
            try:
                intron_c_count = junc_list[intron_c_values[index]]
            except:
                intron_c_count = 0
            # PSI
            if ((intron_a_count + intron_b_count) / 2 + intron_c_count) != 0:
                psi = ((intron_a_count + intron_b_count) / 2) / ((intron_a_count + intron_b_count) / 2 + intron_c_count)
            else:
                psi = None
            psi_list += [psi]
        event_l += [psi_list]
    return(event_l)

def col_mse(sample_id, group_or_not) -> list:
    """
    Returns a list of column names for output files for MSE events.

    Args:
    - sample_id (list): A list of sample IDs.
    - group_or_not (bool): Whether it is a group analysis or not.

    Returns:
    - col (list): A list of column names for output files.
    """

    col = []
    for i in sample_id:
        col += [i + "_junction", i + "_PSI"]
    col = ["event_id", "pos_id", "mse_n", "exon", "intron", "strand", "gene_id", "gene_name", "label"] + col
    return(col)

def mse(junc_dict_all, sample_id, event_df, num_process, minimum_reads, k) -> list:
    """
    Calculate PSI for each sample in the MSE event.

    Args:
    - junc_dict_all (dict): A dictionary containing junction counts for each sample.
    - sample_id (list): A list of sample IDs.
    - event_df (pd.DataFrame): A pandas DataFrame containing information about the MSE event.
    - num_process (int): The number of processes.
    - minimum_reads (int): The minimum number of reads for each junction.
    - k (int): An integer representing the index of the current process.

    Returns:
    - event_l (list): A list of lists containing PSI values for each sample in the MSE event.
    """

    event_l = []
    # Sample ID
    sample_size = len(sample_id)
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    pos_values = event_split_df.pos_id.values
    mse_n_values = event_split_df.mse_n.values
    exon_values = event_split_df.exon.values
    intron_values = event_split_df.intron.values
    strand_values = event_split_df.strand.values
    gene_id_values = event_split_df.gene_id.values
    gene_name_values = event_split_df.gene_name.values
    label_values = event_split_df.label.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index], pos_values[index], mse_n_values[index], exon_values[index], intron_values[index], strand_values[index], gene_id_values[index], gene_name_values[index], label_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            intron_list = intron_values[index].split(";")
            intron_count_list = []
            for j in range(len(intron_list)):
                try:
                    intron_count = junc_list[intron_list[j]]
                except:
                    intron_count = 0
                intron_count_list.append(intron_count)
            inclusion_intron_count_list = intron_count_list[0:-1]
            exclusion_intron_count = intron_count_list[-1]
            # Check minimum junction read count
            if (sum(inclusion_intron_count_list) >= minimum_reads*len(inclusion_intron_count_list)) or (exclusion_intron_count >= minimum_reads):
                # PSI
                if (sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) + exclusion_intron_count != 0:
                    psi = (sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) / ((sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) + exclusion_intron_count)
                else:
                    psi = np.nan
            else:
                psi = np.nan
            intron_count_concat = ";".join([str(x) for x in intron_count_list])
            psi_list += [intron_count_concat, psi]
        event_l += [psi_list]
    return(event_l)

def mse_ind(junc_dict_all, event_df, sample_id, num_process, k) -> list:
    """
    Calculates PSI for each sample in the MSE event.

    Args:
    - junc_dict_all (dict): A dictionary containing junction counts for each sample.
    - event_df (pd.DataFrame): A pandas DataFrame containing information about the MSE event.
    - sample_id (list): A list of sample IDs.
    - num_process (int): The number of processes.
    - k (int): An integer representing the index of the current process.

    Returns:
    - list: A list of lists containing the event ID and PSI values for each sample.
    """

    event_l = []
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    mse_n_values = event_split_df.mse_n.values
    exon_values = event_split_df.exon.values
    intron_values = event_split_df.intron.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            intron_list = intron_values[index].split(";")
            intron_count_list = []
            for j in range(len(intron_list)):
                try:
                    intron_count = junc_list[intron_list[j]]
                except:
                    intron_count = 0
                intron_count_list.append(intron_count)
            inclusion_intron_count_list = intron_count_list[0:-1]
            exclusion_intron_count = intron_count_list[-1]
            # PSI
            if (sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) + exclusion_intron_count != 0:
                psi = (sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) / ((sum(inclusion_intron_count_list) / len(inclusion_intron_count_list)) + exclusion_intron_count)
            else:
                psi = np.nan
            psi_list += [psi]
        event_l += [psi_list]
    return(event_l)

def col_five_three_afe_ale(sample_id, group_or_not) -> list:
    """
    Returns a list of column names for output files.

    Args:
    - sample_id: a list of sample IDs
    - group_or_not (bool): Whether it is a group analysis or not.

    Returns:
    - col: a list of column names for output files
    """

    col = []
    for i in sample_id:
        col += [i + "_junction_a", i + "_junction_b", i + "_PSI"]
    col = ["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name", "label"] + col
    return(col)

def five_three_afe_ale(junc_dict_all, sample_id, event_df, num_process, minimum_reads, k) -> list:
    """
    Calculates PSI for each sample.

    Args:
    - junc_dict_all: a dictionary containing junction reads for each sample
    - sample_id: a list of sample IDs
    - event_df: a pandas DataFrame containing AS event information
    - num_process (int): The number of processes.
    - minimum_reads (int): The minimum number of reads for each junction.
    - k: an integer representing the index of the current process

    Returns:
    - event_l: a list of lists containing PSI values for each sample
    """

    event_l = []
    # Sample ID
    sample_size = len(sample_id)
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    pos_values = event_split_df.pos_id.values
    exon_a_values = event_split_df.exon_a.values
    exon_b_values = event_split_df.exon_b.values
    intron_a_values = event_split_df.intron_a.values
    intron_b_values = event_split_df.intron_b.values
    strand_values = event_split_df.strand.values
    gene_id_values = event_split_df.gene_id.values
    gene_name_values = event_split_df.gene_name.values
    label_values = event_split_df.label.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index], pos_values[index], exon_a_values[index], exon_b_values[index], intron_a_values[index], intron_b_values[index], strand_values[index], gene_id_values[index], gene_name_values[index], label_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_b_count = junc_list[intron_b_values[index]]
            except:
                intron_b_count = 0
            # Check minimum junction read count
            if (intron_a_count >= minimum_reads) or (intron_b_count >= minimum_reads):
                # PSI
                if intron_a_count + intron_b_count != 0:
                    psi = intron_a_count / (intron_a_count + intron_b_count)
                else:
                    psi = np.nan
            else:
                psi = np.nan
            psi_list += [intron_a_count, intron_b_count, psi]
        event_l += [psi_list]
    return(event_l)

def five_three_afe_ale_ind(junc_dict_all, event_df, sample_id, num_process, k) -> list:
    """
    Calculates PSI for each sample in a given event DataFrame for a specific sample ID.

    Args:
    - junc_dict_all (dict): A dictionary containing junction counts for each sample.
    - event_df (pd.DataFrame): A DataFrame containing event information.
    - sample_id (list): A list of sample IDs to calculate PSI for.
    - num_process (int): The number of processes.
    - k (int): The index of the current process.

    Returns:
    - list: A list of lists containing the event ID and PSI values for each sample.
    """

    event_l = []
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    intron_a_values = event_split_df.intron_a.values
    intron_b_values = event_split_df.intron_b.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_b_count = junc_list[intron_b_values[index]]
            except:
                intron_b_count = 0
            # PSI
            if intron_a_count + intron_b_count != 0:
                psi = intron_a_count / (intron_a_count + intron_b_count)
            else:
                psi = None
            psi_list += [psi]
        event_l += [psi_list]
    return(event_l)

def col_mxe(sample_id, group_or_not) -> list:
    """
    Returns a list of column names for mxe events.

    Args:
    - sample_id (list): A list of sample IDs.
    - group_or_not (bool): Whether it is a group analysis or not.

    Returns:
    - col (list): A list of column names for mxe events.
    """

    col = []
    for i in sample_id:
        col += [i + "_junction_a1", i + "_junction_a2", i + "_junction_b1", i + "_junction_b2", i + "_PSI"]
    col = ["event_id", "pos_id", "exon_a", "exon_b", "intron_a1", "intron_a2", "intron_b1", "intron_b2", "strand", "gene_id", "gene_name", "label"] + col
    return(col)

def mxe(junc_dict_all, sample_id, event_df, num_process, minimum_reads, k) -> list:
    """
    Calculates PSI for each sample in the mxe event.

    Args:
    - junc_dict_all (dict): A dictionary containing junction reads for each sample.
    - sample_id (list): A list of sample IDs.
    - event_df (pandas.DataFrame): A pandas DataFrame containing information about the mxe event.
    - num_process (int): The number of processes.
    - minimum_reads (int): The minimum number of reads for each junction.
    - k (int): An integer representing the index of the current process.

    Returns:
    - event_l (list): A list of lists containing PSI values for each sample in the mxe event.
    """

    event_l = []
    # Sample ID
    sample_size = len(sample_id)
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    pos_values = event_split_df.pos_id.values
    exon_a_values = event_split_df.exon_a.values
    exon_b_values = event_split_df.exon_b.values
    intron_a1_values = event_split_df.intron_a1.values
    intron_a2_values = event_split_df.intron_a2.values
    intron_b1_values = event_split_df.intron_b1.values
    intron_b2_values = event_split_df.intron_b2.values
    strand_values = event_split_df.strand.values
    gene_id_values = event_split_df.gene_id.values
    gene_name_values = event_split_df.gene_name.values
    label_values = event_split_df.label.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index], pos_values[index], exon_a_values[index], exon_b_values[index], intron_a1_values[index], intron_a2_values[index], intron_b1_values[index], intron_b2_values[index], strand_values[index], gene_id_values[index], gene_name_values[index], label_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a1_count = junc_list[intron_a1_values[index]]
            except:
                intron_a1_count = 0
            try:
                intron_a2_count = junc_list[intron_a2_values[index]]
            except:
                intron_a2_count = 0
            try:
                intron_b1_count = junc_list[intron_b1_values[index]]
            except:
                intron_b1_count = 0
            try:
                intron_b2_count = junc_list[intron_b2_values[index]]
            except:
                intron_b2_count = 0
            # Check minimum junction read count
            if (intron_a1_count + intron_a2_count >= minimum_reads*2) or (intron_b1_count + intron_b2_count >= minimum_reads*2):
                # PSI
                if intron_a1_count + intron_a2_count + intron_b1_count + intron_b2_count != 0:
                    psi = (intron_a1_count + intron_a2_count) / ((intron_a1_count + intron_a2_count) + (intron_b1_count + intron_b2_count))
                else:
                    psi = np.nan
            else:
                psi = np.nan
            psi_list += [intron_a1_count, intron_a2_count, intron_b1_count, intron_b2_count, psi]
        event_l += [psi_list]
    return(event_l)

def mxe_ind(junc_dict_all, event_df, sample_id, num_process, k) -> list:
    """
    Calculate PSI of MXE events for each sample.

    Parameters:
    - junc_dict_all (dict): A dictionary containing junction counts for each sample.
    - event_df (pd.DataFrame): A DataFrame containing information about each MXE event.
    - sample_id (list): A list of sample IDs.
    - num_process (int): The number of processes.
    - k (int): The index of the current process.

    Returns:
    - list: A list of lists, where each inner list contains the event ID and PSI values for each sample.
    """

    event_l = []
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    intron_a1_values = event_split_df.intron_a1.values
    intron_a2_values = event_split_df.intron_a2.values
    intron_b1_values = event_split_df.intron_b1.values
    intron_b2_values = event_split_df.intron_b2.values
    for index in range(event_split_df.shape[0]):
        psi_list = [event_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a1_count = junc_list[intron_a1_values[index]]
            except:
                intron_a1_count = 0
            try:
                intron_a2_count = junc_list[intron_a2_values[index]]
            except:
                intron_a2_count = 0
            try:
                intron_b1_count = junc_list[intron_b1_values[index]]
            except:
                intron_b1_count = 0
            try:
                intron_b2_count = junc_list[intron_b2_values[index]]
            except:
                intron_b2_count = 0
            # PSI
            if intron_a1_count + intron_a2_count + intron_b1_count + intron_b2_count != 0:
                psi = (intron_a1_count + intron_a2_count) / ((intron_a1_count + intron_a2_count) + (intron_b1_count + intron_b2_count))
            else:
                psi = None
            psi_list += [psi]
        event_l += [psi_list]
    return(event_l)

def col_ri(sample_id, group_or_not) -> list:
    """
    Returns a list of column names for RI events.

    Args:
    - sample_id (list): A list of sample IDs.
    - group_or_not (bool): Whether it is a group analysis or not.

    Returns:
    - col (list): A list of column names for RI events.
    """

    col = []
    for i in sample_id:
        col += [i + "_junction_a_start", i + "_junction_a_end", i + "_junction_a", i + "_PSI"]
    col = ["event_id", "pos_id", "exon_a", "exon_b", "exon_c", "intron_a", "strand", "gene_id", "gene_name", "label"] + col
    return(col)

def ri(junc_dict_all, sample_id, event_df, num_process, minimum_reads, k) -> list:
    """
    Calculates PSI for each sample.

    Args:
    - junc_dict_all (dict): A dictionary containing junction reads for each sample.
    - sample_id (list): A list of sample IDs.
    - event_df (pd.DataFrame): A pandas DataFrame containing information about the events.
    - num_process (int): The number of processes.
    - minimum_reads (int): The minimum number of reads for each junction.
    - k (int): An integer representing the index of the current process.

    Returns:
    - event_l (list): A list of lists containing PSI values for each event.
    """

    event_l = []
    sample_size = len(sample_id)
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    pos_values = event_split_df.pos_id.values
    exon_a_values = event_split_df.exon_a.values
    exon_b_values = event_split_df.exon_b.values
    exon_c_values = event_split_df.exon_c.values
    intron_a_values = event_split_df.intron_a.values
    strand_values = event_split_df.strand.values
    gene_id_values = event_split_df.gene_id.values
    gene_name_values = event_split_df.gene_name.values
    label_values = event_split_df.label.values
    for index in range(event_split_df.shape[0]):
        chr = str(intron_a_values[index].split(":")[0])
        intron_a_start = int(intron_a_values[index].split(":")[1].split("-")[0])
        intron_a_end = int(intron_a_values[index].split(":")[1].split("-")[1])
        intron_a_start_junc = chr + ":" + str(intron_a_start) + "-" + str(intron_a_start + 1)
        intron_a_end_junc = chr + ":" + str(intron_a_end - 1) + "-" + str(intron_a_end)
        psi_list = [event_values[index], pos_values[index], exon_a_values[index], exon_b_values[index], exon_c_values[index], intron_a_values[index], strand_values[index], gene_id_values[index], gene_name_values[index], label_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_a_start_junc_count = junc_list[intron_a_start_junc]
            except:
                intron_a_start_junc_count = 0
            try:
                intron_a_end_junc_count = junc_list[intron_a_end_junc]
            except:
                intron_a_end_junc_count = 0
            # Check minimum junction read count
            if (intron_a_start_junc_count + intron_a_end_junc_count >= minimum_reads*2) or (intron_a_count >= minimum_reads):
                # PSI
                if ((intron_a_start_junc_count + intron_a_end_junc_count) / 2 + intron_a_count) != 0:
                    psi = ((intron_a_start_junc_count + intron_a_end_junc_count) / 2) / ((intron_a_start_junc_count + intron_a_end_junc_count) / 2 + intron_a_count)
                else:
                    psi = np.nan
            else:
                psi = np.nan
            psi_list += [intron_a_start_junc_count, intron_a_end_junc_count, intron_a_count, psi]
        event_l += [psi_list]
    return(event_l)

def ri_ind(junc_dict_all, event_df, sample_id, num_process, k) -> list:
    """
    Calculate PSI for each sample for RI events.

    Args:
    - junc_dict_all (dict): A dictionary containing junction reads for each sample.
    - event_df (pd.DataFrame): A pandas DataFrame containing information about the events.
    - sample_id (list): A list of sample IDs.
    - num_process (int): The number of processes.
    - k (int): An integer representing the index of the current process.

    Returns:
    - list: A list of lists containing the PSI values for each sample for each event.
    """

    event_l = []
    # event list
    AS_event_l = list(set(event_df["event_id"]))
    # split
    event_num_split = np.array_split(AS_event_l, num_process)
    event_split_df = event_df[event_df["event_id"].isin(event_num_split[k])]
    event_split_df = event_split_df.reset_index()
    event_values = event_split_df.event_id.values
    intron_a_values = event_split_df.intron_a.values
    for index in range(event_split_df.shape[0]):
        chr = str(intron_a_values[index].split(":")[0])
        intron_a_start = int(intron_a_values[index].split(":")[1].split("-")[0])
        intron_a_end = int(intron_a_values[index].split(":")[1].split("-")[1])
        intron_a_start_junc = chr + ":" + str(intron_a_start) + "-" + str(intron_a_start + 1)
        intron_a_end_junc = chr + ":" + str(intron_a_end - 1) + "-" + str(intron_a_end)
        psi_list = [event_values[index]]
        for i in sample_id:
            junc_list = junc_dict_all[i]
            try:
                intron_a_count = junc_list[intron_a_values[index]]
            except:
                intron_a_count = 0
            try:
                intron_a_start_junc_count = junc_list[intron_a_start_junc]
            except:
                intron_a_start_junc_count = 0
            try:
                intron_a_end_junc_count = junc_list[intron_a_end_junc]
            except:
                intron_a_end_junc_count = 0
            # PSI
            if ((intron_a_start_junc_count + intron_a_end_junc_count) / 2 + intron_a_count) != 0:
                psi = ((intron_a_start_junc_count + intron_a_end_junc_count) / 2) / ((intron_a_start_junc_count + intron_a_end_junc_count) / 2 + intron_a_count)
            else:
                psi = None
            psi_list += [psi]
        event_l += [psi_list]
    return(event_l)

def diff_se(df, group_list, FDR, dPSI) -> pd.DataFrame:
    """
    Differential splicing analysis for SE.

    Args:
    - df: pandas DataFrame containing junction counts and PSI values for each sample.
    - group_list: list of two strings representing the two groups to compare.
    - FDR (float): False discovery rate.
    - dPSI (float): Minimum delta PSI.

    Returns:
    - result_df: pandas DataFrame containing differential splicing analysis results for SE events.
    """

    result_df = df
    # Drop rows with nan values
    result_df = result_df.dropna()
    result_df = result_df.reset_index()
    group1 = group_list[0]
    group2 = group_list[1]
    group1_junction_a = group1 + "_junction_a"
    group1_junction_b = group1 + "_junction_b"
    group1_junction_c = group1 + "_junction_c"
    group1_PSI = group1 + "_PSI"
    group2_junction_a = group2 + "_junction_a"
    group2_junction_b = group2 + "_junction_b"
    group2_junction_c = group2 + "_junction_c"
    group2_PSI = group2 + "_PSI"
    # Keep columns of junctions and PSI of the two groups
    result_df = result_df[
        ["event_id", "pos_id", "exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name", "label"] + \
        [group1_junction_a, group1_junction_b, group1_junction_c, group1_PSI, group2_junction_a, group2_junction_b, group2_junction_c, group2_PSI]
    ]
    if result_df.shape[0] != 0:
        # delta PSI
        result_df.loc[:, 'dPSI'] = result_df.loc[:, group2_PSI] - result_df.loc[:, group1_PSI]
        # Fisher's exact test, Odds ratio
        result_df = result_df.reset_index()
        result_df = result_df.drop(columns = "index")
        group1_junction_a_values = result_df[group1_junction_a].values
        group1_junction_b_values = result_df[group1_junction_b].values
        group1_junction_c_values = result_df[group1_junction_c].values
        group2_junction_a_values = result_df[group2_junction_a].values
        group2_junction_b_values = result_df[group2_junction_b].values
        group2_junction_c_values = result_df[group2_junction_c].values
        oddsr_junction_a_col = []
        p_junction_a_col = []
        oddsr_junction_b_col = []
        p_junction_b_col = []
        p_maximum_col = []
        for index in range(result_df.shape[0]):
            # Make 2x2 table
            # inc1 - exc
            table_2x2_junction_a = [
                [group1_junction_a_values[index], group1_junction_c_values[index]],
                [group2_junction_a_values[index], group2_junction_c_values[index]]
            ]
            # inc2 - exc
            table_2x2_junction_b = [
                [group1_junction_b_values[index], group1_junction_c_values[index]],
                [group2_junction_b_values[index], group2_junction_c_values[index]]
            ]
            # Fisher's exact test
            oddsr_junction_a, p_junction_a = stats.fisher_exact(table_2x2_junction_a, alternative = 'two-sided')
            oddsr_junction_b, p_junction_b = stats.fisher_exact(table_2x2_junction_b, alternative = 'two-sided')
            p_maximum = max([p_junction_a, p_junction_b])
            oddsr_junction_a_col.append(oddsr_junction_a)
            p_junction_a_col.append(p_junction_a)
            oddsr_junction_b_col.append(oddsr_junction_b)
            p_junction_b_col.append(p_junction_b)
            p_maximum_col.append(p_maximum)
        result_df["OR_junction_a"] = oddsr_junction_a_col
        result_df["p_junction_a"] = p_junction_a_col
        result_df["OR_junction_b"] = oddsr_junction_b_col
        result_df["p_junction_b"] = p_junction_b_col
        result_df["p_maximum"] = p_maximum_col
        # FDR correction
        result_df.loc[:, "q"] = multitest.multipletests(result_df.loc[:, "p_maximum"], method = "fdr_bh")[1]
        result_df.loc[
            (((result_df.loc[:, "OR_junction_a"] >= 3/2) & (result_df.loc[:, "OR_junction_b"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a"] <= 2/3) & (result_df.loc[:, "OR_junction_b"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI),
            "Diff events"
        ] = "Yes"
        result_df.loc[
            ~((((result_df.loc[:, "OR_junction_a"] >= 3/2) & (result_df.loc[:, "OR_junction_b"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a"] <= 2/3) & (result_df.loc[:, "OR_junction_b"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI)),
            "Diff events"
        ] = "No"
        result_df.loc[result_df["q"] == 0, "q"] = 1e-323
    # Rename columns
    result_df = result_df.rename(columns = {
        group1 + "_junction_a": "ref_junction_a",
        group1 + "_junction_b": "ref_junction_b",
        group1 + "_junction_c": "ref_junction_c",
        group1 + "_PSI": "ref_PSI",
        group2 + "_junction_a": "alt_junction_a",
        group2 + "_junction_b": "alt_junction_b",
        group2 + "_junction_c": "alt_junction_c",
        group2 + "_PSI": "alt_PSI"
    })
    return(result_df)

def diff_mse(df, group_list, FDR, dPSI) -> pd.DataFrame:
    """
    Differential splicing analysis for MSE events.

    Args:
    - df (pd.DataFrame): DataFrame containing MSE junction and PSI information.
    - group_list (list): List of two group names to compare
    - FDR (float): False discovery rate
    - dPSI (float): Minimum delta PSI

    Returns:
    - pd.DataFrame: DataFrame containing differential splicing analysis results for MSE events.
    """

    result_df = df
    # Drop rows with nan values
    result_df = result_df.dropna()
    result_df = result_df.reset_index()
    group1 = group_list[0]
    group2 = group_list[1]
    group1_junction = group1 + "_junction"
    group1_PSI = group1 + "_PSI"
    group2_junction = group2 + "_junction"
    group2_PSI = group2 + "_PSI"
    # Keep columns of junctions and PSI of the two groups
    result_df = result_df[
        ["event_id", "pos_id", "mse_n", "exon", "intron", "strand", "gene_id", "gene_name", "label"] + \
        [group1_junction, group1_PSI, group2_junction, group2_PSI]
    ]
    if result_df.shape[0] != 0:
        # delta PSI
        result_df["dPSI"] = result_df[group2_PSI] - result_df[group1_PSI]
        # Fisher's exact test, Odds ratio
        result_df = result_df.reset_index()
        result_df = result_df.drop(columns = "index")
        group1_junction_values = result_df[group1_junction].values
        group2_junction_values = result_df[group2_junction].values
        oddsr_junction_col = []
        oddsr_diff_up_col = []
        oddsr_diff_down_col = []
        p_junction_col = []
        p_maximum_col = []
        for index in range(result_df.shape[0]):
            oddsr_junction_list = []
            oddsr_diff_up_list = []
            oddsr_diff_down_list = []
            p_junction_list = []
            group1_junction_count_list = [int(x) for x in group1_junction_values[index].split(";")]
            group1_inclusion_junction_count_list = group1_junction_count_list[0:-1]
            group1_exclusion_junction_count = group1_junction_count_list[-1]
            group2_junction_count_list = [int(x) for x in group2_junction_values[index].split(";")]
            group2_inclusion_junction_count_list = group2_junction_count_list[0:-1]
            group2_exclusion_junction_count = group2_junction_count_list[-1]
            for j in range(len(group1_inclusion_junction_count_list)):
                group1_inclusion_junction_count = group1_inclusion_junction_count_list[j]
                group2_inclusion_junction_count = group2_inclusion_junction_count_list[j]
                # Make 2x2 table
                table_2x2_junction = [
                    [group1_inclusion_junction_count, group1_exclusion_junction_count],
                    [group2_inclusion_junction_count, group2_exclusion_junction_count]
                ]
                # Fisher's exact test
                oddsr, p = stats.fisher_exact(table_2x2_junction, alternative = 'two-sided')
                oddsr_diff_up_list.append(oddsr >= 3/2)
                oddsr_diff_down_list.append(oddsr <= 2/3)
                oddsr_junction_list.append(oddsr)
                p_junction_list.append(p)
            oddsr_diff_up = all(oddsr_diff_up_list)
            oddsr_diff_down = all(oddsr_diff_down_list)
            oddsr_junction_concat = ";".join([str(x) for x in oddsr_junction_list])
            p_junction_concat = ";".join([str(x) for x in p_junction_list])
            oddsr_junction_col.append(oddsr_junction_concat)
            oddsr_diff_up_col.append(oddsr_diff_up)
            oddsr_diff_down_col.append(oddsr_diff_down)
            p_junction_col.append(p_junction_concat)
            p_maximum_col.append(max(p_junction_list))
        result_df["OR_junction"] = oddsr_junction_col
        result_df["OR_diff_up"] = oddsr_diff_up_col
        result_df["OR_diff_down"] = oddsr_diff_down_col
        result_df["p_junction"] = p_junction_col
        result_df["p_maximum"] = p_maximum_col
        # FDR correction
        result_df.loc[:, "q"] = multitest.multipletests(result_df.loc[:, "p_maximum"], method = "fdr_bh")[1]
        result_df["Diff events"] = result_df.apply(
            lambda x: "Yes" if (x["OR_diff_up"] or x["OR_diff_down"]) and (x["q"] < FDR) and (abs(x["dPSI"]) >= dPSI) else "No",
            axis = 1
        )
        result_df["q"] = result_df["q"].apply(lambda x: 1e-323 if x == 0 else x)
        # Drop columns
        result_df = result_df.drop(columns = ["OR_diff_up", "OR_diff_down"])
    # Rename columns
    result_df = result_df.rename(columns = {
        group1 + "_junction": "ref_junction",
        group1 + "_PSI": "ref_PSI",
        group2 + "_junction": "alt_junction",
        group2 + "_PSI": "alt_PSI"
    })
    return(result_df)

def diff_five_three_afe_ale(df, group_list, FDR, dPSI) -> pd.DataFrame:
    """
    Differential splicing analysis for FIVE, THREE, AFE, and ALE events.

    Args:
    - df: pandas DataFrame containing the PSI values for each sample and each event
    - group_list: list of two strings containing the names of the two groups to compare
    - FDR (float): False discovery rate
    - dPSI (float): Minimum delta PSI

    Returns:
    - result_df: pandas DataFrame containing the differential splicing analysis results for FIVE, THREE, AFE, and ALE events
    """

    result_df = df
    # Drop rows with nan values
    result_df = result_df.dropna()
    result_df = result_df.reset_index()
    group1 = group_list[0]
    group2 = group_list[1]
    group1_junction_a = group1 + "_junction_a"
    group1_junction_b = group1 + "_junction_b"
    group1_PSI = group1 + "_PSI"
    group2_junction_a = group2 + "_junction_a"
    group2_junction_b = group2 + "_junction_b"
    group2_PSI = group2 + "_PSI"
    # Keep columns of junctions and PSI of the two groups
    result_df = result_df[
        ["event_id", "pos_id", "exon_a", "exon_b", "intron_a", "intron_b", "strand", "gene_id", "gene_name", "label"] + \
        [group1_junction_a, group1_junction_b, group1_PSI, group2_junction_a, group2_junction_b, group2_PSI]
    ]
    if result_df.shape[0] != 0:
        # delta PSI
        result_df.loc[:, 'dPSI'] = result_df.loc[:, group2_PSI] - result_df.loc[:, group1_PSI]
        # Fisher's exact test, Odds ratio
        result_df = result_df.reset_index()
        result_df = result_df.drop(columns = "index")
        group1_junction_a_values = result_df[group1_junction_a].values
        group1_junction_b_values = result_df[group1_junction_b].values
        group2_junction_a_values = result_df[group2_junction_a].values
        group2_junction_b_values = result_df[group2_junction_b].values
        oddsr_col = []
        p_col = []
        for index in range(result_df.shape[0]):
            # Make 2x2 table
            table_2x2_intron = [
                [group1_junction_a_values[index], group1_junction_b_values[index]],
                [group2_junction_a_values[index], group2_junction_b_values[index]]
            ]
            # Fisher's exact test
            oddsr, p = stats.fisher_exact(table_2x2_intron, alternative = 'two-sided')
            oddsr_col.append(oddsr)
            p_col.append(p)
        result_df["OR"] = oddsr_col
        result_df["p"] = p_col
        # FDR correction
        result_df.loc[:, "q"] = multitest.multipletests(result_df.loc[:, "p"], method = "fdr_bh")[1]
        result_df.loc[
            ((result_df.loc[:, "OR"] >= 3/2) |
            (result_df.loc[:, "OR"] <= 2/3)) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI),
            "Diff events"
        ] = "Yes"
        result_df.loc[
            ~(((result_df.loc[:, "OR"] >= 3/2) |
            (result_df.loc[:, "OR"] <= 2/3)) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI)),
            "Diff events"
        ] = "No"
        result_df.loc[result_df["q"] == 0, "q"] = 1e-323
    # Rename columns
    result_df = result_df.rename(columns = {
        group1 + "_junction_a": "ref_junction_a",
        group1 + "_junction_b": "ref_junction_b",
        group1 + "_PSI": "ref_PSI",
        group2 + "_junction_a": "alt_junction_a",
        group2 + "_junction_b": "alt_junction_b",
        group2 + "_PSI": "alt_PSI"
    })
    return(result_df)

def diff_mxe(df, group_list, FDR, dPSI) -> pd.DataFrame:
    """
    Differential splicing analysis for MXE.

    Args:
    - df (pd.DataFrame): DataFrame containing MXE junction and PSI information.
    - group_list (list): List of two group names to compare
    - FDR (float): False discovery rate
    - dPSI (float): Minimum delta PSI

    Returns:
    - pd.DataFrame: DataFrame containing differential splicing analysis results for MXE.
    """

    result_df = df
    # Drop rows with nan values
    result_df = result_df.dropna()
    result_df = result_df.reset_index()
    group1 = group_list[0]
    group2 = group_list[1]
    group1_junction_a1 = group1 + "_junction_a1"
    group1_junction_a2 = group1 + "_junction_a2"
    group1_junction_b1 = group1 + "_junction_b1"
    group1_junction_b2 = group1 + "_junction_b2"
    group1_PSI = group1 + "_PSI"
    group2_junction_a1 = group2 + "_junction_a1"
    group2_junction_a2 = group2 + "_junction_a2"
    group2_junction_b1 = group2 + "_junction_b1"
    group2_junction_b2 = group2 + "_junction_b2"
    group2_PSI = group2 + "_PSI"
    # Keep columns of junctions and PSI of the two groups
    result_df = result_df[
        ["event_id", "pos_id", "exon_a", "exon_b", "intron_a1", "intron_a2", "intron_b1", "intron_b2", "strand", "gene_id", "gene_name", "label"] + \
        [group1_junction_a1, group1_junction_a2, group1_junction_b1, group1_junction_b2, group1_PSI, group2_junction_a1, group2_junction_a2, group2_junction_b1, group2_junction_b2, group2_PSI]
    ]
    if result_df.shape[0] != 0:
        # delta PSI
        result_df.loc[:, 'dPSI'] = result_df.loc[:, group2_PSI] - result_df.loc[:, group1_PSI]
        # Fisher's exact test, Odds ratio
        result_df = result_df.reset_index()
        result_df = result_df.drop(columns = "index")
        group1_junction_a1_values = result_df[group1_junction_a1].values
        group1_junction_a2_values = result_df[group1_junction_a2].values
        group1_junction_b1_values = result_df[group1_junction_b1].values
        group1_junction_b2_values = result_df[group1_junction_b2].values
        group2_junction_a1_values = result_df[group2_junction_a1].values
        group2_junction_a2_values = result_df[group2_junction_a2].values
        group2_junction_b1_values = result_df[group2_junction_b1].values
        group2_junction_b2_values = result_df[group2_junction_b2].values
        oddsr_junction_a1b1_col = []
        p_junction_a1b1_col = []
        oddsr_junction_a1b2_col = []
        p_junction_a1b2_col = []
        oddsr_junction_a2b1_col = []
        p_junction_a2b1_col = []
        oddsr_junction_a2b2_col = []
        p_junction_a2b2_col = []
        p_maximum_col = []
        for index in range(result_df.shape[0]):
            # Make 2x2 table
            table_2x2_junction_a1b1 = [
                [group1_junction_a1_values[index], group1_junction_b1_values[index]],
                [group2_junction_a1_values[index], group2_junction_b1_values[index]]
            ]
            table_2x2_junction_a1b2 = [
                [group1_junction_a1_values[index], group1_junction_b2_values[index]],
                [group2_junction_a1_values[index], group2_junction_b2_values[index]]
            ]
            table_2x2_junction_a2b1 = [
                [group1_junction_a2_values[index], group1_junction_b1_values[index]],
                [group2_junction_a2_values[index], group2_junction_b1_values[index]]
            ]
            table_2x2_junction_a2b2 = [
                [group1_junction_a2_values[index], group1_junction_b2_values[index]],
                [group2_junction_a2_values[index], group2_junction_b2_values[index]]
            ]
            # Fisher's exact test
            oddsr_junction_a1b1, p_junction_a1b1 = stats.fisher_exact(table_2x2_junction_a1b1, alternative = 'two-sided')
            oddsr_junction_a1b2, p_junction_a1b2 = stats.fisher_exact(table_2x2_junction_a1b2, alternative = 'two-sided')
            oddsr_junction_a2b1, p_junction_a2b1 = stats.fisher_exact(table_2x2_junction_a2b1, alternative = 'two-sided')
            oddsr_junction_a2b2, p_junction_a2b2 = stats.fisher_exact(table_2x2_junction_a2b2, alternative = 'two-sided')
            p_maximum = max([p_junction_a1b1, p_junction_a1b2, p_junction_a2b1, p_junction_a2b2])
            oddsr_junction_a1b1_col.append(oddsr_junction_a1b1)
            p_junction_a1b1_col.append(p_junction_a1b1)
            oddsr_junction_a1b2_col.append(oddsr_junction_a1b2)
            p_junction_a1b2_col.append(p_junction_a1b2)
            oddsr_junction_a2b1_col.append(oddsr_junction_a2b1)
            p_junction_a2b1_col.append(p_junction_a2b1)
            oddsr_junction_a2b2_col.append(oddsr_junction_a2b2)
            p_junction_a2b2_col.append(p_junction_a2b2)
            p_maximum_col.append(p_maximum)
        result_df["OR_junction_a1b1"] = oddsr_junction_a1b1_col
        result_df["p_junction_a1b1"] = p_junction_a1b1_col
        result_df["OR_junction_a1b2"] = oddsr_junction_a1b2_col
        result_df["p_junction_a1b2"] = p_junction_a1b2_col
        result_df["OR_junction_a2b1"] = oddsr_junction_a2b1_col
        result_df["p_junction_a2b1"] = p_junction_a2b1_col
        result_df["OR_junction_a2b2"] = oddsr_junction_a2b2_col
        result_df["p_junction_a2b2"] = p_junction_a2b2_col
        result_df["p_maximum"] = p_maximum_col
        # FDR correction
        result_df.loc[:, "q"] = multitest.multipletests(result_df.loc[:, "p_maximum"], method = "fdr_bh")[1]
        result_df.loc[
            (((result_df.loc[:, "OR_junction_a1b1"] >= 3/2) & (result_df.loc[:, "OR_junction_a1b2"] >= 3/2) & (result_df.loc[:, "OR_junction_a2b1"] >= 3/2) & (result_df.loc[:, "OR_junction_a2b2"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a1b1"] <= 2/3) & (result_df.loc[:, "OR_junction_a1b2"] <= 2/3) & (result_df.loc[:, "OR_junction_a2b1"] <= 2/3) & (result_df.loc[:, "OR_junction_a2b2"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI),
            "Diff events"
        ] = "Yes"
        result_df.loc[
            ~((((result_df.loc[:, "OR_junction_a1b1"] >= 3/2) & (result_df.loc[:, "OR_junction_a1b2"] >= 3/2) & (result_df.loc[:, "OR_junction_a2b1"] >= 3/2) & (result_df.loc[:, "OR_junction_a2b2"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a1b1"] <= 2/3) & (result_df.loc[:, "OR_junction_a1b2"] <= 2/3) & (result_df.loc[:, "OR_junction_a2b1"] <= 2/3) & (result_df.loc[:, "OR_junction_a2b2"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI)),
            "Diff events"
        ] = "No"
        result_df.loc[result_df["q"] == 0, "q"] = 1e-323
    # Rename columns
    result_df = result_df.rename(columns = {
        group1 + "_junction_a1": "ref_junction_a1",
        group1 + "_junction_a2": "ref_junction_a2",
        group1 + "_junction_b1": "ref_junction_b1",
        group1 + "_junction_b2": "ref_junction_b2",
        group1 + "_PSI": "ref_PSI",
        group2 + "_junction_a1": "alt_junction_a1",
        group2 + "_junction_a2": "alt_junction_a2",
        group2 + "_junction_b1": "alt_junction_b1",
        group2 + "_junction_b2": "alt_junction_b2",
        group2 + "_PSI": "alt_PSI"
    })
    return(result_df)

def diff_ri(df, group_list, FDR, dPSI) -> pd.DataFrame:
    """
    Differential splicing analysis for RI.

    Args:
    - df (pd.DataFrame): Dataframe containing the splicing events.
    - group_list (list): List of two groups to compare.
    - FDR (float): False discovery rate.
    - dPSI (float): Minimum delta PSI.

    Returns:
    pd.DataFrame: Dataframe containing the differential splicing analysis results.
    """

    result_df = df
    # Drop rows with nan values
    result_df = result_df.dropna()
    result_df = result_df.reset_index()
    group1 = group_list[0]
    group2 = group_list[1]
    group1_junction_a = group1 + "_junction_a"
    group1_junction_a_start = group1 + "_junction_a_start"
    group1_junction_a_end = group1 + "_junction_a_end"
    group1_PSI = group1 + "_PSI"
    group2_junction_a = group2 + "_junction_a"
    group2_junction_a_start = group2 + "_junction_a_start"
    group2_junction_a_end = group2 + "_junction_a_end"
    group2_PSI = group2 + "_PSI"
    # Keep columns of junctions and PSI of the two groups
    result_df = result_df[
        ["event_id", "pos_id", "exon_a", "exon_b", "exon_c", "intron_a", "strand", "gene_id", "gene_name", "label"] + \
        [group1_junction_a, group1_junction_a_start, group1_junction_a_end, group1_PSI, group2_junction_a, group2_junction_a_start, group2_junction_a_end, group2_PSI]
    ]
    if result_df.shape[0] != 0:
        # delta PSI
        result_df.loc[:, 'dPSI'] = result_df.loc[:, group2_PSI] - result_df.loc[:, group1_PSI]
        # Fisher's exact test, Odds ratio
        result_df = result_df.reset_index()
        result_df = result_df.drop(columns = "index")
        group1_junction_a_values = result_df[group1_junction_a].values
        group1_junction_a_start_values = result_df[group1_junction_a_start].values
        group1_junction_a_end_values = result_df[group1_junction_a_end].values
        group2_junction_a_values = result_df[group2_junction_a].values
        group2_junction_a_start_values = result_df[group2_junction_a_start].values
        group2_junction_a_end_values = result_df[group2_junction_a_end].values
        oddsr_junction_a_start_col = []
        p_junction_a_start_col = []
        oddsr_junction_a_end_col = []
        p_junction_a_end_col = []
        p_maximum_col = []
        for index in range(result_df.shape[0]):
            # Make 2x2 table
            # inc1 - exc
            table_2x2_junction_a_start = [
                [group1_junction_a_values[index], group1_junction_a_start_values[index]],
                [group2_junction_a_values[index], group2_junction_a_start_values[index]]
            ]
            # inc2 - exc
            table_2x2_junction_a_end = [
                [group1_junction_a_values[index], group1_junction_a_end_values[index]],
                [group2_junction_a_values[index], group2_junction_a_end_values[index]]
            ]
            # Fisher's exact test
            oddsr_junction_a_start, p_junction_a_start = stats.fisher_exact(table_2x2_junction_a_start, alternative = 'two-sided')
            oddsr_junction_a_end, p_junction_a_end = stats.fisher_exact(table_2x2_junction_a_end, alternative = 'two-sided')
            p_maximum = max([p_junction_a_start, p_junction_a_end])
            oddsr_junction_a_start_col.append(oddsr_junction_a_start)
            p_junction_a_start_col.append(p_junction_a_start)
            oddsr_junction_a_end_col.append(oddsr_junction_a_end)
            p_junction_a_end_col.append(p_junction_a_end)
            p_maximum_col.append(p_maximum)
        result_df["OR_junction_a_start"] = oddsr_junction_a_start_col
        result_df["p_junction_a_start"] = p_junction_a_start_col
        result_df["OR_junction_a_end"] = oddsr_junction_a_end_col
        result_df["p_junction_a_end"] = p_junction_a_end_col
        result_df["p_maximum"] = p_maximum_col
        # FDR correction
        result_df.loc[:, "q"] = multitest.multipletests(result_df.loc[:, "p_maximum"], method = "fdr_bh")[1]
        result_df.loc[
            (((result_df.loc[:, "OR_junction_a_start"] >= 3/2) & (result_df.loc[:, "OR_junction_a_end"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a_start"] <= 2/3) & (result_df.loc[:, "OR_junction_a_end"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI),
            "Diff events"
        ] = "Yes"
        result_df.loc[
            ~((((result_df.loc[:, "OR_junction_a_start"] >= 3/2) & (result_df.loc[:, "OR_junction_a_end"] >= 3/2)) |
            ((result_df.loc[:, "OR_junction_a_start"] <= 2/3) & (result_df.loc[:, "OR_junction_a_end"] <= 2/3))) &
            (result_df.loc[:, 'q'] < FDR) &
            (result_df.loc[:, 'dPSI'].abs() >= dPSI)),
            "Diff events"
        ] = "No"
        result_df.loc[result_df["q"] == 0, "q"] = 1e-323
    # Rename columns
    result_df = result_df.rename(columns = {
        group1 + "_junction_a": "ref_junction_a",
        group1 + "_junction_a_start": "ref_junction_a_start",
        group1 + "_junction_a_end": "ref_junction_a_end",
        group1 + "_PSI": "ref_PSI",
        group2 + "_junction_a": "alt_junction_a",
        group2 + "_junction_a_start": "alt_junction_a_start",
        group2 + "_junction_a_end": "alt_junction_a_end",
        group2 + "_PSI": "alt_PSI"
    })
    return(result_df)

def ttest(output_ind_df, group_df, group_list) -> pd.DataFrame:
    """
    Performs a t-test on the PSI values of two groups and adds a column with the resulting p-values to the output dataframe.

    Args:
    - output_ind_df (pd.DataFrame): The dataframe containing the PSI values for each sample.
    - group_df (pd.DataFrame): The dataframe containing the group assignments for each sample.
    - group_list (list): A list of two strings representing the names of the two groups being compared.

    Returns:
    - pd.DataFrame: The input dataframe with an additional column containing the p-values resulting from the t-test.
    """

    output_ind_df = output_ind_df.reset_index()
    output_ind_df = output_ind_df.drop(columns = "index")
    group1 = group_list[0]
    group2 = group_list[1]
    sample_group1 = list(group_df[group_df['group'] == group1]["sample"])
    sample_group1 = [i + "_PSI" for i in sample_group1]
    try:
        sample_group1 = [output_ind_df[i].values for i in sample_group1]
    except:
        logger.debug(f"Sample names do not match.")
        logger.debug(f"Columns: {output_ind_df.columns}")
        logger.debug(f"Sample names: {sample_group1}")
        raise ValueError("Error: Sample names do not match.")
    sample_group2 = list(group_df[group_df['group'] == group2]["sample"])
    sample_group2 = [i + "_PSI" for i in sample_group2]
    sample_group2 = [output_ind_df[i].values for i in sample_group2]
    p_col = []
    for index in range(output_ind_df.shape[0]):
        PSI_group1 = [i[index] for i in sample_group1]
        PSI_group2 = [i[index] for i in sample_group2]
        # t-test
        t, p = stats.ttest_ind(
            PSI_group1,
            PSI_group2,
            equal_var = False,
            nan_policy = "omit"
        )
        p_col.append(p)
    output_ind_df["p_ttest"] = p_col
    return(output_ind_df)

def make_psi_table_sample(sample_list, event_for_analysis_df, junc_dict_all, func_psi, func_col, num_process, minimum_reads) -> pd.DataFrame:
    """
    Make PSI table for each sample.

    Args:
    - sample_list (list): List of sample names.
    - event_for_analysis_df (pd.DataFrame): DataFrame containing the splicing events to be analyzed.
    - junc_dict_all (dict): Dictionary containing the junction information for each sample.
    - func_psi (function): Function to calculate PSI values.
    - func_col (function): Function to make column names.
    - num_process (int): Number of processes to use.
    - minimum_reads (int): Minimum number of reads to be considered.

    Returns:
    - pd.DataFrame: DataFrame containing the PSI values for each sample and each event.

    """

    columns = func_col(sample_list, False)
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_process) as executor:
        futures = [executor.submit(func_psi, junc_dict_all, sample_list, event_for_analysis_df, num_process, minimum_reads, i) for i in range(num_process)]
    output_l = []
    for future in concurrent.futures.as_completed(futures):
        output_l += future.result()
    psi_table_df = pd.DataFrame(
		output_l,
		columns = columns
	)
    return(psi_table_df)

def make_psi_table_group(group_list, event_for_analysis_df, junc_dict_group, func_psi, func_col, num_process, minimum_reads) -> pd.DataFrame:
    """
    Make PSI table for each group.

    Args:
    - group_list (list): List of group names.
    - event_for_analysis_df (pd.DataFrame): DataFrame containing the splicing events to be analyzed.
    - junc_dict_group (dict): Dictionary containing the junction information for each group.
    - func_psi (function): Function to calculate PSI values.
    - func_col (function): Function to make column names.
    - num_process (int): Number of processes to use.
    - minimum_reads (int): Minimum number of reads to be considered.

    Returns:
    - pd.DataFrame: DataFrame containing the PSI values for each group and each event.

    """

    columns = func_col(group_list, True)
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_process) as executor:
        futures = [executor.submit(func_psi, junc_dict_group, group_list, event_for_analysis_df, num_process, minimum_reads, i) for i in range(num_process)]
    output_l = []
    for future in concurrent.futures.as_completed(futures):
        output_l += future.result()
    psi_table_df = pd.DataFrame(
		output_l,
		columns = columns
	)
    return(psi_table_df)

def diff_event(event_for_analysis_df, psi_table_df, junc_dict_all, group_df, group_list, sample_list, func_diff, func_ind, num_process, FDR, dPSI, individual_psi, ttest_bool) -> pd.DataFrame:
    """
    Differential splicing analysis for each splicing event.

    Args:
    - event_for_analysis_df (pd.DataFrame): DataFrame containing the splicing events to be analyzed.
    - psi_table_df (pd.DataFrame): DataFrame containing the PSI values for each sample and each event.
    - junc_dict_all (dict): Dictionary containing the junction information for each sample.
    - group_df (pd.DataFrame): DataFrame containing the group assignments for each sample.
    - group_list (list): List of group names.
    - sample_list (list): List of sample names.
    - func_diff (function): Function to perform differential splicing analysis.
    - func_ind (function): Function to perform individual PSI analysis.
    - num_process (int): Number of processes to use.
    - FDR (float): False discovery rate.
    - dPSI (float): Minimum delta PSI.
    - individual_psi (bool): Whether to perform individual PSI analysis.
    - ttest (bool): Whether to perform t-test.

    Returns:
    - pd.DataFrame: DataFrame containing the differential splicing analysis results for each splicing event.

    """

    output_df = func_diff(psi_table_df, group_list, FDR, dPSI)
    if (output_df.shape[0]) != 0:
        if individual_psi:
            event_for_analysis_df = event_for_analysis_df[event_for_analysis_df["event_id"].isin(output_df["event_id"])]
            with concurrent.futures.ProcessPoolExecutor(max_workers = num_process) as executor:
                futures = [executor.submit(func_ind, junc_dict_all, event_for_analysis_df, sample_list, num_process, i) for i in range(num_process)]
            output_l = []
            for future in concurrent.futures.as_completed(futures):
                output_l += future.result()
            columns_ind = col_ind(sample_list)
            output_ind_df = pd.DataFrame(
                output_l,
                columns = columns_ind
            )
            if ttest_bool:
                output_ind_df = ttest(output_ind_df, group_df, group_list)
            output_df = pd.merge(
                output_df,
                output_ind_df,
                on = "event_id"
            )
        if "q" in output_df.columns:
            output_df = output_df.sort_values(["Diff events", "q"], ascending = [False, True])
    return(output_df)

def make_psi_mtx(psi_table_df) -> pd.DataFrame:
    """
    Make PSI matrix.

    Args:
    - psi_table_df (pd.DataFrame): DataFrame containing the PSI values for each sample and each event.

    Returns:
    - pd.DataFrame: DataFrame containing the PSI values for each event.

    """

    psi_table_df["key"] = psi_table_df["event_id"].str.split("_", expand = True)[1].astype(int)
    psi_table_df = psi_table_df.sort_values("key")
    psi_table_df = psi_table_df.drop(columns = ["key"])
    # Simple PSI matrix
    ind_psi_col = ["event_id", "pos_id"] + [i for i in list(psi_table_df.columns) if i.endswith("_PSI")]
    output_mtx_df = psi_table_df[ind_psi_col]
    output_mtx_df.columns = [i.rstrip("_PSI") for i in list(output_mtx_df.columns)]
    return(psi_table_df, output_mtx_df)

class EventCounter:
    """
    Class to count events.
    """

    def __init__(self, df, threshold):
        self.df = df
        self.labels = ["annotated", "unannotated"]
        self.threshold = threshold

    def count_events(self, diff, direction):
        try:
            return self.df[(self.df["Diff events"] == "Yes") & (self.df["dPSI"] * direction > self.threshold) & (self.df["label"] == diff)].shape[0]
        except:
            return 0

    def count_all_events(self):
        return {
            "up_annotated_num": self.count_events("annotated", 1),
            "up_unannotated_num": self.count_events("unannotated", 1),
            "down_annotated_num": self.count_events("annotated", -1),
            "down_unannotated_num": self.count_events("unannotated", -1),
        }

def save_excel(output_path, SE_df, FIVE_df, THREE_df, MXE_df, RI_df, MSE_df, AFE_df, ALE_df):
    """
    Save excel file.

    Args:
    - output_path (str): Output path.
    - SE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for SE events.
    - FIVE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for FIVE events.
    - THREE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for THREE events.
    - MXE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for MXE events.
    - RI_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for RI events.
    - MSE_df (pd.DataFrame): DataFrame containing the differential splicing events for MSE events.
    - AFE_df (pd.DataFrame): DataFrame containing the differential splicing events for AFE events.
    - ALE_df (pd.DataFrame): DataFrame containing the differential splicing events for ALE events.

    """

    from styleframe import StyleFrame, Styler, utils

    # Style
    style = Styler(
        horizontal_alignment = utils.horizontal_alignments.left,
        border_type = utils.borders.default_grid,
        wrap_text = False
    )
    with StyleFrame.ExcelWriter(output_path + "/results.xlsx") as writer:
        # SE
        SE_sf = StyleFrame(SE_df)
        SE_sf.set_column_width(columns = SE_df.columns, width = 20)
        SE_sf.apply_column_style(cols_to_style = SE_df.columns, styler_obj = style, style_header = True)
        SE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "SE")
        # FIVE
        FIVE_sf = StyleFrame(FIVE_df)
        FIVE_sf.set_column_width(columns = FIVE_df.columns, width = 20)
        FIVE_sf.apply_column_style(cols_to_style = FIVE_df.columns, styler_obj = style, style_header = True)
        try:
            FIVE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "FIVE")
        except:
            pass
        # THREE
        THREE_sf = StyleFrame(THREE_df)
        THREE_sf.set_column_width(columns = THREE_df.columns, width = 20)
        THREE_sf.apply_column_style(cols_to_style = THREE_df.columns, styler_obj = style, style_header = True)
        try:
            THREE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "THREE")
        except:
            pass
        # MXE
        MXE_sf = StyleFrame(MXE_df)
        MXE_sf.set_column_width(columns = MXE_df.columns, width = 20)
        MXE_sf.apply_column_style(cols_to_style = MXE_df.columns, styler_obj = style, style_header = True)
        try:
            MXE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "MXE")
        except:
            pass
        # RI
        RI_sf = StyleFrame(RI_df)
        RI_sf.set_column_width(columns = RI_df.columns, width = 20)
        RI_sf.apply_column_style(cols_to_style = RI_df.columns, styler_obj = style, style_header = True)
        try:
            RI_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "RI")
        except:
            pass
        # MSE
        MSE_sf = StyleFrame(MSE_df)
        MSE_sf.set_column_width(columns = MSE_df.columns, width = 20)
        MSE_sf.apply_column_style(cols_to_style = MSE_df.columns, styler_obj = style, style_header = True)
        try:
            MSE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "MSE")
        except:
            pass
        # AFE
        AFE_sf = StyleFrame(AFE_df)
        AFE_sf.set_column_width(columns = AFE_df.columns, width = 20)
        AFE_sf.apply_column_style(cols_to_style = AFE_df.columns, styler_obj = style, style_header = True)
        try:
            AFE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "AFE")
        except:
            pass
        # ALE
        ALE_sf = StyleFrame(ALE_df)
        ALE_sf.set_column_width(columns = ALE_df.columns, width = 20)
        ALE_sf.apply_column_style(cols_to_style = ALE_df.columns, styler_obj = style, style_header = True)
        try:
            ALE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "ALE")
        except:
            pass

def save_excel_sc(output_path, SE_df, FIVE_df, THREE_df, MXE_df, MSE_df, AFE_df, ALE_df):
    """
    Save excel file.

    Args:
    - output_path (str): Output path.
    - SE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for SE events.
    - FIVE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for FIVE events.
    - THREE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for THREE events.
    - MXE_df (pd.DataFrame): DataFrame containing the differential splicing analysis results for MXE events.
    - MSE_df (pd.DataFrame): DataFrame containing the differential splicing events for MSE events.
    - AFE_df (pd.DataFrame): DataFrame containing the differential splicing events for AFE events.
    - ALE_df (pd.DataFrame): DataFrame containing the differential splicing events for ALE events.

    """

    from styleframe import StyleFrame, Styler, utils

    # Style
    style = Styler(
        horizontal_alignment = utils.horizontal_alignments.left,
        border_type = utils.borders.default_grid,
        wrap_text = False
    )
    with StyleFrame.ExcelWriter(output_path + "/results.xlsx") as writer:
        # SE
        SE_sf = StyleFrame(SE_df)
        SE_sf.set_column_width(columns = SE_df.columns, width = 20)
        SE_sf.apply_column_style(cols_to_style = SE_df.columns, styler_obj = style, style_header = True)
        SE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "SE")
        # FIVE
        FIVE_sf = StyleFrame(FIVE_df)
        FIVE_sf.set_column_width(columns = FIVE_df.columns, width = 20)
        FIVE_sf.apply_column_style(cols_to_style = FIVE_df.columns, styler_obj = style, style_header = True)
        try:
            FIVE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "FIVE")
        except:
            pass
        # THREE
        THREE_sf = StyleFrame(THREE_df)
        THREE_sf.set_column_width(columns = THREE_df.columns, width = 20)
        THREE_sf.apply_column_style(cols_to_style = THREE_df.columns, styler_obj = style, style_header = True)
        try:
            THREE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "THREE")
        except:
            pass
        # MXE
        MXE_sf = StyleFrame(MXE_df)
        MXE_sf.set_column_width(columns = MXE_df.columns, width = 20)
        MXE_sf.apply_column_style(cols_to_style = MXE_df.columns, styler_obj = style, style_header = True)
        try:
            MXE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "MXE")
        except:
            pass
        # MSE
        MSE_sf = StyleFrame(MSE_df)
        MSE_sf.set_column_width(columns = MSE_df.columns, width = 20)
        MSE_sf.apply_column_style(cols_to_style = MSE_df.columns, styler_obj = style, style_header = True)
        try:
            MSE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "MSE")
        except:
            pass
        # AFE
        AFE_sf = StyleFrame(AFE_df)
        AFE_sf.set_column_width(columns = AFE_df.columns, width = 20)
        AFE_sf.apply_column_style(cols_to_style = AFE_df.columns, styler_obj = style, style_header = True)
        try:
            AFE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "AFE")
        except:
            pass
        # ALE
        ALE_sf = StyleFrame(ALE_df)
        ALE_sf.set_column_width(columns = ALE_df.columns, width = 20)
        ALE_sf.apply_column_style(cols_to_style = ALE_df.columns, styler_obj = style, style_header = True)
        try:
            ALE_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "ALE")
        except:
            pass
