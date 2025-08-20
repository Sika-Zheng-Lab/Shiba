import warnings
warnings.simplefilter('ignore')
import argparse
import sys
import os
import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns
import html
import logging

# Configure logging
logger = logging.getLogger(__name__)

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Make plots for alternative splicing events"
	)
	parser.add_argument("-i", "--input", type = str, help = "Directory that contains result files")
	parser.add_argument("-e", "--experiment-table", type = str, help = "Experiment table file")
	parser.add_argument("-s", "--shiba-command", type = str, help = "Shiba command")
	parser.add_argument("-o", "--output", type = str, help = "Directory for output files")
	parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
	args = parser.parse_args()
	return(args)

def load_experiment_table(experiment_table: str):

	# Load experiment table
	experiment_table_df = pd.read_csv(
		experiment_table,
		sep = "\t",
		usecols = ["sample", "group"]
	)
	# Add a column of number for each group, starting from the first group to the last group
	group_order = experiment_table_df["group"].unique().tolist()
	group_order_dict = {group: i for i, group in enumerate(group_order)}
	experiment_table_df["group_order"] = experiment_table_df["group"].map(group_order_dict)
	return experiment_table_df

def load_tpm_pca_table(input_dir: str, experiment_table_df: pd.DataFrame, output_dir: str):

	# Load PCA matrix for TPM
	pca_tpm_df = pd.read_csv(
		os.path.join(input_dir, "pca", "tpm_pca.tsv"),
		sep = "\t"
	)
	pca_tpm_df = pca_tpm_df.rename(columns = {"Unnamed: 0": "sample"})
	pca_tpm_df = pd.merge(pca_tpm_df, experiment_table_df, on = "sample")
	pca_tpm_df = pca_tpm_df.sort_values("group_order")
	# Load contribution table for TPM
	contribution_tpm_df = pd.read_csv(
		os.path.join(input_dir, "pca", "tpm_contribution.tsv"),
		sep = "\t",
		names = ["PC", "contribution"]
	)
	contribution_tpm_PC1 = str((contribution_tpm_df.iloc[0][1]*100).round(2))
	contribution_tpm_PC2 = str((contribution_tpm_df.iloc[1][1]*100).round(2))
	return pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2

def load_psi_pca_table(input_dir: str, experiment_table_df: pd.DataFrame, output_dir: str):

	# Load PCA matrix for PSI
	pca_psi_df = pd.read_csv(
		os.path.join(input_dir, "pca", "psi_pca.tsv"),
		sep = "\t"
	)
	pca_psi_df = pca_psi_df.rename(columns = {"Unnamed: 0": "sample"})
	pca_psi_df = pd.merge(pca_psi_df, experiment_table_df, on = "sample")
	pca_psi_df = pca_psi_df.sort_values("group_order")
	# Load contribution table for PSI
	contribution_psi_df = pd.read_csv(
		os.path.join(input_dir, "pca", "psi_contribution.tsv"),
		sep = "\t",
		names = ["PC", "contribution"]
	)
	contribution_psi_PC1 = str((contribution_psi_df.iloc[0][1]*100).round(2))
	contribution_psi_PC2 = str((contribution_psi_df.iloc[1][1]*100).round(2))
	return pca_psi_df, contribution_psi_PC1, contribution_psi_PC2

def plots_pca(name: str, pca_df: pd.DataFrame, contribution_PC1: str, contribution_PC2: str, output_dir: str):

	# Round PC1 and PC2
	pca_df["PC1"] = pca_df["PC1"].round(2)
	pca_df["PC2"] = pca_df["PC2"].round(2)
	# Color palette
	number_of_groups = pca_df["group"].nunique()
	if number_of_groups <= 8:
		color_palette = px.colors.qualitative.G10
	else:
		color_palette = px.colors.sequential.Viridis
	fig = px.scatter(
		pca_df,
		x = "PC1",
		y = "PC2",
		color = "group",
		color_discrete_sequence = color_palette,
		opacity = 0.5,
		hover_data = ["sample"]
	)
	fig.update_traces(
		marker = dict(
			size = 16,
			line = dict(width = 1, color = 'Black')),
			selector = dict(mode = 'markers')
	)
	fig.update_layout(
		title = dict(
			text = "PCA for " + name,
			font = dict(size = 26, color = 'black'),
			xref = 'paper',
			x = 0.5,
			y = 0.95,
			xanchor = 'center',
		)
	)
	fig.update_layout(
		width = 550,
		height = 400,
		font ={
			"family": "Arial",
			"size": 18
		},
		xaxis_title = "PC1 ({}%)".format(contribution_PC1),
		yaxis_title = "PC2 ({}%)".format(contribution_PC2),
		legend_title = "Group",
	)
	fig.write_html(os.path.join(output_dir, "data", "pca_" + name + ".html"), include_plotlyjs = "cdn")

def plots(AS: str, input_dir: str, output_dir: str):

	# load data
	df = pd.read_csv(
		os.path.join(input_dir, "splicing", "PSI_" + AS + ".txt"),
		sep = "\t"
	)
	if not df.empty:
		# Round dPSI and others
		df["dPSI"] = df["dPSI"].round(2)
		df['-log10(q)'] = -np.log10(df["q"])
		df['-log10(q)'] = df['-log10(q)'].round(2)
		df["ref_PSI"] = df["ref_PSI"].round(2)
		df["alt_PSI"] = df["alt_PSI"].round(2)
	# Volcano plot
	if not df.empty:
		df.loc[(df["Diff events"] == "Yes") & (df["dPSI"] > 0.1), "group"] = "up"
		df.loc[(df["Diff events"] == "Yes") & (df["dPSI"] < -0.1), "group"] = "down"
		df = df.fillna({"group": "others"})
		fig = px.scatter(
			df,
			x = "dPSI",
			y = "-log10(q)",
			color = "group",
			symbol="label",
			opacity = 0.5,
			category_orders = {"group": ["up", "down", "others"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["salmon", "steelblue", "lightgrey"],
			hover_data = ["gene_name", "event_id", "ref_PSI", "alt_PSI", "q"]
		)
		fig.update_traces(
			marker = dict(
				size = 8,
				line = dict(width = 0, color = 'DarkSlateGrey')),
				selector = dict(mode = 'markers')
		)
		fig.update_layout(
			title = dict(
				text = "Volcano plot",
				font = dict(size = 26, color = 'black'),
				xref = 'paper',
				x = 0.5,
				y = 0.95,
				xanchor = 'center'
			)
		)
	else:
		# If there is no data, make empty plot
		fig = px.scatter(
			x = [None],
			y = [None]
		)
		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)
		# Add "No data" text to the empty plot
		fig.add_annotation(
			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')
		)
	fig.update_layout(
		width = 550,
		height = 400,
		font ={
			"family": "Arial",
			"size": 16
		},
		legend_title = "Group, Label",
	)
	fig.write_html(os.path.join(output_dir, "data", "volcano_" + AS + ".html"), include_plotlyjs = "cdn")

	# Scatter
	if not df.empty:
		fig = px.scatter(
			df,
			x = "ref_PSI",
			y = "alt_PSI",
			color = "group",
			symbol="label",
			opacity = 0.5,
			category_orders = {"group": ["up", "down", "others"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["salmon", "steelblue", "lightgrey"],
			hover_data = ["gene_name", "event_id", "dPSI", "q"]
		)
		fig.update_traces(
			marker = dict(size=8,
							line=dict(width=0, color='DarkSlateGrey')),
							selector=dict(mode='markers')
		)
		fig.update_layout(title=dict(text = "Scatter plot",
										font=dict(size=26,
												color='black'),
										xref='paper',
										x=0.5,
										y=0.95,
										xanchor='center'
									)
							)
	else:
		# If there is no data, make empty plot
		fig = px.scatter(
			x = [None],
			y = [None]
		)
		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)
		# Add "No data" text to the empty plot
		fig.add_annotation(
			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')
		)
	fig.update_layout(
		width = 550,
		height= 400,
		font ={
			"family": "Arial",
			"size": 16
		},
		xaxis_title = "PSI (Reference)",
		yaxis_title = "PSI (Alternative)",
		legend_title = "Group, Label",
	)
	fig.write_html(os.path.join(output_dir, "data", "scatter_" + AS + ".html"), include_plotlyjs = "cdn")

	# Bar
	if not df.empty:
		count_df = df[df["group"] != "others"].groupby(["group", "label"], as_index = False).count()[["group", "label", "event_id"]]
		fig = px.bar(
			count_df,
			x = "group",
			y = "event_id",
			color = "label",
			labels = {"event_id": "Count"},
			category_orders = {"group": ["up", "down"], "label": ["annotated", "unannotated"]},
			color_discrete_sequence = ["#9ebcda", "#810f7c"],
			barmode = "relative"
		)
		fig.update_layout(title=dict(text = "Number of DSEs",
										font=dict(size=26,
												color='black'),
										xref='paper',
										x=0.5,
										y=0.95,
										xanchor='center'
									)
							)
	else:
		# if there is no data, make empty plot
		fig = px.bar(
			x = [0],
			y = [0]
		)
		# Remove xlabel and ylabel
		fig.update_xaxes(title = None)
		fig.update_yaxes(title = None)
		# Add "No data" text to the empty plot
		fig.add_annotation(
			text = "No data",
			xref = "paper",
			yref = "paper",
			x = 0.5,
			y = 0.5,
			showarrow = False,
			font = dict(size = 26, color = 'black')
		)
	fig.update_layout(
		width = 550,
		height = 400,
		font ={
			"family": "Arial",
			"size": 16
		},
		xaxis_title = None,
		legend_title = "Group, Label",
	)
	fig.write_html(os.path.join(output_dir, "data", "bar_" + AS + ".html"), include_plotlyjs = "cdn")

def write_summary_html(shiba_command: str, output_dir: str):

	# PCA
	lines_strip_pca_dict = {}
	for pca in ["TPM", "PSI"]:
		with open(os.path.join(output_dir, "data", "pca_" + pca + ".html"), 'r') as f:
			lines = f.readlines()
		lines_strip = [html.escape(line.strip()) for line in lines[2:6]]
		lines_strip_pca_dict[pca] = lines_strip
	# Splicing events
	events = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]
	plottypes = ["volcano", "scatter", "bar"]
	lines_strip_dict = {}
	for event in events:
		lines_strip_dict[event] = {}
		for plottype in plottypes:
			with open(os.path.join(output_dir, "data", plottype + "_" + event + ".html"), 'r') as f:
				lines = f.readlines()
			lines_strip = [html.escape(line.strip()) for line in lines[2:6]]
			lines_strip_dict[event][plottype] = lines_strip

	summary_html = '''
	<!DOCTYPE html>
	<html lang="en">
	<head>
		<meta charset="UTF-8">
		<meta name="viewport" content="width=device-width, initial-scale=1.0">
		<title>Shiba Results Summary</title>
		<style>
			html {{
				scroll-behavior: smooth;
			}}

			body {{
				font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
				margin: 0;
				display: flex;
				flex-direction: column;
				min-height: 100vh;
				width: 95vw;
				background: linear-gradient(135deg, #000000, #3b003b, #7b0067);
				color: #ffffff;
			}}

			header {{
				width: 100vw;
				box-sizing: border-box;
				padding: 20px;
				margin: 0;
				background-color: #1e1e1e;
				text-align: center;
				box-shadow: 0 2px 5px rgba(0, 0, 0, 0.5);
			}}

			header h1 {{
				margin: 0;
				color: #b300b3;
				font-size: 52px;
			}}

			.sidebar {{
				height: calc(100vh - 240px); /* Adjusted to account for header and footer height, and text size */
				width: 250px;
				position: fixed;
				top: 120px; /* Adjusted to account for header height */
				left: 0;
				background-color: #1e1e1e;
				padding-top: 20px;
				display: flex;
				flex-direction: column;
				box-shadow: 2px 0 5px rgba(0, 0, 0, 0.5);
				align-items: center;
			}}

			.sidebar a {{
				padding: 15px 20px;
				text-decoration: none;
				font-size: 18px;
				color: white;
				display: block;
				transition: 0.3s;
				width: 80%;
				text-align: left;
			}}

			.sidebar a:hover {{
				background-color: #7b0067;
			}}

			.sidebar a.active {{
				background-color: #b300b3;
				color: white;
			}}

			.content {{
				margin-left: 250px;
				padding: 20px;
				width: calc(100% - 250px);
				overflow-y: auto; /* Enable vertical scrolling for the content */
				overflow-x: hidden; /* Disable horizontal scrolling */
				flex: 1;
			}}

			.plot {{
				padding: 12px;
				background: rgba(255, 255, 255, 0.1);
				width: calc(100%);
				margin-bottom: 20px;
				border-radius: 12px;
				box-shadow: 0 0 15px rgba(0, 0, 0, 0.3);
			}}

			.plot h2 {{
				margin-top: 0;
				color: #b300b3;
			}}

			.plot-container {{
				display: flex;
				gap: 20px;
				justify-content: flex-start; /* Align items to the start of the container */
				flex-wrap: nowrap; /* Prevent wrapping */
				overflow-x: auto; /* Enable horizontal scrolling if necessary */
			}}

			.plot-container iframe {{
				flex: 0 0 32%; /* Set a flex width for iframes */
				height: 420px; /* Set a fixed height for iframes */
				width: 100%; /* Ensure iframes fill the container */
				border: none;
				border-radius: 12px;
				box-shadow: 0 0 15px rgba(0, 0, 0, 0.3);
			}}

			footer {{
				width: 100vw;
				box-sizing: border-box;
				padding: 20px;
				margin: 0;
				background-color: #1e1e1e;
				text-align: center;
				box-shadow: 0 -2px 5px rgba(0, 0, 0, 0.5);
				position: relative;
			}}
		</style>
	</head>
	<body>
		<header>
			<h1>Shiba Results Summary</h1>
		</header>
		<div class="sidebar">
			<a href="#Shiba command">Shiba command</a>
			<a href="#PCA">PCA</a>
			<a href="#SE">Skipped exon (SE)</a>
			<a href="#FIVE">Alternative 5'ss (FIVE)</a>
			<a href="#THREE">Alternative 3'ss (THREE)</a>
			<a href="#MXE">Mutually exclusive exons (MXE)</a>
			<a href="#RI">Retained intron (RI)</a>
			<a href="#MSE">Multiple skipped exons (MSE)</a>
			<a href="#AFE">Alternative first exon (AFE)</a>
			<a href="#ALE">Alternative last exon (ALE)</a>
		</div>

		<div class="content">
			<div id="Shiba command" class="plot">
				<h2>Shiba command</h2>
				<p class="intro-text">{shiba_command}</p>
			</div>
			<div id="PCA" class="plot">
				<h2>Principal component analysis (PCA)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{pca_tpm_l1}
						{pca_tpm_l2}
						{pca_tpm_l3}
						{pca_tpm_l4}
					"></iframe>
					<iframe srcdoc="
						{pca_psi_l1}
						{pca_psi_l2}
						{pca_psi_l3}
						{pca_psi_l4}
					"></iframe>
				</div>
			</div>
			<div id="SE" class="plot">
				<h2>Skipped exon (SE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_se_l1}
						{volcano_se_l2}
						{volcano_se_l3}
						{volcano_se_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_se_l1}
						{scatter_se_l2}
						{scatter_se_l3}
						{scatter_se_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_se_l1}
						{bar_se_l2}
						{bar_se_l3}
						{bar_se_l4}
					"></iframe>
				</div>
			</div>
			<div id="FIVE" class="plot">
				<h2>Alternative 5'ss (FIVE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_five_l1}
						{volcano_five_l2}
						{volcano_five_l3}
						{volcano_five_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_five_l1}
						{scatter_five_l2}
						{scatter_five_l3}
						{scatter_five_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_five_l1}
						{bar_five_l2}
						{bar_five_l3}
						{bar_five_l4}
					"></iframe>
				</div>
			</div>
			<div id="THREE" class="plot">
				<h2>Alternative 3'ss (THREE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_three_l1}
						{volcano_three_l2}
						{volcano_three_l3}
						{volcano_three_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_three_l1}
						{scatter_three_l2}
						{scatter_three_l3}
						{scatter_three_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_three_l1}
						{bar_three_l2}
						{bar_three_l3}
						{bar_three_l4}
					"></iframe>
				</div>
			</div>
			<div id="MXE" class="plot">
				<h2>Mutually exclusive exons (MXE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_mxe_l1}
						{volcano_mxe_l2}
						{volcano_mxe_l3}
						{volcano_mxe_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_mxe_l1}
						{scatter_mxe_l2}
						{scatter_mxe_l3}
						{scatter_mxe_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_mxe_l1}
						{bar_mxe_l2}
						{bar_mxe_l3}
						{bar_mxe_l4}
					"></iframe>
				</div>
			</div>
			<div id="RI" class="plot">
				<h2>Retained intron (RI)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_ri_l1}
						{volcano_ri_l2}
						{volcano_ri_l3}
						{volcano_ri_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_ri_l1}
						{scatter_ri_l2}
						{scatter_ri_l3}
						{scatter_ri_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_ri_l1}
						{bar_ri_l2}
						{bar_ri_l3}
						{bar_ri_l4}
					"></iframe>
				</div>
			</div>
			<div id="MSE" class="plot">
				<h2>Multiple skipped exons (MSE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_mse_l1}
						{volcano_mse_l2}
						{volcano_mse_l3}
						{volcano_mse_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_mse_l1}
						{scatter_mse_l2}
						{scatter_mse_l3}
						{scatter_mse_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_mse_l1}
						{bar_mse_l2}
						{bar_mse_l3}
						{bar_mse_l4}
					"></iframe>
				</div>
			</div>
			<div id="AFE" class="plot">
				<h2>Alternative first exon (AFE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_afe_l1}
						{volcano_afe_l2}
						{volcano_afe_l3}
						{volcano_afe_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_afe_l1}
						{scatter_afe_l2}
						{scatter_afe_l3}
						{scatter_afe_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_afe_l1}
						{bar_afe_l2}
						{bar_afe_l3}
						{bar_afe_l4}
					"></iframe>
				</div>
			</div>
			<div id="ALE" class="plot">
				<h2>Alternative last exon (ALE)</h2>
				<div class="plot-container">
					<iframe srcdoc="
						{volcano_ale_l1}
						{volcano_ale_l2}
						{volcano_ale_l3}
						{volcano_ale_l4}
					"></iframe>
					<iframe srcdoc="
						{scatter_ale_l1}
						{scatter_ale_l2}
						{scatter_ale_l3}
						{scatter_ale_l4}
					"></iframe>
					<iframe srcdoc="
						{bar_ale_l1}
						{bar_ale_l2}
						{bar_ale_l3}
						{bar_ale_l4}
					"></iframe>
				</div>
			</div>
		</div>
		<footer>
			<p>© 2024 Naoto Kubota</p>
		</footer>
	</body>
	</html>
	'''.format(
		shiba_command = shiba_command,
		pca_tpm_l1 = lines_strip_pca_dict["TPM"][0],
		pca_tpm_l2 = lines_strip_pca_dict["TPM"][1],
		pca_tpm_l3 = lines_strip_pca_dict["TPM"][2],
		pca_tpm_l4 = lines_strip_pca_dict["TPM"][3],
		pca_psi_l1 = lines_strip_pca_dict["PSI"][0],
		pca_psi_l2 = lines_strip_pca_dict["PSI"][1],
		pca_psi_l3 = lines_strip_pca_dict["PSI"][2],
		pca_psi_l4 = lines_strip_pca_dict["PSI"][3],
		volcano_se_l1 = lines_strip_dict["SE"]["volcano"][0],
		volcano_se_l2 = lines_strip_dict["SE"]["volcano"][1],
		volcano_se_l3 = lines_strip_dict["SE"]["volcano"][2],
		volcano_se_l4 = lines_strip_dict["SE"]["volcano"][3],
		scatter_se_l1 = lines_strip_dict["SE"]["scatter"][0],
		scatter_se_l2 = lines_strip_dict["SE"]["scatter"][1],
		scatter_se_l3 = lines_strip_dict["SE"]["scatter"][2],
		scatter_se_l4 = lines_strip_dict["SE"]["scatter"][3],
		bar_se_l1 = lines_strip_dict["SE"]["bar"][0],
		bar_se_l2 = lines_strip_dict["SE"]["bar"][1],
		bar_se_l3 = lines_strip_dict["SE"]["bar"][2],
		bar_se_l4 = lines_strip_dict["SE"]["bar"][3],
		volcano_five_l1 = lines_strip_dict["FIVE"]["volcano"][0],
		volcano_five_l2 = lines_strip_dict["FIVE"]["volcano"][1],
		volcano_five_l3 = lines_strip_dict["FIVE"]["volcano"][2],
		volcano_five_l4 = lines_strip_dict["FIVE"]["volcano"][3],
		scatter_five_l1 = lines_strip_dict["FIVE"]["scatter"][0],
		scatter_five_l2 = lines_strip_dict["FIVE"]["scatter"][1],
		scatter_five_l3 = lines_strip_dict["FIVE"]["scatter"][2],
		scatter_five_l4 = lines_strip_dict["FIVE"]["scatter"][3],
		bar_five_l1 = lines_strip_dict["FIVE"]["bar"][0],
		bar_five_l2 = lines_strip_dict["FIVE"]["bar"][1],
		bar_five_l3 = lines_strip_dict["FIVE"]["bar"][2],
		bar_five_l4 = lines_strip_dict["FIVE"]["bar"][3],
		volcano_three_l1 = lines_strip_dict["THREE"]["volcano"][0],
		volcano_three_l2 = lines_strip_dict["THREE"]["volcano"][1],
		volcano_three_l3 = lines_strip_dict["THREE"]["volcano"][2],
		volcano_three_l4 = lines_strip_dict["THREE"]["volcano"][3],
		scatter_three_l1 = lines_strip_dict["THREE"]["scatter"][0],
		scatter_three_l2 = lines_strip_dict["THREE"]["scatter"][1],
		scatter_three_l3 = lines_strip_dict["THREE"]["scatter"][2],
		scatter_three_l4 = lines_strip_dict["THREE"]["scatter"][3],
		bar_three_l1 = lines_strip_dict["THREE"]["bar"][0],
		bar_three_l2 = lines_strip_dict["THREE"]["bar"][1],
		bar_three_l3 = lines_strip_dict["THREE"]["bar"][2],
		bar_three_l4 = lines_strip_dict["THREE"]["bar"][3],
		volcano_mxe_l1 = lines_strip_dict["MXE"]["volcano"][0],
		volcano_mxe_l2 = lines_strip_dict["MXE"]["volcano"][1],
		volcano_mxe_l3 = lines_strip_dict["MXE"]["volcano"][2],
		volcano_mxe_l4 = lines_strip_dict["MXE"]["volcano"][3],
		scatter_mxe_l1 = lines_strip_dict["MXE"]["scatter"][0],
		scatter_mxe_l2 = lines_strip_dict["MXE"]["scatter"][1],
		scatter_mxe_l3 = lines_strip_dict["MXE"]["scatter"][2],
		scatter_mxe_l4 = lines_strip_dict["MXE"]["scatter"][3],
		bar_mxe_l1 = lines_strip_dict["MXE"]["bar"][0],
		bar_mxe_l2 = lines_strip_dict["MXE"]["bar"][1],
		bar_mxe_l3 = lines_strip_dict["MXE"]["bar"][2],
		bar_mxe_l4 = lines_strip_dict["MXE"]["bar"][3],
		volcano_ri_l1 = lines_strip_dict["RI"]["volcano"][0],
		volcano_ri_l2 = lines_strip_dict["RI"]["volcano"][1],
		volcano_ri_l3 = lines_strip_dict["RI"]["volcano"][2],
		volcano_ri_l4 = lines_strip_dict["RI"]["volcano"][3],
		scatter_ri_l1 = lines_strip_dict["RI"]["scatter"][0],
		scatter_ri_l2 = lines_strip_dict["RI"]["scatter"][1],
		scatter_ri_l3 = lines_strip_dict["RI"]["scatter"][2],
		scatter_ri_l4 = lines_strip_dict["RI"]["scatter"][3],
		bar_ri_l1 = lines_strip_dict["RI"]["bar"][0],
		bar_ri_l2 = lines_strip_dict["RI"]["bar"][1],
		bar_ri_l3 = lines_strip_dict["RI"]["bar"][2],
		bar_ri_l4 = lines_strip_dict["RI"]["bar"][3],
		volcano_mse_l1 = lines_strip_dict["MSE"]["volcano"][0],
		volcano_mse_l2 = lines_strip_dict["MSE"]["volcano"][1],
		volcano_mse_l3 = lines_strip_dict["MSE"]["volcano"][2],
		volcano_mse_l4 = lines_strip_dict["MSE"]["volcano"][3],
		scatter_mse_l1 = lines_strip_dict["MSE"]["scatter"][0],
		scatter_mse_l2 = lines_strip_dict["MSE"]["scatter"][1],
		scatter_mse_l3 = lines_strip_dict["MSE"]["scatter"][2],
		scatter_mse_l4 = lines_strip_dict["MSE"]["scatter"][3],
		bar_mse_l1 = lines_strip_dict["MSE"]["bar"][0],
		bar_mse_l2 = lines_strip_dict["MSE"]["bar"][1],
		bar_mse_l3 = lines_strip_dict["MSE"]["bar"][2],
		bar_mse_l4 = lines_strip_dict["MSE"]["bar"][3],
		volcano_afe_l1 = lines_strip_dict["AFE"]["volcano"][0],
		volcano_afe_l2 = lines_strip_dict["AFE"]["volcano"][1],
		volcano_afe_l3 = lines_strip_dict["AFE"]["volcano"][2],
		volcano_afe_l4 = lines_strip_dict["AFE"]["volcano"][3],
		scatter_afe_l1 = lines_strip_dict["AFE"]["scatter"][0],
		scatter_afe_l2 = lines_strip_dict["AFE"]["scatter"][1],
		scatter_afe_l3 = lines_strip_dict["AFE"]["scatter"][2],
		scatter_afe_l4 = lines_strip_dict["AFE"]["scatter"][3],
		bar_afe_l1 = lines_strip_dict["AFE"]["bar"][0],
		bar_afe_l2 = lines_strip_dict["AFE"]["bar"][1],
		bar_afe_l3 = lines_strip_dict["AFE"]["bar"][2],
		bar_afe_l4 = lines_strip_dict["AFE"]["bar"][3],
		volcano_ale_l1 = lines_strip_dict["ALE"]["volcano"][0],
		volcano_ale_l2 = lines_strip_dict["ALE"]["volcano"][1],
		volcano_ale_l3 = lines_strip_dict["ALE"]["volcano"][2],
		volcano_ale_l4 = lines_strip_dict["ALE"]["volcano"][3],
		scatter_ale_l1 = lines_strip_dict["ALE"]["scatter"][0],
		scatter_ale_l2 = lines_strip_dict["ALE"]["scatter"][1],
		scatter_ale_l3 = lines_strip_dict["ALE"]["scatter"][2],
		scatter_ale_l4 = lines_strip_dict["ALE"]["scatter"][3],
		bar_ale_l1 = lines_strip_dict["ALE"]["bar"][0],
		bar_ale_l2 = lines_strip_dict["ALE"]["bar"][1],
		bar_ale_l3 = lines_strip_dict["ALE"]["bar"][2],
		bar_ale_l4 = lines_strip_dict["ALE"]["bar"][3]
	)
	with open(os.path.join(output_dir, "summary.html"), 'w') as f:
		f.write(summary_html)
	return 0

def load_splicing_summary_table(input_dir: str) -> pd.DataFrame:
	summary_df = pd.read_csv(os.path.join(input_dir, "splicing", "summary.txt"), sep = "\t")
	return summary_df

def barplot_splicing(summary_df: pd.DataFrame, fig_path: str):
	
	g = sns.catplot(
		data = summary_df,
		x = "Number",
		y = "AS",
		order = ["SE", "FIVE", "THREE", "MXE", "RI", "AFE", "ALE", "MSE"],
		hue = "Direction",
		hue_order = ["Up", "Down"],
		col = "Label",
		col_order = ["annotated", "unannotated"],
		palette = ["#d01c8b", "#4dac26"],
		kind = "bar",
		linewidth = 1,
		edgecolor = "black",
		height = 4, aspect = 0.8,
	)
	g.set_axis_labels("# DSEs", "")
	g.set_titles(col_template="{col_name}")
	sns.move_legend(g, "upper right", bbox_to_anchor=(1.05, 0.95), title = "Direction")
	
	# Set x-axis limit to 50 if max value is less than 40
	max_value = summary_df["Number"].max()
	if max_value < 40:
		for ax in g.axes.flatten():
			ax.set_xlim(0, 50)
	
	# Add bar labels for all facets and all containers
	for ax in g.axes.flatten():
		for c in ax.containers:
			# Only add labels for bars with non-zero values
			labels = [str(int(v)) if v != 0 else '' for v in c.datavalues]
			ax.bar_label(c, labels=labels, label_type='edge', padding=2)
	# Adjust spacing between subplots
	plt.subplots_adjust(wspace=0.3)  # Increase horizontal spacing between plots
	plt.savefig(fig_path, dpi = 400, bbox_inches = "tight")

def main():

	# Get arguments
	args = get_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Starting making plots....")

	# Set variables
	input_dir = args.input
	output_dir = args.output
	# Make directory
	os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)
	# Load experiment table
	logger.info("Loading experiment table....")
	experiment_table_df = load_experiment_table(args.experiment_table)
	# TPM
	logger.info("Making plots for PCA....")
	logger.debug("TPM....")
	pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2 = load_tpm_pca_table(input_dir, experiment_table_df, output_dir)
	plots_pca("TPM", pca_tpm_df, contribution_tpm_PC1, contribution_tpm_PC2, output_dir)
	# PSI
	logger.debug("PSI....")
	pca_psi_df, contribution_psi_PC1, contribution_psi_PC2 = load_psi_pca_table(input_dir, experiment_table_df, output_dir)
	plots_pca("PSI", pca_psi_df, contribution_psi_PC1, contribution_psi_PC2, output_dir)
	# Splicing events
	logger.info("Making plots for splicing events....")
	AS_list = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]
	for AS in AS_list:
		plots(AS, input_dir, output_dir)
	# Write summary html
	logger.info("Writing summary html....")
	write_summary_html(args.shiba_command, output_dir)
	# Splicing summary
	logger.info("Making barplot for splicing summary....")
	png_dir = os.path.join(output_dir, "png")
	os.makedirs(png_dir, exist_ok=True)
	fig_path = os.path.join(png_dir, "barplot_splicing_summary.png")
	barplot_splicing(load_splicing_summary_table(input_dir), fig_path)

	logger.info("Making plots completed!")

if __name__ == '__main__':
    main()
