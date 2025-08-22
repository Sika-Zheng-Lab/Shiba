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
from template_renderer import HTMLTemplateRenderer, get_splicing_event_config

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
	"""Write summary HTML using modern template system with individual event files."""
	
	# Initialize template renderer
	template_dir = os.path.join(os.path.dirname(__file__), "templates")
	renderer = HTMLTemplateRenderer(template_dir)
	
	# Copy static files (CSS, JS) to output directory
	renderer.copy_static_files(output_dir)
	
	# Load PCA content
	pca_tpm_content = load_plot_content(output_dir, "pca_TPM.html")
	pca_psi_content = load_plot_content(output_dir, "pca_PSI.html")
	
	# Prepare splicing events data
	splicing_events = []
	event_configs = get_splicing_event_config()
	
	for config in event_configs:
		event_data = {
			'id': config['id'],
			'icon': config['icon'],
			'title': config['title'],
			'description': config['description'],
			'volcano_content': load_plot_content(output_dir, f"volcano_{config['code']}.html"),
			'scatter_content': load_plot_content(output_dir, f"scatter_{config['code']}.html"),
			'bar_content': load_plot_content(output_dir, f"bar_{config['code']}.html")
		}
		splicing_events.append(event_data)
		
		# Generate standalone individual HTML file for this event (main file)
		standalone_individual_html = renderer.render_individual_event_standalone_html(event_data)
		individual_filename = f"{config['id']}.html"
		with open(os.path.join(output_dir, individual_filename), 'w', encoding='utf-8') as f:
			f.write(standalone_individual_html)
		logger.info(f"Generated standalone individual HTML: {individual_filename}")
	
	# Prepare template data for index page
	template_data = {
		'shiba_command': shiba_command,
		'pca_tpm_content': pca_tpm_content,
		'pca_psi_content': pca_psi_content,
		'splicing_events': splicing_events
	}
	
	# Generate standalone index/overview HTML (main entry point)
	standalone_index_content = renderer.render_index_standalone_html(template_data)
	with open(os.path.join(output_dir, "index.html"), 'w', encoding='utf-8') as f:
		f.write(standalone_index_content)
	
	# Generate standalone summary HTML (traditional single-page version)
	standalone_html_content = renderer.render_standalone_summary_html(template_data)
	with open(os.path.join(output_dir, "summary.html"), 'w', encoding='utf-8') as f:
		f.write(standalone_html_content)
	
	logger.info("HTML files generated successfully!")
	logger.info("  - index.html (main overview with embedded CSS/JS)")
	logger.info("  - Individual event files (se.html, five.html, etc.) - completely self-contained")
	logger.info("  - summary.html (traditional single-page version with embedded CSS/JS)")
	logger.info("All files are standalone and can be moved anywhere!")
	return 0

def load_plot_content(output_dir: str, filename: str) -> str:
	"""Load plot content from an HTML file."""
	file_path = os.path.join(output_dir, "data", filename)
	
	try:
		with open(file_path, 'r', encoding='utf-8') as f:
			content = f.read()
			
		# For PCA files, return the full content for direct embedding
		if filename.startswith('pca_'):
			return content
		
		# For other plot files, extract just the essential lines for iframe embedding
		lines = content.split('\n')
		essential_lines = [html.escape(line.strip()) for line in lines[2:6]]
		return '\n'.join(essential_lines)
		
	except FileNotFoundError:
		logger.warning(f"Plot file not found: {filename}")
		return create_empty_plot_content("Plot not available")
	except Exception as e:
		logger.warning(f"Error loading plot file {filename}: {e}")
		return create_empty_plot_content("Plot not available")

def create_empty_plot_content(message: str) -> str:
	"""Create content for empty/missing plots."""
	return f'''
		<div style="display: flex; align-items: center; justify-content: center; height: 100%; 
		           font-family: Arial, sans-serif; color: #64748b; text-align: center;">
			<div>
				<i class="fas fa-chart-line" style="font-size: 3rem; margin-bottom: 1rem; opacity: 0.5;"></i>
				<p style="font-size: 1.1rem; margin: 0;">{html.escape(message)}</p>
			</div>
		</div>
	'''

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
