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
		# Remove fixed width to allow responsive sizing
		height = 500,  # Increased from 400 to 500 for taller plots
		font ={
			"family": "Arial",
			"size": 18
		},
		xaxis_title = "PC1 ({}%)".format(contribution_PC1),
		yaxis_title = "PC2 ({}%)".format(contribution_PC2),
		legend_title = "Group",
		# Add responsive layout settings
		autosize = True,
		margin = dict(l=60, r=60, t=80, b=60)
	)
	# Use responsive configuration when writing HTML
	fig.write_html(
		os.path.join(output_dir, "data", "pca_" + name + ".html"), 
		include_plotlyjs = "cdn",
		config = {
			'responsive': True,
			'displayModeBar': True,
			'displaylogo': False,
			'modeBarButtonsToRemove': ['pan2d', 'lasso2d']
		}
	)

def plots(AS: str, input_dir: str, output_dir: str):

	# load data
	df = pd.read_csv(
		os.path.join(input_dir, "splicing", "PSI_" + AS + ".txt"),
		sep = "\t"
	)
	if not df.empty:
		# Round dPSI and others
		df["dPSI"] = df["dPSI"].round(4)
		df['-log10(q)'] = -np.log10(df["q"])
		df['-log10(q)'] = df['-log10(q)'].round(4)
		df["ref_PSI"] = df["ref_PSI"].round(4)
		df["alt_PSI"] = df["alt_PSI"].round(4)
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
			color_discrete_sequence = ["#d01c8b", "#4dac26", "lightgrey"],
			hover_data = ["gene_name", "event_id", "ref_PSI", "alt_PSI", "q"]
		)
		fig.update_traces(
			marker = dict(
				size = 8,
				line = dict(width = 0, color = 'DarkSlateGrey')),
				selector = dict(mode = 'markers')
		)
		# Remove the title from the volcano plot
		fig.update_layout(title=None)
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
		# Remove fixed width for responsive sizing
		height = 400,
		font ={
			"family": "Arial",
			"size": 16
		},
		legend_title = "Group, Label",
		# Add responsive layout settings
		autosize = True,
		margin = dict(l=60, r=60, t=80, b=60)
	)
	# Use responsive configuration when writing HTML
	fig.write_html(
		os.path.join(output_dir, "data", "volcano_" + AS + ".html"), 
		include_plotlyjs = "cdn",
		config = {
			'responsive': True,
			'displayModeBar': True,
			'displaylogo': False,
			'modeBarButtonsToRemove': ['pan2d', 'lasso2d']
		}
	)

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
			color_discrete_sequence = ["#d01c8b", "#4dac26", "lightgrey"],
			hover_data = ["gene_name", "event_id", "dPSI", "q"]
		)
		fig.update_traces(
			marker = dict(size=8,
							line=dict(width=0, color='DarkSlateGrey')),
							selector=dict(mode='markers')
		)
		fig.update_layout(title=None)
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
		# Remove fixed width for responsive sizing
		height = 400,
		font ={
			"family": "Arial",
			"size": 16
		},
		xaxis_title = "PSI (Reference)",
		yaxis_title = "PSI (Alternative)",
		legend_title = "Group, Label",
		# Add responsive layout settings
		autosize = True,
		margin = dict(l=60, r=60, t=80, b=60)
	)
	# Use responsive configuration when writing HTML
	fig.write_html(
		os.path.join(output_dir, "data", "scatter_" + AS + ".html"), 
		include_plotlyjs = "cdn",
		config = {
			'responsive': True,
			'displayModeBar': True,
			'displaylogo': False,
			'modeBarButtonsToRemove': ['pan2d', 'lasso2d']
		}
	)

def calculate_event_count(input_dir: str, AS: str) -> int:
	"""Calculate the number of differential splicing events for a given AS type."""
	try:
		df = pd.read_csv(
			os.path.join(input_dir, "splicing", "PSI_" + AS + ".txt"),
			sep = "\t"
		)
		if df.empty:
			return 0
		
		# Count all events
		all_events_count = len(df)
		return all_events_count

	except FileNotFoundError:
		logger.warning(f"PSI file not found for {AS}")
		return 0
	except Exception as e:
		logger.warning(f"Error calculating event count for {AS}: {e}")
		return 0

def get_shiba_version() -> str:
	"""Get Shiba version from VERSION file."""
	try:
		# Get the directory where this script is located
		script_dir = os.path.dirname(os.path.abspath(__file__))
		# Go up one level to the main Shiba directory
		version_path = os.path.join(os.path.dirname(script_dir), "VERSION")
		
		with open(version_path, 'r', encoding='utf-8') as f:
			version = f.read().strip()
			return version
	except FileNotFoundError:
		logger.warning("VERSION file not found")
		return "unknown"
	except Exception as e:
		logger.warning(f"Error reading VERSION file: {e}")
		return "unknown"

def load_splicing_summary_image(output_dir: str) -> str:
	"""Load the splicing summary PNG and convert to base64 for embedding."""
	import base64
	
	png_path = os.path.join(output_dir, "png", "barplot_splicing_summary.png")
	
	try:
		with open(png_path, 'rb') as img_file:
			img_data = img_file.read()
			img_base64 = base64.b64encode(img_data).decode('utf-8')
			return f'<img src="data:image/png;base64,{img_base64}" alt="Differential Splicing Events Summary" class="splicing-summary-image">'
	except FileNotFoundError:
		logger.warning(f"Splicing summary image not found: {png_path}")
		return '<p style="text-align: center; color: #64748b;">Summary chart not available</p>'
	except Exception as e:
		logger.warning(f"Error loading splicing summary image: {e}")
		return '<p style="text-align: center; color: #64748b;">Summary chart not available</p>'

def write_summary_html(shiba_command: str, input_dir: str, output_dir: str):
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
		event_count = calculate_event_count(input_dir, config['code'])
		event_data = {
			'id': config['id'],
			'icon': config['icon'],
			'title': config['title'],
			'description': config['description'],
			'event_count': event_count,
			'volcano_content': load_plot_content(output_dir, f"volcano_{config['code']}.html"),
			'scatter_content': load_plot_content(output_dir, f"scatter_{config['code']}.html")
		}
		splicing_events.append(event_data)
		
		# Generate standalone individual HTML file for this event (main file)
		standalone_individual_html = renderer.render_individual_event_standalone_html(event_data)
		individual_filename = f"{config['id']}.html"
		individual_filepath = os.path.join(output_dir, "data", individual_filename)
		with open(individual_filepath, 'w', encoding='utf-8') as f:
			f.write(standalone_individual_html)
		logger.info(f"Generated standalone individual HTML: data/{individual_filename}")
	
	# Load splicing summary image
	splicing_summary_content = load_splicing_summary_image(output_dir)
	
	# Get Shiba version
	shiba_version = get_shiba_version()
	
	# Prepare template data for index page
	template_data = {
		'shiba_command': shiba_command,
		'shiba_version': shiba_version,
		'pca_tpm_content': pca_tpm_content,
		'pca_psi_content': pca_psi_content,
		'splicing_summary_content': splicing_summary_content,
		'splicing_events': splicing_events
	}
	
	# Generate standalone index/overview HTML (main entry point)
	standalone_index_content = renderer.render_index_standalone_html(template_data)
	with open(os.path.join(output_dir, "summary.html"), 'w', encoding='utf-8') as f:
		f.write(standalone_index_content)
	
	logger.info("HTML files generated successfully!")
	logger.info("  - summary.html (main overview page with links to individual event files)")
	logger.info("  - Individual event files in data/ (data/se.html, data/five.html, etc.) - completely self-contained")
	logger.info("Note: summary.html requires the data/ directory with individual event files to function properly.")
	return 0

def load_plot_content(output_dir: str, filename: str) -> str:
	"""Load plot content from an HTML file for direct embedding."""
	file_path = os.path.join(output_dir, "data", filename)
	
	try:
		with open(file_path, 'r', encoding='utf-8') as f:
			content = f.read()
		
		# For direct embedding, we need to extract the plot div and script
		# but keep them as interactive Plotly content
		import re
		
		# Look for the main plot div with plotly-graph-div class
		div_pattern = r'<div[^>]*class="plotly-graph-div"[^>]*>.*?</div>'
		div_match = re.search(div_pattern, content, re.DOTALL)
		
		# Look for the script that initializes the plot
		script_pattern = r'<script type="text/javascript">.*?Plotly\.newPlot.*?</script>'
		script_match = re.search(script_pattern, content, re.DOTALL)
		
		if div_match and script_match:
			# Return the div and script together for direct embedding
			return f"{div_match.group()}\n{script_match.group()}"
		
		# Fallback: try to find any div and script tags
		div_fallback = re.search(r'<div[^>]*id="[^"]*"[^>]*>.*?</div>', content, re.DOTALL)
		script_fallback = re.search(r'<script[^>]*>.*?Plotly.*?</script>', content, re.DOTALL)
		
		if div_fallback and script_fallback:
			return f"{div_fallback.group()}\n{script_fallback.group()}"
		
		# If regex fails, try a simpler line-by-line approach
		lines = content.split('\n')
		plot_lines = []
		in_body = False
		
		for line in lines:
			# Skip head content, only get body content
			if '<body>' in line:
				in_body = True
				continue
			if '</body>' in line:
				break
			if in_body and line.strip():
				# Skip doctype, html, head tags
				if not any(tag in line.lower() for tag in ['<!doctype', '<html', '<head', '</head', '<meta', '<title', '<link']):
					plot_lines.append(line)
		
		if plot_lines:
			return '\n'.join(plot_lines)
		
		# Final fallback: return full content
		return content
		
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
	g.set_axis_labels("# Differentially spliced events", "")
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
	
	# Splicing summary
	logger.info("Making barplot for splicing summary....")
	png_dir = os.path.join(output_dir, "png")
	os.makedirs(png_dir, exist_ok=True)
	fig_path = os.path.join(png_dir, "barplot_splicing_summary.png")
	barplot_splicing(load_splicing_summary_table(input_dir), fig_path)

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
	write_summary_html(args.shiba_command, input_dir, output_dir)

	logger.info("Making plots completed!")

if __name__ == '__main__':
    main()
