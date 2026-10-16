import warnings
warnings.simplefilter('ignore')
import argparse
import sys
import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import html
import logging
from template_renderer import HTMLTemplateRenderer, get_splicing_event_config
from lib import figure_style as style
from pathlib import Path

# Configure logging
logger = logging.getLogger(__name__)

def get_args():

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "Make plots for alternative splicing events"
	)
	parser.add_argument("--single-cell", action="store_true", help="Read flat scShiba results without bulk PCA")
	parser.add_argument("--only-psi", action="store_true", help="Single-cell PSI-only report without differential tests")
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
	contribution_tpm_PC1 = str((contribution_tpm_df.iloc[0, 1]*100).round(2))
	contribution_tpm_PC2 = str((contribution_tpm_df.iloc[1, 1]*100).round(2))
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
	contribution_psi_PC1 = str((contribution_psi_df.iloc[0, 1]*100).round(2))
	contribution_psi_PC2 = str((contribution_psi_df.iloc[1, 1]*100).round(2))
	return pca_psi_df, contribution_psi_PC1, contribution_psi_PC2

def plots_pca(name, pca_df, contribution_PC1, contribution_PC2, output_dir):
    style.configure()
    fig, ax = plt.subplots(figsize=(4.8, 3.6), layout='constrained')
    interactive = go.Figure()
    symbols = ['circle', 'triangle-up', 'square', 'diamond', 'triangle-down', 'cross', 'x', 'triangle-left']
    markers = ['o', '^', 's', 'D', 'v', 'P', 'X', '<']
    for index, (group, frame) in enumerate(pca_df.groupby('group', sort=False)):
        if 'group_order' in frame:
            index = int(frame.group_order.iloc[0])
        interactive.add_trace(go.Scatter(x=frame.PC1, y=frame.PC2, mode='markers', name=str(group),
            text=frame['sample'] if 'sample' in frame else frame.index.astype(str),
            marker=dict(size=10, color=style.PALETTE[index % 8], symbol=symbols[index % 8]),
            hovertemplate='%{text}<br>PC1=%{x:.4g}<br>PC2=%{y:.4g}<extra>%{fullData.name}</extra>'))
        ax.scatter(frame.PC1, frame.PC2, s=18, color=style.PALETTE[index % 8],
                   marker=markers[index % 8], label=str(group), linewidths=.3, edgecolors='white')
    ax.set_xlabel(f'PC1 ({contribution_PC1}%)')
    ax.set_ylabel(f'PC2 ({contribution_PC2}%)')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    interactive.update_layout(xaxis_title=f'PC1 ({contribution_PC1}%)', yaxis_title=f'PC2 ({contribution_PC2}%)')
    content = style.export(fig, output_dir, f'pca_{name}', f'{name} PCA; each point represents a sample.', interactive=interactive)
    directory = Path(output_dir)/'data'; directory.mkdir(parents=True, exist_ok=True)
    (directory/f'pca_{name}.html').write_text(style.document(name+' PCA', [content]))


def plots(AS, input_dir, output_dir, single_cell=False):
    style.configure()
    frame = pd.read_csv((Path(input_dir) if single_cell else Path(input_dir)/'splicing')/f'PSI_{AS}.txt', sep='\t')
    for mode, x, y, xlabel, ylabel in [
        ('volcano', 'dPSI', 'q', 'Observed ΔPSI', '−log10(FDR-adjusted P)'),
        ('scatter', 'ref_PSI', 'alt_PSI', 'Reference PSI', 'Alternative PSI')]:
        fig, ax = plt.subplots(figsize=(4.5, 3.4), layout='constrained')
        interactive = go.Figure()
        if len(frame):
            valid = np.isfinite(frame[x]) & np.isfinite(frame[y])
            if mode == 'volcano': valid &= frame[y].between(0, 1)
            data = frame.loc[valid]
            significant = data['Diff events'].eq('Yes')
            for mask, label, color, marker in [(~significant, 'Not called', style.PALETTE[0], 'o'),
                    (significant & data.dPSI.ge(0), 'Increased PSI', style.PALETTE[2], '^'),
                    (significant & data.dPSI.lt(0), 'Decreased PSI', style.PALETTE[1], 'v')]:
                subset = data.loc[mask]
                if len(subset):
                    values = -np.log10(subset[y].clip(lower=np.finfo(float).tiny)) if mode == 'volcano' else subset[y]
                    custom = subset.reindex(columns=['event_id', 'gene_name', 'q']).fillna('').to_numpy()
                    interactive.add_trace(go.Scattergl(x=subset[x], y=values, mode='markers', name=label, customdata=custom,
                        marker=dict(size=9, opacity=.65, color=color, symbol={'o':'circle','^':'triangle-up','v':'triangle-down'}[marker]),
                        hovertemplate='%{customdata[0]} · %{customdata[1]}<br>x=%{x:.4g}<br>y=%{y:.4g}<br>q=%{customdata[2]}<extra>%{fullData.name}</extra>'))
                    ax.scatter(subset[x], values, color=color, marker=marker, s=10, linewidths=0, alpha=.65, label=label, rasterized=True)
            if len(data): ax.legend()
        else: ax.text(.5, .5, 'No event results', ha='center', transform=ax.transAxes)
        if mode == 'scatter':
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.plot([0, 1], [0, 1], color='#aaaaaa', lw=.7, zorder=0)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        interactive.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)
        if mode == 'scatter':
            interactive.update_xaxes(range=[0,1]); interactive.update_yaxes(range=[0,1])
            interactive.add_shape(type='line', x0=0, x1=1, y0=0, y1=1, line=dict(color='#aaaaaa',width=1))
        if not interactive.data:
            interactive.add_annotation(text='No finite event results', x=.5, y=.5, xref='paper', yref='paper', showarrow=False)
        content = style.export(fig, output_dir, f'{mode}_{AS}', f'{AS}: existing differential-event calls; missing coordinates omitted. Zero q-values clipped to machine minimum for log display.', interactive=interactive)
        directory = Path(output_dir)/'data'; directory.mkdir(parents=True, exist_ok=True)
        (directory/f'{mode}_{AS}.html').write_text(style.document(AS, [content]))

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

def write_summary_html(shiba_command, input_dir, output_dir, single_cell=False):
    sections = [style.command_section(input_dir, shiba_command),
                '<p>Figures use existing Fisher results without recalculating tests. Increased/decreased PSI refers to alternative minus reference. '
                'Missing coordinates are omitted; zero q-values are clipped to machine minimum for logarithmic display. ' +
                ('PSI is aggregated by cell group; cells are not independent biological replicates.</p>' if single_cell else
                 'PCA includes the samples in the supplied PCA files, which may contain groups beyond the pairwise comparison.</p>'),
                '<h2>Differential splicing</h2>', load_plot_content(output_dir, 'barplot_splicing_summary.html')]
    for name in (() if single_cell else ('TPM', 'PSI')):
        sections += ['<h2>'+name+' PCA</h2>', load_plot_content(output_dir, f'pca_{name}.html')]
    for event in get_splicing_event_config():
        code = event['code']
        sections += ['<h2>'+html.escape(event['title'])+'</h2>',
                     load_plot_content(output_dir, f'volcano_{code}.html'),
                     load_plot_content(output_dir, f'scatter_{code}.html')]
    (Path(output_dir)/'summary.html').write_text(style.document('scShiba Fisher analysis' if single_cell else 'Shiba Fisher analysis', sections), encoding='utf-8')
    return 0

def write_single_cell_report(input_dir, output_dir, command=None, only_psi=False):
    """Render aggregated single-cell results using the shared Fisher figures."""
    source, output = Path(input_dir), Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    if only_psi:
        sections = [style.command_section(source, command),
                    '<p>PSI-only analysis: differential testing was not performed.</p>']
        for event in get_splicing_event_config():
            frame = pd.read_csv(source / f"PSI_{event['code']}.txt", sep='\t')
            sections += ['<h2>'+html.escape(event['title'])+'</h2>',
                         f'<p>{len(frame)} events; first 200 rows shown.</p>',
                         frame.head(200).to_html(index=False, escape=True)]
        (output/'summary.html').write_text(style.document('scShiba PSI analysis', sections), encoding='utf-8')
        return
    summary = pd.read_csv(source/'summary.txt', sep='\t')
    barplot_splicing(summary, str(output/'png/barplot_splicing_summary.png'),
                    str(output/'pdf/barplot_splicing_summary.pdf'))
    for event in get_splicing_event_config():
        plots(event['code'], source, output, single_cell=True)
    write_summary_html(command, source, output, single_cell=True)


def load_plot_content(output_dir, filename):
    interactive = Path(output_dir)/'interactive'/(Path(filename).stem+'.html')
    if interactive.exists(): return interactive.read_text(encoding='utf-8')
    path = Path(output_dir)/'png'/(Path(filename).stem+'.png')
    if path.exists(): return style.image_tag(path, Path(filename).stem.replace('_', ' '))
    return '<p>Plot not available.</p>'

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

def barplot_splicing(summary_df, png_path, pdf_path):
    style.configure()
    types = ['SE', 'FIVE', 'THREE', 'MXE', 'RI', 'AFE', 'ALE', 'MSE']
    from plotly.subplots import make_subplots
    interactive = make_subplots(rows=1, cols=2, shared_xaxes=True, subplot_titles=['Annotated', 'Unannotated'])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharex=True, layout='constrained')
    for index, (ax, label) in enumerate(zip(axes, ['annotated', 'unannotated'])):
        data = summary_df[summary_df.Label.eq(label)]
        for direction, offset, color in [('Up', -.18, style.PALETTE[2]), ('Down', .18, style.PALETTE[1])]:
            values = data[data.Direction.eq(direction)].groupby('AS').Number.sum().reindex(types, fill_value=0)
            interactive.add_trace(go.Bar(x=values.to_list(), y=types, orientation='h', name=direction, legendgroup=direction, showlegend=index==0, marker_color=color), row=1, col=index+1)
            ax.barh(np.arange(len(types))+offset, values, height=.33, color=color, label=direction)
        ax.set_yticks(np.arange(len(types)), types); ax.invert_yaxis()
        ax.set_xlabel('Differentially spliced events')
        ax.set_title(label.capitalize(), fontsize=8)
        ax.text(-.12, 1.08, chr(65+index), transform=ax.transAxes, weight='bold', fontsize=10)
    from matplotlib.ticker import MaxNLocator
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axes[1].legend()
    interactive.update_layout(barmode='group')
    maximum = float(summary_df.Number.max()) if len(summary_df) else 0
    interactive.update_xaxes(title_text='Differentially spliced events', range=[0,max(1,maximum*1.1)])
    interactive.update_yaxes(autorange='reversed')
    output = Path(png_path).parent.parent
    style.export(fig, output, 'barplot_splicing_summary', 'Counts of existing significant events, by annotation and direction.', interactive=interactive)

def main():

	# Get arguments
	args = get_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Starting making plots....")
	if args.single_cell:
		write_single_cell_report(args.input, args.output, args.shiba_command, args.only_psi)
		return
	manifest_path = os.path.join(args.input, "splicing", "analysis.json")
	if os.path.isfile(manifest_path):
		import json
		with open(manifest_path) as handle:
			regression = json.load(handle).get("stat_method") in ("beta_binomial", "beta-binomial")
		if regression:
			from lib.factorial_report import write_report
			write_report(args.input, args.output, args.experiment_table, args.shiba_command)
			return

	# Set variables
	input_dir = args.input
	output_dir = args.output
	# Make directory
	os.makedirs(os.path.join(output_dir, "data"), exist_ok=True)
	
	# Splicing summary
	logger.info("Making barplot for splicing summary....")
	png_dir = os.path.join(output_dir, "png")
	pdf_dir = os.path.join(output_dir, "pdf")
	os.makedirs(png_dir, exist_ok=True)
	os.makedirs(pdf_dir, exist_ok=True)
	png_path = os.path.join(png_dir, "barplot_splicing_summary.png")
	pdf_path = os.path.join(pdf_dir, "barplot_splicing_summary.pdf")
	barplot_splicing(load_splicing_summary_table(input_dir), png_path, pdf_path)

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
