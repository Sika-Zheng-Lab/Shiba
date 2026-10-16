"""Self-contained, contrast-aware report with publication-style interactive figures."""
import html
import json
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from . import figure_style as style


def write_report(input_dir, output_dir, experiment_table=None, shiba_command=None):
    fonts = style.configure()
    import matplotlib.pyplot as plt
    source = Path(input_dir) / 'splicing'
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source / 'analysis.json').read_text())
    results = pd.read_csv(source / 'event_statistics.tsv', sep='\t')
    prediction_path = source / 'model_predictions.tsv'
    predictions = pd.read_csv(prediction_path, sep='\t') if prediction_path.exists() else pd.DataFrame()
    sections = [style.command_section(input_dir, shiba_command), '<p>Formula: <code>' + html.escape(manifest['formula']) + '</code></p>',
                '<p>Tests concern logit-scale contrasts. FDR correction is performed separately per contrast across coverage-passing events. Filtered events are excluded; failed fits are not discoveries. Model PSI is an equal-weight component mean, distinct from observed Shiba PSI.</p>']
    settings = {label: manifest[key] for key, label in [
        ('n_samples', 'Samples'), ('n_events', 'Input events'),
        ('n_events_passed', 'Coverage-passing events'), ('fdr', 'FDR threshold'),
        ('p_adjust', 'Adjustment'), ('minimum_reads', 'Minimum reads'),
        ('min_sample_fraction', 'Minimum sample fraction'),
        ('filter_group_columns', 'Coverage grouping'), ('min_effect', 'Minimum PSI effect')]
        if key in manifest}
    if settings:
        sections.append(pd.DataFrame({'Setting': list(settings), 'Value': [str(v) for v in settings.values()]}).to_html(index=False, escape=True))
    contrasts = list(results.contrast.drop_duplicates())
    types = list(results.event_type.drop_duplicates())
    fig, ax = plt.subplots(figsize=(7.2, max(2.6, len(contrasts)*.5)), layout='constrained')
    counts = results[results['Diff events'].eq('Yes')].groupby(['contrast', 'event_type']).size()
    interactive = go.Figure()
    left = np.zeros(len(contrasts))
    for i, kind in enumerate(types):
        values = np.array([counts.get((c, kind), 0) for c in contrasts])
        interactive.add_trace(go.Bar(x=values, y=contrasts, orientation='h', name=kind, marker_color=style.PALETTE[i % 8]))
        ax.barh(np.arange(len(contrasts)), values, left=left, height=.55, label=kind, color=style.PALETTE[i % 8])
        left += values
    ax.set_yticks(np.arange(len(contrasts)), [str(c) for c in contrasts])
    ax.set_xlabel('Significant events'); ax.invert_yaxis()
    if types: ax.legend(ncol=4, loc='upper center', bbox_to_anchor=(.5, 1.3))
    if not left.any(): ax.text(.5, .5, 'No significant events', transform=ax.transAxes, ha='center')
    interactive.update_layout(barmode='stack', xaxis_title='Significant events', yaxis=dict(autorange='reversed'))
    sections += ['<h2>Discovery summary</h2>', style.export(fig, output, 'barplot_splicing_summary', 'Counts use the existing differential-event calls; no additional tests are performed.', interactive=interactive)]
    count_table = pd.DataFrame({kind: [counts.get((c, kind), 0) for c in contrasts] for kind in types}, index=contrasts)
    count_table.index.name = 'Contrast'
    sections.append(count_table.to_html(escape=True))
    if results.empty: sections.append('<p>No event results available.</p>')
    for index, (contrast, table) in enumerate(results.groupby('contrast', sort=False)):
        sections.append('<h2>' + html.escape(str(contrast)) + '</h2>')
        status = table.groupby('status', dropna=False).size().rename('Events').reset_index()
        sections.append(status.to_html(index=False, escape=True))
        valid = np.isfinite(table.conservative_component_effect) & np.isfinite(table.q_event) & table.q_event.between(0, 1)
        plotted = table.loc[valid]
        interactive = go.Figure()
        fig, ax = plt.subplots(figsize=(4.8, 3.4), layout='constrained')
        for significant, label, color, marker in [(False, 'Not called', style.PALETTE[0], 'o'), (True, 'Significant', style.PALETTE[2], '^')]:
            part = plotted[plotted['Diff events'].eq('Yes') == significant]
            if len(part):
                interactive.add_trace(go.Scattergl(x=part.conservative_component_effect, y=-np.log10(part.q_event.clip(lower=np.finfo(float).tiny)),
                    mode='markers', name=label, customdata=part[['event_id','event_type','q_event']].to_numpy(),
                    marker=dict(size=9, opacity=.65, color=color, symbol='triangle-up' if significant else 'circle'),
                    hovertemplate='%{customdata[0]} (%{customdata[1]})<br>Effect=%{x:.4g}<br>q=%{customdata[2]}<extra>%{fullData.name}</extra>'))
                ax.scatter(part.conservative_component_effect, -np.log10(part.q_event.clip(lower=np.finfo(float).tiny)), s=10, alpha=.65, linewidths=0, color=color, marker=marker, label=label, rasterized=True)
        ax.set_xlabel('Conservative component effect (logit contrast)')
        ax.set_ylabel('−log10(FDR-adjusted P)')
        if len(plotted): ax.legend()
        else: ax.text(.5, .5, 'No finite test results', ha='center', transform=ax.transAxes)
        interactive.update_layout(xaxis_title='Conservative component effect (logit contrast)', yaxis_title='−log10(FDR-adjusted P)')
        if not len(plotted): interactive.add_annotation(text='No finite test results', x=.5, y=.5, xref='paper', yref='paper', showarrow=False)
        sections.append(style.export(fig, output, f'contrast_{index+1}', f'{contrast}: {len(plotted):,} plotted; {len(table)-len(plotted):,} rows with missing/invalid coordinates omitted. Zero q-values are clipped to machine minimum only for the logarithm.', interactive=interactive))
    if len(predictions) and {'kind', 'name', 'value', 'event_id', 'event_type'}.issubset(predictions):
        effects = predictions[predictions.kind.eq('effect') & np.isfinite(predictions.value)]
        sections.append('<h2>Model PSI effects</h2>')
        for index, (name, table) in enumerate(effects.groupby('name', sort=False)):
            chosen = table.assign(magnitude=table.value.abs()).sort_values(['magnitude','event_type','event_id'], ascending=[False,True,True]).head(20).sort_values('value')
            fig, ax = plt.subplots(figsize=(6.4, max(2.5, .19*len(chosen)+.8)), layout='constrained')
            interactive = go.Figure(go.Scatter(x=chosen.value, y=chosen.event_type.astype(str)+': '+chosen.event_id.astype(str), mode='markers', marker=dict(color=style.PALETTE[1], size=10)))
            interactive.update_layout(xaxis_title=str(name)+' (PSI difference)')
            ax.scatter(chosen.value, np.arange(len(chosen)), s=14, color=style.PALETTE[1])
            ax.axvline(0, color='#999999', lw=.7)
            ax.set_yticks(np.arange(len(chosen)), chosen.event_type.astype(str)+': '+chosen.event_id.astype(str))
            ax.set_xlabel(str(name) + ' (PSI difference)')
            sections.append(style.export(fig, output, f'psi_effect_{index+1}', 'Largest absolute model PSI effects: at most 20 events, selected by magnitude without a significance filter. Point estimates only; no confidence intervals were computed.', interactive=interactive))
    else: sections.append('<p>No model predictions available.</p>')
    if experiment_table:
        import plots
        metadata = plots.load_experiment_table(experiment_table)
        for name, loader, file in [('TPM', plots.load_tpm_pca_table, 'tpm_pca.tsv'), ('PSI', plots.load_psi_pca_table, 'psi_pca.tsv')]:
            if (Path(input_dir)/'pca'/file).exists():
                frame, pc1, pc2 = loader(input_dir, metadata, output_dir)
                plots.plots_pca(name, frame, pc1, pc2, output_dir)
                sections += ['<h2>'+name+' PCA</h2>', plots.load_plot_content(output_dir, f'pca_{name}.html')]
    for title, table in [('Event results', results), ('Model predictions', predictions)]:
        sections += ['<h2>'+title+'</h2><p>First 200 rows embedded for inspection; complete TSVs remain in the analysis output directory.</p>', '<div class="table">'+table.head(200).to_html(index=False, escape=True)+'</div>']
    (output/'figure_style.json').write_text(json.dumps({'font_files':fonts, 'dpi':300, 'formats':['png','pdf','svg']}, indent=2))
    (output/'summary.html').write_text(style.document('Shiba beta-binomial regression', sections), encoding='utf-8')
