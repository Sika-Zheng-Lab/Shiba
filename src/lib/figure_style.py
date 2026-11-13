"""Shared publication figure style; Arial preferred, with default-font fallback."""
import base64
import html
import json
import os
import uuid
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager, ft2font

PALETTE = ['#737373', '#487C91', '#B66A55', '#7E8D62', '#8C789B', '#B59B63', '#608F88', '#A67F8D']


def configure():
    font_directory = os.environ.get('SHIBA_FONT_DIR')
    if font_directory:
        for path in font_manager.findSystemFonts(fontpaths=[font_directory]):
            font_manager.fontManager.addfont(path)
    paths = {}
    family = 'Arial'
    try:
        for weight in ('normal', 'bold'):
            paths[weight] = font_manager.findfont(font_manager.FontProperties(family='Arial', weight=weight), fallback_to_default=False)
    except ValueError:
        family = 'DejaVu Sans'
        paths = {weight: font_manager.findfont(font_manager.FontProperties(family=family, weight=weight))
                 for weight in ('normal', 'bold')}
    matplotlib.rcParams.update({
        'font.family': family, 'font.size': 8, 'axes.labelsize': 8,
        'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.linewidth': .7, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': False, 'xtick.direction': 'out', 'ytick.direction': 'out',
        'xtick.major.width': .7, 'ytick.major.width': .7,
        'xtick.major.size': 3, 'ytick.major.size': 3, 'lines.linewidth': 1.1,
        'lines.markersize': 4, 'legend.frameon': False, 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'savefig.dpi': 300,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'axes.prop_cycle': matplotlib.cycler(color=PALETTE),
    })
    return paths


def image_tag(path, caption):
    data = base64.b64encode(Path(path).read_bytes()).decode('ascii')
    return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption, quote=True)}"><figcaption>{html.escape(caption)}</figcaption></figure>'


def export(fig, output, name, caption, interactive=None):
    import matplotlib.pyplot as plt
    output = Path(output)
    fig.canvas.draw()
    for extension in ('png', 'pdf', 'svg'):
        directory = output / extension
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(directory / (name + '.' + extension), bbox_inches='tight', dpi=300)
    plt.close(fig)
    if interactive is not None:
        content = plot_html(interactive, caption)
        directory = output / 'interactive'; directory.mkdir(parents=True, exist_ok=True)
        (directory / (name + '.html')).write_text(content, encoding='utf-8')
        return content
    return image_tag(output / 'png' / (name + '.png'), caption)


def document(title, sections):
    from plotly.offline import get_plotlyjs
    script = '<script>' + get_plotlyjs() + '</script>'
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>' + html.escape(title) + '</title><style>'
            'body{font:17px Arial,sans-serif;color:#222;max-width:1100px;margin:32px auto;padding:0 20px;line-height:1.5}'
            'h1{font-size:28px}h2{font-size:23px;margin-top:32px}figure{margin:20px 0}'
            'img{max-width:100%;height:auto}figcaption{font-size:15px;color:#555}'
            'table{border-collapse:collapse;font-size:14px}th,td{padding:5px 8px;border-bottom:1px solid #ddd;text-align:left}'
            '.table{overflow:auto;max-height:480px}code{overflow-wrap:anywhere}pre.command{white-space:pre-wrap;overflow-wrap:anywhere;padding:14px;border:1px solid #ddd;font-size:15px}'
            '</style>' + script + '</head><body><h1>' + html.escape(title) + '</h1>' + ''.join(sections) + INTERACTIVE_LOADER + '</body></html>')


def plot_html(fig, caption):
    """Shared browser theme; CSS font stack falls back on the viewer's machine."""
    fig.update_layout(template='none', font=dict(family='Arial, sans-serif', size=16, color='#222'),
                      paper_bgcolor='white', plot_bgcolor='white', height=520,
                      margin=dict(l=125, r=45, t=60, b=90),
                      legend=dict(font=dict(size=14), bgcolor='rgba(0,0,0,0)'),
                      hoverlabel=dict(font=dict(family='Arial, sans-serif', size=15)))
    fig.update_xaxes(showgrid=False, showline=True, linewidth=1, linecolor='#444',
                     ticks='outside', tickwidth=1, ticklen=4, zeroline=False,
                     title_font=dict(size=18), tickfont=dict(size=14))
    fig.update_yaxes(showgrid=False, showline=True, linewidth=1, linecolor='#444',
                     ticks='outside', tickwidth=1, ticklen=4, zeroline=False,
                     title_font=dict(size=18), tickfont=dict(size=14))
    fig.update_annotations(font=dict(size=16))
    identifier = 'plot-' + uuid.uuid4().hex
    # JSON is data, not executable code; escape angle brackets for HTML safety.
    payload = fig.to_json().replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    lazy = 'true' if any(trace.type == 'scattergl' for trace in fig.data) else 'false'
    return ('<figure><div class="shiba-plot" data-lazy="'+lazy+'" id="'+identifier+'" style="height:520px"></div>'
            '<script type="application/json" id="'+identifier+'-data">'+payload+'</script>'
            '<figcaption>'+html.escape(caption)+'</figcaption></figure>')


INTERACTIVE_LOADER = r"""<script>
const options = {responsive:true, displaylogo:false, scrollZoom:false,
                 toImageButtonOptions:{format:'svg'}};
const observer = new IntersectionObserver(entries => {
  for (const entry of entries) {
    const el = entry.target;
    if (entry.isIntersecting && !el._rendering && !el._active) {
      el._rendering = true;
      const spec = el._saved || JSON.parse(document.getElementById(el.id+'-data').textContent);
      Plotly.newPlot(el, spec.data, spec.layout, options).then(() => {
        el._rendering = false; el._active = true;
        if (!el._visible) {el._saved = {data:el.data, layout:el.layout}; Plotly.purge(el); el._active = false;}
      });
    } else if (!entry.isIntersecting && el._active) {
      el._saved = {data:el.data, layout:el.layout};
      Plotly.purge(el); el._active = false;
    }
    el._visible = entry.isIntersecting;
  }
}, {rootMargin:'200px'});
document.querySelectorAll('.shiba-plot').forEach(el => {
  if (el.dataset.lazy === 'true') observer.observe(el);
  else {
    const spec = JSON.parse(document.getElementById(el.id+'-data').textContent);
    Plotly.newPlot(el, spec.data, spec.layout, options).then(() => {el._active = true;});
  }
});
</script>"""


def command_section(input_dir, command=None):
    """Prefer the launcher's command; otherwise use recorded pipeline provenance."""
    if not command:
        for path in (Path(input_dir).parent / 'report.json', Path(input_dir) / 'report.json'):
            if path.exists():
                try:
                    command = json.loads(path.read_text()).get('run', {}).get('command')
                except (ValueError, OSError, AttributeError):
                    continue
                if command:
                    break
    if not command:
        return '<h2>Analysis command</h2><p>Not recorded. Supply --shiba-command when regenerating this report.</p>'
    return '<h2>Analysis command</h2><pre class="command">' + html.escape(str(command)) + '</pre>'
