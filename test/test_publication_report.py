"""Publication exports, self-contained HTML, and missing-data behavior."""
import base64
import json
import sys
import tempfile
import unittest
from pathlib import Path
from html.parser import HTMLParser
from unittest.mock import patch
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from lib import figure_style
from lib.factorial_report import write_report
import plots

class Images(HTMLParser):
    def __init__(self):
        super().__init__(); self.sources=[]; self.scripts=[]
    def handle_starttag(self, tag, attrs):
        attrs=dict(attrs)
        if tag=='img': self.sources.append(attrs['src'])
        if tag=='script': self.scripts.append(attrs)

class TestPublicationReport(unittest.TestCase):
    def test_font_fallback(self):
        original = figure_style.font_manager.findfont
        def without_arial(properties, *args, **kwargs):
            if 'Arial' in properties.get_family():
                raise ValueError('missing Arial')
            return original(properties, *args, **kwargs)
        with patch.object(figure_style.font_manager, 'findfont', side_effect=without_arial):
            paths = figure_style.configure()
            self.assertIn('DejaVuSans', paths['normal'])
            import matplotlib.pyplot as plt
            with tempfile.TemporaryDirectory() as tmp:
                fig, ax = plt.subplots();ax.set_xlabel('Fallback font')
                figure_style.export(fig,tmp,'fallback','Fallback')
                self.assertTrue((Path(tmp)/'pdf/fallback.pdf').exists())

    def test_interactive_payload_round_trip(self):
        import re
        import plotly.graph_objects as go
        label = '</script><test>&'
        content = figure_style.plot_html(go.Figure(go.Scatter(x=[1], y=[2], text=[label])), 'Caption')
        payload = re.search(r'<script type="application/json"[^>]*>(.*?)</script>', content).group(1)
        self.assertNotIn('<', payload)
        spec = json.loads(payload)
        self.assertEqual(spec['data'][0]['text'], [label])
        self.assertEqual(spec['layout']['font']['family'], 'Arial, sans-serif')
        self.assertFalse(spec['layout']['xaxis']['showgrid'])

    def test_command_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); results=root/'results';results.mkdir()
            self.assertIn('Not recorded', figure_style.command_section(results))
            (root/'report.json').write_text(json.dumps({'run':{'command':'python shiba.py "config file.yaml"'}}))
            self.assertIn('config file.yaml', figure_style.command_section(results))
            self.assertIn('explicit &lt;command&gt;', figure_style.command_section(results,'explicit <command>'))

    def test_factorial_missing_predictions_and_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'input'/'splicing'; source.mkdir(parents=True)
            (source/'analysis.json').write_text(json.dumps({'formula':'genotype * treatment <test>'}))
            frame=pd.DataFrame({'contrast':['a','a','b'], 'event_type':['SE']*3,
                'event_id':['SE_1','SE_2','SE_3'], 'q_event':[0,np.nan,.2],
                'conservative_component_effect':[1,np.nan,-1],
                'status':['ok','filtered_low_coverage','ok'], 'Diff events':['Yes','No','No']})
            frame.to_csv(source/'event_statistics.tsv',sep='\t',index=False)
            for empty in (False, True):
                if empty: frame.iloc[:0].to_csv(source/'event_statistics.tsv',sep='\t',index=False)
                out=root/str(empty);write_report(root/'input',out,shiba_command='python shiba.py config.yaml --tag <demo>')
                content=(out/'summary.html').read_text(); parsed=Images();parsed.feed(content)
                self.assertIn('class="shiba-plot"', content)
                self.assertTrue(all('src' not in script for script in parsed.scripts))
                self.assertIn('&lt;test&gt;',content)
                self.assertIn('python shiba.py config.yaml --tag &lt;demo&gt;', content)
                self.assertIn('No model predictions',content)
                for src in parsed.sources:
                    self.assertTrue(src.startswith('data:image/png;base64,'))
                    self.assertTrue(base64.b64decode(src.split(',')[1]).startswith(b'\x89PNG'))
                self.assertTrue((out/'svg/barplot_splicing_summary.svg').exists())
                self.assertIn('Arial',(out/'svg/barplot_splicing_summary.svg').read_text())

    def test_fisher_and_pca_interactive_with_static_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'input/splicing').mkdir(parents=True)
            frame=pd.DataFrame({'dPSI':[.3,-.4], 'q':[.01,.02], 'ref_PSI':[.2,.8],
                                'alt_PSI':[.5,.4], 'Diff events':['Yes','Yes']})
            frame.to_csv(root/'input/splicing/PSI_SE.txt',sep='\t',index=False)
            plots.plots('SE',root/'input',root/'out')
            pca=pd.DataFrame({'group':['ctrl','treated'], 'PC1':[.123456,2.5], 'PC2':[1.,-1.]})
            original=pca.copy(deep=True)
            plots.plots_pca('PSI',pca,'60','20',root/'out')
            pd.testing.assert_frame_equal(pca,original)
            summary = pd.DataFrame({'AS':['SE','SE'], 'Direction':['Up','Down'], 'Label':['annotated','unannotated'], 'Number':[1,1]})
            plots.barplot_splicing(summary,root/'out/png/barplot_splicing_summary.png',root/'out/pdf/barplot_splicing_summary.pdf')
            plots.write_summary_html('test <command>',root/'input',root/'out')
            text=(root/'out/summary.html').read_text();parsed=Images();parsed.feed(text)
            self.assertEqual(text.count('class="shiba-plot"'),4)
            self.assertTrue(all('src' not in script for script in parsed.scripts))
            self.assertIn('Plotly.newPlot', text)
            for name in ['volcano_SE','scatter_SE','pca_PSI']:
                self.assertTrue((root/'out/pdf'/f'{name}.pdf').exists())
