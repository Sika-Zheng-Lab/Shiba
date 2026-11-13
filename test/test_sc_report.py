"""Single-cell reports share the bulk interactive theme without requiring PCA."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT))
import plots
import scshiba


def fixture(root):
    source = root/'results'
    source.mkdir(parents=True)
    rows = []
    for event in plots.get_splicing_event_config():
        code = event['code']
        pd.DataFrame({
            'event_id':[code+'_1',code+'_2'], 'gene_name':['Gene A','Gene B'],
            'ref_PSI':[.2,.6], 'alt_PSI':[.7,.5], 'dPSI':[.5,-.1],
            'q':[.001,.8], 'Diff events':['Yes','No']
        }).to_csv(source/f'PSI_{code}.txt', sep='\t', index=False)
        rows.append([code,'Up','annotated',1])
    pd.DataFrame(rows,columns=['AS','Direction','Label','Number']).to_csv(source/'summary.txt',sep='\t',index=False)
    (root/'report.json').write_text(json.dumps({'run':{'command':'scshiba.py "my config.yaml"'}}))
    return source


class TestSingleCellReport(unittest.TestCase):
    def test_interactive_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=fixture(root); output=root/'plots'
            plots.write_single_cell_report(source,output)
            content=(output/'summary.html').read_text()
            self.assertIn('scShiba Fisher analysis',content)
            self.assertIn('my config.yaml',content)
            self.assertEqual(content.count('class="shiba-plot"'),17)
            self.assertNotIn('TPM PCA',content)
            self.assertNotIn('Plot not available',content)
            self.assertEqual(len(list((output/'svg').glob('*.svg'))),17)
            plots.write_single_cell_report(source,output,only_psi=True)
            content=(output/'summary.html').read_text()
            self.assertIn('differential testing was not performed',content)
            self.assertNotIn('class="shiba-plot"',content)

    def test_launcher(self):
        from argparse import Namespace
        with tempfile.TemporaryDirectory() as tmp:
            config=dict(workdir=tmp,experiment_table='experiment.tsv',gtf='annotation.gtf',
                        reference_group='Ref',alternative_group='Alt',fdr=.05,
                        delta_psi=.1,minimum_reads=10,only_psi=False,excel=False)
            for only in (False,True):
                config['only_psi']=only
                with patch.object(scshiba,'parse_args',return_value=Namespace(config='config.yaml',process=1,start_step=4,verbose=False)), \
                     patch.object(scshiba.general,'load_config',return_value=config), \
                     patch.object(scshiba.general,'validate_config',return_value=[]), \
                     patch.object(scshiba.general,'execute_command',return_value=0) as execute, \
                     patch.object(scshiba.general,'generate_report'), patch.object(scshiba.time,'sleep'):
                    scshiba.main()
                command=execute.call_args.args[0]
                self.assertIn('--single-cell',command)
                self.assertIn('-s',command)
                self.assertEqual('--only-psi' in command,only)

if __name__ == '__main__':
    unittest.main()
