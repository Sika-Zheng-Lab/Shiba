"""Config-to-pipeline routing and real factorial execution through shiba.py."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import shutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import shiba
from test.test_factorial import make_fixture, has_backend


class TestLauncher(unittest.TestCase):
    def fixture(self, root):
        experiment = root / 'experiment.tsv'
        experiment.write_text('sample\tbam\tgroup\na\t/a.bam\tA\nb\t/b.bam\tB\n')
        gtf = root / 'annotation.gtf'
        gtf.write_text('# fixture\n')
        return dict(workdir=str(root / 'out'), experiment_table=str(experiment), gtf=str(gtf))

    def launch(self, root, config):
        path = root / 'config.yaml'
        path.write_text(yaml.safe_dump(config))
        with patch.object(sys, 'argv', ['shiba.py', str(path), '-p', '4']), \
             patch.object(shiba.time, 'sleep'), \
             patch.object(shiba.general, 'generate_report'), \
             patch.object(shiba.general, 'execute_command', return_value=0) as execute:
            shiba.main()
        return [call.args[0] for call in execute.call_args_list]

    def test_full_pipeline_routing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for method in (None, "fisher", "beta_binomial"):
                regression = method == "beta_binomial"
                config = self.fixture(root)
                if regression:
                    config.update(stat_method='beta_binomial', formula='group', coef=['groupB'],
                                  ttest=True, filter_group=['group'])
                else:
                    config.update(reference_group='A', alternative_group='B')
                if method is not None:
                    config["stat_method"] = method
                commands = self.launch(root, config)
                self.assertEqual(len(commands), 7)
                self.assertTrue(all(c[0] == sys.executable for c in commands))
                psi = commands[3]
                self.assertNotIn('-b', psi)
                self.assertEqual('--stat-method' in psi, regression)
                if regression:
                    self.assertEqual(psi[psi.index('--formula') + 1], 'group')
                    self.assertEqual(psi[psi.index('-t') + 1], 'False')
                    self.assertEqual(psi[psi.index('--filter-group') + 1], 'group')
                    self.assertEqual(commands[4][commands[4].index('-r') + 1], 'NA')

    def test_removed_cli_flags(self):
        for script in ('psi.py', 'scpsi.py'):
            result = subprocess.run([sys.executable, str(ROOT / 'src' / script),
                                     'junctions', 'events', 'output', '--beta-binomial'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('unrecognized arguments: --beta-binomial', result.stderr)

    def test_removed_config_stops_before_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.fixture(root)
            config.update(beta_binomial=False, reference_group='A', alternative_group='B')
            path = root / 'config.yaml'
            path.write_text(yaml.safe_dump(config))
            with patch.object(sys, 'argv', ['shiba.py', str(path)]), \
                 patch.object(shiba.time, 'sleep'), \
                 patch.object(shiba.general, 'execute_command') as execute:
                with self.assertRaises(SystemExit):
                    shiba.main()
                execute.assert_not_called()

    @unittest.skipUnless(has_backend(), 'R packages glmmTMB/jsonlite unavailable')
    def test_real_factorial_from_yaml(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = make_fixture(root)
            config = self.fixture(root)
            out = root / 'out'
            (out / 'junctions').mkdir(parents=True)
            shutil.copy(root / 'junctions.bed', out / 'junctions' / 'junctions.bed')
            shutil.copytree(root / 'events', out / 'events')
            config.update(stat_method='beta_binomial', sample_metadata=str(root / 'samples.tsv'),
                          formula='genotype * treatment + batch', categorical=['batch'],
                          coef=['genotypeKO:treatmentDrug'],
                          reference_levels={'genotype': 'WT', 'treatment': 'Control'})
            metadata['bam'] = '/unused.bam'
            metadata['group'] = 'all'
            metadata.to_csv(config['experiment_table'], sep='\t', index=False)
            path = root / 'config.yaml'
            path.write_text(yaml.safe_dump(config))
            result = subprocess.run([sys.executable, str(ROOT / 'shiba.py'), str(path),
                                     '--start-step', '4', '--mame', '-p', '2'],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            import pandas as pd
            frame = pd.read_csv(out / 'results/splicing/event_statistics.tsv', sep='\t')
            self.assertEqual(frame.set_index('event_id').loc['SE_1', 'Diff events'], 'Yes')
            self.assertNotIn('p_beta', frame.columns)
            import json
            manifest = json.loads((out / 'results/splicing/analysis.json').read_text())
            self.assertEqual(manifest['stat_method'], 'beta_binomial')
