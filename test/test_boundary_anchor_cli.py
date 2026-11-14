"""Boundary counting and merging must use matching defaults in every CLI mode."""
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import bam2junc
import sc2junc
from src.lib.general import apply_config_defaults


class TestBoundaryAnchorCLI(unittest.TestCase):
    def test_counting_and_merging_defaults_and_overrides(self):
        for module, run_args in (
            (bam2junc, ['--bam', 'sample.bam']),
            (sc2junc, ['--alignment', 'sample.bam', '--barcode', 'barcodes.tsv']),
        ):
            commands = [
                ['-i', 'samples.tsv', '-g', 'annotation.gtf', '-r', 'EVENT_RI.txt', '-o', 'counts.bed'],
                ['run', *run_args, '-g', 'annotation.gtf', '--prefix', 'sample'],
                ['merge', '-i', 'samples.tsv', '-d', 'counts', '-r', 'EVENT_RI.txt', '-o', 'counts.bed'],
            ]
            for command in commands:
                for option, expected in (([], 6), (['-b', '1'], 1), (['--boundary-anchor', '10'], 10)):
                    with self.subTest(module=module.__name__, command=command, option=option):
                        args = module.build_parser().parse_args(command + option)
                        self.assertEqual(args.boundary_anchor, expected)

    def test_pipeline_defaults_preserve_explicit_values(self):
        for mode in ('bulk', 'sc'):
            for supplied, expected in (({}, 6), ({'boundary_anchor_length': 1}, 1),
                                       ({'boundary_anchor_length': 10}, 10)):
                with self.subTest(mode=mode, supplied=supplied):
                    config = supplied.copy()
                    apply_config_defaults(config, mode=mode)
                    self.assertEqual(config['boundary_anchor_length'], expected)
