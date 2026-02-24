import unittest
import tempfile
import os
import sys

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
from lib.general import (
    validate_file_exists,
    validate_experiment_table_columns,
    validate_groups_bulk,
    validate_groups_sc,
    validate_config_types,
    validate_config,
)


class TestValidateFileExists(unittest.TestCase):
    def test_existing_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(b"test")
            tmp_path = f.name
        try:
            errors = validate_file_exists(tmp_path, "Test file")
            self.assertEqual(errors, [])
        finally:
            os.unlink(tmp_path)

    def test_nonexistent_file(self):
        errors = validate_file_exists("/nonexistent/path/file.txt", "Test file")
        self.assertEqual(len(errors), 1)
        self.assertIn("Test file not found", errors[0])


class TestValidateExperimentTableColumns(unittest.TestCase):
    def _write_tsv(self, header, rows=None):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write(header + "\n")
        if rows:
            for row in rows:
                f.write(row + "\n")
        f.close()
        return f.name

    def test_bulk_all_columns_present(self):
        path = self._write_tsv("sample\tbam\tgroup", ["s1\tb1.bam\tRef"])
        try:
            errors = validate_experiment_table_columns(path, ["sample", "bam", "group"])
            self.assertEqual(errors, [])
        finally:
            os.unlink(path)

    def test_bulk_missing_column(self):
        path = self._write_tsv("sample\tbam", ["s1\tb1.bam"])
        try:
            errors = validate_experiment_table_columns(path, ["sample", "bam", "group"])
            self.assertEqual(len(errors), 1)
            self.assertIn("group", errors[0])
        finally:
            os.unlink(path)

    def test_sc_all_columns_present(self):
        path = self._write_tsv("barcode\tSJ", ["/path/to/bc.tsv\t/path/to/SJ"])
        try:
            errors = validate_experiment_table_columns(path, ["barcode", "SJ"])
            self.assertEqual(errors, [])
        finally:
            os.unlink(path)

    def test_empty_file(self):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.close()
        try:
            errors = validate_experiment_table_columns(f.name, ["sample", "bam", "group"])
            self.assertEqual(len(errors), 1)
            self.assertIn("empty", errors[0])
        finally:
            os.unlink(f.name)

    def test_file_not_found(self):
        errors = validate_experiment_table_columns("/nonexistent.tsv", ["sample"])
        self.assertEqual(len(errors), 1)
        self.assertIn("not found", errors[0])


class TestValidateGroupsBulk(unittest.TestCase):
    def _write_experiment(self, rows):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("sample\tbam\tgroup\n")
        for row in rows:
            f.write(row + "\n")
        f.close()
        return f.name

    def test_valid_groups(self):
        path = self._write_experiment([
            "s1\tb1.bam\tRef",
            "s2\tb2.bam\tAlt",
        ])
        try:
            errors = validate_groups_bulk(path, "Ref", "Alt")
            self.assertEqual(errors, [])
        finally:
            os.unlink(path)

    def test_reference_not_found(self):
        path = self._write_experiment([
            "s1\tb1.bam\tGroupA",
            "s2\tb2.bam\tGroupB",
        ])
        try:
            errors = validate_groups_bulk(path, "Ref", "GroupB")
            self.assertEqual(len(errors), 1)
            self.assertIn('reference_group "Ref" not found', errors[0])
            self.assertIn("GroupA", errors[0])
        finally:
            os.unlink(path)

    def test_alternative_not_found(self):
        path = self._write_experiment([
            "s1\tb1.bam\tRef",
            "s2\tb2.bam\tRef",
        ])
        try:
            errors = validate_groups_bulk(path, "Ref", "Alt")
            self.assertEqual(len(errors), 1)
            self.assertIn('alternative_group "Alt" not found', errors[0])
        finally:
            os.unlink(path)

    def test_both_not_found(self):
        path = self._write_experiment([
            "s1\tb1.bam\tX",
            "s2\tb2.bam\tY",
        ])
        try:
            errors = validate_groups_bulk(path, "Ref", "Alt")
            self.assertEqual(len(errors), 2)
        finally:
            os.unlink(path)

    def test_same_reference_and_alternative(self):
        path = self._write_experiment([
            "s1\tb1.bam\tRef",
            "s2\tb2.bam\tRef",
        ])
        try:
            errors = validate_groups_bulk(path, "Ref", "Ref")
            self.assertTrue(any("must be different" in e for e in errors))
        finally:
            os.unlink(path)

    def test_file_not_found(self):
        errors = validate_groups_bulk("/nonexistent.tsv", "Ref", "Alt")
        self.assertTrue(any("not found" in e for e in errors))


class TestValidateGroupsSc(unittest.TestCase):
    def _write_barcode(self, rows):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("barcode\tgroup\n")
        for row in rows:
            f.write(row + "\n")
        f.close()
        return f.name

    def _write_sc_experiment(self, barcode_paths):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("barcode\tSJ\n")
        for bp in barcode_paths:
            f.write(f"{bp}\t/path/to/SJ\n")
        f.close()
        return f.name

    def test_valid_groups(self):
        bc = self._write_barcode(["ACGT\tCluster-1", "TGCA\tCluster-2"])
        exp = self._write_sc_experiment([bc])
        try:
            errors = validate_groups_sc(exp, "Cluster-1", "Cluster-2")
            self.assertEqual(errors, [])
        finally:
            os.unlink(bc)
            os.unlink(exp)

    def test_reference_not_found(self):
        bc = self._write_barcode(["ACGT\tCluster-1", "TGCA\tCluster-2"])
        exp = self._write_sc_experiment([bc])
        try:
            errors = validate_groups_sc(exp, "Missing", "Cluster-2")
            self.assertEqual(len(errors), 1)
            self.assertIn('reference_group "Missing" not found', errors[0])
        finally:
            os.unlink(bc)
            os.unlink(exp)

    def test_same_reference_and_alternative(self):
        bc = self._write_barcode(["ACGT\tCluster-1"])
        exp = self._write_sc_experiment([bc])
        try:
            errors = validate_groups_sc(exp, "Cluster-1", "Cluster-1")
            self.assertTrue(any("must be different" in e for e in errors))
        finally:
            os.unlink(bc)
            os.unlink(exp)

    def test_barcode_file_not_found(self):
        exp = self._write_sc_experiment(["/nonexistent/barcode.tsv"])
        try:
            errors = validate_groups_sc(exp, "Ref", "Alt")
            self.assertTrue(any("Barcode file not found" in e for e in errors))
        finally:
            os.unlink(exp)


class TestValidateConfigTypes(unittest.TestCase):
    def test_valid_bulk_config(self):
        config = {
            "fdr": 0.05,
            "delta_psi": 0.1,
            "minimum_reads": 10,
            "minimum_anchor_length": 6,
            "minimum_intron_length": 70,
            "maximum_intron_length": 500000,
            "strand": "XS",
        }
        errors = validate_config_types(config, mode="bulk")
        self.assertEqual(errors, [])

    def test_valid_sc_config(self):
        config = {
            "fdr": 0.05,
            "delta_psi": 0.1,
            "minimum_reads": 10,
        }
        errors = validate_config_types(config, mode="sc")
        self.assertEqual(errors, [])

    def test_fdr_out_of_range(self):
        errors = validate_config_types({"fdr": 1.5}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("fdr", errors[0])

    def test_fdr_zero(self):
        errors = validate_config_types({"fdr": 0}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("fdr", errors[0])

    def test_fdr_not_number(self):
        errors = validate_config_types({"fdr": "abc"}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("fdr must be a number", errors[0])

    def test_delta_psi_out_of_range(self):
        errors = validate_config_types({"delta_psi": -0.1}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("delta_psi", errors[0])

    def test_minimum_reads_negative(self):
        errors = validate_config_types({"minimum_reads": -1}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("minimum_reads", errors[0])

    def test_minimum_reads_not_int(self):
        errors = validate_config_types({"minimum_reads": "abc"}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("minimum_reads must be an integer", errors[0])

    def test_invalid_strand(self):
        errors = validate_config_types({"strand": "invalid"}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("strand", errors[0])

    def test_strand_valid_values(self):
        for strand in ["XS", "0", "1", "2"]:
            errors = validate_config_types({"strand": strand}, mode="bulk")
            self.assertEqual(errors, [], f"strand={strand} should be valid")

    def test_max_intron_less_than_min(self):
        config = {
            "minimum_intron_length": 100,
            "maximum_intron_length": 50,
        }
        errors = validate_config_types(config, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("maximum_intron_length", errors[0])

    def test_negative_anchor_length(self):
        errors = validate_config_types({"minimum_anchor_length": -1}, mode="bulk")
        self.assertEqual(len(errors), 1)
        self.assertIn("minimum_anchor_length", errors[0])

    def test_sc_mode_ignores_bulk_params(self):
        config = {
            "minimum_anchor_length": -999,
            "strand": "INVALID",
        }
        errors = validate_config_types(config, mode="sc")
        self.assertEqual(errors, [])


class TestValidateConfig(unittest.TestCase):
    """Integration tests for validate_config facade function."""

    def _write_tsv(self, header, rows):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write(header + "\n")
        for row in rows:
            f.write(row + "\n")
        f.close()
        return f.name

    def test_valid_bulk_config(self):
        gtf = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False)
        gtf.write(b"dummy gtf")
        gtf.close()
        exp = self._write_tsv("sample\tbam\tgroup", [
            "s1\tb1.bam\tRef",
            "s2\tb2.bam\tAlt",
        ])
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": gtf.name,
            "experiment_table": exp,
            "reference_group": "Ref",
            "alternative_group": "Alt",
            "fdr": 0.05,
            "delta_psi": 0.1,
            "minimum_reads": 10,
            "minimum_anchor_length": 6,
            "minimum_intron_length": 70,
            "maximum_intron_length": 500000,
            "strand": "XS",
        }
        try:
            errors = validate_config(config, mode="bulk")
            self.assertEqual(errors, [])
        finally:
            os.unlink(gtf.name)
            os.unlink(exp)

    def test_bulk_config_with_wrong_groups(self):
        gtf = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False)
        gtf.write(b"dummy gtf")
        gtf.close()
        exp = self._write_tsv("sample\tbam\tgroup", [
            "s1\tb1.bam\tGroupA",
            "s2\tb2.bam\tGroupB",
        ])
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": gtf.name,
            "experiment_table": exp,
            "reference_group": "Ref",
            "alternative_group": "Alt",
            "fdr": 0.05,
            "delta_psi": 0.1,
            "minimum_reads": 10,
        }
        try:
            errors = validate_config(config, mode="bulk")
            group_errors = [e for e in errors if "not found in experiment table" in e]
            self.assertEqual(len(group_errors), 2)
        finally:
            os.unlink(gtf.name)
            os.unlink(exp)

    def test_bulk_config_missing_ref_and_alt(self):
        gtf = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False)
        gtf.write(b"dummy gtf")
        gtf.close()
        exp = self._write_tsv("sample\tbam\tgroup", ["s1\tb1.bam\tRef"])
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": gtf.name,
            "experiment_table": exp,
            "fdr": 0.05,
        }
        try:
            errors = validate_config(config, mode="bulk")
            self.assertTrue(any("reference_group is not specified" in e for e in errors))
            self.assertTrue(any("alternative_group is not specified" in e for e in errors))
        finally:
            os.unlink(gtf.name)
            os.unlink(exp)

    def test_bulk_config_only_psi_skips_group_check(self):
        gtf = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False)
        gtf.write(b"dummy gtf")
        gtf.close()
        exp = self._write_tsv("sample\tbam\tgroup", ["s1\tb1.bam\tRef"])
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": gtf.name,
            "experiment_table": exp,
            "only_psi": True,
            "only_psi_group": True,
            "fdr": 0.05,
        }
        try:
            errors = validate_config(config, mode="bulk")
            group_errors = [e for e in errors if "reference_group" in e or "alternative_group" in e]
            self.assertEqual(len(group_errors), 0)
        finally:
            os.unlink(gtf.name)
            os.unlink(exp)

    def test_multiple_errors_collected(self):
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": "/nonexistent.gtf",
            "experiment_table": "/nonexistent.tsv",
            "reference_group": "Ref",
            "alternative_group": "Ref",  # same as reference
            "fdr": 5.0,  # out of range
            "delta_psi": -1,  # out of range
        }
        errors = validate_config(config, mode="bulk")
        # Should have: gtf not found, experiment not found, ref==alt, fdr, delta_psi
        self.assertGreaterEqual(len(errors), 4)

    def test_valid_sc_config(self):
        gtf = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False)
        gtf.write(b"dummy gtf")
        gtf.close()
        bc = self._write_tsv("barcode\tgroup", [
            "ACGT\tCluster-1",
            "TGCA\tCluster-2",
        ])
        exp = self._write_tsv("barcode\tSJ", [
            f"{bc}\t/path/to/SJ",
        ])
        config = {
            "workdir": "/tmp/test_workdir",
            "gtf": gtf.name,
            "experiment_table": exp,
            "reference_group": "Cluster-1",
            "alternative_group": "Cluster-2",
            "fdr": 0.05,
            "delta_psi": 0.1,
            "minimum_reads": 10,
        }
        try:
            errors = validate_config(config, mode="sc")
            self.assertEqual(errors, [])
        finally:
            os.unlink(gtf.name)
            os.unlink(bc)
            os.unlink(exp)


if __name__ == "__main__":
    unittest.main()
