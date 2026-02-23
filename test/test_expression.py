import unittest
import pandas as pd
import os
import sys
# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
try:
    from lib.expression import ExpressionProcessor
    HAS_EXPRESSION = True
except ImportError:
    HAS_EXPRESSION = False

@unittest.skipUnless(HAS_EXPRESSION, "pysam not installed")
class TestExpressionProcessor(unittest.TestCase):
    def setUp(self):
        # Sample dataframe to test ExpressionProcessor
        self.df = pd.DataFrame({
            "Gene": ["Gene1", "Gene2", "Gene3"],
            "Length": [1000, 2000, 1500],
            "Sample1": [50, 100, 75],
            "Sample2": [60, 120, 90]
        })
        self.processor = ExpressionProcessor(self.df.copy())

    def test_TPM(self):
        processed_df = self.processor.TPM()
        # Expected TPM values based on calculation
        expected_df = pd.DataFrame({
            "Gene": ["Gene1", "Gene2", "Gene3"],
            "Sample1": [333333.33, 333333.33, 333333.33],
            "Sample2": [333333.33, 333333.33, 333333.33]
        })
        pd.testing.assert_frame_equal(processed_df.set_index("Gene"), expected_df.set_index("Gene"), atol=0.01)

    def test_CPM(self):
        processed_df = self.processor.CPM()
        # Corrected expected CPM values
        expected_df = pd.DataFrame({
            "Gene": ["Gene1", "Gene2", "Gene3"],
            "Sample1": [200000.0, 400000.0, 400000.0],
            "Sample2": [200000.0, 400000.0, 400000.0]
        })
        total_sample1 = sum(self.df["Sample1"])
        total_sample2 = sum(self.df["Sample2"])
        expected_df["Sample1"] = [
            (50 / total_sample1) * 1e6,
            (100 / total_sample1) * 1e6,
            (75 / total_sample1) * 1e6,
        ]
        expected_df["Sample2"] = [
            (60 / total_sample2) * 1e6,
            (120 / total_sample2) * 1e6,
            (90 / total_sample2) * 1e6,
        ]
        pd.testing.assert_frame_equal(processed_df.set_index("Gene"), expected_df.set_index("Gene"), atol=0.01)

if __name__ == "__main__":
    unittest.main()
