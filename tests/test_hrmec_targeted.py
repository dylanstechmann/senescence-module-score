import unittest

from senescore import hrmec_targeted as ht
from senescore.genes import GENE_SETS


class TargetedPlanTests(unittest.TestCase):
    def test_removals_are_the_frozen_lists_and_sit_inside_senmayo(self):
        self.assertEqual(set(ht.REMOVALS), {"recurring3", "late_top10", "etoposide_top10", "union_top"})
        self.assertEqual(len(ht.LATE_TOP10), 10)
        self.assertEqual(len(ht.ETOPOSIDE_TOP10), 10)
        senmayo = set(GENE_SETS["senmayo"]["symbols"])
        for label, genes in ht.REMOVALS.items():
            self.assertTrue(set(genes) <= senmayo, label)
            self.assertEqual(len(genes), len(set(genes)), label)
        self.assertEqual(set(ht.REMOVALS["union_top"]),
                         set(ht.RECURRING3) | set(ht.LATE_TOP10) | set(ht.ETOPOSIDE_TOP10))


if __name__ == "__main__":
    unittest.main()
