"""
Tests for miscellaneous utilities.
"""
from __future__ import print_function, division, absolute_import
try:
    import cPickle as pickle
except ImportError:
    import pickle
import unittest

from rdkit import Chem

from rdkit_utils import PicklableMol


class TestPicklableMol(unittest.TestCase):
    def setUp(self):
        """
        Set up for tests.
        """
        self.mol = Chem.MolFromSmiles('CC(C)(C)NC[C@@H](C1=CC(=C(C=C1)O)CO)O')
        self.mol.SetProp('_Name', 'levalbuterol')

    def test_picklable_mol(self):
        """Test PicklableMol."""
        mol = pickle.loads(pickle.dumps(PicklableMol(self.mol),
                                        pickle.HIGHEST_PROTOCOL))
        assert mol.HasProp('_Name')
        assert mol.GetProp('_Name') == self.mol.GetProp('_Name')

        # make sure stereochemistry is preserved
        for a, b in zip(mol.GetAtoms(), self.mol.GetAtoms()):
            assert a.GetChiralTag() == b.GetChiralTag()

        # full comparison
        assert mol.ToBinary() == self.mol.ToBinary()
