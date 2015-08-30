import vcf
import SDST.strelka
import os
from nose.tools import assert_almost_equal

TESTFILE='strelkaRaw.vcf.gz'

def _dpath(path=''):
    """get path to a data file (relative to the directory this
	test lives in)"""
    return os.path.realpath(os.path.join(os.path.dirname(__file__), path))


def test_addHeader():
    """Make sure that addStrelkaHeaders adds at least one header to vcf reader"""
    reader = vcf.Reader(open(os.path.join(_dpath(),TESTFILE)))
    SDST.strelka.addStrelkaHeaders(reader)
    assert reader.infos['TUMREF']

def test_modifyStrelkaRow():
    """Make sure that modifyStrelkaRow adds appropriate info to model record"""
    reader = vcf.Reader(open(os.path.join(_dpath(),TESTFILE)))
    row = SDST.strelka.modifyStrelkaRow(next(reader))
    assert row.INFO['TUMREF']==2
    assert row.INFO['TUMALT']==5
    assert row.INFO['NORMREF']==10
    assert row.INFO['NORMALT']==0
    assert_almost_equal(row.INFO['TUMVAF'],5/7)
    assert row.INFO['TUMVARFRACTION']==1

