from src.util import convert_AA_1to3, convert_AA_3to1

###############################
# Test amino acid conversions #
###############################


def test_convert_AA_1to3_short():
    aa = ['M', 'H']
    assert(convert_AA_1to3(aa) == ['MET', 'HIS'])

def test_convert_AA_1to3_one():
    aa = ['L']
    assert(convert_AA_1to3(aa) == ['LEU'])

def test_convert_AA_1to3_zero_length():
    aa = []
    assert(convert_AA_1to3(aa) == [])

def test_convert_AA_3to1_short():
    aa = ['MET', 'HIS']
    assert(convert_AA_3to1(aa) == ['M', 'H'])

def test_convert_AA_3to1_one():
    aa = ['LEU']
    assert(convert_AA_3to1(aa) == ['L'])

def test_convert_AA_3to1_zero_length():
    aa = []
    assert(convert_AA_3to1(aa) == [])


