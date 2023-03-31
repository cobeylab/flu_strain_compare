from src.classes.flu_compare import compare_seq_no_reference, FluPngs

###################################################
# Test comparison with identical numbering scheme #
###################################################

def test_no_reference_short():
    fixture = [FluPngs("1", "PNGS1", 1/3),
            FluPngs("2", "PNGS2", 2/3),
            FluPngs("3", "PNGS3", 1/3)]
    comparison = compare_seq_no_reference([set(["1","2"]), set(["2", "3"])], lambda p: p+1)
    assert(all([x in fixture for x in comparison]) and all([x in comparison for x in fixture]))

def test_no_reference_no_common():
    fixture = [FluPngs("1", "PNGS1", 1/4),
            FluPngs("2", "PNGS2", 1/4),
            FluPngs("3", "PNGS3", 1/4),
            FluPngs("4", "PNGS4", 1/4)]
    comparison = compare_seq_no_reference([set(["1","2"]), set(["3", "4"])], lambda p: p+1)
    assert(all([x in fixture for x in comparison]) and all([x in comparison for x in fixture]))

def test_no_reference_identical():
    fixture = [FluPngs("1", "PNGS1", 1),
            FluPngs("2", "PNGS2", 1)]
    comparison = compare_seq_no_reference([set(["1","2"]), set(["1", "2"])], lambda p: p+1)
    assert(all([x in fixture for x in comparison]) and all([x in comparison for x in fixture]))

def test_no_reference_empty():
    comparison = compare_seq_no_reference([set([]), set([])], lambda p: p+1)
    assert(comparison == [])

##########################
# Test equality of PNGSs #
##########################


def test_equal_PNGS():
    assert(FluPngs("1", "PNGS1", 0) == FluPngs("1", "PNGS1", 0))

def test_not_equal_PNGS():
    assert(FluPngs("1", "PNGS1", 1) != FluPngs("1", "PNGS1", 0))
    assert(FluPngs("1", "PNGS1", 1) != FluPngs("1", "PNGS0", 1))
    assert(FluPngs("1", "PNGS1", 1) != FluPngs("0", "PNGS1", 1))
