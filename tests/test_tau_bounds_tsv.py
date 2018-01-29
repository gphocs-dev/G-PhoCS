import csv

try:
    import pytest
except:
    pass

TAU_BOUNDS_TSV_FILE = 'out/sample-tau-bounds.tsv'

tau_bounds_file = open(TAU_BOUNDS_TSV_FILE)


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def test_taus_lower_then_bounds():
    tau_bounds = csv.reader(tau_bounds_file, delimiter='\t')
    next(tau_bounds)  # skip header

    for i, row in enumerate(tau_bounds):
        for j, (bound, tau) in enumerate(chunks(row[1:], 2)):
            assert bound >= tau, "bad bound in row %d in pair %d" % (i + 2, j)





if __name__ == '__main__':
    pytest.main()
