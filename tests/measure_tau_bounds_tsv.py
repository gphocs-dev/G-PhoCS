import csv
import sys
from collections import defaultdict


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def measure_average_diff_between_bound_and_tau(tau_bounds_f):
    with open(tau_bounds_f) as f:
        tau_bounds = csv.reader(f, delimiter='\t')
        header = next(tau_bounds)
        names = {j: (bound_col, tau_col) for j, (bound_col, tau_col) in enumerate(chunks(header[1:], 2))}
        sum_diffs = defaultdict(float)

        num_rows = 0
        for row in tau_bounds:
            for j, (bound, tau) in enumerate(chunks(row[1:], 2)):
                diff = float(bound) - float(tau)
                assert diff > 0.
                sum_diffs[j] += diff
            num_rows += 1

        print("for file '%s':" % tau_bounds_f)
        for k, v in sum_diffs.items():
            print("%.40f is the average diff of cols %s" % (v / num_rows, names[k]))


if __name__ == '__main__':
    tau_bounds_file = sys.argv[1]
    measure_average_diff_between_bound_and_tau(tau_bounds_file)
