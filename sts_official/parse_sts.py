"""
parse_sts.py -- summarises the official STS finalAnalysisReport.txt files in sts_official/reports/
into results/nist_official_sts.json (per test: number of sub-tests, how many fail the proportion or
uniformity criterion as flagged by STS with '*', minimum proportion observed).
"""
import json, os, re, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ROW = re.compile(r"^\s*((?:\d+\s+){10})(\S+)\s*(\*?)\s+(\d+)/(\d+)\s*(\*?)\s+(\S+)\s*$")

def parse(path):
    tests = {}
    for line in open(path, errors="replace"):
        m = ROW.match(line)
        if not m:
            continue
        p = m.group(2); pflag = m.group(3) == "*"; k, n = int(m.group(4)), int(m.group(5))
        prop_flag = m.group(6) == "*"; name = m.group(7)
        t = tests.setdefault(name, {"subtests": 0, "failed": 0, "min_proportion": 1.0, "n_min": n})
        t["subtests"] += 1
        t["failed"] += int(pflag or prop_flag)
        t["min_proportion"] = min(t["min_proportion"], k / n if n else 0.0)
        t["n_min"] = min(t["n_min"], n)
    return tests

if __name__ == "__main__":
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, "reports", "*.txt"))):
        cfg = os.path.basename(f)[:-4]
        tests = parse(f)
        out[cfg] = {"tests": tests, "tests_passed": sum(t["failed"] == 0 for t in tests.values()),
                    "tests_total": len(tests)}
        print(cfg, f"{out[cfg]['tests_passed']}/{out[cfg]['tests_total']} tests with no flagged sub-test")
    json.dump(out, open(os.path.join(ROOT, "results", "nist_official_sts.json"), "w"), indent=1)
