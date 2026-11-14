"""Compare event-level PSI and differential calls from two Shiba output directories."""

import argparse
import csv
import json
import math
from pathlib import Path


EVENTS = ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE")


def read_events(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != len({row["event_id"] for row in rows}):
        raise ValueError(f"Duplicate event IDs in {path}")
    return {row["event_id"]: row for row in rows}


def different(first, second, tolerance):
    if first == second:
        return False
    try:
        left, right = float(first), float(second)
    except (TypeError, ValueError):
        return True
    if math.isnan(left) and math.isnan(right):
        return False
    return not math.isclose(left, right, rel_tol=0, abs_tol=tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument("--report", required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    report = {}

    for event in EVENTS:
        before = read_events(Path(args.baseline) / f"PSI_{event}.txt")
        after = read_events(Path(args.candidate) / f"PSI_{event}.txt")
        shared = before.keys() & after.keys()
        columns = before[next(iter(before))].keys() if before else ()
        psi_columns = [name for name in columns if name.endswith("_PSI") or name == "dPSI"]
        tracked = psi_columns + [name for name in ("Diff events", "q", "q_beta") if name in columns]
        changed = {
            name: sum(different(before[key].get(name), after[key].get(name), args.tolerance)
                      for key in shared)
            for name in tracked
        }
        changed_events = sum(
            any(different(before[key].get(name), after[key].get(name), args.tolerance)
                for name in tracked)
            for key in shared
        )
        report[event] = {
            "baseline_rows": len(before),
            "candidate_rows": len(after),
            "shared_rows": len(shared),
            "baseline_only_rows": len(before.keys() - after.keys()),
            "candidate_only_rows": len(after.keys() - before.keys()),
            "changed_events": changed_events,
            "changed_columns": changed,
        }

    with open(args.report, "w") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
