#!/usr/bin/env python3
"""Reproducible BB calibration and runtime experiment; not a claim of universal FDR control.

python test/simulate_factorial.py --output /tmp/calibration --events 200 --replicates 8
The two SE inclusion junctions deliberately share the same counts: components
are perfectly dependent, so this checks that shared evidence is not doubled.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from test_factorial import make_fixture, ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=200)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--depth", type=int, default=100)
    parser.add_argument("--rho", type=float, default=.05)
    parser.add_argument("--seed", type=int, default=731)
    parser.add_argument("--processes", type=int, default=2)
    args = parser.parse_args()
    if args.events < 4 or args.replicates < 2 or args.depth < 10 or not 0 < args.rho < 1:
        parser.error("Require events >= 4, replicates >= 2, depth >= 10 and 0 < rho < 1")
    args.output.mkdir(parents=True, exist_ok=False)
    meta = make_fixture(args.output, args.replicates)
    original = pd.read_csv(args.output / "events" / "EVENT_SE.txt", sep="\t").iloc[0]
    rng = np.random.default_rng(args.seed)
    g = (meta.genotype == "KO").astype(int).to_numpy()
    t = (meta.treatment == "Drug").astype(int).to_numpy()
    concentration = 1 / args.rho - 1
    events, counts, truth = [], [], []
    for i in range(args.events):
        interaction = 0. if i < args.events // 2 else 1.8
        # Nonzero main effects under the interaction null.
        mu = 1 / (1 + np.exp(-(-1.5 + .5*g + .3*t + interaction*g*t)))
        p = rng.beta(mu * concentration, (1-mu)*concentration)
        k = rng.binomial(args.depth, p); e = args.depth - k
        row = original.copy(); row['event_id'] = f"SE_{i+1}"; row['pos_id'] = f"SE@{i+1}"
        base = (i+1)*1000
        for name, start, end, values in [('a',base,base+100,k),('b',base+200,base+300,k),('c',base,base+300,e)]:
            key = f"chr1:{start}-{end}"; row[f'intron_{name}'] = key
            counts.append(['chr1', start, end, key, *values])
        events.append(row); truth.append(interaction)
    pd.DataFrame(events).to_csv(args.output / 'events' / 'EVENT_SE.txt', sep='\t', index=False)
    pd.DataFrame(counts, columns=['chr','start','end','ID',*meta['sample']]).to_csv(args.output/'junctions.bed', sep='\t', index=False)
    start = time.monotonic()
    subprocess.run([sys.executable, str(ROOT/'src'/'psi.py'), str(args.output/'junctions.bed'),
                    str(args.output/'events'), str(args.output/'results'), '--stat-method','beta-binomial',
                    '--sample-metadata',str(args.output/'samples.tsv'),'--formula','genotype * treatment',
                    '--reference-level','genotype=WT','--reference-level','treatment=Control',
                    '--coef','genotypeKO:treatmentDrug','-p',str(args.processes)], check=True)
    elapsed = time.monotonic() - start
    result = pd.read_csv(args.output/'results'/'event_statistics.tsv', sep='\t')
    result['truth'] = result.event_id.map({f'SE_{i+1}': value for i,value in enumerate(truth)})
    result.to_csv(args.output/'truth_and_results.tsv', sep='\t', index=False)
    null = result[result.truth == 0]; alternative = result[result.truth != 0]
    false_positive = int(null.p_event.lt(.05).sum())
    interval = binomtest(false_positive, len(null)).proportion_ci()
    called = result['Diff events'].eq('Yes')
    stats = pd.read_csv(args.output/'results'/'component_statistics.tsv', sep='\t')
    report = dict(seed=args.seed, events=args.events, replicates_per_cell=args.replicates,
                  depth=args.depth, rho=args.rho, elapsed_seconds=elapsed,
                  untestable=int(result.status.eq('untestable').sum()),
                  null_rejection_rate=false_positive/len(null),
                  null_rejection_95pct_interval=[interval.low,interval.high],
                  realized_false_discovery_proportion=float((called & result.truth.eq(0)).sum()/max(1,called.sum())),
                  power=float(alternative['Diff events'].eq('Yes').mean()),
                  median_estimated_rho=float(stats.loc[stats.status.eq('ok'),'rho'].median()),
                  note='One simulated run; false discovery proportion is not an estimate of expected FDR across simulations.')
    (args.output/'calibration.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
