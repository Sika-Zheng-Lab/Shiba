"""Integer junction comparisons for Shiba's event-centric regression."""
from itertools import product

import numpy as np
import pandas as pd

EVENT_TYPES = ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE")
COUNT_COLUMNS = ["event_type", "event_id", "component_id", "sample", "inclusion", "exclusion"]
COMPONENT_COLUMNS = ["event_type", "event_id", "component_id", "inclusion_junction", "exclusion_junction"]


def component_pairs(event_type, row):
    """Return ordered (inclusion ID, exclusion ID) pairs; no pooled trials."""
    if event_type == "SE":
        pairs = [(row["intron_a"], row["intron_c"]), (row["intron_b"], row["intron_c"])]
    elif event_type == "MSE":
        introns = row["intron"].split(";")
        if len(introns) < 2:
            raise ValueError("MSE requires inclusion and skipping junctions")
        pairs = [(j, introns[-1]) for j in introns[:-1]]
    elif event_type in ("FIVE", "THREE"):
        pairs = [(row["intron_a"], row["intron_b"])]
    elif event_type == "MXE":
        pairs = list(product([row["intron_a1"], row["intron_a2"]],
                             [row["intron_b1"], row["intron_b2"]]))
    elif event_type in ("AFE", "ALE"):
        pairs = list(product(row["intron_a"].split(";"), row["intron_b"].split(";")))
    elif event_type == "RI":
        intron = row["intron_a"]
        chrom, positions = intron.rsplit(":", 1)
        start, end = map(int, positions.split("-"))
        pairs = [(f"{chrom}:{start}-{start + 1}", intron),
                 (f"{chrom}:{end - 1}-{end}", intron)]
    else:
        raise ValueError(f"Unknown event type: {event_type}")
    if any(not a or not b or a == b for a, b in pairs) or len(set(pairs)) != len(pairs):
        raise ValueError(f"Invalid or duplicate components for {row['event_id']}")
    return pairs


def extract_components(event_type, events, junctions, samples):
    """Extract a bounded batch of events from JunctionData, treating absent IDs as zero."""
    counts, definitions = [], []
    columns = [junctions.sample_to_col[s] for s in samples]
    zero = np.zeros(len(samples), dtype=np.int64)
    for row in events.to_dict("records"):
        for index, (inc, exc) in enumerate(component_pairs(event_type, row), 1):
            component = f"c{index}"
            definitions.append([event_type, row["event_id"], component, inc, exc])
            a = junctions.junction_to_idx.get(inc)
            b = junctions.junction_to_idx.get(exc)
            k = zero if a is None else junctions.array[a, columns]
            e = zero if b is None else junctions.array[b, columns]
            if np.any(k + e >= 2**53):
                raise ValueError(f"Component total exceeds exact R integer precision: {row['event_id']}")
            for sample, inclusion, exclusion in zip(samples, k, e):
                counts.append([event_type, row["event_id"], component, sample,
                               int(inclusion), int(exclusion)])
    return pd.DataFrame(counts, columns=COUNT_COLUMNS), pd.DataFrame(definitions, columns=COMPONENT_COLUMNS)
