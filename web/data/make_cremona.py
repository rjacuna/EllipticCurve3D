#!/usr/bin/env python3
"""Build data/cremona/ from John Cremona's tables (https://github.com/JohnCremona/ecdata).

  python3 make_cremona.py /path/to/ecdata        # a checkout, or any directory holding allcurves/ and alllabels/

Writes one JSON shard per range of 10,000 conductors, <lo>-<hi>.json, plus index.json listing them. A shard
maps a conductor (as a string) to its isogeny classes in Cremona's order,

    [cremona_letter, lmfdb_letter, [[lmfdb_number, a1, a2, a3, a4, a6], ...]]

with the curves in Cremona's order inside the class, so that Cremona's number is the position + 1 and the
LMFDB number is stored. Both label systems are therefore available offline for every curve of conductor
up to 500,000. Sage's bundled database stops at 9,999; these files come straight from the tables.
"""
import glob, json, os, sys, time

if len(sys.argv) != 2:
    sys.exit(__doc__)
src = sys.argv[1]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cremona")
os.makedirs(out, exist_ok=True)
files = sorted(glob.glob(os.path.join(src, "allcurves", "allcurves.*")), key=lambda f: int(f.rsplit(".", 1)[1].split("-")[0]))
if not files:
    sys.exit(f"no allcurves/allcurves.* files under {src}")
index = {"source": "https://github.com/JohnCremona/ecdata (allcurves, alllabels)", "built": time.strftime("%Y-%m-%d"),
         "shards": [], "curves": 0, "max_conductor": 0}
for f in files:
    rng = f.rsplit(".", 1)[1]
    lo, hi = (int(v) for v in rng.split("-"))
    labels = {}                                              # (N, cremona letter, cremona number) -> (lmfdb letter, lmfdb number)
    with open(os.path.join(src, "alllabels", "alllabels." + rng)) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 6:
                labels[(p[0], p[1], p[2])] = (p[4], int(p[5]))
    data = {}                                                # N -> {cremona letter: [lmfdb letter, [(num, [lmfdb num, a1..a6]), ...]]}
    n = 0
    with open(f) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 4:
                continue
            N, cl, num, ainvs = p[0], p[1], int(p[2]), json.loads(p[3])
            ll, ln = labels[(N, cl, p[2])]
            entry = data.setdefault(N, {}).setdefault(cl, [ll, []])
            entry[1].append((num, [ln] + ainvs))
            n += 1
    shard = {N: [[cl, e[0], [c for _, c in sorted(e[1])]] for cl, e in classes.items()] for N, classes in data.items()}
    name = f"{lo:06d}-{hi:06d}.json"
    with open(os.path.join(out, name), "w") as fh:
        json.dump(shard, fh, separators=(",", ":"))
    size = os.path.getsize(os.path.join(out, name))
    index["shards"].append({"lo": lo, "hi": hi, "file": name, "curves": n, "bytes": size})
    index["curves"] += n
    index["max_conductor"] = max(index["max_conductor"], hi)
    print(f"{name}: {n:7d} curves, {size / 1e6:5.1f} MB")
with open(os.path.join(out, "index.json"), "w") as fh:
    json.dump(index, fh, indent=1)
print(f"total: {index['curves']} curves up to conductor {index['max_conductor']}, {sum(s['bytes'] for s in index['shards']) / 1e6:.0f} MB")
