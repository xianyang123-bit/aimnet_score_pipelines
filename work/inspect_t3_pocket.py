import json
import pickle
import sys

import lmdb


path = sys.argv[1]
env = lmdb.open(path, subdir=False, readonly=True, lock=False, readahead=False)
with env.begin() as txn:
    rows = list(txn.cursor())
env.close()

print(json.dumps({"path": path, "records": len(rows), "keys": [k.decode(errors="replace") for k, _ in rows[:5]]}, indent=2))
if rows:
    row = pickle.loads(rows[0][1])
    summary = {}
    for key, value in row.items():
        item = {"type": type(value).__name__}
        if hasattr(value, "shape"):
            item["shape"] = list(value.shape)
        elif hasattr(value, "__len__") and not isinstance(value, (str, bytes, dict)):
            item["length"] = len(value)
        if key in {"pocket", "smi"}:
            item["value"] = value
        if key in {"pocket_atoms", "atoms"}:
            item["sample"] = list(value[:20])
        summary[key] = item
    print(json.dumps(summary, indent=2, default=str))
