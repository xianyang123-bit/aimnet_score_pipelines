#!/usr/bin/env python3
import json
from collections import Counter, defaultdict

import pandas as pd

casf_path = "/home/xianyang/LigUnity-official/test_datasets/casf_label_seq.json"
t3_path = "/data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset/targets.csv.gz"

casf = json.load(open(casf_path))
t3 = pd.read_csv(t3_path)
layer_by_up = dict(zip(t3["uniprot"], t3["layer"]))

complex_counts = Counter()
pdbs = defaultdict(list)
for record in casf:
    up = record.get("uniprot")
    layer = layer_by_up.get(up, "not_in_T3")
    complex_counts[layer] += 1
    pdbs[layer].append(record["pockets"][0])

unique_casf = {record.get("uniprot") for record in casf}
target_counts = Counter(layer_by_up.get(up, "not_in_T3") for up in unique_casf)
result = {
    "casf_complexes": len(casf),
    "casf_unique_uniprot": len(unique_casf),
    "complex_counts_by_t3_layer": dict(complex_counts),
    "target_counts_by_t3_layer": dict(target_counts),
    "pdbs_by_t3_layer": dict(pdbs),
    "t3_target_counts": t3.groupby("layer").size().astype(int).to_dict(),
    "t3_active_counts": t3.groupby("layer")["n_actives"].sum().astype(int).to_dict(),
}
print(json.dumps(result, indent=2))
with open("/home/xianyang/aimnet2_score_pipelines/work/casf_t3_overlap.json", "w") as handle:
    json.dump(result, handle, indent=2)
