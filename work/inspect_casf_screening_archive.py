import sys
import tarfile
from collections import Counter, defaultdict

archive = sys.argv[1]
prefix = "CASF-2016/decoys_screening/"
targets = Counter()
examples = defaultdict(list)

with tarfile.open(archive, "r:gz") as tf:
    for member in tf:
        name = member.name
        if not name.startswith(prefix) or not member.isfile():
            continue
        rest = name[len(prefix):]
        parts = rest.split("/")
        if len(parts) < 2:
            continue
        target = parts[0]
        targets[target] += 1
        if len(examples[target]) < 15:
            examples[target].append(rest)

print("targets", len(targets))
for target in sorted(targets)[:12]:
    print(target, targets[target])
    for name in examples[target]:
        print(" ", name)
