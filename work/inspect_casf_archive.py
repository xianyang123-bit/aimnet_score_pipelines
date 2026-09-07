#!/usr/bin/env python3
import collections
import pathlib
import sys
import tarfile


archive = pathlib.Path(sys.argv[1])
prefix_counts = {depth: collections.Counter() for depth in range(1, 5)}
coordinate_counts = collections.Counter()
coordinate_sizes = collections.Counter()
coordinate_examples = []
all_examples = []
file_count = 0
total_bytes = 0

coordinate_suffixes = (".mol2", ".mol2.gz", ".sdf", ".sdf.gz", ".pdbqt")

with tarfile.open(archive, "r|gz") as handle:
    for member in handle:
        if not member.isfile():
            continue
        file_count += 1
        total_bytes += member.size
        parts = pathlib.PurePosixPath(member.name).parts
        for depth in prefix_counts:
            prefix_counts[depth]["/".join(parts[:depth])] += 1
        if len(all_examples) < 100:
            all_examples.append(member.name)
        lower = member.name.lower()
        if lower.endswith(coordinate_suffixes):
            key = "/".join(parts[: min(4, len(parts))])
            coordinate_counts[key] += 1
            coordinate_sizes[key] += member.size
            if len(coordinate_examples) < 200:
                coordinate_examples.append(f"{member.name}|{member.size}")

print(f"FILES|{file_count}")
print(f"UNCOMPRESSED_BYTES|{total_bytes}")
for depth in prefix_counts:
    print(f"PREFIX_DEPTH_{depth}")
    for name, count in prefix_counts[depth].most_common(100):
        print(f"{count}|{name}")
print("COORDINATE_GROUPS")
for name, count in coordinate_counts.most_common():
    print(f"{count}|{coordinate_sizes[name]}|{name}")
print("COORDINATE_EXAMPLES")
print("\n".join(coordinate_examples))
print("FIRST_FILES")
print("\n".join(all_examples))
