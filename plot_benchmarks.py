# Plots benchmarks of Solver1 in C, Fortran and Python

import sys
import matplotlib.pyplot as plt
import numpy as np

with open("benchmarks.csv") as f:
  for i in range(6):
    f.readline()
  while True:
    line = f.readline()
    if not line.startswith("#"):
      break
  langs = [l.strip('"') for l in line.strip().split(",")]
  data = {l:[] for l in langs}
  for line in f:
    tokens = [float(x) for x in line.strip().split(",")]
    for i in range(len(tokens)):
      data[langs[i]].append(tokens[i])

# Compute averages
num_langs = len(langs)
num_runs = len(data[langs[0]])
results = []
for l in langs:
  mean = sum(data[l])/len(data[l])
  minv = min(data[l])
  maxv = max(data[l])
  results.append((l, mean, minv, maxv))
results.sort(key=lambda x: x[1])
ref_idx = 1
ref_value = results[ref_idx][1]

use_log_scale = "--log" in sys.argv

plt.figure(figsize=(9,5))

xs = range(len(langs))
labels, means, mins, maxs = zip(*results)
labels = [x.replace("+numpy", "\n+numpy") for x in labels]
colors = ["C0", "C1", "C2", "C4", "C5", "C6", "C7", "C8", "C9"]
yerrs = np.abs(np.array([mins, maxs]) - means)
plt.bar(xs, means, yerr=yerrs, capsize=3, width=0.5, tick_label=labels, log=False, color=colors)
for i in xs:
  x = xs[i]
  value = means[i]
  ratio = value/ref_value
  if ratio >= 100:
    text = f"{ratio:.1f}"
  else:
    text = f"{ratio:.2f}"
  y = value + yerrs[1][i]
  plt.annotate(text, xy=(x, y), xytext=(0, 5), textcoords="offset pixels", color="k", ha="center", fontsize=14)

plt.ylabel("Execution time [s]")
plt.grid(ls=":", axis="y")
title = "Euler1D, Sod shock tube, NX=5000, Lax-Friedrichs solver, {} runs avg".format(num_runs)
if use_log_scale:
  plt.yscale("log")
  title +=" (log scale)"
  y1,y2 = plt.ylim()
  plt.ylim(y1,y2*1.2)
else:
  title +=" (lin scale)"
plt.title(title)

notes = [
  "Numbers above bars are ratios relative to C/C++",
  "Erorr bars show min/max of 10 runs"
]
plt.annotate("\n".join(notes), xy=(0.01, 0.99), ha="left", va="top", xycoords="axes fraction", fontsize=9)

plt.tight_layout()

if "--save" in sys.argv:
  if use_log_scale:
    out_fname = "benchmark_log.png"
  else:
    out_fname = "benchmark_lin.png"
  plt.savefig(out_fname)
else:
  plt.show()