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

# Compute averages and ratios
num_langs = len(langs)
num_runs = len(data[langs[0]])
results = []
for l in langs:
  mean = sum(data[l])/len(data[l])
  minv = min(data[l])
  maxv = max(data[l])
  results.append([l, 0, mean, minv, maxv])
results.sort(key=lambda x: x[2])
ref_idx = 1
ref_value = results[ref_idx][2]
for i in range(len(results)):
  results[i][1] = results[i][2] / ref_value

colors = ["C0", "C1", "C2", "C4", "C5", "C6", "C7", "C8", "C9"]
title = "Euler1D, Sod shock tube, NX=5000, Lax-Friedrichs solver"

# ----------------------------
# Figure 1: linear scale, Python separated

results1 = [r for r in results if not "Python" in r[0] or "numpy" in r[0]]
results2 = [r for r in results if r not in results1]
len1 = len(results1)
len2 = len(results2)

fig, axs = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [len1, len2]})

for subnum, _results in enumerate([results1, results2]):

  plt.sca(axs[subnum])

  xs = range(len(_results))
  labels, ratios, means, mins, maxs = zip(*_results)
  labels = [x.replace("+numpy", "\n+numpy") for x in labels]
  yerrs = np.abs(np.array([mins, maxs]) - means)
  plt.bar(xs, ratios, yerr=yerrs, capsize=3, width=0.5, tick_label=labels, log=False, color=colors[subnum*len1:(subnum+1)*len1])
  for i in xs:
    x = xs[i]
    ratio = ratios[i]
    if ratio >= 100:
      text = f"{ratio:.1f}"
    else:
      text = f"{ratio:.2f}"
    y = ratio + yerrs[1][i]
    plt.annotate(text, xy=(x, y), xytext=(0, 5), textcoords="offset pixels", color="k", ha="center", fontsize=14)

  plt.grid(ls=":", axis="y") 
  plt.xlim(-0.5, len(_results)-0.5)
  ylims = plt.ylim()
  plt.ylim(ylims[0], ylims[1]*1.02)
  
  if subnum == 0:
    plt.ylabel("Execution time ratio relative to C/C++")
    notes = [
      "Erorr bars show min/max of 10 runs"
    ]
    plt.annotate("\n".join(notes), xy=(0.01, 0.99), ha="left", va="top", xycoords="axes fraction", fontsize=9)
  else:
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()

plt.suptitle(title + " (lin scale)", y=0.955)

plt.tight_layout()
# plt.subplots_adjust(wspace=0)

if "--save" in sys.argv:
  out_fname = "benchmark_lin.png"
  plt.savefig(out_fname)
else:
  plt.show()

# ----------------------------
# Figure 2: logarithmic scale

plt.figure(figsize=(9,5))

xs = range(len(results))
labels, ratios, means, mins, maxs = zip(*results)
print(labels)
labels = [x.replace("+numpy", "\n+numpy") for x in labels]
yerrs = np.abs(np.array([mins, maxs]) - means)
plt.bar(xs, ratios, yerr=yerrs, capsize=3, width=0.5, tick_label=labels, log=False, color=colors[:len(results)])
for i in xs:
  x = xs[i]
  ratio = ratios[i]
  if ratio >= 100:
    text = f"{ratio:.1f}"
  else:
    text = f"{ratio:.2f}"
  y = ratio + yerrs[1][i]
  plt.annotate(text, xy=(x, y), xytext=(0, 5), textcoords="offset pixels", color="k", ha="center", fontsize=14)

plt.grid(ls=":", axis="y") 
plt.xlim(-0.5, len(results)-0.5)
plt.yscale("log")
ylims = plt.ylim()
plt.ylim(ylims[0], ylims[1]*1.2)

plt.ylabel("Execution time ratio relative to C/C++")
notes = [
  "Erorr bars show min/max of 10 runs"
]
plt.annotate("\n".join(notes), xy=(0.01, 0.99), ha="left", va="top", xycoords="axes fraction", fontsize=9)

plt.title(title + " (log scale)")

plt.tight_layout()
# plt.subplots_adjust(wspace=0)

if "--save" in sys.argv:
  out_fname = "benchmark_log.png"
  plt.savefig(out_fname)
else:
  plt.show()