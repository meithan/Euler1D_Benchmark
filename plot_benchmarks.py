# Plots benchmarks of Solver1 in C, Fortran and Python

import sys
import matplotlib.pyplot as plt

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

num_runs = len(data[langs[0]])

use_log_scale = "--log" in sys.argv

plt.figure(figsize=(8,5))

xs = range(len(langs))
avgs = [sum(data[l])/len(data[l]) for l in langs]
colors = ["C0", "C1", "C2", "C4", "C5", "C6", "C7", "C8"]
labels = [l.replace('+n', "\n+n") for l in langs]
plt.bar(xs, avgs, width=0.4, tick_label=labels, log=False, color=colors)
ref_value = avgs[0]
for x, value, color in zip(xs, avgs, colors):
  text = f"{value/ref_value:.2f}x"
  plt.annotate(text, xy=(x, value), xytext=(0, 5), textcoords="offset pixels", color="k", ha="center", fontsize=14)

plt.ylabel("Execution time [s]")
plt.grid(ls=":", axis="y")
title = "Euler1D, Sod shock tube, NX=5000, Lax-Friedrichs, {} runs avg".format(num_runs)
if use_log_scale:
  plt.yscale("log")
  title +=" (log scale)"
  y1,y2 = plt.ylim()
  plt.ylim(y1,y2*1.2)
else:
  title +=" (lin scale)"
plt.title(title)

plt.annotate("Performance ratios vs C", xy=(0.99, 0.99), ha="right", va="top", xycoords="axes fraction", fontsize=9)

plt.tight_layout()

if "--save" in sys.argv:
  if use_log_scale:
    out_fname = "benchmark_log.png"
  else:
    out_fname = "benchmark_lin.png"
  plt.savefig(out_fname)
else:
  plt.show()