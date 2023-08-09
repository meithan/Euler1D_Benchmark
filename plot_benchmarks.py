# Plots benchmarks of Solver1 in C, Fortran and Python

import sys
import matplotlib.pyplot as plt
import numpy as np

data = {}
with open("benchmarks.csv") as f:
  for line in f:
    if line.startswith("#"): continue
    tokens =line.strip().split(",")
    lang = tokens[0]
    data[lang] = [float(x) for x in tokens[1:]]
langs = list(data.keys())

# Compute averages and ratios
num_langs = len(langs)
num_runs = len(data[langs[0]])
results = []
for l in langs:
  mean = sum(data[l])/len(data[l])
  minv = min(data[l])
  maxv = max(data[l])
  results.append([l, 0, mean, minv, maxv, 0, 0])
results.sort(key=lambda x: x[2])
ref_lang = "Fortran"
ref_idx = [x[0] for x in results].index(ref_lang)
ref_value = results[ref_idx][2]
for i in range(len(results)):
  results[i][1] = results[i][2] / ref_value
  results[i][5] = results[i][3] / ref_value
  results[i][6] = results[i][4] / ref_value

colors = {
  "Fortran": "C0",
  "C/C++": "C1",
  "Rust": "C2",
  "Java": "C4",
  "Python3.11+Numpy": "C5",
  "Python3.8+Numpy": "C6",
  "Julia": "C7",
  "Julia+LoopVec": "C3",
  "Python3.11": "C8",
  "Python3.8": "C9",  
}

title = "Euler1D, Sod shock tube, NX=5000, Lax-Friedrichs solver"

# ----------------------------
# Figure 1: linear scale, Python separated

results1 = [r for r in results if not "Python" in r[0] or "Numpy" in r[0]]
results2 = [r for r in results if r not in results1]
len1 = len(results1)
len2 = len(results2)

fig, axs = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [len1, len2]})

for subnum, _results in enumerate([results1, results2]):

  plt.sca(axs[subnum])

  xs = range(len(_results))
  langs, ratios, means, minvs, maxvs, minrs, maxrs = zip(*_results)
  labels = langs
  labels = [x.replace("+Numpy", "\n+Numpy") for x in labels]
  labels = [x.replace("Julia+", "Julia\n+") for x in labels]
  yerrs = np.abs(np.array([minrs, maxrs]) - ratios)
  plt.bar(xs, ratios, yerr=yerrs, capsize=3, width=0.5, tick_label=labels, log=False, color=[colors[l] for l in langs])
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
    plt.ylabel(f"Execution time ratio relative to {ref_lang}")
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
langs, ratios, means, minvs, maxvs, minrs, maxrs = zip(*results)
labels = langs
labels = [x.replace("+Numpy", "\n+Numpy") for x in labels]
labels = [x.replace("Julia+", "Julia\n+") for x in labels]
yerrs = np.abs(np.array([minrs, maxrs]) - ratios)
plt.bar(xs, ratios, yerr=yerrs, capsize=3, width=0.5, tick_label=labels, log=False, color=[colors[l] for l in langs])
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

plt.ylabel(f"Execution time ratio relative to {ref_lang}")
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