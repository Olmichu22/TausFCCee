import ROOT 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pandas as pd
import os
import textwrap

import matplotlib.pyplot as plt
PI0 = "π⁰"
PI = "π"
MU = "μ"
E = "e"
N = "n"
NEUTRINO = "ν"
TAU = "τ"
GAMMA = "γ"
RHO = "ρ"
A1 = r"$a_1$"

def recodify_key(key, photon=False):
  if photon:
    if key < 0:
      if key == -13:
        new_key = f"{MU}"
      elif key == -11:
        new_key = f"{E}"
      elif key <= -20:
        new_key = f"h{N}"
      else:
        new_key = "Unknown"
    elif key == 0:
      new_key = f"h"
    elif key < 10:
      new_key = f"h{2*key}{GAMMA}"
    elif key == 10:
      new_key = f"3h"
    else:
      new_key = f"3h{2*(key-10)}{GAMMA}"
  else:
    if key < 0:
      if key == -13:
        new_key = f"{TAU} \\rightarrow {MU}2{NEUTRINO}"
      elif key == -11:
        new_key = f"{TAU} \\rightarrow {E}2{NEUTRINO}"
      elif key <= -20:
        new_key = f"{PI}{N}"
      else:
        new_key = "Unknown"
    elif key == 0:
      new_key = f"{PI}"
    elif key < 10:
        new_key = f"{PI}{key}{PI0}"
    elif key == 10:
      new_key = f"3{PI}"
    else:
      new_key = f"{3}{PI}{key-10}{PI0}"
  return new_key

def plot_absolute(graphs, absolute_keys, colors, xaxis,
                             min_absolute_value, max_absolute_value,
                             title, outputpath, mig_str):
    fig, ax = plt.subplots()
    # fig.set_size_inches(8, 6)
    markers = ['v', 'x', '+', '*', "o"]

    prec_keys = [key.split("->")[0] for key in absolute_keys]
    consec_keys = [key.split("->")[1] for key in absolute_keys]

    prec_keys_set = list(set(prec_keys))
    consec_keys_set = list(set(consec_keys))
    prec_keys_set.sort()
    consec_keys_set.sort()
    # Color for prec_keys_set
    prec_colors = {key: colors[i % len(colors)] for i, key in enumerate(prec_keys_set)}
    # Linestyle for consec_keys_set
    consec_markers = {key: markers[i % len(markers)] for i, key in enumerate(consec_keys_set)}

    for i, key in enumerate(absolute_keys):
        # print(graphs[key])
        x, y = graphs[key]['x'], graphs[key]['y']

        id_0, id_1 = key.split("->")
        id_0_rec = recodify_key(int(id_0))
        id_1_rec = recodify_key(int(id_1), photon=True)
        label = f"{id_0_rec} → {id_1_rec}"
        color = colors[i % len(colors)]

        ax.plot(
            x,
            y,
            marker=consec_markers[id_1],
            linestyle="-",
            color=prec_colors[id_0],
            label=label,
            linewidth=0.8,
            markersize=5,
        )
        if i == 0:
            ax.set_xlabel(xaxis)
            ax.set_ylabel("Counts")
            ax.set_ylim(0.7 * min_absolute_value if min_absolute_value!=0 else -200, 1.1 * max_absolute_value)
            
            ax.set_yticks(np.arange(0, max_absolute_value + 250, 250))
    # Horizontal grid lines
    ax.minorticks_on()
    ax.yaxis.grid(True, linestyle='-', alpha=0.5)

    ax.legend(ncol=2, loc='upper left', bbox_to_anchor=(0.6, 1.15))
    # ax.legend(ncol=2, loc="best")
    wrapped_title = "\n".join(textwrap.wrap(title, width=20))
    ax.set_title(wrapped_title, fontsize=11, loc='left', pad=20, fontweight='bold',)

    fig.tight_layout()
    fig.savefig(outputpath + f"graphs_plot_{mig_str}.png", bbox_inches='tight')
    plt.close(fig)


def plot_metric(graphs, metric_keys, colors, xaxis,
                           title, outputpath, metric, mig_str):
    fig, ax = plt.subplots()
    markers = ['v', 'x', '+', '*', "o"]
    line_markers = ["solid", "dotted", "dashed", "dashdot", (0, (1, 10)), (0, (3, 5, 1, 5))]
    
    prec_keys = [key.split("->")[0] for key in metric_keys]
    consec_keys = [key.split("->")[1].split(" ")[0] for key in metric_keys]

    prec_keys_set = list(set(prec_keys))
    consec_keys_set = list(set(consec_keys))
    prec_keys_set.sort()
    consec_keys_set.sort()
    # Color for prec_keys_set
    prec_colors = {key: colors[i % len(colors)] for i, key in enumerate(prec_keys_set)}
    # Linestyle for consec_keys_set
    consec_markers = {key: markers[i % len(markers)] for i, key in enumerate(consec_keys_set)}

    consec_line_markers = {key: line_markers[i % len(line_markers)] for i, key in enumerate(consec_keys_set)}
    
    # Check number of events
    key_0 = list(graphs.keys())[0]
    n_els = len(graphs[key_0]['y'])

    for i, key in enumerate(metric_keys):

        # color = colors[i % len(colors)]
        x, y = graphs[key]['x'], graphs[key]['y']
        id_0, id_1 = key.split("->")
        id_0_rec = recodify_key(int(id_0))
        if metric in id_1 or "recall" in id_1:
            id_1 = id_1.split(" ")[0]
        id_1_rec = recodify_key(int(id_1), photon=True)
        key = f"{id_0_rec} → {id_1_rec}"
        if n_els <= 25:
          ax.plot(
              x,
              y,
              marker=consec_markers[id_1],
              linestyle="-",
              color=prec_colors[id_0],
              label=key,
              linewidth=0.8,
              markersize=5,
          )
        else:
          ax.plot(
              x,
              y,
              linestyle=consec_line_markers[id_1],
              color=prec_colors[id_0],
              label=key,
              linewidth=2,
              markersize=5,
          )
          
        if i == 0:
            ax.set_xlabel(xaxis, fontsize=12)
            ax.set_ylabel("Cuenta normalizada", fontsize=12)
            ax.set_ylim(-0.1, 1.1)
            ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.minorticks_on()
    ax.yaxis.grid(True, linestyle='-', alpha=0.5)
    
    # fig.subplots_adjust(top=5)  # Ajusta el margen superior para dar espacio a la leyenda
    ax.legend(ncol=2, loc='upper left', bbox_to_anchor=(0.55, 1.15), fontsize=12)
    # ax.legend(ncol=2, loc="best")
    # title = title + f"\n({metric.capitalize()})"
    wrapped_title = "\n".join(textwrap.wrap(title, width=30))
    ax.set_title(wrapped_title, fontweight='bold',fontsize=15, loc='left', pad=20)

    #aumentamos tamaño de figura
    fig.set_size_inches(8, 5)
    fig.tight_layout()
    fig.savefig(outputpath + f"graphs_plot_{metric}_{mig_str}.png", bbox_inches='tight')
    plt.close(fig)


def plot_accuracy(graphs, colors, xaxis,
                           title, outputpath):
  fig, ax = plt.subplots()

  for i, key in enumerate(graphs.keys()):
    x = graphs[key]['x']
    y = graphs[key]['y']
    ax.plot(x, y, marker='o', linestyle='-', color=colors[i],
            label=key, linewidth=0.8, markersize=5)
  
  ax.set_xlabel(xaxis)
  ax.set_ylabel("Accuracy")
  ax.set_ylim(0, 1.1)
  ax.set_yticks(np.arange(0, 1.1, 0.1))
  ax.minorticks_on()
  ax.yaxis.grid(True, linestyle='-', alpha=0.5)
  
  # Set title
  ax.set_title(title, fontsize=11, loc='left', pad=20, fontweight='bold')
  ax.legend(ncol=2, loc='upper left', bbox_to_anchor=(0.6, 1.15))
  fig.tight_layout()
  fig.savefig(outputpath + f"graphs_plot_accuracy.png", bbox_inches='tight')
  plt.close(fig)