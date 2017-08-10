#! /usr/bin/env python

import sys
import os
import math

import matplotlib as mpl

# Use TrueType (42) fonts rather than Type 3 fonts
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
tex_font_settings = {
        "text.usetex": True,
        "font.family": "sans-serif",
        "text.latex.preamble" : [
                "\\usepackage[T1]{fontenc}",
                "\\usepackage[cm]{sfmath}",
                ]
}

mpl.rcParams.update(tex_font_settings)

import matplotlib.pyplot as plt
from matplotlib import gridspec

import project_util

def spreadsheet_iter(spreadsheet_path, sep = '\t'):
    with open(spreadsheet_path, "r") as file_stream:
        header = next(file_stream).strip().split(sep)
        for row_idx, row in enumerate(file_stream):
            if row.strip() == '':
                continue
            r = [el.strip() for el in row.strip().split(sep)]
            if len(r) != len(header):
                raise Exception('row {0} of spreadsheet {1} has {2} columns, '
                        'header has {3}'.format(row_idx + 1, spreadsheet_path,
                                len(r), len(header)))
            yield dict(zip(header, r))

def get_dict_from_spreadsheet(spreadsheet_path, sep = '\t'):
    ss_iter = spreadsheet_iter(spreadsheet_path,
            sep = sep)
    row_dict = next(ss_iter)
    d = dict(zip(row_dict.keys(),
            [[row_dict[k]] for k in row_dict.keys()]))
    for row_dict in ss_iter:
        for k in row_dict.keys():
            d[k].append(row_dict[k])
    return d

def get_results():
    path = os.path.join(project_util.RESULTS_DIR,
            "results.txt")
    return get_dict_from_spreadsheet(path)

def plot_marginal_likelihood_results(mls, glm_mls):
    plt.close('all')
    fig = plt.figure(figsize = (5.0, 4.0))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(mls, glm_mls)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 4.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    xlabel_text = ax.set_xlabel("Quadrature log marginal likelihood", size = 14)
    ylabel_text = ax.set_ylabel("ABC-GLM log marginal likelihood", size = 14)

    gs.update(left = 0.13, right = 0.95, bottom = 0.17, top = 0.95)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "marginal-likelihood-results.pdf")
    plt.savefig(plot_path)

def plot_scatter(x, y, x_label, y_label, plot_path):
    mn = min(x + y)
    mx = max(x + y)
    axis_buffer = math.fabs(mx - mn) * 0.05
    axis_min = mn - axis_buffer
    axis_max = mx + axis_buffer
    plt.close('all')
    fig = plt.figure(figsize = (5.0, 4.0))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(x, y)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 4.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    identity_line, = ax.plot(
            [axis_min, axis_max],
            [axis_min, axis_max])
    plt.setp(identity_line,
            color = '0.7',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    xlabel_text = ax.set_xlabel(x_label, size = 14)
    ylabel_text = ax.set_ylabel(y_label, size = 14)

    gs.update(left = 0.13, right = 0.95, bottom = 0.17, top = 0.95)

    plt.savefig(plot_path)

def main_cli(argv = sys.argv):
    d = get_results()
    for k in d:
        d[k] = [float(x) for x in d[k]]
    plot_marginal_likelihood_results(d["trapezoidal_10000_ln_ml"], d["glm_ln_ml"])
    plot_scatter(
            x = d["true_edge_length"],
            y = d["mcmc_mean_edge_length"],
            x_label = "True branch length",
            y_label = "MCMC estimated branch length",
            plot_path = os.path.join(project_util.RESULTS_DIR, "mcmc-branch-length-plot.pdf"))
    plot_scatter(
            x = d["true_edge_length"],
            y = d["glm_mode_edge_length"],
            x_label = "True branch length",
            y_label = "ABC-GLM estimated branch length",
            plot_path = os.path.join(project_util.RESULTS_DIR, "glm-branch-length-plot.pdf"))


if __name__ == "__main__":
    main_cli()
