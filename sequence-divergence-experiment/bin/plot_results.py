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


def plot_marginal_likelihood_results(
        correct_model_ml,
        correct_model_glm_ml,
        vague_model_ml,
        vague_model_glm_ml):
    plt.close('all')
    fig = plt.figure(figsize = (6.0, 4.75))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(correct_model_ml, correct_model_glm_ml)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 5.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    xlabel_text = ax.set_xlabel("Quadrature log marginal likelihood of correct model", size = 12)
    ylabel_text = ax.set_ylabel("ABC-GLM log marginal likelihood of correct model", size = 12)

    gs.update(left = 0.1, right = 0.98, bottom = 0.1, top = 0.98)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "correct-model-marginal-likelihoods.pdf")
    plt.savefig(plot_path)

    plt.close('all')
    fig = plt.figure(figsize = (6.0, 4.5))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(vague_model_ml, vague_model_glm_ml)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 5.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    xlabel_text = ax.set_xlabel("Quadrature log marginal likelihood of vague model", size = 12)
    ylabel_text = ax.set_ylabel("ABC-GLM log marginal likelihood of vague model", size = 12)

    gs.update(left = 0.1, right = 0.98, bottom = 0.1, top = 0.98)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "vague-model-marginal-likelihoods.pdf")
    plt.savefig(plot_path)


def plot_bayes_factor_results(quadrature_ln_bf, glm_ln_bf):
    mn = min(quadrature_ln_bf + glm_ln_bf)
    mx = max(quadrature_ln_bf + glm_ln_bf)
    axis_buffer = math.fabs(mx - mn) * 0.05
    axis_min = mn - axis_buffer
    axis_max = mx + axis_buffer

    plt.close('all')
    fig = plt.figure(figsize = (5.25, 4.0))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(quadrature_ln_bf, glm_ln_bf)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 5.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    identity_line, = ax.plot(
            [-1.0, 1.0],
            [-1.0, 1.0])
    plt.setp(identity_line,
            color = '0.7',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    # Draw line associated with failing to penalize for essentially zero
    # likelihood across the second (almost) half of vague prior which is
    # roughly log(1/2); i.e., the upper half of the vague prior (0.1-0.2)
    # should drag down marginal likelihood by about 1/2
    # More precisely it should reduce ML by (0.2 - 0.0001) / (0.1 - 0.0001)
    prior_lower = 1.0 / 10000.0
    vague_penalty = math.log((0.2 - prior_lower) / (0.1 - prior_lower))
    penalty_line, = ax.plot(
            [-1.0,                 1.0],
            [-1.0 - vague_penalty, 1.0 - vague_penalty])
    plt.setp(penalty_line,
            color = '0.7',
            linestyle = '--',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    xlabel_text = ax.set_xlabel("Quadrature log Bayes factor", size = 14)
    ylabel_text = ax.set_ylabel("ABC-GLM log Bayes factor", size = 14)

    gs.update(left = 0.13, right = 0.99, bottom = 0.12, top = 0.99)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "ln-bayes-factors.pdf")
    plt.savefig(plot_path)

    quadrature_bf = [math.exp(lnbf) for lnbf in quadrature_ln_bf]
    glm_bf = [math.exp(lnbf) for lnbf in glm_ln_bf]
    mn = min(quadrature_bf + glm_bf)
    mx = max(quadrature_bf + glm_bf)
    axis_buffer = math.fabs(mx - mn) * 0.05
    axis_min = mn - axis_buffer
    axis_max = mx + axis_buffer

    plt.close('all')
    fig = plt.figure(figsize = (5.25, 4.0))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(quadrature_bf, glm_bf)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 5.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    identity_line, = ax.plot(
            [0.0, 3.0],
            [0.0, 3.0])
    plt.setp(identity_line,
            color = '0.7',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    # Draw line associated with failing to penalize for essentially zero
    # likelihood across the second (almost) half of vague prior which is
    # roughly log(1/2); i.e., the upper half of the vague prior (0.1-0.2)
    # should drag down marginal likelihood by about 1/2
    # More precisely it should reduce ML by (0.2 - 0.0001) / (0.1 - 0.0001)
    prior_lower = 1.0 / 10000.0
    vague_penalty = (0.2 - prior_lower) / (0.1 - prior_lower)
    penalty_line, = ax.plot(
            [0.0,                 3.0],
            [0.0, 3.0 / vague_penalty])
    plt.setp(penalty_line,
            color = '0.7',
            linestyle = '--',
            linewidth = 1.0,
            marker = '',
            zorder = 0)
    xlabel_text = ax.set_xlabel("Quadrature Bayes factor", size = 14)
    ylabel_text = ax.set_ylabel("ABC-GLM Bayes factor", size = 14)

    gs.update(left = 0.13, right = 0.99, bottom = 0.12, top = 0.99)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "bayes-factors.pdf")
    plt.savefig(plot_path)


def plot_edge_estimates(
        true_edge_length,
        correct_model_glm_estimate,
        correct_model_mcmc_estimate,
        vague_model_glm_estimate,
        vague_model_mcmc_estimate):
    mn = min(true_edge_length + correct_model_glm_estimate + \
            correct_model_mcmc_estimate + vague_model_glm_estimate + \
            vague_model_mcmc_estimate)
    mx = max(true_edge_length + correct_model_glm_estimate + \
            correct_model_mcmc_estimate + vague_model_glm_estimate + \
            vague_model_mcmc_estimate)
    axis_buffer = math.fabs(mx - mn) * 0.05
    axis_min = mn - axis_buffer
    axis_max = mx + axis_buffer
    plt.close('all')
    fig = plt.figure(figsize = (7.0, 5.5))
    gs = gridspec.GridSpec(2, 2,
            wspace = 0.0,
            hspace = 0.0)

    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(true_edge_length, correct_model_glm_estimate)
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
    ax.text(0.5, 1.0,
            "ABC-GLM",
            size = 14.0,
            horizontalalignment = "center",
            verticalalignment = "bottom",
            transform = ax.transAxes)

    ax = plt.subplot(gs[0, 1])
    line, = ax.plot(true_edge_length, correct_model_mcmc_estimate)
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
    ax.text(0.5, 1.0,
            "MCMC",
            size = 14.0,
            horizontalalignment = "center",
            verticalalignment = "bottom",
            transform = ax.transAxes)
    ax.text(1.0, 0.5,
            "Correct model",
            size = 14.0,
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = 270.0,
            transform = ax.transAxes)

    ax = plt.subplot(gs[1, 0])
    line, = ax.plot(true_edge_length, vague_model_glm_estimate)
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
    ax.text(-0.15, 1.0,
            "Estimated branch length",
            size = 16.0,
            horizontalalignment = "right",
            verticalalignment = "center",
            rotation = 90.0,
            transform = ax.transAxes)
    ax.text(1.0, -0.13,
            "True branch length",
            size = 16.0,
            horizontalalignment = "center",
            verticalalignment = "top",
            transform = ax.transAxes)

    ax = plt.subplot(gs[1, 1])
    line, = ax.plot(true_edge_length, vague_model_mcmc_estimate)
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
    ax.text(1.0, 0.5,
            "Vague model",
            size = 14.0,
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = 270.0,
            transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    gs.update(left = 0.10, right = 0.97, bottom = 0.10, top = 0.96)

    plot_path = os.path.join(project_util.RESULTS_DIR,
            "branch-length-estimates.pdf")
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
    plot_marginal_likelihood_results(
            correct_model_ml = d["correct_model_trapezoidal_10000_ln_ml"],
            correct_model_glm_ml = d["correct_model_glm_ln_ml"],
            vague_model_ml = d["vague_model_trapezoidal_10000_ln_ml"],
            vague_model_glm_ml = d["vague_model_glm_ln_ml"])
    plot_bayes_factor_results(
            quadrature_ln_bf = d["trapezoidal_10000_ln_bf"],
            glm_ln_bf = d["glm_ln_bf"])
    plot_edge_estimates(
            true_edge_length = d["true_edge_length"],
            correct_model_glm_estimate = d["correct_model_glm_mode_edge_length"],
            correct_model_mcmc_estimate = d["correct_model_mcmc_mean_edge_length"],
            vague_model_glm_estimate = d["vague_model_glm_mode_edge_length"],
            vague_model_mcmc_estimate = d["vague_model_mcmc_mean_edge_length"])


if __name__ == "__main__":
    main_cli()
