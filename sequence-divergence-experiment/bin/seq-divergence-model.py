#! /usr/bin/env python

import sys
import os
import math
import random
import datetime
import argparse
import itertools
import re
import logging
import unittest
import subprocess

import dendropy

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
_LOG = logging.getLogger(os.path.basename(__file__))


number_pattern_str = r"[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?"
glm_density_pattern_str = r"marginal\s+density:\s+(?P<ml>" + number_pattern_str + ")"
glm_density_pattern = re.compile(glm_density_pattern_str, re.VERBOSE)


def count_differences(seq1, seq2):
    n_diffs = 0
    for i in range(len(seq1)):
        if seq1[i].index != seq2[i].index:
            n_diffs += 1
    return n_diffs


class WindowOperator(object):
    def __init__(self,
            window_size = 0.1,
            auto_optimize = True,
            auto_optimize_delay = 1000):
        self.window_size = window_size
        self.auto_optimize_delay = auto_optimize_delay
        self.number_of_proposals = 0
        self.number_of_acceptions = 0
        self.number_of_proposals_for_tuning = 0
        self.target_acceptance_probability = 0.44
        self.auto_optimize = auto_optimize

    def get_acceptance_rate(self):
        if self.number_of_proposals < 1:
            return None
        return self.number_of_acceptions / float(self.number_of_proposals)

    def set_window_size(self, window_size):
        assert window_size > 0.0
        self.window_size = window_size

    def get_window_addend(self, rng):
        return (rng.random() * 2 * self.window_size) - self.window_size

    def propose(self, value, rng):
        value += self.get_window_addend(rng)
        return value, 0.0

    def reject(self):
        self.number_of_proposals += 1
        if self.number_of_proposals >= self.auto_optimize_delay:
            self.number_of_proposals_for_tuning += 1

    def accept(self):
        self.number_of_proposals += 1
        self.number_of_acceptions += 1
        if self.number_of_proposals >= self.auto_optimize_delay:
            self.number_of_proposals_for_tuning += 1

    def optimize(self, log_alpha):
        if (not self.auto_optimize) or (
                self.number_of_proposals < self.auto_optimize_delay):
            return
        delta = self.calc_delta(log_alpha)
        delta += math.log(self.window_size)
        self.set_window_size(math.exp(delta))

    def calc_delta(self, log_alpha):
        if (not self.auto_optimize) or (
                self.number_of_proposals < self.auto_optimize_delay):
            return 0.0
        target = self.target_acceptance_probability
        count = self.number_of_proposals_for_tuning + 1.0
        delta_p = ((1.0 / count) * (math.exp(min(log_alpha, 0.0)) - target))
        mx = sys.float_info.max
        if ((delta_p > -mx) and (delta_p < mx)):
            return delta_p
        return 0.0

class ScaleOperator(object):
    def __init__(self,
            scale = 0.1,
            auto_optimize = True,
            auto_optimize_delay = 1000):
        self.scale = scale
        self.auto_optimize_delay = auto_optimize_delay
        self.number_of_proposals = 0
        self.number_of_acceptions = 0
        self.number_of_proposals_for_tuning = 0
        self.target_acceptance_probability = 0.44
        self.auto_optimize = auto_optimize

    def get_acceptance_rate(self):
        if self.number_of_proposals < 1:
            return None
        return self.number_of_acceptions / float(self.number_of_proposals)

    def set_scale(self, scale):
        assert scale > 0.0
        self.scale = scale

    def get_multiplier(self, rng):
        return math.exp(self.scale * ((2.0 * rng.random()) - 1.0))

    def propose(self, value, rng):
        m = self.get_multiplier(rng)
        value *= m
        return value, math.log(m)

    def reject(self):
        self.number_of_proposals += 1
        if self.number_of_proposals >= self.auto_optimize_delay:
            self.number_of_proposals_for_tuning += 1

    def accept(self):
        self.number_of_proposals += 1
        self.number_of_acceptions += 1
        if self.number_of_proposals >= self.auto_optimize_delay:
            self.number_of_proposals_for_tuning += 1

    def optimize(self, log_alpha):
        if (not self.auto_optimize) or (
                self.number_of_proposals < self.auto_optimize_delay):
            return
        delta = self.calc_delta(log_alpha)
        delta += math.log(self.scale)
        self.set_scale(math.exp(delta))

    def calc_delta(self, log_alpha):
        if (not self.auto_optimize) or (
                self.number_of_proposals < self.auto_optimize_delay):
            return 0.0
        target = self.target_acceptance_probability
        count = self.number_of_proposals_for_tuning + 1.0
        delta_p = ((1.0 / count) * (math.exp(min(log_alpha, 0.0)) - target))
        mx = sys.float_info.max
        if ((delta_p > -mx) and (delta_p < mx)):
            return delta_p
        return 0.0


class EdgeModel(object):
    """
    A class for handling a model with an edge between two sequences that evolve
    via a Jukes-Cantor (1969) model.

    >>> m = EdgeModel(seed = 111, prior_lower = 0.1, prior_upper = 0.2, seq_length = 100)
    >>> m.power = 0.0
    >>> mcmc_sample = m.mcmc(1000000, 10)
    >>> len(mcmc_sample)
    100001
    >>> edge_sample = [e for g, l, e in mcmc_sample[1:]]
    >>> len(edge_sample)
    100000
    >>> prior_mean = sum(edge_sample) / len(edge_sample)
    >>> abs(prior_mean - 0.15) < 0.0005
    True
    >>> m = EdgeModel(seed = 111, edge_length = 0.03, seq_length = 10000)
    >>> mcmc_sample = m.mcmc(10000, 10)
    >>> len(mcmc_sample)
    1001
    >>> edge_sample = [e for g, l, e in mcmc_sample[1:]]
    >>> len(edge_sample)
    1000
    >>> post_mean = sum(edge_sample) / len(edge_sample)
    >>> abs(post_mean - 0.03) < 0.005
    True
    """

    def __init__(self,
            seed,
            seq1 = None,
            seq2 = None,
            prior_lower = 0.000001,
            prior_upper = 0.1,
            edge_length = None,
            seq_length = 1000):
        self.ctmc = dendropy.model.discrete.Jc69()
        self.mutation_rate = 1.0
        self.power = 1.0
        self.seed = seed
        self.rng = random.Random()
        self.rng.seed(self.seed)
        self.prior_lower = prior_lower
        self.prior_upper = prior_upper
        self.edge_length = edge_length
        if self.edge_length is None:
            self.edge_length = self.draw_edge_from_prior()
        if not seq1:
            self.seq_length = seq_length
            self.seq1, self.seq2 = self.simulate_sequences()
        else:
            assert len(seq1) == len(seq2)
            self.seq_length = len(seq1)
            self.seq1 = [self.ctmc.state_alphabet[x] for x in seq1]
            self.seq2 = [self.ctmc.state_alphabet[x] for x in seq2]
        self.stored_edge_length = None
        self.lnl = None
        self.power_lnl = None
        self.stored_lnl = None
        self.stored_power_lnl = None
        self.window_operator = WindowOperator(
                window_size = (self.prior_upper - self.prior_lower) / 10.0,
                auto_optimize = True)
        self.scale_operator = ScaleOperator(
                scale = 0.2,
                auto_optimize = True)
        self.calc_ln_likelihood()
    
    def draw_edge_from_prior(self):
        return self.rng.uniform(self.prior_lower, self.prior_upper)

    def draw_sequences_from_prior(self):
        edge_len = self.draw_edge_from_prior()
        return self.simulate_sequences(edge_len)

    def get_edge_prior_probability(self):
        return 1.0 / (self.prior_upper - self.prior_lower)

    def get_edge_ln_prior_probability(self):
        return math.log(self.get_edge_prior_probability())

    def simulate_sequences(self, edge_length = None):
        edge_len = edge_length
        if edge_len is None:
            edge_len = self.edge_length
        sequence1 = self.ctmc.stationary_sample(
                seq_len = self.seq_length,
                rng = self.rng)
        sequence2 = self.ctmc.simulate_descendant_states(
                ancestral_states = sequence1,
                edge_length = edge_len,
                mutation_rate = self.mutation_rate,
                rng = self.rng)
        return sequence1, sequence2

    def get_proportion_of_variable_sites(self):
        ndiffs = count_differences(self.seq1, self.seq2)
        return ndiffs / float(self.seq_length)

    def calc_ln_likelihood(self, edge_length = None):
        edge_len = edge_length
        if edge_len is None:
            edge_len = self.edge_length
        self.lnl = 0.0
        # pm = self.ctmc.pmatrix(tlen = self.edge_length,
        #         rate = self.mutation_rate)
        lnp_same = math.log(self.ctmc.pij(
                state_i = 0,
                state_j = 0,
                tlen = edge_len,
                rate = self.mutation_rate))
        lnp_diff = math.log(self.ctmc.pij(
                state_i = 0,
                state_j = 1,
                tlen = edge_len,
                rate = self.mutation_rate))
        for i in range(self.seq_length):
            if self.seq1[i].index == self.seq2[i].index:
                self.lnl += lnp_same
            else:
                self.lnl += lnp_diff
            # lnl += math.log(pm[self.seq1[i].index][self.seq2[i].index])
            # lnl += math.log(self.ctmc.pij(
            #         state_i = self.seq1[i].index,
            #         state_j = self.seq2[i].index,
            #         tlen = self.edge_length,
            #         rate = self.mutation_rate))
        self.power_lnl = self.lnl * self.power
        return self.lnl

    def calc_power_ln_likelihood(self):
        if self.power == 0.0:
            self.power_lnl = 0.0
            return self.power_lnl
        self.calc_ln_likelihood()
        self.power_lnl *= self.power
        return self.power_lnl

    def store_model_state(self):
        self.stored_lnl = self.lnl
        self.stored_power_lnl = self.power_lnl
        self.stored_edge_length = self.edge_length

    def restore_model_state(self):
        self.lnl = self.stored_lnl
        self.power_lnl = self.stored_power_lnl
        self.edge_length = self.stored_edge_length

    def perform_window_update(self):
        e = self.edge_length
        new_e, hastings_ratio = self.window_operator.propose(e, self.rng)
        if (new_e < self.prior_lower) or (new_e > self.prior_upper):
            self.window_operator.reject()
            return
        self.store_model_state()
        self.edge_length = new_e
        self.calc_power_ln_likelihood()
        lnl_ratio = self.power_lnl - self.stored_power_lnl
        acceptance_prob = lnl_ratio + hastings_ratio
        u = self.rng.random()
        if ((acceptance_prob >= 0.0) or (u < math.exp(acceptance_prob))):
            self.window_operator.accept()
        else:
            self.window_operator.reject()
            self.restore_model_state()
        self.window_operator.optimize(acceptance_prob)

    def perform_scale_update(self):
        e = self.edge_length
        new_e, hastings_ratio = self.scale_operator.propose(e, self.rng)
        if (new_e < self.prior_lower) or (new_e > self.prior_upper):
            self.scale_operator.reject()
            return
        self.store_model_state()
        self.edge_length = new_e
        self.calc_power_ln_likelihood()
        lnl_ratio = self.power_lnl - self.stored_power_lnl
        acceptance_prob = lnl_ratio + hastings_ratio
        u = self.rng.random()
        if u < math.exp(acceptance_prob):
            self.scale_operator.accept()
        else:
            self.scale_operator.reject()
            self.restore_model_state()
        self.scale_operator.optimize(acceptance_prob)

    def mcmc(self, number_of_generations, sample_frequency,
            screen_frequency = None):
        window_size = (self.prior_upper - self.prior_lower) / 10.0
        scale = 0.2
        if self.seq_length < 2000:
            scale = 0.5
            window_size = (self.prior_upper - self.prior_lower) / 5.0
        auto_optimize = True
        if self.power < 0.1:
            window_size = (self.prior_upper - self.prior_lower) / 1.5 
            scale = 1.0
            auto_optimize = False
        self.window_operator = WindowOperator(
                window_size = window_size,
                auto_optimize = auto_optimize)
        self.scale_operator = ScaleOperator(
                scale = scale,
                auto_optimize = auto_optimize)
        if screen_frequency is None:
            screen_frequency = sample_frequency * 10
        self.edge_length = self.draw_edge_from_prior()
        self.calc_power_ln_likelihood()
        if screen_frequency > 0:
            _LOG.info("gen\tlnl\tedge_length")
            _LOG.info("{gen}\t{lnl}\t{edge_length}".format(
                    gen = 0,
                    lnl = self.power_lnl,
                    edge_length = self.edge_length))
        mcmc_samples = [(0, self.power_lnl, self.edge_length)]
        for i in range(number_of_generations):
            self.perform_scale_update()
            # self.perform_window_update()
            if ((screen_frequency > 0) and ((i + 1) % screen_frequency == 0)):
                _LOG.info("{gen}\t{lnl}\t{edge_length}".format(
                        gen = i + 1,
                        lnl = self.power_lnl,
                        edge_length = self.edge_length))
            if (i + 1) % sample_frequency == 0:
                mcmc_samples.append((i + 1, self.power_lnl, self.edge_length))
        if (screen_frequency > 0):
            _LOG.info("Mean edge length: {0}".format(
                    sum([e for g, l, e in mcmc_samples[1:]]) / (len(
                            mcmc_samples) - 1)))
            _LOG.info("Window operator acceptance rate: {0}".format(
                    self.window_operator.get_acceptance_rate()))
            _LOG.info("Scale operator acceptance rate: {0}".format(
                    self.scale_operator.get_acceptance_rate()))
        return mcmc_samples

    def get_lnl_samples_from_prior(self, number_of_samples, screen_frequency = 1000):
        lnls = []
        for i in range(number_of_samples):
            if ((screen_frequency > 0) and ((i + 1) % screen_frequency == 0)):
                _LOG.info("Importance sample {0} of {1}".format(
                        i + 1, number_of_samples))
            edge_len = self.draw_edge_from_prior()
            lnls.append(self.calc_ln_likelihood(edge_length = edge_len))
        return lnls

    def get_monte_carlo_marginal_likelihood_approximation(self, number_of_samples):
        lnls = self.get_lnl_samples_from_prior(number_of_samples)
        assert len(lnls) == number_of_samples
        mean_lnl = sum(lnls) / len(lnls)
        scaled_ml = sum(math.exp(l - mean_lnl) for l in lnls) / len(lnls)
        ln_ml = math.log(scaled_ml) + mean_lnl
        return ln_ml

    def rectangular_integration_of_posterior(self, number_of_steps):
        ln_prior = self.get_edge_ln_prior_probability()
        step_size = (self.prior_upper - self.prior_lower) / number_of_steps
        position = self.prior_lower
        ln_areas = []
        total_ln_density = 0.0
        while position < self.prior_upper:
            if (position + step_size) > self.prior_upper:
                step_size = self.prior_upper - position
            midpoint = position + (step_size / 2.0)
            ln_density = self.calc_ln_likelihood(
                    edge_length = midpoint) + ln_prior
            total_ln_density += ln_density
            ln_slice_area = math.log(step_size) + ln_density 
            ln_areas.append(ln_slice_area)
            position += step_size
        assert ((len(ln_areas) == number_of_steps) or
                (len(ln_areas) == number_of_steps + 1)), (
                "unexpected number of steps {0}".format(len(ln_areas)))
        mean_ln_area = sum(ln_areas) / len(ln_areas)
        mean_ln_density = total_ln_density / len(ln_areas)
        # Reduce underflow by shifting all log areas before exponentiating
        scaled_areas = [math.exp(lna - mean_ln_area) for lna in ln_areas]
        scaled_ml = sum(scaled_areas)
        # Shift the log marginal likelihood back
        ln_ml = math.log(scaled_ml) + mean_ln_area
        return ln_ml, mean_ln_density

    def trapezoidal_integration_of_posterior(self, number_of_steps,
            mean_ln_posterior_density):
        ln_prior = self.get_edge_ln_prior_probability()
        step_size = (self.prior_upper - self.prior_lower) / number_of_steps
        position = self.prior_lower
        scaled_ml = 0.0
        n = 0
        while position < self.prior_upper:
            if (position + step_size) > self.prior_upper:
                step_size = self.prior_upper - position
            midpoint = position + (step_size / 2.0)
            left_ln_density = self.calc_ln_likelihood(
                    edge_length = position) + ln_prior
            right_ln_density = self.calc_ln_likelihood(
                    edge_length = position + step_size) + ln_prior
            # Shift densities to reduce underflow
            scaled_left_density = math.exp(left_ln_density - mean_ln_posterior_density)
            scaled_right_density = math.exp(right_ln_density - mean_ln_posterior_density)
            scaled_mean_density = (scaled_left_density + scaled_right_density) / 2.0
            scaled_ml += scaled_mean_density * step_size
            n += 1
            position += step_size
        assert ((n == number_of_steps) or
                (n == number_of_steps + 1)), (
                "unexpected number of steps {0}".format(n))
        # Shift the log marginal likelihood back
        ln_ml = math.log(scaled_ml) + mean_ln_posterior_density
        return ln_ml


    def abc_sample(self):
        edge_len = self.draw_edge_from_prior()
        seq_1, seq_2 = self.simulate_sequences(edge_len)
        div = count_differences(seq_1, seq_2) / float(self.seq_length)
        return edge_len, div

    def abc_rejection(self,
            number_of_prior_samples,
            number_of_posterior_samples,
            screen_frequency = 1000):
        observed_div = count_differences(self.seq1, self.seq2) / float(self.seq_length)
        posterior_samples = []
        for i in range(number_of_posterior_samples):
            edge_len, sim_div = self.abc_sample()
            distance = math.fabs(observed_div - sim_div)
            posterior_samples.append((distance, edge_len, sim_div))
        posterior_samples.sort()
        for i in range(number_of_prior_samples):
            if ((i + 1) % screen_frequency == 0):
                _LOG.info("ABC prior sample {0} of {1}".format(
                        i + 1, number_of_prior_samples))
            edge_len, sim_div = self.abc_sample()
            distance = math.fabs(observed_div - sim_div)
            if distance < posterior_samples[-1][0]:
                posterior_samples.pop()
                posterior_samples.append((distance, edge_len, sim_div))
                posterior_samples.sort()
        return posterior_samples

def arg_is_positive_int(i):
    try:
        if int(i) < 1:
            raise
    except:
        msg = '{0!r} is not a positive integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_float(i):
    try:
        if float(i) <= 0.0:
            raise
    except:
        msg = '{0!r} is not a positive real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return float(i)

def arg_is_nonnegative_float(i):
    try:
        if float(i) < 0:
            raise
    except:
        msg = '{0!r} is not a non-negative real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)



class ABCestimatorRegressWorker(object):
    def __init__(self,
            observed_path,
            posterior_path,
            exe_path = 'ABCestimator',
            output_prefix = None,
            parameter_indices = [1],
            num_posterior_samples = 1000,
            bandwidth = 0.002,
            num_posterior_quantiles = 500,
            ):
        self.process = None
        self.finished = False
        self.exit_code = None
        self.stdout = None
        self.stderr = None
        self.exe_path = exe_path
        if not self.exe_path:
            self.exe_path = 'ABCestimator'
        self.observed_path = observed_path
        self.posterior_path = posterior_path
        self.parameter_indices = parameter_indices
        self.num_posterior_samples = num_posterior_samples
        self.num_posterior_quantiles = int(num_posterior_quantiles)
        self.bandwidth = bandwidth
        self.failed = False
        self.output_prefix = output_prefix
        if self.output_prefix is None:
            self.output_prefix = os.path.join(os.path.curdir, "abc-glm")
        self.cfg_path = self.output_prefix + ".cfg"
        self.cmd = [self.exe_path, self.cfg_path]
        self._write_config()
        self.marginal_likelihood = None
        self.mode_edge_length = None

    def start(self):
        self._start()

    def _start(self):
        _LOG.debug('Starting process with following command:\n\t'
                '{0}'.format(' '.join([str(x) for x in self.cmd])))
        sout = subprocess.PIPE
        serr = subprocess.PIPE
        self.process = subprocess.Popen(self.cmd,
                stdout = sout,
                stderr = serr,
                shell = False)
        self.stdout, self.stderr = self.process.communicate()
        self.exit_code = self.process.wait()
        if self.exit_code != 0:
            _LOG.error('execution failed')
            raise Exception('ABCestimator failed.\ninvocation:\n{0}\n'
                    'stderr:\n{1}'.format(
                    ' '.join([str(x) for x in self.cmd]),
                    self.get_stderr()))
        try:
            self._post_process()
        except:
            _LOG.error("Error during post-processing")
            raise
        self.finished = True

    def _post_process(self):
        self._parse_glm_estimates()
        self._parse_abctoolbox_stdout()

    def _parse_glm_estimates(self):
        path = self.output_prefix + "-PosteriorCharacteristics_Obs0.txt"
        with open(path, "r") as stream:
            header = stream.next().strip().split()
            assert header == ["what", "edge_length"]
            modes = stream.next().strip().split()
            assert modes[0] == "mode"
            assert len(modes) == 2
            self.mode_edge_length = float(modes[1])
    
    def _parse_abctoolbox_stdout(self):
        glm_ml_matches = glm_density_pattern.findall(self.stdout)
        assert len(glm_ml_matches) == 1
        self.marginal_likelihood = float(glm_ml_matches[0])
        
    def _write_config(self):
        with open(self.cfg_path, "w") as cfg:
            cfg.write('estimationType standard\n')
            cfg.write('simName {0}\n'.format(self.posterior_path))
            cfg.write('obsName {0}\n'.format(self.observed_path))
            cfg.write('params {0}\n'.format(','.join(
                    [str(i+1) for i in self.parameter_indices])))
            cfg.write('diracPeakWidth {0}\n'.format(self.bandwidth))
            cfg.write('posteriorDensityPoints {0}\n'.format(
                    self.num_posterior_quantiles))
            cfg.write('stadardizeStats 0\n')
            cfg.write('writeRetained 0\n')
            cfg.write('maxReadSims {0}\n'.format(
                    self.num_posterior_samples + 100))
            cfg.write('outputPrefix {0}\n'.format(self.output_prefix + "-"))


def main_cli(argv = sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('--seq-length',
            action = 'store',
            type = arg_is_positive_int,
            default = 10000,
            help = 'Length of sequences.')
    # parser.add_argument('--mc-samples',
    #         action = 'store',
    #         type = arg_is_positive_int,
    #         default = 500000,
    #         help = ('Number of Monte Carlo samples to use to approximate'
    #                 'the marginal likelihood.'))
    # parser.add_argument('--mc-reps',
    #         action = 'store',
    #         type = arg_is_positive_int,
    #         default = 3,
    #         help = ('Number of replicate Monte Carlo approximations of the'
    #                 'marginal likelihood.'))
    parser.add_argument('--abc-prior-samples',
            action = 'store',
            type = arg_is_positive_int,
            default = 100000,
            help = ('Number of prior samples for ABC analyses.'))
    parser.add_argument('--abc-posterior-samples',
            action = 'store',
            type = arg_is_positive_int,
            default = 1000,
            help = ('Number of posterior samples for ABC analyses.'))
    parser.add_argument('--mcmc-generations',
            action = 'store',
            type = arg_is_positive_int,
            default = 10100,
            help = ('Number of MCMC generations to run.'))
    parser.add_argument('--mcmc-sample-frequency',
            action = 'store',
            type = arg_is_positive_int,
            default = 10,
            help = ('Frequency with which to sample MCMC chain.'))
    parser.add_argument('--prior-lower',
            action = 'store',
            type = arg_is_nonnegative_float,
            default = 0.000001,
            help = ('Lower limit on uniform edge length prior.'))
    parser.add_argument('--prior-upper',
            action = 'store',
            type = arg_is_nonnegative_float,
            default = 0.1,
            help = ('Upper limit on uniform edge length prior.'))
    parser.add_argument('--vague-prior-lower',
            action = 'store',
            type = arg_is_nonnegative_float,
            default = 0.0,
            help = ('Lower limit on vague uniform edge length prior.'))
    parser.add_argument('--vague-prior-upper',
            action = 'store',
            type = arg_is_nonnegative_float,
            default = 0.15,
            help = ('Upper limit on vague uniform edge length prior.'))
    parser.add_argument('--output-prefix',
            action = 'store',
            type = str,
            default = os.path.join(os.path.curdir, "rep-1"),
            help = ('Prefix for output files.'))
    parser.add_argument('seed',
            metavar='SEED',
            type = arg_is_positive_int,
            help = 'Seed for random number generator.')

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    m = EdgeModel(seed = args.seed,
            prior_lower = args.prior_lower,
            prior_upper = args.prior_upper,
            seq_length = args.seq_length)

    true_parameter_path = args.output_prefix + "-true-parameter-value.txt"
    with open(true_parameter_path, "w") as out:
        out.write("edge_length\n{0}\n".format(m.edge_length))

    observed_div = m.get_proportion_of_variable_sites() 
    observed_stat_path = args.output_prefix + "-observed-sum-stat.txt"
    with open(observed_stat_path, "w") as out:
        out.write("proportion_of_variable_sites\n{0}\n".format(
                observed_div))

    start_time = datetime.datetime.now()

    for i in range(2):
        # First pass, inference under true model
        # Second pass, inference under model with vague prior
        if i == 1:
            m.prior_lower = args.vague_prior_lower
            m.prior_upper = args.vague_prior_upper
            m.output_prefix += "-vague-model"

        # ml_estimates = []
        # for i in range(args.mc_reps):
        #     _LOG.info("MC replicate {0} of {1}".format(i + 1,
        #             args.mc_reps))
        #     ml = m.get_monte_carlo_marginal_likelihood_approximation(
        #             args.mc_samples)
        #     ml_estimates.append(ml)
        # ml_estimates.sort()
        # ml_path = args.output_prefix + "-mc-ml-estimates.txt"
        # with open(ml_path, "w") as out:
        #     for e in ml_estimates:
        #         out.write("{0}\n".format(e))
        for number_of_steps in [100, 1000, 10000]:
            rect_ml, mean_density = m.rectangular_integration_of_posterior(
                    number_of_steps = number_of_steps)
            trap_ml = m.trapezoidal_integration_of_posterior(
                    number_of_steps = number_of_steps,
                    mean_ln_posterior_density = mean_density)
            rect_path = args.output_prefix + "-rectangular-ml-estimate-{0}.txt".format(
                    number_of_steps)
            trap_path = args.output_prefix + "-trapezoidal-ml-estimate-{0}.txt".format(
                    number_of_steps)
            with open(rect_path, "w") as out:
                out.write("{0}\n".format(rect_ml))
            with open(trap_path, "w") as out:
                out.write("{0}\n".format(trap_ml))

        abc_sample = m.abc_rejection(
                number_of_prior_samples = args.abc_prior_samples,
                number_of_posterior_samples = args.abc_posterior_samples)
        abc_path = args.output_prefix + "-abc-posterior-sample.txt"
        with open(abc_path, "w") as out:
            out.write("distance\tedge_length\tproportion_of_variable_sites\n")
            for d, e, p in abc_sample:
                out.write("{d}\t{e}\t{p}\n".format(d = d, e = e, p = p))

        glm_worker = ABCestimatorRegressWorker(
                observed_path = observed_stat_path,
                posterior_path = abc_path,
                exe_path = "ABCestimator",
                output_prefix = args.output_prefix + "-abctoolbox-glm",
                parameter_indices = [1],
                num_posterior_samples = args.abc_posterior_samples,
                bandwidth = 2.0 / args.abc_posterior_samples,
                num_posterior_quantiles = args.abc_posterior_samples / 2)
        glm_worker.start()

        glm_ml_path = args.output_prefix + "-abctoolbox-glm-ml-estimate.txt"
        with open(glm_ml_path, "w") as out:
            out.write("{0}\n".format(glm_worker.marginal_likelihood))
        glm_edge_path = args.output_prefix + "-abctoolbox-glm-edge-mode.txt"
        with open(glm_edge_path, "w") as out:
            out.write("{0}\n".format(glm_worker.mode_edge_length))

        mcmc_sample = m.mcmc(
                number_of_generations = args.mcmc_generations,
                sample_frequency = args.mcmc_sample_frequency)
        mcmc_path = args.output_prefix + "-mcmc-sample.txt"
        with open(mcmc_path, "w") as out:
            out.write("generation\tln_likelihood\tedge_length\n")
            for g, l, e in mcmc_sample:
                out.write("{g}\t{l}\t{e}\n".format(g = g, l = l, e = e))
        mcmc_edge_samples = [e for g, l, e in mcmc_sample[1:]]
        mcmc_mean_edge = sum(mcmc_edge_samples) / len(mcmc_edge_samples)
        mcmc_edge_path = args.output_prefix + "-mcmc-edge-mean.txt"
        with open(mcmc_edge_path, "w") as out:
            out.write("{0}\n".format(mcmc_mean_edge))

        m.power = 0.0

        mcmc_sample = m.mcmc(
                number_of_generations = args.mcmc_generations,
                sample_frequency = args.mcmc_sample_frequency)
        mcmc_path = args.output_prefix + "-mcmc-prior-sample.txt"
        with open(mcmc_path, "w") as out:
            out.write("generation\tln_likelihood\tedge_length\n")
            for g, l, e in mcmc_sample:
                out.write("{g}\t{l}\t{e}\n".format(g = g, l = l, e = e))

    stop_time = datetime.datetime.now()
    run_time = stop_time - start_time
    _LOG.info("Run time: {0}".format(str(run_time)))

if __name__ == "__main__":
    if "--run-tests" in sys.argv:

        sys.stderr.write("""
*********************************************************************
Running test suite using the following Python executable and version:
{0}
{1}
*********************************************************************
\n""".format(sys.executable, sys.version))

        import doctest

        # doctest.testmod(verbose = True)
        suite = unittest.TestSuite()
        suite.addTest(doctest.DocTestSuite())

        tests = unittest.defaultTestLoader.loadTestsFromName(
               os.path.splitext(os.path.basename(__file__))[0])
        suite.addTests(tests)

        runner = unittest.TextTestRunner(verbosity = 2)
        runner.run(suite)

        sys.exit(0)

    main_cli()
