# effect_size_simulation.py
"""
Will simulation population-level data for a case-control study. Will then
simulate many experiments at different sample sizes. For each it will compute
estimates of sample effect size as well as effect size inflation. Can produce a
plot of effect sizes as pdf and can produce a csv file with statistics about the
degree of effect size inflation.

Example usage:
python effect_size_simulation.py --es 0.5 --n_exp 10000 --espdf2save es_subplot.pdf --csv2save test.csv

python effect_size_simulation.py --n_exp 10000 --esinfpdf2save es_inf_plot.pdf
"""

# import libraries
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
from optparse import OptionParser
import random


# set seed for reproducibility
np.random.seed(1)


# function to parse input arguments
def parse_args():
    """
    Parse arguments.
    """
    parser=OptionParser()
    parser.add_option('--es',"",dest='es',help="Population effect size to use ex: --es 0.5",default=0.5)
    parser.add_option('--n_exp',"",dest='n_exp',help="Number of experiments to simulat ex: --n_exp 10000",default=10000)
    parser.add_option('--espdf2save',"",dest='espdf2save',help="PDF filename for effect size figure ex: --espdf2save plot.pdf",default=None)
    parser.add_option('--csv2save',"",dest='csv2save',help="csv file to save results of simulation over range of effect sizes ex: --csv2save res.csv",default=None)
    parser.add_option('--esinfpdf2save',"",dest='esinfpdf2save',help="PDF filename for effect size inflation figure ex: --esinfpdf2save es_inf.pdf",default=None)
    (options,args) = parser.parse_args()
    return(options)


# function to compute Cohen's d
def cohens_d(x, y, DIM='rows', SIGN = True):
    """
    Compute Cohen's d.
    """

    if DIM == 'rows':
        dim = 0
    elif DIM == 'columns':
        dim = 1

    # n-1 for x and y
    lx = x.shape[dim]-1
    ly = y.shape[dim]-1

    # mean difference
    if SIGN:
        md = np.mean(x, axis = dim, dtype = np.float64) - np.mean(y, axis = dim, dtype = np.float64)
    else:
        md = np.abs(np.mean(x, axis = dim, dtype = np.float64) - np.mean(y, axis = dim, dtype = np.float64))

    # pooled variance
    csd = (lx * np.var(x, axis = dim, dtype = np.float64)) + (ly * np.var(y, axis = dim, dtype = np.float64))
    csd = np.sqrt(csd/(lx + ly))

    # compute cohen's d
    d  = md/csd

    return(d)


# function to generate population data
def generate_population(pop1_params, pop2_params, pop_size):
    """
    Generate population-level data.

    pop1_params, pop2_params
        First element of pop1_params and pop2_params is the population mean.
        Second element is population standard deviation.

    pop_size
        Size of the population.

    Output will be two arrays of normally distributed population-level data.
    """

    # grab population parameters
    pop1_mu = pop1_params[0]
    pop1_sd = pop1_params[1]
    pop2_mu = pop2_params[0]
    pop2_sd = pop2_params[1]

    # generate random normally distributed data at the population-level
    pop1_data = np.random.normal(pop1_mu, pop1_sd, pop_size)
    pop2_data = np.random.normal(pop2_mu, pop2_sd, pop_size)

    return([pop1_data,pop2_data])


def randomly_select_sample(pop_data, sample_size):
    """
    Randomly sample data from population without replacement.
    """

    sample_data = np.array(random.sample(pop_data, sample_size))
    # sample_data = np.random.choice(pop_data, sample_size, replace = False)

    return(sample_data)


# function to evaluate results of experiment
def evaluate_experiment(samp1_data, samp2_data, alpha_thresh = 0.05):
    """
    Calculate effect size and hypothesis test stats on experiment.
    """

    # compute effect size
    d = cohens_d(samp1_data, samp2_data)

    # run hypothesis test
    hp_res = stats.ttest_ind(samp1_data, samp2_data)

    # grab t-stat
    t = hp_res.statistic

    # grab p-value
    p = hp_res.pvalue

    # grab indicator of whether H0 is rejected
    alpha_thresh = np.array(alpha_thresh, dtype = float)
    reject_h0 = p < alpha_thresh

    return((d, t, p, reject_h0))


# function to simulate an experiment
def simulate_experiment(pop1_data,pop2_data,sample_size):
    """
    Simulate random sampling from population.
    """

    # randomly sample group1
    # [samp1_data, samp1_idx] = randomly_select_sample(pop1_data, sample_size)
    samp1_data = randomly_select_sample(pop1_data, sample_size)

    # randomly sample group2
    # [samp2_data, samp2_idx] = randomly_select_sample(pop2_data, sample_size)
    samp2_data = randomly_select_sample(pop2_data, sample_size)

    return((samp1_data, samp2_data))


# function to plot all simulated effect sizes
def plot_histogram(effect_size, nbins = 100, face_color = "gray"):
    """
    Plot histogram of all simulated experiments
    """

    [n, bins, patches] = plt.hist(effect_size, nbins, facecolor = face_color)


# function to plot just statistically significant effect sizes
def plot_sigexp_histogram(effect_size, mask, nbins = 100, face_color = "red"):
    """
    Plot histogram of just statistically significant experiments.
    """

    es2plot = effect_size[mask]
    [n, bins, patches] = plt.hist(es2plot, nbins, facecolor = face_color)

# function to make an effect size plot
def make_effectsize_plot(effect_size, mask, popES, sample_size):

    # plot all simulated effect sizes
    plot_histogram(effect_size, nbins = 100, face_color = "gray")

    # plot just the statistically significant effect sizes
    plot_sigexp_histogram(effect_size, mask, nbins = 100, face_color = "red")

    # plot the true population effect size
    plt.axvline(x = popES, color = "green")

    # insert x and y-axis labels
    plt.xlabel('Effect Size')
    plt.ylabel("Count")

    # insert grid
    plt.grid()

    # insert plot title
    plt.title(sample_size)

# function for making specific subplots
def make_subplot(axarr, x, y, effect_size, mask, popES, ci, sample_size,
    nbins = 100, gridline_width = 0.5, xlimits = [-1, 2], font_name = "Arial"):
    """
    Make a subplot.
    """

    # plot all simulated effect sizes
    axarr[x, y].hist(effect_size, nbins, facecolor = "gray")

    # plot just the statistically significant effect sizes
    axarr[x, y].hist(effect_size[mask], nbins, facecolor = "red")

    # plot the true population effect size
    axarr[x, y].axvline(x = popES, color = (0,1,0))

    # plot confidence intervals for effect size
    axarr[x, y].axvline(x = ci[0], color = "black")
    axarr[x, y].axvline(x = ci[1], color = "black")

    # insert x and y-axis labels
    axarr[x, y].set_xlabel('Effect Size', fontname=font_name)
    axarr[x, y].set_ylabel("Count", fontname=font_name)

    # insert subplot title
    title_str = "n = %s" % (str(sample_size))
    axarr[x, y].set_title(title_str, fontname=font_name, fontweight='bold')

    # insert grid
    axarr[x, y].grid(linewidth = gridline_width)

    # set x-axis limits
    axarr[x, y].set_xlim(xlimits)

# function to put together multiple subplots
def make_effectsize_subplots(nrows, ncols, effect_size, mask, popES, ci,
    sample_sizes, nbins = 100, gridline_width = 0.5, xlimits = [-1, 2],
    font_name = "Arial"):
    """
    Put all subplots together.
    """

    # Four axes, returned as a 2-d array
    f, axarr = plt.subplots(nrows, ncols)

    # matrix with indices for each subplot
    ss_mat = np.array(range(0,len(sample_sizes))).reshape([nrows,ncols])

    # loop over rows and columns in subplot
    for irow in range(0,nrows):
        for icol in range(0,ncols):

            # grab input for making specific subplot
            ss_idx = ss_mat[irow,icol]
            sample_size2use = sample_sizes[ss_idx]
            d2use = effect_size[:, ss_idx]
            mask2use = mask[:,ss_idx].astype(bool)
            ci2use = ci[ss_idx,:]

            # draw specific subplot
            make_subplot(axarr = axarr, x = irow, y = icol, effect_size = d2use,
                mask = mask2use, popES = popES, ci = ci2use,
                sample_size = sample_size2use, nbins = nbins,
                gridline_width = gridline_width,
                xlimits = xlimits)

    # set tight layout to eliminate overlap of subplots
    plt.tight_layout()
    plt.suptitle("Effect Size d = %0.02f" % popES, fontname=font_name,
        fontweight='bold')

# function to find percentiles
def find_percentile(effect_size, ci_interval):
    """
    Find effect size percentile scores.
    """

    low_ci = stats.scoreatpercentile(effect_size, ci_interval[0])
    high_ci = stats.scoreatpercentile(effect_size, ci_interval[1])
    return([low_ci, high_ci])

# function for computing effect size inflation stats
def compute_inflation_stats(effect_size, h0, sample_sizes, popES,
    ci_interval = [0.5, 99.5]):
    """
    Calculate statistics for effect size inflation.
    """

    # make h0 boolean mask
    mask = h0.astype(bool)

    # pre-allocate variables to use in loop over sample sizes
    avg_inflate = np.zeros([len(sample_sizes)])
    median_inflate = np.zeros([len(sample_sizes)])
    min_inflate = np.zeros([len(sample_sizes)])
    max_inflate = np.zeros([len(sample_sizes)])
    low_ci_inflate = np.zeros([len(sample_sizes)])
    high_ci_inflate = np.zeros([len(sample_sizes)])

    # loop over sample sizes
    for ss_index, sample_size in enumerate(sample_sizes):
        # grab mask of statistically significant experiments
        mask2use = mask[:,ss_index]

        # grab effect sizes from statistically significant experiments
        data2use = effect_size[mask2use, ss_index]

        # average effect size for rejected H0
        avg_inflate[ss_index] = data2use.mean()

        # median effect size for rejected H0
        median_inflate[ss_index] = np.median(data2use)

        # min effect size for rejected H0
        min_inflate[ss_index] = data2use.min()

        # max effect size for rejected H0
        max_inflate[ss_index] = data2use.max()

        # ci of rejected H0
        low_ci_inflate[ss_index] = stats.scoreatpercentile(data2use, ci_interval[0])
        high_ci_inflate[ss_index] = stats.scoreatpercentile(data2use, ci_interval[1])

    # average effect size percent increase
    avg_es_inflate_increase = ((avg_inflate-popES)/popES)*100

    # dictionary with effect size inflation stats
    es_inflation_stats = {"n":sample_sizes, "mean_d":avg_inflate,
        "low_ci":low_ci_inflate, "high_ci":high_ci_inflate,
        "median_d":median_inflate,
        "mean_d_percent_increase":avg_es_inflate_increase,
        "popD":np.ones([len(sample_sizes)])*popES}

    # data frame with effect size inflation stats to output
    es_inf_stats_df = pd.DataFrame(data = es_inflation_stats)
    return(es_inf_stats_df)


# function for running main part of the simulation
def run_main_simulation(pop_mean1, pop_sd1, pop_mean2, pop_sd2, sample_sizes,
    pop_size = 10000000, n_exp = 10000, ci_interval = [0.5, 99.5]):
    """
    Main function for effect size simulation.
    """

    # population parameters
    pop1_params = [pop_mean1, pop_sd1]
    pop2_params = [pop_mean2, pop_sd2]

    # generate population level data
    [pop1_data, pop2_data] = generate_population(pop1_params, pop2_params, pop_size)

    # compute population effect size
    D  = cohens_d(pop1_data, pop2_data)

    # # pre-allocate variables
    d = np.zeros([n_exp,len(sample_sizes)])
    t = np.zeros([n_exp,len(sample_sizes)])
    p = np.zeros([n_exp,len(sample_sizes)])
    h0 = np.zeros([n_exp,len(sample_sizes)])

    # Simulate n experiments
    for iexp in np.arange(n_exp):

        # loop over different sample sizes
        for ss_index, sample_size in enumerate(sample_sizes):

            print("Running experiment #%d at n=%d" % (iexp+1, sample_size))

            # simulate random samples from population
            [samp1_data, samp2_data] = simulate_experiment(pop1_data,pop2_data,sample_size)
            # grab statistics to evaluate results of experiment
            [es, tstat, pval, reject_h0] = evaluate_experiment(samp1_data, samp2_data)
            d[iexp,ss_index] = es
            t[iexp,ss_index] = tstat
            p[iexp,ss_index] = pval
            h0[iexp,ss_index] = reject_h0

    # find percentile scores for effect sizes
    ci = np.zeros([len(sample_sizes),len(ci_interval)])
    for ss_index, sample_size in enumerate(sample_sizes):
        ci[ss_index,:] = find_percentile(effect_size = d[:,ss_index],
            ci_interval = ci_interval)

    # Return sample effect sizes (d), tstats (t), pvals (p), rejected null
    # hypotheses (h0), true population ES (D), and confidence intervals for
    # effect sizes (ci)
    return([d, t, p, h0, D, ci])


# function to run simulation over range of effect sizes
def effect_size_inflation_sim(es_range, pop_sd1, pop_mean2, pop_sd2,
    sample_sizes, n_exp = 10000, pop_size = 10000000):
    """
    Run simulation over a range of effect sizes and sample sizes
    """

    # pre-allocate results array
    es_inf_res = np.zeros([len(sample_sizes),es_range.shape[0]])

    # loop over range of effect sizes
    for es_idx, es in enumerate(es_range):

        # run main simulation
        [d, t, p, h0, D, ci] = run_main_simulation(pop_mean1 = es,
            pop_sd1 = pop_sd1, pop_mean2 = pop_mean2, pop_sd2 = pop_sd2,
            sample_sizes = sample_sizes, pop_size = pop_size, n_exp = n_exp)

        # calcular effect size inflation stats
        es_inflation_stats = compute_inflation_stats(effect_size = d, h0 = h0,
            sample_sizes = sample_sizes, popES = D, ci_interval = [0.5, 99.5])

        # grab result and put into results array
        es_inf_res[:, es_idx] = es_inflation_stats["mean_d_percent_increase"]

    es_inf_res = es_inf_res.T
    return(es_inf_res)


# function to plot effect size inflation over range of effect sizes
def plot_es_inflation(es_inf_res, es_range, sample_sizes, gridline_width = 0.5,
    fig_size = (12,12), font_name = "Arial", font_size = 18):
    """
    Plot average effect size inflation over range of effect sizes and sample
    sizes.
    """

    # make figure a specific size
    plt.figure(figsize = fig_size)

    # make plot
    plt.plot(es_inf_res[1:,:])

    # add grid lines
    plt.grid(linewidth = gridline_width)

    # add appropriate xticks
    plt.xticks(range(0,len(es_range)-1), es_range[1:])

    # make legend
    ss_legend = []
    for ss in sample_sizes:
        ss_legend.append("n = %d" % ss)
    plt.legend(ss_legend)

    # add x and y-axis labels
    plt.ylabel("Average Effect Size Inflation (Percent Increase)",
        fontname = font_name, fontsize = font_size)
    plt.xlabel("Population Effect Size", fontname = font_name,
        fontsize = font_size)

    # tight layout
    plt.tight_layout()

    # add plot title
    plt.title("Average Effect Size Inflation", fontname = font_name,
        fontweight='bold', fontsize = font_size)







# boilerplate code to call main code for executing
if __name__ == '__main__':

    # parse arguments
    opts = parse_args()

    # run simulation
    pop_mean1 = np.array(opts.es, dtype=float)
    pop_sd1 = np.array(1, dtype = float)
    pop_mean2 = np.array(0, dtype = float)
    pop_sd2 = np.array(1, dtype = float)
    n_exp = np.array(opts.n_exp, dtype = int)
    sample_sizes = [20, 50, 100, 200, 1000, 2000]
    pop_size = 10000000
    [d, t, p, h0, D, ci] = run_main_simulation(pop_mean1 = pop_mean1,
        pop_sd1 = pop_sd1, pop_mean2 = pop_mean2, pop_sd2 = pop_sd2,
        sample_sizes = sample_sizes, pop_size = pop_size, n_exp = n_exp)

    # make effect size plot
    if opts.espdf2save is not None:
        n_subplot_rows = 3
        n_subplot_cols = 2
        make_effectsize_subplots(nrows = n_subplot_rows, ncols = n_subplot_cols,
            effect_size = d, mask = h0, popES = D, ci = ci,
            sample_sizes = sample_sizes, nbins = 100, gridline_width = 0.5,
            xlimits = [-1, 2])

        # save effect size plot as pdf
        plt.savefig(opts.espdf2save)

    # compute effect size inflation stats
    if opts.csv2save is not None:
        es_inflation_stats = compute_inflation_stats(effect_size = d, h0 = h0,
            sample_sizes = sample_sizes, popES = D, ci_interval = [0.5, 99.5])
        print(es_inflation_stats)

        # write effect size inflation stats to csv file
        es_inflation_stats.to_csv(opts.csv2save, index = False)


    if opts.esinfpdf2save is not None:
        # run simulation over a range of effect sizes
        es_range = np.linspace(0,2,21)
        es_inf_res = effect_size_inflation_sim(es_range = es_range,
            pop_sd1 = pop_sd1, pop_mean2 = pop_mean2, pop_sd2 = pop_sd2,
            sample_sizes =  sample_sizes, n_exp = n_exp, pop_size = pop_size)

        # plot effect size inflation over a range of effect sizes and sample sizes
        plot_es_inflation(es_inf_res, es_range, sample_sizes, gridline_width = 0.5,
            fig_size = (10,8))
        # save plot as pdf
        plt.savefig(opts.esinfpdf2save)
