# heterogeneity_simulation.py
"""
Here we will simulate ASD population with 5 subtypes. Each subtype has a 20%
prevalence in the population. Both ASD and TD populations will be set to the
same mean and sd, so there will be no overall case-control difference. We
then simulate experiments with different sample sizes, and compute what is
the sample prevalence for the 5 ASD subtypes.

Example usage:
python heterogeneity_simulation.py --n_exp 10000 --n_subgrp 5 --mu_subgrp '-1,-0.5,0,0.5,1' --sample_sizes '20,50,100,200,1000,2000'--pop_sd 1 --subgrp_pop_n 200000 --kds2save kdsplot.pdf --pdf2save sample_prevalence_plot.pdf
"""

# import libraries
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import random
import pickle
from optparse import OptionParser

# set seed for reproducibility
np.random.seed(1)



# function to parse input arguments
def parse_args():
    """
    Parse arguments.
    """
    parser=OptionParser()
    parser.add_option('--n_exp',"",dest='n_exp',help="Number of experiments to simulation ex: --n_exp 10000",default=10000)
    parser.add_option('--n_subgrp',"",dest='n_subgrp',help="Number of ASD subgroups to simulation ex: --n_subgrp 5",default=5)
    parser.add_option('--mu_subgrp',"",dest='mu_subgrp',help="Effect size for each subgroup ex: --mu_subgrp '-1,-0.5,0,0.5,1'",default='-1,-0.5,0,0.5,1')
    parser.add_option('--sample_sizes',"",dest='sample_sizes',help="Sample sizes to simulate ex: --sample_sizes '20,50,100,200,1000,2000'",default='20,50,100,200,1000,2000')
    parser.add_option('--pop_sd',"",dest='pop_sd',help="Population standard deviation ex: --pop_sd 1",default=1)
    parser.add_option('--subgrp_pop_n',"",dest='subgrp_pop_n',help="Size of each subgroup in the population ex: --subgrp_pop_n 200000",default=200000)
    parser.add_option('--ksd2save',"",dest='ksd2save',help="PDF filename of ksdensity plot figure to save ex: --ksd2save ksd_plot.pdf",default=None)
    parser.add_option('--pdf2save',"",dest='pdf2save',help="PDF filename of figure to save ex: --pdf2save plot.pdf",default=None)
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



# function for generating population-level ASD data with subgrps
def generate_asd_data(mu_subgrp, sd, subgrp_pop_n, n_subgrp):
    """
    Generate simulated ASD populaton data with n subgrps.
    """

    # population n accumulated across all subgrps
    pop_n = n_subgrp * subgrp_pop_n

    # pre-allocate data array to be [n, n_subgrp]
    data = np.zeros([subgrp_pop_n,n_subgrp])

    # generate gaussian data
    for mu_idx, mu in enumerate(mu_subgrp):
        data[:,mu_idx] = np.random.normal(mu, sd, subgrp_pop_n)

    # stack subgrp data on top of each other
    data_stacked = data.reshape((pop_n,1))
    return((data, data_stacked))



# function for making ASD subgrp labels
def make_subgrp_labels(data):
    """
    Make ASD subgrp labels.
    """

    # subgrp and population parameters from data
    subgrp_n = data.shape[0]
    n_subgrp = data.shape[1]
    pop_n = n_subgrp * subgrp_n

    # make a matrix of subgrp labels
    labels2use = range(1,n_subgrp+1)
    subgrp_mat = np.zeros([subgrp_n, n_subgrp])
    for subgrp_idx, subgrp_label in enumerate(labels2use):
        subgrp_mat[:,subgrp_idx] = np.matlib.repmat(1*(subgrp_idx+1),1,subgrp_n)

    # reshape subgrp matrix into a stacked array
    subgrp_mat = subgrp_mat.T
    subgrp_labels = subgrp_mat.reshape((pop_n,1))

    n_subgrps = np.unique(subgrp_labels)

    return([subgrp_labels, n_subgrps])



# function for generating population-level non-ASD data
def generate_nonasd_data(asd_mu, asd_sd, pop_n):
    """
    Generate simulated non-ASD populaton data.
    """

    # non-ASD standard deviation is half that of ASD
    sd = asd_sd/2
    data = np.random.normal(asd_mu, sd, pop_n)

    return(data)



# function to randomly select sample from population without replacement
def randomly_select_sample(pop_data, sample_size):
    """
    Randomly sample data from population without replacement.
    """

    # make indices 1:len(pop_data)
    idx = range(len(pop_data))

    # randomly permute indices
    rand_idx = np.random.permutation(idx)

    # grab the first n randomly permuted indices
    sample_idx = rand_idx[range(sample_size)]

    # sample data
    sample_data  = pop_data[sample_idx]

    return([sample_data, sample_idx])



# function to count frequency of each subgrp in a sample
def subgrp_sample_freq(sample_subgrp_labels, n_subgrp, sample_size):
    """
    Count how many individuals from each subgrp are present in a sample.
    """

    subgrps = np.arange(n_subgrp)+1

    # pre-allocate
    subgrp_freq = np.arange(n_subgrp, dtype = float)

    for subgrp_idx, subgrp in enumerate(subgrps):
        mask = sample_subgrp_labels==subgrp
        subgrp_freq[subgrp_idx] = sum(mask)

    # calculate prevalence of subgrp in the sample
    subgrp_prevalence = subgrp_freq/sample_size

    return([subgrp_freq, subgrp_prevalence])



# function to simulate an experiment
def simulate_experiment(pop1_data_stacked, pop2_data, sample_size, pop1_data):
    """
    Simulate random sampling from population.
    """

    # make subgrp on pop1_data without stacking
    [subgrp_labels, n_subgrp] = make_subgrp_labels(pop1_data)

    # randomly sample group1
    [samp1_data, samp1_idx] = randomly_select_sample(pop1_data_stacked,
        sample_size)
    samp1_subgrp_labels = subgrp_labels[samp1_idx]

    # randomly sample group2
    [samp2_data, samp2_idx] = randomly_select_sample(pop2_data, sample_size)

    # calculate sample effect sizes
    d = cohens_d(samp1_data, samp2_data)

    # calculate subgrp prevalence in sample
    [subgrp_freq, subgrp_prevalence] = subgrp_sample_freq(sample_subgrp_labels = samp1_subgrp_labels,
        n_subgrp = len(n_subgrp), sample_size = sample_size)

    # make dictionary with results to output
    results = {"effectsize":d, "subgrp_prevalence":subgrp_prevalence,
        "samp1_data":samp1_data, "samp1_idx":samp1_idx, "samp2_data":samp1_data,
        "samp2_idx":samp2_idx, "sample_size":sample_size}

    return(results)



# function for running simulation across many experiments
def run_main_simulation(pop1_data_stacked, pop2_data, n_exp, sample_size,
    pop1_data):
    """
    Run main simulation at specific sample size.
    """

    # make vector of experiments
    experiments = np.arange(n_exp)+1

    # number of subgrps
    n_subgrp = pop1_data.shape[1]

    # pre-allocate
    eff_size = np.zeros([n_exp, 1])
    subgrp_prevalence = np.zeros([n_exp, n_subgrp])

    # loop over experiments
    for idx, i_exp in enumerate(experiments):

        print("Running experiment #%d at n=%d" % (i_exp, sample_size))

        # simulate experiment
        results = simulate_experiment(pop1_data_stacked = pop1_data_stacked,
            pop2_data = pop2_data, sample_size = sample_size,
            pop1_data = pop1_data)

        # grab results of experiment
        eff_size[idx,:] = results["effectsize"]
        subgrp_prevalence[idx,:] = results["subgrp_prevalence"]

    # make dictionary for output
    out_res = {"effectsize":eff_size, "subgrp_prevalence":subgrp_prevalence}

    return(out_res)



# function for running simulation over a range of sample sizes
def run_simulation_over_sample_sizes(pop1_data_stacked, pop2_data, n_exp,
    sample_sizes, pop1_data):
    """
    Run the main simulation over a range of sample sizes.
    """

    # build results dictionary
    res_keys = []
    for sample_size in sample_sizes:
        res_keys.append("n%d_effectsize" % sample_size)
        res_keys.append("n%d_subgrp_prevalence" % sample_size)
    results = {key: [] for key in res_keys}

    for sample_size in sample_sizes:
        # run main simulation
        res = run_main_simulation(pop1_data_stacked, pop2_data, n_exp,
            sample_size, pop1_data)

        # grab effect size and subgrp sample prevalence
        key2use = "n%d_effectsize" % sample_size
        results[key2use] = res["effectsize"]

        key2use = "n%d_subgrp_prevalence" % sample_size
        results[key2use] = res["subgrp_prevalence"]

    return(results)



# function for plotting histograms
def plot_subgrp_histograms(data, fig_size = [12,12], nbins = 100,
    alpha_level = 0.7, gridline_width = 0.5):
    """
    Make histograms for each ASD subgrp.
    """

    # create legend labels
    label2use = []
    for i in range(1,data.shape[1]+1):
        label2use.append("ASD%d" % i)

    # make figure specific size
    plt.figure(figsize = fig_size)

    # plot histograms
    plt.hist(data, nbins, alpha = alpha_level, label = label2use)

    # add legend
    plt.legend()

    # add grid lines
    plt.grid(linewidth = gridline_width)

    # add x and y-axis labels
    plt.xlabel("DV")
    plt.ylabel("Count")

    # show plot
    plt.show()



# generate subsample of data from population for ksdensity plot
def generate_data_for_ksdensity_plot(grp_type, mu, sd, n4plot = 10000):

    if grp_type == "ASDsubgrps":
        # number of subgrps
        ngrp = 5
        # pre-allocate
        data4plot = np.zeros([n4plot,ngrp])

        # loop over each subgrp
        for m_idx, m in enumerate(mu):
            # generate gaussian data
            data4plot[:,m_idx] = np.random.normal(m, sd, n4plot)

        data4plot_stacked = data4plot.reshape((ngrp*n4plot,1))

    elif grp_type == "non-ASD" or grp_type == "ASD":
        # number of subgrps
        ngrp = 1
        data4plot = np.random.normal(mu, sd, n4plot)
        data4plot_stacked = data4plot.reshape((n4plot,1))

    return(data4plot, data4plot_stacked)



# function for making kernel density plot
def plot_subgrp_ksdensity(mu_subgrp, sd, n4plot = 10000,
    xlimits = [-6,6], gridline_width = 0.5, fig_size = [12,12],
    grp_type = "ASDsubgrps", axarr=None, sp_idx = None):
    """
    Make kernel density plots for each ASD subgrp.
    """

    # generate data for plot
    [data4plot, data_stacked] = generate_data_for_ksdensity_plot(grp_type = grp_type,
        mu = mu_subgrp, sd = sd, n4plot = n4plot)

    # loop over each subgrp
    for mu_idx, mu in enumerate(mu_subgrp):
        data2use = data4plot[:,mu_idx]
        # get kernel density
        density = stats.gaussian_kde(data2use)
        # make xs for plotting
        xs = np.linspace(xlimits[0], xlimits[1], density.n)
        # create legend label
        label2use = "ASD%d" % (mu_idx+1)
        if axarr is None:
            plt.figure(figsize = fig_size)

            # make plot
            plt.plot(xs, density(xs), label = label2use)

            # add grid lines
            plt.grid(linewidth = gridline_width)

            # add legend
            plt.legend()

            # add x and y-axis labels
            plt.xlabel("DV")
            plt.ylabel("Count")

        else:
            axarr[sp_idx].plot(xs, density(xs), label = label2use)

            # add grid lines
            axarr[sp_idx].grid(linewidth = gridline_width)

            # add legend
            axarr[sp_idx].legend()

            # add x and y-axis labels
            axarr[sp_idx].set_xlabel("DV")
            axarr[sp_idx].set_ylabel("Count")



# function for making kernel density plot
def plot_pop_ksdensity(grp, mu, sd, n4plot = 10000, xlimits = [-6,6],
    gridline_width = 0.5, fig_size = [12,12], axarr = None, sp_idx = None):
    """
    Make kernel density plots for each population.
    """

    # generate data for plot
    [data4plot, data_stacked] = generate_data_for_ksdensity_plot(grp_type = grp,
        mu = mu, sd = sd, n4plot = n4plot)

    if grp=="ASDsubgrps":
        data2use = data_stacked
    elif grp == "non-ASD" or grp == "ASD":
        data2use = data4plot

    # get kernel density
    density = stats.gaussian_kde(data2use)

    # make xs for plotting
    xs = np.linspace(xlimits[0], xlimits[1], density.n)

    if axarr is None:
        plt.figure(figsize = fig_size)

        # make plot
        plt.plot(xs,density(xs), label = grp)

        # add grid lines
        plt.grid(linewidth = gridline_width)

        # add legend
        plt.legend()

        # add x and y-axis labels
        plt.xlabel("DV")
        plt.ylabel("Count")

    else:
        axarr[sp_idx].plot(xs,density(xs), label = grp)

        # add grid lines
        axarr[sp_idx].grid(linewidth = gridline_width)

        # add legend
        axarr[sp_idx].legend()

        # add x and y-axis labels
        axarr[sp_idx].set_xlabel("DV")
        axarr[sp_idx].set_ylabel("Count")



# make all ksdensity plots
def make_ksdensity_subplots(grand_mu, grand_sd, mu_subgrp, sd, n4plot = 10000,
    xlimits = [-5,5], gridline_width = 0.5, fig_size = [12,12],
    nrows = 3, ncols = 1):
    """
    Put together all ksdensity subplots into one figure.
    """

    # Four axes, returned as a 2-d array
    f, axarr = plt.subplots(nrows, ncols)
    f.set_size_inches(fig_size[0],fig_size[1])

    # non-ASD subplot
    subplt_idx = 0
    mu2use = grand_mu
    sd2use = grand_sd/2
    plot_pop_ksdensity(grp = "non-ASD", mu = mu2use,
        sd = sd2use, n4plot = n4plot, xlimits = xlimits,
        gridline_width = gridline_width, fig_size = fig_size,
        axarr = axarr, sp_idx = subplt_idx)

    # ASD subplot
    subplt_idx = 1
    mu2use = grand_mu
    sd2use = grand_sd
    plot_pop_ksdensity(grp = "ASD", mu = mu2use,
        sd = sd2use, n4plot = n4plot, xlimits = xlimits,
        gridline_width = gridline_width, fig_size = fig_size,
        axarr = axarr, sp_idx = subplt_idx)

    # ASD subgrps subplot
    subplt_idx = 2
    plot_subgrp_ksdensity(mu_subgrp = mu_subgrp, sd = sd, n4plot = n4plot,
        xlimits = xlimits, gridline_width = gridline_width,
        grp_type = "ASDsubgrps", fig_size = fig_size, axarr = axarr,
        sp_idx = subplt_idx)



# function to plot sample subgrp prevalences as histogram
def make_sample_subgrp_prevalence_plot(results, sample_size, axarr = None,
    sp_idx = None, gridline_width = 0.5, nbins = 20, xlimits = [-0.05, 0.6],
    fig_size = [12,12],
    legend_labels = ["ASD1","ASD2","ASD3", "ASD4","ASD5"]):
    """
    Plot sample subgrp prevalences as histogram.
    """

    # key for retrieving data from results dictionary
    key = "n%d_subgrp_prevalence" % sample_size

    if axarr is None:
        # make figure specific size
        plt.figure(figsize = fig_size)
        plt.hist(results[key], bins = nbins)

        # add legend
        plt.legend(legend_labels)

        # add grid lines
        plt.grid(linewidth = gridline_width)

        # set x-axis limits
        plt.xlim(xlimits)

        # add x and y-axis labels
        plt.xlabel("Sample Prevalence")
        plt.ylabel("Count")

    else:
        x = sp_idx[0]
        y = sp_idx[1]

        data2plot = results[key]
        axarr[x, y].hist(data2plot, bins = nbins)

        # insert subplot title
        title_str = "n = %s" % (str(sample_size))
        axarr[x, y].set_title(title_str)

        # add legend
        axarr[x, y].legend(legend_labels)

        # add grid lines
        axarr[x, y].grid(linewidth = gridline_width)

        # set x-axis limits
        axarr[x, y].set_xlim(xlimits)

        # add x and y-axis labels
        axarr[x, y].set_xlabel("Sample Prevalence")
        axarr[x, y].set_ylabel("Count")

    # set tight layout to eliminate overlap of subplots
    plt.tight_layout()



def make_all_sample_prevalence_subplots(results, sample_sizes,
    gridline_width = 0.5, nbins = 20, xlimits = [-0.05, 0.6],
    fig_size = [12,12], nrows = 3, ncols = 2,
    legend_labels = ["ASD1","ASD2","ASD3", "ASD4","ASD5"]):

    # Four axes, returned as a 2-d array
    f, axarr = plt.subplots(nrows, ncols)
    f.set_size_inches(fig_size[0],fig_size[1])

    sp_mat = np.arange(len(sample_sizes)).reshape([nrows,ncols])

    # loop over rows and columns in subplot
    for irow in range(0,nrows):
        for icol in range(0,ncols):

            # grab input for making specific subplot
            ss_idx = sp_mat[irow,icol]
            sample_size = sample_sizes[ss_idx]
            sp_idx = [irow,icol]
            make_sample_subgrp_prevalence_plot(results = results,
                sample_size = sample_size, axarr = axarr, sp_idx = sp_idx,
                gridline_width = 0.5, nbins = 20, xlimits = [-0.05, 0.6],
                legend_labels = legend_labels)



# boilerplate code to call main code for executing
if __name__ == '__main__':

    # parse arguments
    opts = parse_args()

    # n experiments to simulate
    n_exp = np.array(opts.n_exp, dtype = int)

    # n subgrps
    n_subgrp = np.array(opts.n_subgrp, dtype = int)

    # sample size of each subgroup
    subgrp_pop_n = np.array(opts.subgrp_pop_n, dtype = int)

    # mean for ASD subgrps from -1 to 1 in intervals of 0.5
    # mu_subgrp = np.linspace(-1,1,n_subgrp)
    mu_subgrp = opts.mu_subgrp
    mu_subgrp = np.array(mu_subgrp.split(','),dtype = float)

    # sample_sizes
    sample_sizes = opts.sample_sizes
    sample_sizes = np.array(sample_sizes.split(','), dtype = int)

    # population SD
    sd = np.array(opts.pop_sd, dtype = int)

    # total population n
    pop_n = n_subgrp*subgrp_pop_n


    # generate ASD data
    [asd_data, asd_data_stacked] = generate_asd_data(mu_subgrp, sd,
        subgrp_pop_n, n_subgrp)

    # generate non-ASD data
    nonasd_data = generate_nonasd_data(asd_data_stacked.mean(),
        asd_data_stacked.std(), pop_n)

    # make ksdensity plots
    make_ksdensity_subplots(grand_mu = asd_data_stacked.mean(),
        grand_sd = asd_data_stacked.std(), mu_subgrp = mu_subgrp, sd = sd,
        n4plot = 10000, xlimits = [-5,5], gridline_width = 0.5,
        fig_size = [12,12], nrows = 3, ncols = 1)
    # save figure
    if opts.ksd2save  is not None:
        plt.savefig(opts.ksd2save)

    # run main simulation over range of sample sizes
    # n_exp = 10000
    # sample_sizes = [20, 50, 100, 200, 1000, 2000]
    results = run_simulation_over_sample_sizes(pop1_data_stacked = asd_data_stacked,
        pop2_data = nonasd_data, n_exp = n_exp, sample_sizes = sample_sizes,
        pop1_data = asd_data)

    # # Saving the objects:
    # fname2save = "/Users/mvlombardo/Dropbox/Jupyter/objs.pickle"
    # with open(fname2save, 'w') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([results], f)

    # # Getting back the objects:
    # with open(fname2save) as f:  # Python 3: open(..., 'rb')
    #     results = pickle.load(f)

    # make sample prevalence plot
    # make_sample_subgrp_prevalence_plot(results, sample_size = 20, subplt = False,
    #     sp_idx = None, gridline_width = 0.5, nbins = 20, xlimits = [-0.05, 0.6],
    #     legend_labels = ["ASD1","ASD2","ASD3", "ASD4","ASD5"])

    make_all_sample_prevalence_subplots(results = results,
        sample_sizes = sample_sizes, gridline_width = 0.5, nbins = 20,
        xlimits = [-0.05, 0.6], fig_size = [12,12], nrows = 3, ncols = 2,
        legend_labels = ["ASD1","ASD2","ASD3", "ASD4","ASD5"])
    # save figure
    if opts.pdf2save is not None:
        plt.savefig(opts.pdf2save)
