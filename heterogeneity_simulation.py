# heterogeneity_simulation.py
"""
Here we will simulate ASD population with 5 subtypes. Each subtype has a 20%
prevalence in the population. Both ASD and TD populations will be set to the
same mean and sd, so there will be no overall case-control difference. We
then simulate experiments with different sample sizes, and compute what is
the sample prevalence for the 5 ASD subtypes.
"""

# import libraries
import sys
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import random

# set seed for reproducibility
np.random.seed(1)



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



# function for generating population-level non-ASD data
def generate_nonasd_data(asd_mu, asd_sd, pop_n):
    """
    Generate simulated non-ASD populaton data.
    """

    # non-ASD standard deviation is half that of ASD
    sd = asd_sd/2
    data = np.random.normal(asd_mu, sd, pop_n)

    return(data)


def randomly_select_sample(pop_data, sample_size):
    """
    Randomly sample data from population without replacement.
    """

    # make indices 1:len(pop_data)
    idx = range(len(pop_data))

    # randomly permute indices
    rand_idx = np.random.permute(idx)

    # grab the first n randomly permuted indices
    sample_idx = rand_idx[range(sample_size)]

    # sample data
    sample_data  = pop_data[idx2use]

    return([sample_data, sample_idx])


# function to simulate an experiment
def simulate_experiment(pop1_data,pop2_data,sample_size):
    """
    Simulate random sampling from population.
    """

    # randomly sample group1
    [samp1_data, samp1_idx] = randomly_select_sample(pop1_data, sample_size)

    # randomly sample group2
    [samp2_data, samp2_idx] = randomly_select_sample(pop2_data, sample_size)

    return((samp1_data, samp2_data))



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
    grp_type = "ASDsubgrps", subplt = False, sp_idx = None):
    """
    Make kernel density plots for each ASD subgrp.
    """

    # # make figure specific size
    plt.figure(figsize = fig_size)

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
        if not subplt:
            # make plot
            plt.plot(xs, density(xs), label = label2use)
        else:
            plt.subplot(sp_idx[0],sp_idx[1],sp_idx[2]).plot(xs, density(xs),
            label = label2use)
        # plt.plot(xs, density(xs), label = label2use)

    # add grid lines
    plt.grid(linewidth = gridline_width)

    # add legend
    plt.legend()

    # add x and y-axis labels
    plt.xlabel("DV")
    plt.ylabel("Count")

    # show plot
    plt.show()



# function for making kernel density plot
def plot_pop_ksdensity(grp, mu, sd, n4plot = 10000, xlimits = [-6,6],
    gridline_width = 0.5, fig_size = [12,12], subplt = False, sp_idx = None):
    """
    Make kernel density plots for each population.
    """

    # # make figure specific size
    plt.figure(figsize = fig_size)

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

    if not subplt:
        # make plot
        plt.plot(xs,density(xs), label = grp)
    else:
        plt.subplot(sp_idx[0],sp_idx[1],sp_idx[2]).plot(xs,density(xs),
        label = grp)
    # plt.plot(xs,density(xs), label = grp)

    # add grid lines
    plt.grid(linewidth = gridline_width)

    # add legend
    plt.legend()

    # add x and y-axis labels
    plt.xlabel("DV")
    plt.ylabel("Count")

    # show plot
    plt.show()


# make all ksdensity plots
def make_ksdensity_subplots(grand_mu, grand_sd, mu_subgrp, sd, n4plot = 10000,
    xlimits = [-5,5], gridline_width = 0.5, fig_size = [12,12],
    nrows = 3, ncols = 1):
    """
    Put together all ksdensity subplots into one figure.
    """

    # non-ASD subplot
    subplt_idx = 1
    mu2use = grand_mu
    sd2use = grand_sd/2
    plot_pop_ksdensity(grp = "non-ASD", mu = mu2use,
        sd = sd2use, n4plot = n4plot, xlimits = xlimits,
        gridline_width = gridline_width, fig_size = fig_size,
        subplt = True, sp_idx = [nrows,ncols,subplt_idx])

    # ASD subplot
    subplt_idx = 2
    mu2use = grand_mu
    sd2use = grand_sd
    plot_pop_ksdensity(grp = "ASD", mu = mu2use,
        sd = sd2use, n4plot = n4plot, xlimits = xlimits,
        gridline_width = gridline_width, fig_size = fig_size,
        subplt = True, sp_idx = [nrows,ncols,subplt_idx])

    # ASD subgrps subplot
    subplt_idx = 3
    plot_subgrp_ksdensity(mu_subgrp, sd, n4plot = n4plot, xlimits = xlimits,
        gridline_width = gridline_width, grp_type = "ASDsubgrps",
        fig_size = fig_size, subplt = True, sp_idx = [nrows,ncols,subplt_idx])




# boilerplate code to call main code for executing
if __name__ == '__main__':

    # n subgrps
    n_subgrp = 5

    # sample size of each subgroup
    subgrp_pop_n = 200000

    # total population n
    pop_n = n_subgrp*subgrp_pop_n

    # mean for ASD subgrps from -1 to 1 in intervals of 0.5
    mu_subgrp = np.linspace(-1,1,n_subgrp)

    # population SD
    sd = 1

    # generate ASD data
    [asd_data, asd_data_stacked] = generate_asd_data(mu_subgrp, sd,
        subgrp_pop_n, n_subgrp)

    # generate non-ASD data
    nonasd_data = generate_nonasd_data(asd_data_stacked.mean(),
        asd_data_stacked.std(), pop_n)

    # make ksdensity plots
    make_ksdensity_subplots(grand_mu = asd_data_stacked.mean(),
        grand_sd = asd_data_stacked.std(), mu_subgrp = mu_subgrp, sd = sd,
        n4plot = 10000, xlimits = [-5,5], gridline_width = 0.5, fig_size = [12,12],
        nrows = 3, ncols = 1)
