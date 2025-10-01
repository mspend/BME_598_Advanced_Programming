#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 11:36:00 2025

@author: maurispendlove
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('InClassActivity5.2.pdf') as pdf:

    ###################################
    # Normal Distribution
    ###################################
    
    # Stats.norm is a normal continuous random variable
    # The location (``loc``) keyword specifies the mean.
    # The scale (``scale``) keyword specifies the standard deviation.
    
    sample = stats.norm.rvs(loc=0, scale=1, size=1000)
    
    print(sample.size)
    
    print("The sample mean is " + str(sample.mean()))
    
    print('The sample standard deviation is ' , str(sample.std()))
    
    # initialize a figureand set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Create and plot a curve based on PDF
    x = np.linspace(stats.norm.ppf(0.01), stats.norm.ppf(0.99), 100)
    ax.plot(x, stats.norm.pdf(x), 'r-', lw=5, alpha=0.6, label='empirical')
    
    # Create the histogram:
    ax.hist(sample, density=True, bins='auto', histtype='stepfilled')
    ax.set_xlim([x[0], x[-1]])
    ax.legend(loc='best', frameon=False)
    
    # Add a title and axis labels
    fig.suptitle('A Standard Normal Distribution')
    ax.set_xlabel('Randomly Sampled Numbers')
    ax.set_ylabel('Frequency')
    
    pdf.savefig()
    
    ###################################
    # Binomial distribution
    ###################################
    
    # This distribution appears when you repeat the same kind of an experiment multiple (but fixed number of) times, and you record whether a certain event happened (which we call a success) or didn’t (which we call a failure).
    
    # discrete
    # n = number of trials
    # p = probability of success in each trial
    
    # Generate 1000 samples from a binomial distribution (n=10, p=0.5)
    binom_sample = stats.binom.rvs(n=10, p=0.5, size=1000)
    
    # initialize a figure and set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Add a title and legends
    fig.suptitle('A Binomial Distribution')
    ax.set_xlabel('Number of Successes')
    ax.set_ylabel('Frequency or Probability Density')
    
    # Create the histogram:
    ax.hist(binom_sample, density=True, bins=range(12), histtype='stepfilled')
    pdf.savefig()
    plt.close()
    
    ###################################
    # Poisson distribution
    ###################################
    
    # The Poisson random variable counts the number of successes in 
    #  independent Bernoulli trials in the limit as 
    #  and 
    #  where the probability of success in each trial is 
    #  and 
    #  is a constant. It can be used to approximate the Binomial random variable or in its own right to count the number of events that occur in the interval 
    #  for a process satisfying certain “sparsity” constraints
    
    # lambda or mu is the average number of events in a fixed interval (mean), or the rate parameter
    
    # Generate 1000 samples from a poisson distribution (lambda=5)
    poisson_sample = stats.poisson.rvs(mu=5, size=1000)
    
    # initialize a figure and set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Add a title and legends
    fig.suptitle('A Poisson Distribution')
    ax.set_xlabel('Number of Successes')
    ax.set_ylabel('Frequency or Probability Density')
    
    # Create the histogram:
    ax.hist(poisson_sample, density=True, bins=range(12), histtype='stepfilled')
    pdf.savefig()
    plt.close()
    
    
    # Create a combined histogram
    fig, axs = plt.subplots(1, 2, figsize=(6, 4),layout='constrained')
    
    axs[0].hist(binom_sample, density=True, bins=range(12), histtype='stepfilled',color='skyblue',label='Binomial')
    
    axs[0].set_xlabel('Number of Successes')
    axs[0].set_ylabel('Frequency or Probability Density')
    axs[0].set_title('Binomial Distribution')
    
    
    axs[1].hist(poisson_sample, density=True, bins=range(12), histtype='stepfilled',color='green',label='Poisson')
    
    axs[1].set_xlabel('Number of Successes')
    axs[1].set_ylabel('Frequency or Probability Density')
    axs[1].set_title('Poisson Distribution')
    
    fig.legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
    pdf.savefig()
    plt.close() 
    
    ###################################
    # Uniform Distribution
    ###################################
    
    # In a uniform distribution, the probabiliyt of each outcome is equal
    uniform_sample = stats.uniform.rvs(loc=0, scale=1, size=30)
    
    print(uniform_sample.mean())# Calculate sample mean and standard deviation
    
    sample_mean = uniform_sample.mean()
    sample_std = uniform_sample.std()
    print(f"Mean: {sample_mean:.3f}, Std Dev: {sample_std:.3f}")
    
    # initialize a figure
    fig, ax = plt.subplots(1, 1)
    
    # Plot a histogram
    ax.hist(uniform_sample, density=True, bins=10, histtype='stepfilled', alpha=0.6, label='Uniform Sample')
    
    
    # Overlay normal distribution with same mean and std dev
    x = np.linspace(0, 1, 200)  # x-values from 0 to 1 (range of uniform)
    normal_pdf = stats.norm.pdf(x, loc=sample_mean, scale=sample_std)
    ax.plot(x, normal_pdf, 'r-', label='Normal (matched mean/std)')
    
    # Labels and legend
    ax.set_xlabel('Value')
    ax.set_ylabel('Probability Density')
    ax.set_title('Histogram of Uniform Sample with Normal Overlay')
    ax.legend()
    
    pdf.savefig()
    plt.close()
    
    
    
    ###################################
    # Exponential Distribution
    ###################################
    
    # lambda = the rate at which events occur
    # scale = 1/lambda
    expon_sample = stats.expon.rvs(loc=1, scale=1/2, size=1000)
    
    # initialize a figure and set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Add a title and legends
    fig.suptitle('An Exponential Distribution')
    ax.set_xlabel('Value')
    ax.set_ylabel('Probability Density')
    
    # Create the histogram:
    ax.hist(expon_sample, density=True, bins='auto', histtype='stepfilled')
    pdf.savefig()
    plt.close()
    
    
    
    # Back to normal distribution
    
    mu = 100
    sigma = 15
    normal_sample = stats.norm.rvs(loc=mu, scale=sigma, size=10000)
    
    # probability that a value is less than 120
    prob = stats.norm.cdf(120, loc=mu, scale=sigma)
    
    # 95th percentile
    percentile = stats.norm.ppf(0.95, loc=mu, scale=sigma)
    
    
    # initialize a figure and set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Add a title and legends
    fig.suptitle('A Normal Distribution of Blood Pressure')
    ax.set_xlabel('Blood pressure')
    ax.set_ylabel('Probability Density')
    
    # Create the histogram:
    ax.hist(normal_sample, density=True, bins='auto', histtype='stepfilled')
    
    # Add extra features
    ax.axvline(x=120, color='red', linestyle='--', linewidth=2, label='Probability of blood pressure < 120')
    
    ax.axvline(x=percentile, color='skyblue', linestyle='--', linewidth=2, label='95th Percentile')
    plt.legend()
    pdf.savefig()
    plt.close()
    
    
    # Fit blood pressure data
    patient_bps = [118, 125, 130, 110, 135, 142, 128, 120, 138, 145, 132, 126, 140, 150, 122, 134, 129, 136, 144, 127, 131, 119, 133, 141, 137, 124, 121, 139, 147, 143]
    
    (mu, sigma) = stats.norm.fit(patient_bps)
    
    # probability that a value is greater than 140
    prob_over_140 = 1- stats.norm.cdf(140, loc=mu, scale=sigma)
    
    # initialize a figure and set its axes from 0-1
    fig, ax = plt.subplots(1, 1)
    
    # Create the histogram:
    ax.hist(patient_bps, density=True, bins='auto', histtype='stepfilled')
    
    
    # Add a title and legends
    fig.suptitle('Patient Blood Pressures')
    ax.set_xlabel('Blood pressure')
    ax.set_ylabel('Probability Density')
    
    # Overlay fitted normal PDF
    x = np.linspace(min(patient_bps)-5, max(patient_bps)+5, 200)
    fitted_pdf = stats.norm.pdf(x, loc=mu, scale=sigma)
    ax.plot(x, fitted_pdf, 'r-', label='Fitted Normal PDF')

    
    # Add vertical line at 140 mmHg
    ax.axvline(x=140, color='black', linestyle='--', linewidth=2, label='High Blood Pressure')
    
    # Add legend
    plt.legend()
    
    pdf.savefig()
    plt.close()

