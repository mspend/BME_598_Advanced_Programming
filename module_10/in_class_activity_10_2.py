#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:33:10 2025

@author: maurispendlove
"""

import matplotlib
# matplotlib.use('Agg')
# matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

fev = pd.read_table(filepath_or_buffer='fev_dat.txt', delimiter='\t')
