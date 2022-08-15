#!/usr/bin/env python3
import argparse
import os
import time
import warnings

import hyp_analysis_utils as hau
import numpy as np
import pandas as pd
import xgboost as xgb
import yaml
from analysis_classes import ModelApplication, TrainingAnalysis
from hipe4ml import analysis_utils
from hipe4ml.model_handler import ModelHandler

CENTRALITY_CLASS=[[0, 90]]
CT_BINS = [1, 2, 4, 6, 8, 10, 14, 18, 23, 35]
PT_BINS = [2, 10]
COLUMNS = ['V0CosPA','pt','ProngsDCA','PiProngPvDCAXY','He3ProngPvDCAXY','He3ProngPvDCA','PiProngPvDCA','NpidClustersHe3','TPCnSigmaHe3','TPCnSigmaPi']
APPLICATION_COLUMNS = ['score', 'm', 'ct', 'pt', 'centrality', 'Matter']

prefix = 'var2'
signal_path = '/home/fmazzasc/Hypertriton/HypertritonML/Tables/2Body/'
mc_table_name = signal_path + f'SignalTable_{prefix}.root'
df_applied_mc = hau.get_applied_mc(mc_table_name, CENTRALITY_CLASS, PT_BINS, CT_BINS, COLUMNS, APPLICATION_COLUMNS, 2, '')

print(df_applied_mc.head())
df_applied_mc.to_parquet(os.path.dirname(signal_path) + f'/applied_mc_df_{prefix}.parquet.gzip', compression='gzip')

