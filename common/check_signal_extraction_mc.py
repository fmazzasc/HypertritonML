#!/usr/bin/env python3
import argparse
import os
import time

import hyp_analysis_utils as hau
import hyp_plot_utils as hpu
import numpy as np
import pandas as pd
import ROOT
import yaml

import math

ROOT.gROOT.LoadMacro('RooCustomPdfs/RooDSCBShape.cxx++')
# ROOT.gInterpreter.ProcessLine('RooCustomPdfs/RooDSCBShape.h')

from ROOT import RooDSCBShape

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

ROOT.gROOT.SetBatch()

np.random.seed(42)

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('config', help='Path to the YAML configuration file')
parser.add_argument('-s', '--significance',
                    help='Use the BDTefficiency selection from the significance scan', action='store_true')
parser.add_argument('-dbshape', '--dbshape',
                    help='Fit using DSCBShape', action='store_true')

parser.add_argument('-matter', '--matter', help='Run with matter', action='store_true')
parser.add_argument('-antimatter', '--antimatter', help='Run with antimatter', action='store_true')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
###############################################################################

SPLIT = ''

if args.matter:
    SPLIT = '_matter'

if args.antimatter:
    SPLIT = '_antimatter'

###############################################################################

# define some globals
FILE_PREFIX = params['FILE_PREFIX'] + SPLIT


DATA_PATH = os.path.expandvars(params['DATA_PATH'])
MC_PATH = os.path.expandvars(params['MC_PATH'])


BKG_MODELS = params['BKG_MODELS'] if 'BKG_MODELS' in params else ['expo']

CENT_CLASS = params['CENTRALITY_CLASS'][0]
PT_BINS = params['PT_BINS']
CT_BINS = params['CT_BINS']

EFF_MIN, EFF_MAX, EFF_STEP = params['BDT_EFFICIENCY']
EFF_ARRAY = np.around(np.arange(EFF_MIN, EFF_MAX+EFF_STEP, EFF_STEP), 2)

SIGNIFICANCE_SCAN = args.significance
DBSHAPE = args.dbshape

FIX_EFF = 0.70 if not SIGNIFICANCE_SCAN else 0
###############################################################################

###############################################################################
# input/output files
results_dir = os.environ['HYPERML_RESULTS_{}'.format(params['NBODY'])]
tables_dir = os.path.dirname(DATA_PATH)
efficiency_dir = os.environ['HYPERML_EFFICIENCIES_{}'.format(params['NBODY'])]

# mc file
tables_dir = os.path.dirname(MC_PATH)
file_name = tables_dir + f'/applied_mc_df_{FILE_PREFIX}.parquet.gzip'
file_name_var1 = tables_dir + '/applied_mc_df_var1.parquet.gzip'
file_name_var2 = tables_dir + '/applied_mc_df_var2.parquet.gzip'
## consider matbud variations
mc_dfs = {'nomin': pd.read_parquet(file_name, engine='fastparquet'), 'var1': pd.read_parquet(file_name_var1, engine='fastparquet'), 'var2': pd.read_parquet(file_name_var2, engine='fastparquet')}

# significance scan output
file_name = results_dir + f'/Efficiencies/{FILE_PREFIX}_sigscan.npy'
sigscan_dict = np.load(file_name, allow_pickle=True).item()

# output file
suffix = "_dscb"
file_name = results_dir + f'/{FILE_PREFIX}_signal_extraction_mc_{suffix}.root'
output_file = ROOT.TFile(file_name, 'recreate')
###############################################################################

start_time = time.time()
###############################################################################
# define support globals


RECO_SHIFT_H2 = {}
MC_MASS = 2.99131

# prepare histograms for the analysis
for model in BKG_MODELS:
    RECO_SHIFT_H2[model] = {}
    for var_type in mc_dfs.keys():
            RECO_SHIFT_H2[model][var_type] = ROOT.TH2D(f'reco_shift_{model}_{var_type}', ';#it{c}t cm;BDT efficiency; m (MeV/c^{2})',
                                len(CT_BINS) - 1, np.array(CT_BINS, dtype='double'), len(EFF_ARRAY) - 1, np.array(EFF_ARRAY, dtype='double'))


# useful methods
def get_eff_index(eff):
    idx = (eff - EFF_MIN + EFF_STEP) * 100
    if isinstance(eff, np.ndarray):
        return idx.astype(int)

    return int(idx)


def get_effscore_dict(ctbin):
    info_string = f'090_210_{ctbin[0]}{ctbin[1]}'
    file_name = efficiency_dir + f'/Eff_Score_{info_string}.npy'

    return {round(e[0], 2): e[1] for e in np.load(file_name).T}

##############################################################################


def fill_histo(histo, ctbin, eff, entry, entry_error):
    bin_idx = histo.FindBin((ctbin[0] + ctbin[1]) / 2, round(eff + 0.005, 3))

    histo.SetBinContent(bin_idx, entry)
    histo.SetBinError(bin_idx, entry_error)


def fill_reco_shift(model, var_type, ctbin, eff, shift, shift_error):
    fill_histo(RECO_SHIFT_H2[model][var_type], ctbin, eff, shift, shift_error)

###############################################################################


# significance-scan/fixed efficiencies switch
if not SIGNIFICANCE_SCAN:
    eff_best_array = np.full(len(CT_BINS) - 1, FIX_EFF)
else:
    eff_best_array = [round(sigscan_dict[f'ct{ctbin[0]}{ctbin[1]}pt210'][0], 2) for ctbin in zip(
        CT_BINS[:-1], CT_BINS[1:])]

# efficiency ranges for sampling the systematics
syst_eff_ranges = np.asarray([list(range(int(x * 100) - 10, int(x * 100) + 11)) for x in eff_best_array]) / 100
print("RANGES: ", syst_eff_ranges)

eff_best_it = iter(eff_best_array)
eff_range_it = iter(syst_eff_ranges)

###############################################################################
# actual analysis
for ctbin in zip(CT_BINS[:-1], CT_BINS[1:]):
    eff_best = next(eff_best_it)
    eff_range = next(eff_range_it)
    print("----------------------------------------------------")
    print("ct bin: ", ctbin , "eff_best: ", eff_best, "eff_range: ", eff_range)

    for var_type in mc_dfs.keys():
        print("var_type: ", var_type)
        ct_dir = output_file.mkdir(f'ct{ctbin[0]}{ctbin[1]}_{var_type}')
        ct_dir.cd()
        score_dict = get_effscore_dict(ctbin)
        mc_df = mc_dfs[var_type]
        # get data slice for this ct bin
        mc_slice = mc_df.query('@ctbin[0]<ct<@ctbin[1] and 2.960<m<3.040')

        for eff in eff_range:

            mass = ROOT.RooRealVar('m', 'm_{^{3}He+#pi}', 2.960, 3.040, 'GeV/c^{2}')
            mass.setVal(MC_MASS)

            # get the data slice as a RooDataSet
            tsd = score_dict[eff]
            mc_array = np.array(mc_slice.query('score>@tsd')['m'].values, dtype=np.float64)
            roo_mc_slice = hau.ndarray2roo(mc_array[:10000], mass) if len(mc_array) < 10000 else hau.ndarray2roo(mc_array, mass)
            mu = ROOT.RooRealVar('mu', 'hypertriton mass', 2.989, 2.993, 'GeV/c^{2}')
            sigma = ROOT.RooRealVar('sigma', 'hypertriton width', 0.0001, 0.004, 'GeV/c^{2}')
            a1 = ROOT.RooRealVar('a1', 'a1', 0, 5.)
            a2 = ROOT.RooRealVar('a2', 'a2', 0, 10.)
            n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
            n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
            signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)
            ROOT.RooMsgService.instance().setSilentMode(True)
            fit_results_mc_dscb = signal.fitTo(roo_mc_slice, ROOT.RooFit.Range(2.960, 3.040), ROOT.RooFit.NumCPU(64), ROOT.RooFit.Save())
            reco_shift = MC_MASS - mu.getVal()
            reco_shift_err = mu.getError()

            frame = mass.frame(80)
            frame.SetName(f'mc_eff{eff:.2f}_{model}_{var_type}')

            roo_mc_slice.plotOn(frame, ROOT.RooFit.Name('MC'))
            signal.plotOn(frame, ROOT.RooFit.Name('signal pdf'),
                          ROOT.RooFit.LineColor(ROOT.kBlue))

            # add info to plot
            pinfo = ROOT.TPaveText(0.537, 0.474, 0.937, 0.875, 'NDC')
            pinfo.SetBorderSize(0)
            pinfo.SetFillStyle(0)
            pinfo.SetTextAlign(30+3)
            pinfo.SetTextFont(42)
            frame.addObject(pinfo)
            frame.Write()

            # loop over the possible background models
            for model in BKG_MODELS:
                fill_reco_shift(model, var_type , ctbin, eff, reco_shift *1e3, reco_shift_err*1e3)

output_file.cd()
for var_type in mc_dfs.keys():
    for model in BKG_MODELS:
        RECO_SHIFT_H2[model][var_type].Write()

output_file.Close()
###############################################################################
print(
    f'--- analysis time: {((time.time() - start_time) / 60):.2f} minutes ---')
