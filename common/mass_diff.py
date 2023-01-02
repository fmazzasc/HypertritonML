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

ROOT.gROOT.SetBatch()

np.random.seed(42)

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('config', help='Path to the YAML configuration file')
parser.add_argument('-s', '--significance', help='Use the BDTefficiency selection from the significance scan', action='store_true')
parser.add_argument('-syst', '--systematics', help='Run systematic uncertanties estimation', action='store_true')
parser.add_argument('-dbshape', '--dbshape', help='Fit using DSCBShape', action='store_true')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
###############################################################################

###############################################################################
# define some globals
FILE_PREFIX = params['FILE_PREFIX']
FILE_PREFIX_MAT = params['FILE_PREFIX'] + "_matter"
FILE_PREFIX_ANTIMAT = params['FILE_PREFIX'] + "_antimatter"

BKG_MODELS = params['BKG_MODELS'] if 'BKG_MODELS' in params else ['expo']
CENT_CLASS = params['CENTRALITY_CLASS'][0]
PT_BINS = params['PT_BINS']
CT_BINS = params['CT_BINS']

EFF_MIN, EFF_MAX, EFF_STEP = params['BDT_EFFICIENCY']
EFF_ARRAY = np.around(np.arange(EFF_MIN, EFF_MAX+EFF_STEP, EFF_STEP), 2)

SIGNIFICANCE_SCAN = args.significance
SYSTEMATICS = args.systematics
DBSHAPE = args.dbshape

SYSTEMATICS_COUNTS = 100000
FIX_EFF = 0.70 if not SIGNIFICANCE_SCAN else 0
###############################################################################

###############################################################################
# input/output files
results_dir = os.environ['HYPERML_RESULTS_{}'.format(params['NBODY'])]
efficiency_dir = os.environ['HYPERML_EFFICIENCIES_{}'.format(params['NBODY'])]

# significance scan output
file_name_sigsc = results_dir + f'/Efficiencies/{FILE_PREFIX}/sigscan.npy'
sigscan_dict = np.load(file_name_sigsc, allow_pickle=True).item()


suffix = "" if not DBSHAPE else "_dscb"
# output file
file_name = results_dir + f'/{FILE_PREFIX}_mass_diff{suffix}.root'
output_file = ROOT.TFile(file_name, 'recreate')
###############################################################################

file_name_mat = results_dir + f'/{FILE_PREFIX_MAT}_signal_extraction{suffix}.root'
signal_extr_file_mat = ROOT.TFile(file_name_mat, 'read')

file_name_antimat = results_dir + f'/{FILE_PREFIX_ANTIMAT}_signal_extraction{suffix}.root'
signal_extr_file_antimat = ROOT.TFile(file_name_antimat, 'read')


###############################################################################
start_time = time.time()
###############################################################################

###############################################################################
# define support globals
MASS_H2_MAT = {}
MASS_SHIFT_H2_MAT = {}
RECO_SHIFT_H2_MAT = {}
MASS_SHIFT_BEST_MAT = {}
RECO_SHIFT_BEST_MAT = {}
MASS_BEST_MAT = {}

MASS_H2_ANTIMAT = {}
MASS_SHIFT_H2_ANTIMAT = {}
RECO_SHIFT_H2_ANTIMAT = {}
MASS_SHIFT_BEST_ANTIMAT = {}
RECO_SHIFT_BEST_ANTIMAT = {}
MASS_BEST_ANTIMAT = {}



# support for systematics
MASS_H2_DBSHAPE = {}
RECO_SHIFT_H2_DBSHAPE = {}


MC_MASS = 2.99131
lamdba_mass_shift = 0.036




# prepare histograms for the analysis
for model in BKG_MODELS:
        MASS_H2_MAT[model] = signal_extr_file_mat.Get(f'mass_{model}')
        MASS_BEST_MAT[model] = MASS_H2_MAT[model].ProjectionX(f'mass_best_{model}_mat')
        RECO_SHIFT_H2_MAT[model] = signal_extr_file_mat.Get(f'reco_shift_{model}')
        RECO_SHIFT_BEST_MAT[model] = MASS_BEST_MAT[model].Clone("reco_shift_best_mat")

        MASS_H2_ANTIMAT[model] = signal_extr_file_antimat.Get(f'mass_{model}')
        MASS_BEST_ANTIMAT[model] = MASS_H2_ANTIMAT[model].ProjectionX(f'mass_best_{model}_antimat')
        RECO_SHIFT_H2_ANTIMAT[model] = signal_extr_file_antimat.Get(f'reco_shift_{model}')
        RECO_SHIFT_BEST_ANTIMAT[model] = MASS_BEST_ANTIMAT[model].Clone("reco_shift_best_antimat")



# helper methods
def get_eff_index(eff):
    idx = (eff - EFF_MIN + EFF_STEP) * 100
    if isinstance(eff, np.ndarray):
        return idx.astype(int)
    return int(idx)

def fill_histo_best(histo, ctbin, entry, entry_error):
    bin_idx = histo.FindBin((ctbin[0] + ctbin[1]) / 2)
    histo.SetBinContent(bin_idx, entry)
    histo.SetBinError(bin_idx, entry_error)


def get_measured_h2(h2, bkg, ctbin, eff):
    bin_idx = h2[bkg].FindBin((ctbin[0] + ctbin[1]) / 2, round(eff + 0.005, 3))
    var = h2[bkg].GetBinContent(bin_idx)
    error = h2[bkg].GetBinError(bin_idx)
    return var, error

def get_th1(eff, ctbin, file):
    histo = file.Get(f'ct{ctbin[0]}{ctbin[1]}/histo_{eff}__m')
    histo.SetDirectory(0)
    return histo



# def get_h1_frame()
###############################################################################

# significance-scan/fixed efficiencies switch
# if not SIGNIFICANCE_SCAN:
# eff_best_array = np.full(len(CT_BINS) - 1, FIX_EFF)
# else:
eff_best_array = [round(sigscan_dict[f'ct{ctbin[0]}{ctbin[1]}pt210'][0], 2) for ctbin in zip(CT_BINS[:-1], CT_BINS[1:])]
syst_eff_ranges = np.asarray([list(range(int(x * 100) - 10, int(x * 100) + 11)) for x in eff_best_array]) / 100

cv_list = []

for ctbin in zip(CT_BINS[:-1], CT_BINS[1:]):
    eff = round(sigscan_dict[f'ct{ctbin[0]}{ctbin[1]}pt210'][0], 2)
    print("BDT eff: ", eff, "ctbin: ", ctbin)

    for model in BKG_MODELS:

        mass_antimat, mass_error_antimat = get_measured_h2(MASS_H2_ANTIMAT, model, ctbin, eff)
        reco_shift_antimat, reco_shift_error_antimat = get_measured_h2(RECO_SHIFT_H2_ANTIMAT, model, ctbin, eff)
        fill_histo_best(RECO_SHIFT_BEST_ANTIMAT[model], ctbin, reco_shift_antimat, reco_shift_error_antimat)
        fill_histo_best(MASS_BEST_ANTIMAT[model], ctbin, mass_antimat, mass_error_antimat)
        MASS_BEST_ANTIMAT[model].Add(RECO_SHIFT_BEST_ANTIMAT[model])
        print("mass_antimat: ", mass_antimat)
        print("reco_shift_antimat: ", reco_shift_antimat)

        mass_mat, mass_error_mat = get_measured_h2(MASS_H2_MAT, model, ctbin, eff)
        print("mass_mat: ", mass_mat)
        reco_shift_mat, reco_shift_error_mat = get_measured_h2(RECO_SHIFT_H2_MAT, model, ctbin, eff)
        print("reco_shift_mat: ", reco_shift_mat)
        fill_histo_best(RECO_SHIFT_BEST_MAT[model], ctbin, reco_shift_mat, reco_shift_error_mat)
        fill_histo_best(MASS_BEST_MAT[model], ctbin, mass_mat, mass_error_mat)
        MASS_BEST_MAT[model].Add(RECO_SHIFT_BEST_MAT[model])

        print(MASS_BEST_MAT[model].GetBinContent(1))
        print(MASS_BEST_ANTIMAT[model].GetBinContent(1))


        diff_hist = MASS_BEST_MAT[model].Clone(f'diff_hist_{model}')
        diff_hist.GetYaxis().SetTitle("(m_{ {}^{3}_{#Lambda}H} - m_{{}_{#bar{#Lambda}}^{3}#bar{H}})/m_{ {}^{3}_{#Lambda}H}")
        diff_hist.Add(MASS_BEST_ANTIMAT[model], -1)
        # diff_hist.Scale(1/2991.31)
        diff_hist.Fit('pol0', 'Q')
        




        # reco_shift, reco_shift_error = get_measured_h2(RECO_SHIFT_H2_ANTIMAT, model, ctbin, eff)
        # fill_histo_best(RECO_SHIFT_BEST_ANTIMAT[model], ctbin, reco_shift, reco_shift_error)
        # fill_histo_best(MASS_BEST_ANTIMAT[model], ctbin, mass, mass_error)

        mass_mat = mass_mat + reco_shift_mat
        mass_antimat = mass_antimat + reco_shift_antimat
        mass_antimat = mass_antimat*10**(-3)
        histo_mat = get_th1(eff, ctbin, signal_extr_file_mat)
        histo_antimat = get_th1(eff, ctbin, signal_extr_file_antimat)



        mass_line_antimat = ROOT.TLine(mass_antimat, 0, mass_antimat, histo_antimat.GetMaximum())
        mass_line_antimat.SetLineColor(ROOT.kBlue)
        mass_line_antimat.SetLineStyle(ROOT.kDashed)


        mass_line_mat = ROOT.TLine(mass_mat, 0, mass_mat, histo_mat.GetMaximum())
        mass_line_mat.SetLineColor(ROOT.kRed)
        mass_line_mat.SetLineWidth(2)
        mass_line_mat.SetLineStyle(ROOT.kDashed)


        cv = ROOT.TCanvas(f'cv_{model}_ct{ctbin[0]}{ctbin[1]}', f'cv_{model}_ct{ctbin[0]}{ctbin[1]}', 800, 600)
        histo_mat.SetLineColor(ROOT.kRed)
        histo_antimat.SetLineColor(ROOT.kBlue)
        histo_mat.Draw()
        histo_antimat.Draw('same')

        mass_line_mat.Draw()
        mass_line_antimat.Draw()



        output_file.cd()
        cv.Write()


if SYSTEMATICS:
    isWritten = False
    # systematics histos
    mass_dist_syst = ROOT.TH1D('syst_mass_diff', '; (m_{ {}^{3}_{#Lambda}H} - m_{{}_{#bar{#Lambda}}^{3}#bar{H}})/m_{ {}^{3}_{#Lambda}H} ;counts', 200, -0.6 *10**(-4), 1.6*10**(-4))
    mass_dist_stat = ROOT.TH1D('stat_mass_diff', '; (m_{ {}^{3}_{#Lambda}H} - m_{{}_{#bar{#Lambda}}^{3}#bar{H}})/m_{ {}^{3}_{#Lambda}H} ;counts', 200, -0.6*10**(-4), 1.6*10**(-4))
    tmp_mass_mat = MASS_H2_MAT[BKG_MODELS[0]].ProjectionX('tmp_mass_mat')
    tmp_mass_antimat = MASS_H2_ANTIMAT[BKG_MODELS[0]].ProjectionX('tmp_mass_antimat')
    combinations = set()


    for _ in range(SYSTEMATICS_COUNTS):
        tmp_mass_mat.Reset()
        tmp_mass_antimat.Reset()

        bkg_list = []
        eff_list = []
        bkg_idx_list = []
        eff_idx_list = []

        # loop over ctbins
        for ctbin_idx in range(len(CT_BINS)-1):
            # random bkg model
            bkg_index = np.random.randint(0, len(BKG_MODELS))
            bkg_idx_list.append(bkg_index)
            bkg_list.append(BKG_MODELS[bkg_index])

            # randon BDT efficiency in the defined range
            eff = np.random.choice(syst_eff_ranges[ctbin_idx])
            eff_list.append(eff)
            eff_index = get_eff_index(eff)
            eff_idx_list.append(eff_index)

        # convert indexes into hash and if already sampled skip this combination
        combo = ''.join(map(str, bkg_idx_list + eff_idx_list))
        if combo in combinations:
            continue

        # if indexes are good measure B_{Lambda}
        ctbin_idx = 1
        ct_bin_it = iter(zip(CT_BINS[:-1], CT_BINS[1:]))


        for model, eff in zip(bkg_list, eff_list):
            ctbin = next(ct_bin_it)
            mass_mat, mass_error_mat = get_measured_h2(MASS_H2_MAT, model, ctbin, eff)
            mass_antimat, mass_error_antimat = get_measured_h2(MASS_H2_ANTIMAT, model, ctbin, eff)
            reco_shift_mat, reco_shift_error_mat = get_measured_h2(RECO_SHIFT_H2_MAT, model, ctbin, eff)
            reco_shift_antimat, reco_shift_error_antimat = get_measured_h2(RECO_SHIFT_H2_ANTIMAT, model, ctbin, eff)

            mass_mat += reco_shift_mat
            mass_error_mat = np.sqrt(mass_error_mat**2 + reco_shift_error_mat**2)
            mass_antimat += reco_shift_antimat
            mass_error_antimat = np.sqrt(mass_error_antimat**2 + reco_shift_error_antimat**2)


            tmp_mass_mat.SetBinContent(ctbin_idx, mass_mat)
            tmp_mass_mat.SetBinError(ctbin_idx, mass_error_mat)
            tmp_mass_antimat.SetBinContent(ctbin_idx, mass_antimat)
            tmp_mass_antimat.SetBinError(ctbin_idx, mass_error_antimat)

            ctbin_idx += 1

        tmp_mass_mat.Add(tmp_mass_antimat, -1)
        tmp_mass_mat.Scale(1/2991.31)

        tmp_mass_mat.Fit('pol0', 'Q')
        fit_func = tmp_mass_mat.GetFunction('pol0')
        mass_diff = fit_func.GetParameter(0)
        mass_diff_error = fit_func.GetParError(0)

        if fit_func.GetNDF() == 0:
            chi2red = 100
        else:
            chi2red = fit_func.GetChisquare()/fit_func.GetNDF()

        if chi2red < 2.:
            if chi2red<1.5 and isWritten==False:
                print(eff_list)
                tmp_mass_mat.GetYaxis().SetTitle('(m_{{}^{3}_{#Lambda}H} - m_{{}_{#bar{#Lambda}}^{3}#bar{H}})/m_{{}^{3}_{#Lambda}H}')
                tmp_mass_mat.GetXaxis().SetTitle('#it{c}t (cm)')

                tmp_mass_mat.Write()
                isWritten = True

            mass_dist_syst.Fill(mass_diff)
            mass_dist_stat.Fill(mass_diff_error)
            combinations.add(combo)

    mass_dist_syst.Write()
    mass_dist_stat.Write()



output_file.cd()
for model in BKG_MODELS:
    RECO_SHIFT_BEST_MAT[model].Write()
    MASS_BEST_MAT[model].Write()
    RECO_SHIFT_BEST_ANTIMAT[model].Write()
    MASS_BEST_ANTIMAT[model].Write()
    diff_hist.Write()

for cv in cv_list:
    cv.Write()

output_file.Close()



###############################################################################
print(f'--- analysis time: {((time.time() - start_time) / 60):.2f} minutes ---')
