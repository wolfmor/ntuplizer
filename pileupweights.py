#!/usr/bin/env python

# This plot takes the output for the two different 2016 eras and plots the results.

import os
import sys
import ROOT as r

# convert to /pb
scale_factor = 1e-6
scale_factor_lumiToPU = ((69200) / 11246)

## PU in data
datainfile_16preVFP = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root'
datainfile_16postVFP = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root'
datainfile_16all = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root'
datainfile_17 = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root'
datainfile_18 = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root'

## MC inst. lumi in mu barn
# infile1 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2016BF.root"
# infile2 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2016GH.root"
# infile_17 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2017_shifts.root"
# infile_18 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2018_shifts.root"

# # these are the distributions of the n_trueInteractions variable
# from: https://raw.githubusercontent.com/cms-sw/cmssw/a65c2e1a23f2e7fe036237e2e34cda8af06b8182/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py
mc_profile_16 = [
    1.00402360149e-05, 5.76498797172e-05, 7.37891400294e-05, 0.000110932895295, 0.000158857714773,
    0.000368637432599, 0.000893114107873, 0.00189700774575, 0.00358880167437, 0.00636052573486,
    0.0104173961179, 0.0158122597405, 0.0223785660712, 0.0299186888073, 0.0380275944896,
    0.0454313901624, 0.0511181088317, 0.0547434577348, 0.0567906239028, 0.0577145461461,
    0.0578176902735, 0.0571251566494, 0.0555456541498, 0.053134383488, 0.0501519041462,
    0.0466815838899, 0.0429244592524, 0.0389566776898, 0.0348507152776, 0.0307356862528,
    0.0267712092206, 0.0229720184534, 0.0193388653099, 0.0159602510813, 0.0129310510552,
    0.0102888654183, 0.00798782770975, 0.00606651703058, 0.00447820948367, 0.00321589786478,
    0.0022450422045, 0.00151447388514, 0.000981183695515, 0.000609670479759, 0.000362193408119,
    0.000211572646801, 0.000119152364744, 6.49133515399e-05, 3.57795801581e-05, 1.99043569043e-05,
    1.13639319832e-05, 6.49624103579e-06, 3.96626216416e-06, 2.37910222874e-06, 1.50997403362e-06,
    1.09816650247e-06, 7.31298519122e-07, 6.10398791529e-07, 3.74845774388e-07, 2.65177281359e-07,
    2.01923536742e-07, 1.39347583555e-07, 8.32600052913e-08, 6.04932421298e-08, 6.52536630583e-08,
    5.90574603808e-08, 2.29162474068e-08, 1.97294602668e-08, 1.7731096903e-08, 3.57547932012e-09,
    1.35039815662e-09, 8.50071242076e-09, 5.0279187473e-09, 4.93736669066e-10, 8.13919708923e-10,
    5.62778926097e-09, 5.15140589469e-10, 8.21676746568e-10, 0.0, 1.49166873577e-09,
    8.43517992503e-09, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
]
# from: https://raw.githubusercontent.com/cms-sw/cmssw/435f0b04c0e318c1036a6b95eb169181bbbe8344/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py
mc_profile_17 = [
    1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
    0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
    0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
    0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
    0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
    0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
    0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
    0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
    0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
    0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
    0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
    0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
    0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
    3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
    1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
]
# from: https://raw.githubusercontent.com/cms-sw/cmssw/a65c2e1a23f2e7fe036237e2e34cda8af06b8182/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py
mc_profile_18 = [
    8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
    0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
    0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
    0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
    0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
    0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
    0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
    0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
    0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
    0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
    0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
    0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
    0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
    2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
    3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
    4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
]

output = "pileupweights.root"

fout = r.TFile(output, "recreate")


### 2016

# f1 = r.TFile(infile1)
# h1star = f1.Get("pileup")
# h1star.Scale(scale_factor)
# # h1.Scale(scale_factor_lumiToPU)
#
# f2 = r.TFile(infile2)
# h2star = f2.Get("pileup")
# h2star.Scale(scale_factor)
# # h2.Scale(scale_factor_lumiToPU)

f16preVFP = r.TFile(datainfile_16preVFP)
h16preVFP = f16preVFP.Get("pileup")
h16preVFP.Sumw2()
h16preVFP.Scale(scale_factor)
# h3.Scale(1/scale_factor_lumiToPU)
# h3.Scale(1e-6)

f16postVFP = r.TFile(datainfile_16postVFP)
h16postVFP = f16postVFP.Get("pileup")
h16postVFP.Sumw2()
h16postVFP.Scale(scale_factor)

h16mc_norm = r.TH1F("h16mc_norm", "h16mc_norm", 99, 0, 99)
# h2 = r.TH1F("h2", "h2", 100, 0, 100)
# h3 = r.TH1F("h3", "h3", 100, 0, 100)

for abin in range(99):

    h16mc_norm.SetBinContent(abin+1, mc_profile_16[abin])
    h16mc_norm.SetBinError(abin+1, 0.)

    # print "bin no.", abin, h3star.GetBinContent(abin + 1)
    # bincontent_data = h3star.GetBinContent(abin + 1)
    # errorLow_data = h3star.GetBinErrorLow(abin + 1)
    # errorUp_data = h3star.GetBinErrorUp(abin + 1)
    # h3.Fill(abin, bincontent_data)
    # h3.SetBinError(abin + 1, errorUp_data)
    #
    # bincontent_MC1 = h1star.GetBinContent(abin + 1)
    # errorLow_MC1 = h1star.GetBinErrorLow(abin + 1)
    # errorUp_MC1 = h1star.GetBinErrorUp(abin + 1)
    # h1.Fill(abin, bincontent_MC1)
    # h1.SetBinError(abin + 1, errorUp_MC1)
    #
    # bincontent_MC2 = h2star.GetBinContent(abin + 1)
    # errorLow_MC2 = h2star.GetBinErrorLow(abin + 1)
    # errorUp_MC2 = h2star.GetBinErrorUp(abin + 1)
    # h2.Fill(abin, bincontent_MC2)
    # h2.SetBinError(abin + 1, errorUp_MC2)

r.gStyle.SetOptStat(0)
c1 = r.TCanvas("c1", "c1", 1200, 600)
c1.Divide(2, 1)
c1.cd(1)

h16all = h16preVFP.Clone("sum")
h16all.Add(h16postVFP)
h16all.SetLineColor(r.kBlack)
h16all.Draw("hist")

h16preVFP.SetLineColor(r.kBlue)
h16preVFP.Draw("hist same")
h16postVFP.SetLineColor(r.kRed)
h16postVFP.Draw("hist same")

# h3.SetLineColor(r.kBlack)
# h3.SetMarkerColor(r.kBlack)
# h3.SetMarkerStyle(8)
# h3.Draw("p e same")
# # h3.Draw("p esame")

h16all.GetXaxis().SetTitle("Pileup")
h16all.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")
h16all.GetYaxis().SetTitleOffset(1.5)

l = r.TLegend(0.6, 0.5, 0.9, 0.6)
l.AddEntry(h16all, "2016 all")
l.AddEntry(h16preVFP, "2016 eras B-F")
l.AddEntry(h16postVFP, "2016 eras G-H")
# l.AddEntry(h3, "Run2 data 2016", "p")
l.SetBorderSize(0)
l.Draw()

h16all_norm = h16all.Clone()
h16preVFP_norm = h16preVFP.Clone()
h16postVFP_norm = h16postVFP.Clone()
# h3n = h3.Clone()

h16all_norm.Scale(1.0 / h16all_norm.Integral())
h16preVFP_norm.Scale(1.0 / h16preVFP_norm.Integral())
h16postVFP_norm.Scale(1.0 / h16postVFP_norm.Integral())
# h3n.Scale(1.0/h3n.Integral())
# h3n.Scale(1.0 / h3n.Integral())

c1.cd(2)
h16all_norm.Draw("hist")
h16preVFP_norm.Draw("hist same")
h16postVFP_norm.Draw("hist same")
# h3n.Draw("p same")

h16mc_norm.SetLineColor(r.kGreen+2)
h16mc_norm.Draw("hist same")

h16all_norm.SetTitle("pileup, normalized to 1")
h16all_norm.GetYaxis().SetTitle("Event fraction")
h16all_norm.SetMaximum(h16all_norm.GetMaximum() * 1.2)

l2 = l.Clone()
l2.AddEntry(h16mc_norm, "2016 MC")
l2.Draw()

c1.Print("pileup_2016_eras.C")
c1.Print("pileup_2016_eras.pdf")
c1.Print("pileup_2016_eras.png")

####

c2 = r.TCanvas("c2", "c2", 1200, 600)
c2.Divide(3, 1)
c2.cd(1)
h16preVFP_ratio = h16preVFP_norm.Clone()

# # print h1n.GetXaxis().GetNbins()
# print hs1r.GetXaxis().GetNbins()
# hs1r.GetXaxis().SetRangeUser(0., 100.)

h16preVFP_ratio.Divide(h16mc_norm)
h16preVFP_ratio.Draw()
h16preVFP_ratio.SetTitle("Reweighting 2016 MC #rightarrow Data B-F")
# hs1r.SetLineColor(r.kBlue)
# hs1r.SetMarkerColor(r.kBlue)
h16preVFP_ratio.GetXaxis().SetTitle("Pileup")
h16preVFP_ratio.GetYaxis().SetTitle("Weight factor")
h16preVFP_ratio.GetYaxis().SetRangeUser(-3, 5)

c2.cd(2)
h16postVFP_ratio = h16postVFP_norm.Clone()

# hs2r.GetXaxis().SetRangeUser(0., 100.)
# h3n.GetXaxis().SetRangeUser(0., 100.)

h16postVFP_ratio.Divide(h16mc_norm)
h16postVFP_ratio.Draw()
# hs1r.SetLineColor(r.kRed)
# hs1r.SetMarkerColor(r.kRed)
h16postVFP_ratio.SetTitle("Reweighting 2016 MC #rightarrow Data G-H")
h16postVFP_ratio.GetXaxis().SetTitle("Pileup")
h16postVFP_ratio.GetYaxis().SetTitle("Weight factor")
h16postVFP_ratio.GetYaxis().SetRangeUser(-3, 5)

c2.cd(3)
h16all_ratio = h16all_norm.Clone()

# hs2r.GetXaxis().SetRangeUser(0., 100.)
# h3n.GetXaxis().SetRangeUser(0., 100.)

h16all_ratio.Divide(h16mc_norm)
h16all_ratio.Draw()
# hs1r.SetLineColor(r.kRed)
# hs1r.SetMarkerColor(r.kRed)
h16all_ratio.SetTitle("Reweighting 2016 MC #rightarrow Data B-H")
h16all_ratio.GetXaxis().SetTitle("Pileup")
h16all_ratio.GetYaxis().SetTitle("Weight factor")
h16all_ratio.GetYaxis().SetRangeUser(-3, 5)

c2.Print("ratios_2016.C")
c2.Print("ratios_2016.pdf")
c2.Print("ratios_2016.png")

fout.cd()
h16preVFP_ratio.Write("puweight_2016_HIPM")
h16postVFP_ratio.Write("puweight_2016")
h16all_ratio.Write("puweight_2016_full")

### 2017
if True:

    # f1_17 = r.TFile(infile_17)
    # h1star_17 = f1_17.Get("pileup")
    # h1star_17.Scale(scale_factor)

    f17 = r.TFile(datainfile_17)
    h17 = f17.Get("pileup")
    h17.Sumw2()
    h17.Scale(scale_factor)

    h17mc_norm = r.TH1F("h17mc_norm", "h17mc_norm", 99, 0, 99)
    # h3_17 = r.TH1F("h3_17", "h3_17", 100, 0, 100)

    for abin in range(99):

        h17mc_norm.SetBinContent(abin+1, mc_profile_17[abin])
        h17mc_norm.SetBinError(abin+1, 0.)

        # bincontent_data = h3star_17.GetBinContent(abin + 1)
        # errorLow_data = h3star_17.GetBinErrorLow(abin + 1)
        # errorUp_data = h3star_17.GetBinErrorUp(abin + 1)
        # h3_17.Fill(abin, bincontent_data)
        # h3_17.SetBinError(abin + 1, errorUp_data)
        #
        # bincontent_MC1 = h1star_17.GetBinContent(abin + 1)
        # errorLow_MC1 = h1star_17.GetBinErrorLow(abin + 1)
        # errorUp_MC1 = h1star_17.GetBinErrorUp(abin + 1)
        # h1_17.Fill(abin, bincontent_MC1)
        # h1_17.SetBinError(abin + 1, errorUp_MC1)

    r.gStyle.SetOptStat(0)
    c3 = r.TCanvas("c3", "c3", 1200, 600)
    c3.Divide(2, 1)
    c3.cd(1)

    h17.SetLineColor(r.kBlack)
    h17.Draw("hist")

    # h3_17.SetLineColor(r.kBlack)
    # h3_17.SetMarkerColor(r.kBlack)
    # h3_17.SetMarkerStyle(8)
    # h3_17.Draw("p e same")

    h17.GetXaxis().SetTitle("Pileup")
    h17.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")

    l3 = r.TLegend(0.6, 0.5, 0.9, 0.6)
    l3.AddEntry(h17, "2017")
    # l3.AddEntry(h3_17, "Run2 data 2017", "p")
    l3.SetBorderSize(0)
    l3.Draw()

    h17_norm = h17.Clone()
    # h3n_17 = h3_17.Clone()

    h17_norm.Scale(1.0 / h17_norm.Integral())
    # h3n_17.Scale(1.0 / h3n_17.Integral())

    c3.cd(2)
    h17_norm.Draw("hist")
    # h3n_17.Draw("p same")

    h17mc_norm.SetLineColor(r.kGreen+2)
    h17mc_norm.Draw("hist same")

    h17_norm.SetTitle("pileup, normalized to 1")
    h17_norm.GetYaxis().SetTitle("Event fraction")
    h17_norm.SetMaximum(h17_norm.GetMaximum() * 1.2)

    l4 = l3.Clone()
    l4.AddEntry(h17mc_norm, "2017 MC")
    l4.Draw()

    c3.Print("pileup_2017_eras.C")
    c3.Print("pileup_2017_eras.pdf")
    c3.Print("pileup_2017_eras.png")

    ####

    c4 = r.TCanvas("c4", "c4", 1200, 600)

    h17_ratio = h17_norm.Clone()

    # hs1r_17.GetXaxis().SetRangeUser(0, 100)

    h17_ratio.Divide(h17mc_norm)
    h17_ratio.Draw()
    h17_ratio.SetTitle("Reweighting 2017 MC #rightarrow Data")
    h17_ratio.GetXaxis().SetTitle("Pileup")
    h17_ratio.GetYaxis().SetTitle("Weight factor")
    h17_ratio.GetYaxis().SetRangeUser(-1, 3)

    c4.Print("ratios_2017.C")
    c4.Print("ratios_2017.pdf")
    c4.Print("ratios_2017.png")

    fout.cd()
    # fout.Write()
    h17_ratio.Write("puweight_2017")

### 2018
if True:

    # f1_18 = r.TFile(infile_18)
    # h1star_18 = f1_18.Get("pileup")
    # h1star_18.Scale(scale_factor)

    f18 = r.TFile(datainfile_18)
    h18 = f18.Get("pileup")
    h18.Sumw2()
    h18.Scale(scale_factor)

    h18mc_norm = r.TH1F("h18mc_norm", "h18mc_norm", 99, 0, 99)
    # h3_18 = r.TH1F("h3_18", "h3_18", 100, 0, 100)

    for abin in range(99):

        h18mc_norm.SetBinContent(abin+1, mc_profile_18[abin])
        h18mc_norm.SetBinError(abin+1, 0.)

        # bincontent_data = h3star_18.GetBinContent(abin + 1)
        # errorLow_data = h3star_18.GetBinErrorLow(abin + 1)
        # errorUp_data = h3star_18.GetBinErrorUp(abin + 1)
        # h3_18.Fill(abin, bincontent_data)
        # h3_18.SetBinError(abin + 1, errorUp_data)
        #
        # bincontent_MC1 = h1star_18.GetBinContent(abin + 1)
        # errorLow_MC1 = h1star_18.GetBinErrorLow(abin + 1)
        # errorUp_MC1 = h1star_18.GetBinErrorUp(abin + 1)
        # h1_18.Fill(abin, bincontent_MC1)
        # h1_18.SetBinError(abin + 1, errorUp_MC1)

    r.gStyle.SetOptStat(0)
    c4 = r.TCanvas("c4", "c4", 1200, 600)
    c4.Divide(2, 1)
    c4.cd(1)

    h18.SetLineColor(r.kBlack)
    h18.Draw("hist")

    # h3_18.SetLineColor(r.kBlack)
    # h3_18.SetMarkerColor(r.kBlack)
    # h3_18.SetMarkerStyle(8)
    # h3_18.Draw("p e same")

    h18.GetXaxis().SetTitle("Pileup")
    h18.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")

    l5 = r.TLegend(0.6, 0.5, 0.9, 0.6)
    l5.AddEntry(h18, "2018")
    # l5.AddEntry(h3_18, "Run2 data 2018", "p")
    l5.SetBorderSize(0)
    l5.Draw()

    h18_norm = h18.Clone()
    # h3n_18 = h3_18.Clone()

    h18_norm.Scale(1.0 / h18_norm.Integral())
    # h3n_18.Scale(1.0 / h3n_18.Integral())

    c4.cd(2)
    h18_norm.Draw("hist")
    # h3n_18.Draw("p same")

    h18mc_norm.SetLineColor(r.kGreen+2)
    h18mc_norm.Draw("hist same")

    h18_norm.SetTitle("pileup, normalized to 1")
    h18_norm.GetYaxis().SetTitle("Event fraction")
    h18_norm.SetMaximum(h18_norm.GetMaximum() * 1.2)

    l5 = l4.Clone()
    l5.AddEntry(h18mc_norm, "2018 MC")
    l5.Draw()

    c4.Print("pileup_2018_eras.C")
    c4.Print("pileup_2018_eras.pdf")
    c4.Print("pileup_2018_eras.png")

    ####

    c5 = r.TCanvas("c5", "c5", 1200, 600)

    h18_ratio = h18_norm.Clone()

    # hs1r_18.GetXaxis().SetRangeUser(0., 100.)

    h18_ratio.Divide(h18mc_norm)
    h18_ratio.Draw()
    h18_ratio.SetTitle("Reweighting 2018 MC #rightarrow Data")
    h18_ratio.GetXaxis().SetTitle("Pileup")
    h18_ratio.GetYaxis().SetTitle("Weight factor")
    h18_ratio.GetYaxis().SetRangeUser(-1, 3)

    c5.Print("ratios_2018.C")
    c5.Print("ratios_2018.pdf")
    c5.Print("ratios_2018.png")

    fout.cd()
    # fout.Write()
    h18_ratio.Write("puweight_2018")

fout.Write()
print "just created", fout.GetName()
fout.Close()

print "Press ENTER to exit..."
raw_input()
