#!/usr/bin/env python

# This plot takes the output for the two different 2016 eras and plots the results.

import os
import ROOT as r

# convert to /pb
scale_factor = 1e-6
scale_factor_lumiToPU = ((69200)/11246)

## PU in data
datainfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root'
datainfile_17 = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root'
datainfile_18 = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root'

## MC inst. lumi in mu barn
infile1 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2016BF.root"
infile2 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2016GH.root"

infile_17= "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2017_shifts.root"
infile_18 = "/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/pileup_2018_shifts.root"

output = "pileupweights.root"

fout = r.TFile(output, "recreate")

f1 = r.TFile(infile1)
h1star = f1.Get("pileup")
h1star.Scale(scale_factor)
#h1.Scale(scale_factor_lumiToPU)

f2 = r.TFile(infile2)
h2star = f2.Get("pileup")
h2star.Scale(scale_factor)
#h2.Scale(scale_factor_lumiToPU)

f3 = r.TFile(datainfile)
h3star = f3.Get("pileup")
h3star.Sumw2()
h3star.Scale(scale_factor)
#h3.Scale(1/scale_factor_lumiToPU)
#h3.Scale(1e-6)

h1 = r.TH1F("h1", "h1", 100, 0, 100)
h2 = r.TH1F("h2", "h2", 100, 0, 100)
h3 = r.TH1F("h3", "h3", 100, 0, 100)

for abin in range(0,100):
	
	print "bin no.", abin, h3star.GetBinContent(abin+1)
	bincontent_data = h3star.GetBinContent(abin+1)
	errorLow_data = h3star.GetBinErrorLow(abin+1)
	errorUp_data = h3star.GetBinErrorUp(abin+1)
	h3.Fill(abin, bincontent_data)
	h3.SetBinError(abin+1,errorUp_data)
	
	bincontent_MC1 = h1star.GetBinContent(abin+1)
	errorLow_MC1 = h1star.GetBinErrorLow(abin+1)
	errorUp_MC1 = h1star.GetBinErrorUp(abin+1)
	h1.Fill(abin, bincontent_MC1)
	h1.SetBinError(abin+1,errorUp_MC1)
	
	bincontent_MC2 = h2star.GetBinContent(abin+1)
	errorLow_MC2 = h2star.GetBinErrorLow(abin+1)
	errorUp_MC2 = h2star.GetBinErrorUp(abin+1)
	h2.Fill(abin, bincontent_MC2)
	h2.SetBinError(abin+1,errorUp_MC2)


r.gStyle.SetOptStat(0)
c1 = r.TCanvas("c1", "c1", 1200, 600)
c1.Divide(2,1)
c1.cd(1)

hs = h1.Clone("sum")
hs.Add(h2)
hs.SetLineColor(r.kBlack)
hs.Draw("hist")

h1.SetLineColor(r.kBlue)
h1.Draw("hist same")
h2.SetLineColor(r.kRed)
h2.Draw("hist same")

h3.SetLineColor(r.kBlack)
h3.SetMarkerColor(r.kBlack)
h3.SetMarkerStyle(8)
h3.Draw("p e same")
#h3.Draw("p esame")

hs.GetXaxis().SetTitle("Pileup")
hs.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")
hs.GetYaxis().SetTitleOffset(1.5)

l = r.TLegend(0.6, 0.5, 0.9, 0.6)
l.AddEntry(hs, "all 2016")
l.AddEntry(h1, "2016 eras B-F")
l.AddEntry(h2, "2016 eras G-H")
l.AddEntry(h3, "Run2 data 2016", "p")
l.SetBorderSize(0)
l.Draw()

hsn = hs.Clone()
h1n = h1.Clone()
h2n = h2.Clone()
#h3n = h3.Clone()
h3n = h3.Clone()


h1n.Scale(1.0/h1n.Integral())
hsn.Scale(1.0/hsn.Integral())
h2n.Scale(1.0/h2n.Integral())
#h3n.Scale(1.0/h3n.Integral())
h3n.Scale(1.0/h3n.Integral()) 

c1.cd(2)
hsn.Draw("hist")
h1n.Draw("hist same")
h2n.Draw("hist same")
#h3n.Draw("p same")
h3n.Draw("p same")

hsn.SetTitle("Pileup, normalized to 1")
hsn.GetYaxis().SetTitle("Event fraction")
hsn.SetMaximum(h1n.GetMaximum()*1.1)

l2 = l.Clone()
l2.Draw()

c1.Print("pileup_2016_eras.C")
c1.Print("pileup_2016_eras.pdf")
c1.Print("pileup_2016_eras.png")

####

c2 = r.TCanvas("c2", "c2", 1200, 600)
c2.Divide(2,1)
c2.cd(1)
hs1r = h3n.Clone()

#print h1n.GetXaxis().GetNbins()
print hs1r.GetXaxis().GetNbins()
hs1r.GetXaxis().SetRangeUser(0., 100.)


hs1r.Divide(h1n)
hs1r.Draw()
hs1r.SetTitle("Reweighting overall #rightarrow B-F")
#hs1r.SetLineColor(r.kBlue)
#hs1r.SetMarkerColor(r.kBlue)
hs1r.GetXaxis().SetTitle("Pileup")
hs1r.GetYaxis().SetTitle("Weight factor")
hs1r.GetYaxis().SetRangeUser(-3, 5)

c2.cd(2)
hs2r = h3n.Clone()

hs2r.GetXaxis().SetRangeUser(0., 100.)
h3n.GetXaxis().SetRangeUser(0., 100.)
hs2r.Divide(h2n)
hs2r.Draw()
#hs1r.SetLineColor(r.kRed)
#hs1r.SetMarkerColor(r.kRed)
hs2r.SetTitle("Reweighting overall #rightarrow G-H")
hs2r.GetXaxis().SetTitle("Pileup")
hs2r.GetYaxis().SetTitle("Weight factor")
hs2r.GetYaxis().SetRangeUser(-3, 5)

c2.Print("ratios_2016.C")
c2.Print("ratios_2016.pdf")
c2.Print("ratios_2016.png")

fout.cd()
hs1r.Write("puweight_2016_HIPM")
hs2r.Write("puweight_2016")


	### 2017
if True:
	f1_17 = r.TFile(infile_17)
	h1star_17 = f1_17.Get("pileup")
	h1star_17.Scale(scale_factor)

	f3_17 = r.TFile(datainfile_17)
	h3star_17 = f3_17.Get("pileup")
	h3star_17.Sumw2()
	h3star_17.Scale(scale_factor)

	h1_17 = r.TH1F("h1_17", "h1_17", 100, 0, 100)
	h3_17 = r.TH1F("h3_17", "h3_17", 100, 0, 100)

	for abin in range(0,100):
		
		bincontent_data = h3star_17.GetBinContent(abin+1)
		errorLow_data = h3star_17.GetBinErrorLow(abin+1)
		errorUp_data = h3star_17.GetBinErrorUp(abin+1)
		h3_17.Fill(abin, bincontent_data)
		h3_17.SetBinError(abin+1,errorUp_data)
		
		bincontent_MC1 = h1star_17.GetBinContent(abin+1)
		errorLow_MC1 = h1star_17.GetBinErrorLow(abin+1)
		errorUp_MC1 = h1star_17.GetBinErrorUp(abin+1)
		h1_17.Fill(abin, bincontent_MC1)
		h1_17.SetBinError(abin+1,errorUp_MC1)


	r.gStyle.SetOptStat(0)
	c3 = r.TCanvas("c3", "c3", 1200, 600)
	c3.Divide(2,1)
	c3.cd(1)
	
	h1_17.SetLineColor(r.kBlue)
	h1_17.Draw("hist same")

	h3_17.SetLineColor(r.kBlack)
	h3_17.SetMarkerColor(r.kBlack)
	h3_17.SetMarkerStyle(8)
	h3_17.Draw("p e same")

	h1_17.GetXaxis().SetTitle("Pileup")
	h1_17.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")

	l3 = r.TLegend(0.6, 0.5, 0.9, 0.6)
	l3.AddEntry(h1_17, "2017")
	l3.AddEntry(h3_17, "Run2 data 2017", "p")
	l3.SetBorderSize(0)
	l3.Draw()

	h1n_17 = h1_17.Clone()
	h3n_17 = h3_17.Clone()

	h1n_17.Scale(1.0/h1n_17.Integral())
	h3n_17.Scale(1.0/h3n_17.Integral()) 

	c3.cd(2)
	h1n_17.Draw("hist")
	h3n_17.Draw("p same")

	h1n_17.SetTitle("Pileup, normalized to 1")
	h1n_17.GetYaxis().SetTitle("Event fraction")
	h1n_17.SetMaximum(h1n_17.GetMaximum()*1.1)

	l3.Draw()

	c3.Print("pileup_2017_eras.C")
	c3.Print("pileup_2017_eras.pdf")
	c3.Print("pileup_2017_eras.png")

	####

	c4 = r.TCanvas("c4", "c4", 1200, 600)

	hs1r_17 = h3n_17.Clone()
	hs1r_17.GetXaxis().SetRangeUser(0, 100)

	hs1r_17.Divide(h1n_17)
	hs1r_17.Draw()
	hs1r_17.SetTitle("Reweighting overall")
	hs1r_17.GetXaxis().SetTitle("Pileup")
	hs1r_17.GetYaxis().SetTitle("Weight factor")
	hs1r_17.GetYaxis().SetRangeUser(-1, 3)

	c4.Print("ratios_2017.C")
	c4.Print("ratios_2017.pdf")
	c4.Print("ratios_2017.png")

	fout.cd()
	#fout.Write()
	hs1r_17.Write("puweight_2017")
	
	### 2018
if True:
	f1_18 = r.TFile(infile_18)
	h1star_18 = f1_18.Get("pileup")
	h1star_18.Scale(scale_factor)

	f3_18 = r.TFile(datainfile_18)
	h3star_18 = f3_18.Get("pileup")
	h3star_18.Sumw2()
	h3star_18.Scale(scale_factor)

	h1_18 = r.TH1F("h1_18", "h1_18", 100, 0, 100)
	h3_18 = r.TH1F("h3_18", "h3_18", 100, 0, 100)

	for abin in range(0,100):
		
		bincontent_data = h3star_18.GetBinContent(abin+1)
		errorLow_data = h3star_18.GetBinErrorLow(abin+1)
		errorUp_data = h3star_18.GetBinErrorUp(abin+1)
		h3_18.Fill(abin, bincontent_data)
		h3_18.SetBinError(abin+1,errorUp_data)
		
		bincontent_MC1 = h1star_18.GetBinContent(abin+1)
		errorLow_MC1 = h1star_18.GetBinErrorLow(abin+1)
		errorUp_MC1 = h1star_18.GetBinErrorUp(abin+1)
		h1_18.Fill(abin, bincontent_MC1)
		h1_18.SetBinError(abin+1,errorUp_MC1)


	r.gStyle.SetOptStat(0)
	c3 = r.TCanvas("c3", "c3", 1200, 600)
	c3.Divide(2,1)
	c3.cd(1)
	
	h1_18.SetLineColor(r.kBlue)
	h1_18.Draw("hist same")

	h3_18.SetLineColor(r.kBlack)
	h3_18.SetMarkerColor(r.kBlack)
	h3_18.SetMarkerStyle(8)
	h3_18.Draw("p e same")

	h1_18.GetXaxis().SetTitle("Pileup")
	h1_18.GetYaxis().SetTitle("Recorded lumi (pb^{-1})")

	l3 = r.TLegend(0.6, 0.5, 0.9, 0.6)
	l3.AddEntry(h1_18, "2018")
	l3.AddEntry(h3_18, "Run2 data 2018", "p")
	l3.SetBorderSize(0)
	l3.Draw()

	h1n_18 = h1_18.Clone()
	h3n_18 = h3_18.Clone()

	h1n_18.Scale(1.0/h1n_18.Integral())
	h3n_18.Scale(1.0/h3n_18.Integral()) 

	c3.cd(2)
	h1n_18.Draw("hist")
	h3n_18.Draw("p same")

	h1n_18.SetTitle("Pileup, normalized to 1")
	h1n_18.GetYaxis().SetTitle("Event fraction")
	h1n_18.SetMaximum(h1n_18.GetMaximum()*1.1)

	l3.Draw()

	c3.Print("pileup_2018_eras.C")
	c3.Print("pileup_2018_eras.pdf")
	c3.Print("pileup_2018_eras.png")

	####

	c4 = r.TCanvas("c4", "c4", 1200, 600)

	hs1r_18 = h3n_18.Clone()
	hs1r_18.GetXaxis().SetRangeUser(0., 100.)

	hs1r_18.Divide(h1n_18)
	hs1r_18.Draw()
	hs1r_18.SetTitle("Reweighting overall")
	hs1r_18.GetXaxis().SetTitle("Pileup")
	hs1r_18.GetYaxis().SetTitle("Weight factor")
	hs1r_18.GetYaxis().SetRangeUser(-1, 3)

	c4.Print("ratios_2018.C")
	c4.Print("ratios_2018.pdf")
	c4.Print("ratios_2018.png")

	fout.cd()
	#fout.Write()
	hs1r_18.Write("puweight_2018")

print "just created", fout.GetName()
fout.Close()

print "Press ENTER to exit..."
raw_input()
