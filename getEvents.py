#!/usr/bin/env python
import ROOT as r
from sys import argv, exit

# if len(argv) < 2:
# 	print "balls"
# 	exit()

in_files = ["Zinv", "DY", "Data", "WJets", "DiBoson", "TTbar"]
selection = ["OneMuon_", "DiMuon_", "Photon_", ""][-1]

in_directory = "btag_zero_%s200_275" % selection
in_hist = "AlphaT_2"
# in_directory = "count_total"
# in_hist = "count"

c1 = r.TCanvas()

run_new = [False, True][1]

for process in in_files:
	full_path = "../rootfiles/Root_Files_17Dec_Full2013_noISRRW_noSITV_fixedXS_fixedCode/Had_" + process + ".root"
	if run_new: full_path_new = "../rootfiles/Root_Files_28Jan_Full2013_noISRRW_noSITV_fixedXS_globalAlphaT_fixedCode/backup_had_pho/Had_"+process+".root"
	# print ">>> Opening files:", full_path
	# if run_new: print full_path_new

	file_original = r.TFile.Open(full_path)
	# print file_original
	if run_new: file_new = r.TFile.Open(full_path_new)

	hist = file_original.Get(in_directory+"/"+in_hist)
	# print hist
	if run_new: hist_new = file_new.Get(in_directory+"/"+in_hist)

	higher = 50.
	print "\n>  ", process
	# print "   ", hist.GetEntries(), hist_new.GetEntries()
	# old_val = hist.Integral(int(float(0.55)/0.01)+1,int(float(higher)/0.01))
	old_val = hist.GetEntries()
	print old_val
	if run_new:
		# new_val = hist_new.Integral(int(float(0.55)/0.01)+1,int(float(higher)/0.01))
		new_val = hist_new.GetEntries()
		print new_val
		if old_val>0.: print "ratio:", new_val/old_val

	# div = hist.Clone()
	# div.Divide(hist_new)
	# div.Draw()
	# div.Rebin(10)

	if run_new:
		hist_new.Draw("hist")
		hist_new.Rebin(10)
		hist_new.SetLineColor(r.kBlue)

	hist.Draw("same")
	hist.Rebin(10)
	hist.SetMarkerStyle(3)
	hist.SetLineColor(r.kRed)

	c1.Print("%s%s.pdf" % (selection, process))
