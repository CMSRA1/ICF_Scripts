import plot_grabber as grabr
from itertools import product
import ROOT as r
from time import sleep

def normalise(hist = None):
    entries = hist.Integral()
    if entries:
        hist.Scale(float(1./float(entries)))

def find_max_bin(hist = None):
    max = 0.
    for i in range(1, hist.GetNbinsX()+1):
        val = hist.GetBinContent(i)
        if val > max:
            max = val
    return max

HTbins = ["200_275","275_325","325_375","375_475","475_575","575_675","675_775","775_875","875_975","975_1075","1075"][3:]
ht_scheme = ["incl", "excl"][0]
n_j = ["le3j", "ge4j", "ge2j"][:2]
n_b = ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b", "ge1b"][:2]
selec = ["OneMuon", "DiMuon", "Had"][2]
my_path = "probably_fucked_by_finejet_bug/Root_Files_16Nov_SMS_aT_0p53_v0"
my_title = "ComMinBiasDPhi"

if ht_scheme == "incl":
    my_iter = [HTbins]
    ht_string = "inclHT"
elif ht_scheme == "excl":
    my_iter = HTbins

c1 = r.TCanvas()

r.gStyle.SetOptStat(0)

colours = [r.kBlack, r.kRed, r.kGreen, r.kBlue, r.kOrange, r.kViolet, r.kYellow-2]

lsp_masses = [240, 230, 220, 210, 190, 170]
# lsp_masses = [240, 210, 170]

the_colours = {}
for lsp, col in zip(lsp_masses, colours):
    the_colours[lsp] = col

lsp_masses = [240, 230, 220, 210, 190, 170][1:2]

for htbin in my_iter:
    if ht_scheme == "excl":
        ht_string = htbin
    for nb, nj in product(n_b, n_j):
        hists = {}
        max_hist = None
        max_hist_val = 0.
        if len(lsp_masses) == 6:
            lg = r.TLegend(0.65, 0.65, 0.9, 0.9)
        elif len(lsp_masses) == 1:
            lg = r.TLegend(0.65, 0.83, 0.9, 0.9)
        for n, lspmass in enumerate(lsp_masses):
            hist = grabr.grab_plots(f_path = "/Users/chrislucas/SUSY/AnalysisCode/rootfiles/%s/Had_SMS_T2cc_250_%d.root" % (my_path, lspmass),
                                    sele = selec, h_title = my_title, njet = nj,
                                    btag = nb, ht_bins = htbin)

            hist.RebinX(10)
            hist.SetLineColor(the_colours[lspmass])
            hist.SetLineWidth(2)
            normalise(hist)

            hists[lspmass] = hist #append to list for plotting later

            # find the hist with largest y-range
            this_max = find_max_bin(hist)
            if this_max > max_hist_val:
                max_hist = hist
                max_hist_val = this_max

        #loop through hists for plotting
        max_hist.Draw("hist") #plot tallest hist first
        max_hist.GetXaxis().SetRangeUser(0., 3.2)
        for key in hists:
            print "> (250, %d): %.2f" % (key, float(hists[key].Integral(hists[key].FindBin(0.), hists[key].FindBin(0.3))/hists[key].Integral(1, 100)))
            if hists[key] == max_hist: continue
            hists[key].Draw("histsame")        
        for lsp in lsp_masses:
            lg.AddEntry(hists[lsp], "T2cc (250, %d)" % lsp, "l")

        max_hist.GetYaxis().SetTitle("Unity Normalised Events")
        max_hist.GetXaxis().SetTitleOffset(1.1)
        max_hist.GetYaxis().SetTitleOffset(1.1)
        max_hist.GetXaxis().SetLabelSize(0.03)
        max_hist.GetYaxis().SetLabelSize(0.03)
        max_hist.SetTitle("%s, %s" % (nb, nj))

        lg.SetFillColor(0)
        lg.Draw()

        # line = r.TLine()
        # line.SetLineStyle(2)
        # line.SetLineColor(r.kRed)
        # line.DrawLineNDC(.175, .1, .175, .9)

        c1.Print("dPhiStar_T2cc_250_%d_%s_%s.pdf" % (lsp_masses[0], nb, nj))
        

        # sleep(6)
        # exit()
