import plot_grabber as grabr
import ROOT as r
from itertools import product

def normalise(hist = None):
    entries = hist.Integral()
    if entries:
        hist.Scale(float(1./float(entries)))

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

for htbin in my_iter:
    if ht_scheme == "excl":
        ht_string = htbin
    for nb, nj in product(n_b, n_j):
        lg = r.TLegend(0.7, 0.7, 0.9, 0.9)
        for n, beef in enumerate(["_noDeadECAL", ""]):
        # for beef in ["", "_noISRRW"]:
            hist = grabr.grab_plots(f_path = "/Users/chrislucas/SUSY/AnalysisCode/rootfiles/%s/Had_SMS_T2cc_250_230%s.root" % (my_path, beef),
                                        sele = selec, h_title = my_title, njet = nj,
                                        btag = nb, ht_bins = htbin)
            # normalise(hist)
            if n==0: hist.Draw("histe")
            else:
                hist.Draw("histesame")
                hist.SetLineColor(r.kRed)
            # hist.Scale(0.000957)
            if beef:
                lg.AddEntry(hist, "No DeadECAL", "l")
            else:
                lg.AddEntry(hist, "Nominal", "l")
            # c1.SetLogy(1)
            hist.GetXaxis().SetRangeUser(0., 3.2)
            hist.RebinX(10)
            hist.SetTitle("%s_%s%s" % (nb, nj, beef))
        lg.SetFillColor(0)
        lg.Draw()
        c1.Print("%s_%s%s.pdf" % (nb, nj, beef))