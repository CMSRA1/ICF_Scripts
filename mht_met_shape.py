import ROOT as r
from array import array
from itertools import product
from copy import deepcopy

#==============================#
#   To-do
#   1. import data and MC shape hists 
#   2. plot ratios for each alphaT thresh
#   3. victory
#==============================#

r.gROOT.SetBatch(r.kTRUE)

def set_palette(name="", ncontours=100):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.50, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.81, 1.00, 0.12, 0.00, 0.00]

    stops_ = array('d', stops)
    red_ = array('d', red)
    green_ = array('d', green)
    blue_ = array('d', blue)

    npoints = len(stops_)
    r.TColor.CreateGradientColorTable(npoints, stops_, red_, green_, blue_, ncontours)
    r.gStyle.SetNumberContours(ncontours)

def jet_string(jet = ""):
    if jet == "le3j":
        return "_2"
    elif jet == "ge4j":
        return "_3"
    else:
        return "_all"

def btag_string(btag = ""):
    d = {
            "eq0b": "btag_zero",
            "eq1b": "btag_one",
            "eq2b": "btag_two",
            "eq3b": "btag_three",
            "ge0b": "",
            "ge1b": "btag_morethanzero"
        }

    return d[btag]

def alphat_thresh(bin = None):
    d = {
            1: "0p00",
            2: "0p50",
            3: "0p55",
            4: "0p65",
            5: "0p75",
            6: "0p85",
        }
    return d[bin]

def mc_samples():
    return ["Zinv", "WJets", "DY", "TTbar", "DiBoson", "SingleTop"]

def lumi(sele = "mu"):
    d = {
            "had": 18.493,
            "OneMuon": 19.131,
            "DiMuon": 19.131,
            "ph": 19.12,
    }

    print "> Lumi corr (%s): %.3f (*10.)" % (sele, d[sele])

    return d[sele]

def sb_corr(samp = ""):
    d = {
            "Zinv": .94,
            "WJets": .93,
            "DY": .94,
            "TTbar": 1.18,
            "DiBoson": 1.18,
            "SingleTop": 1.18,
    }

    print "> Sb corr (%s): %.2f" % (samp, d[samp])

    return d[samp]

def trig_eff(sele = "OneMuon", ht = "", njet = ""):

    # if inclusive jet selection, use ge4j effs
    if njet == "all":
        njet = "2"

    if "OneMuon" in sele:
        d = {"150_2": 0.872,"150_3": 0.881,
                "200_2": 0.875,"200_3": 0.881,
                "275_2": 0.878,"275_3": 0.882,
                "325_2": 0.879,"325_3": 0.884,
                "375_2": 0.881,"375_3": 0.886,
                "475_2": 0.882,"475_3": 0.888,
                "575_2": 0.884,"575_3": 0.889,
                "675_2": 0.885,"675_3": 0.890,
                "775_2": 0.886,"775_3": 0.891,
                "875_2": 0.888,"875_3": 0.890,
                "975_2": 0.887,"975_3": 0.890,
                "1075_2":0.884,"1075_3":0.896,}
    if "DiMuon" in sele:
        d = {"150_2": 0.984,"150_3": 0.984,
                "200_2": 0.985,"200_3": 0.984,
                "275_2": 0.985,"275_3": 0.984,
                "325_2": 0.986,"325_3": 0.986,
                "375_2": 0.986,"375_3": 0.985,
                "475_2": 0.986,"475_3": 0.986,
                "575_2": 0.986,"575_3": 0.986,
                "675_2": 0.987,"675_3": 0.986,
                "775_2": 0.986,"775_3": 0.986,
                "875_2": 0.987,"875_3": 0.986,
                "975_2": 0.987,"975_3": 0.988,
                "1075_2":0.987,"1075_3":0.987,}

    print "> Trig corr (%s): %.3f" % (ht+"_"+njet, d[ht+"_"+njet])

    return d[ht+"_"+njet]

def get_dirs(htbins = None, sele = "", btag = "", keyword = ""):

    btag_str = btag_string(btag)

    # if htbins is a string, make it a single length list
    if type(htbins) == type(""):
        htbins = [htbins]

    out = []
    for ht in htbins:
        this_bin = []
        if btag_str: this_bin.append(btag_str)
        this_bin.append(sele)
        this_bin.append(ht)
        if keyword: this_bin.append(keyword)
        out.append("_".join(this_bin))
    return out

def grab_plots(f_path = "", h_title = "", sele = "OneMuon", njet = "", btag = "", ht_bins = []):
    if f_path:
        f = r.TFile.Open(f_path)
    else:
        return

    h_total = None
    for d in get_dirs(htbins = ht_bins, sele = sele, btag = btag):
        h = f.Get("%s/%s_%s" % (d, h_title, jet_string(njet)).Clone()
        if "Data" not in f_path:
            # apply ht bin trig effs
            h.Scale( trig_eff(sele = sele,
                            ht = d.split("_")[-2] if "1075" != d[-4:] else d.split("_")[-1],
                            njet = h_title.split("_")[-1]) )
            h.Scale( sb_corr(f_path.split("_")[-1].split(".root")[0]) )
            h.Scale( lumi(sele)*10. )
        if not h_total:
            h_total = h.Clone()
        else:
            h_total.Add(h)

    h_total_clone = deepcopy(h_total)
    f.Close()

    return h_total_clone

# def combine_MC(hists = {}):
    
#     h_mc_total = None
#     for key in hists:
#         # assert key in mc_samples()

#         # scale by the sb correction
#         hists[key].Scale(sb_corr(key))

#         # scale by luminosity
#         hists[key].Scale(lumi()*10.)

#         if not h_mc_total:
#             h_mc_total = hists[key].Clone()
#         else:
#             h_mc_total.Add(hists[key])

#     return h_mc_total



def draw_ratio(h_d = None, h_mc = None, o_name = "out/test.pdf"):

    if not bin:
        return      

    c1 = r.TCanvas("c1", "c1", 1600, 1200)
    r.gStyle.SetOptStat(0)
    
    lg = r.TLegend(0.75, 0.730, 0.9, 0.89)

    pd1 = r.TPad("pd1", "pd1", 0., 0.3, 1., 1.)
    pd1.SetBottomMargin(0.005)
    pd1.SetRightMargin(0.05)
    pd1.Draw()

    pd2 = r.TPad("pd2", "pd2", 0., 0.02, 1., 0.3)
    pd2.SetTopMargin(0.05)
    pd2.SetBottomMargin(0.26)
    pd2.SetRightMargin(0.05)
    pd2.SetGridx(1)
    pd2.SetGridy(1)
    pd2.Draw()
    pd1.cd()

    xmin = 0.
    xmax = 1.8

    h_d.Draw("hist")
    h_d.RebinX(4)
    h_d.GetXaxis().SetRangeUser(xmin, xmax)
    
    h_mc.Draw("histsame")
    h_mc.SetLineColor(r.kRed+1)
    h_mc.RebinX(4)
    h_mc.GetXaxis().SetRangeUser(xmin, xmax)

    lg.AddEntry(h_d, "Data", "L")
    lg.AddEntry(h_mc, "MC", "L")

    lg.SetFillColor(0)
    lg.SetFillStyle(0)
    lg.SetLineColor(0)
    lg.SetLineStyle(0)
    lg.Draw()

    pd2.cd()

    h1_rat = h_d.Clone()
    h1_rat.Divide(h_mc)

    h1_rat.SetMarkerStyle(4)
    h1_rat.SetMarkerSize(.7)
    h1_rat.SetLineWidth(1)
    h1_rat.SetLineColor(r.kBlack)
    h1_rat.GetYaxis().SetTitle("Ratio")
    h1_rat.GetYaxis().SetRangeUser(0.5, 1.5)

    h1_rat.SetLabelSize(0.12, "X")
    h1_rat.SetLabelSize(0.07, "Y")
    h1_rat.SetTitleSize(0.13, "X")
    h1_rat.SetTitleSize(0.11, "Y")
    h1_rat.SetTitleOffset(0.25, "Y")
    h1_rat.SetTitleOffset(.9, "X")
    h1_rat.Draw("p")

    c1.Print(o_name)

    return h1_rat

def mhtmet_vs_at_analysis(hd = None, hmc = None, file_name = "out/out.pdf"):

    canv = r.TCanvas()

    set_palette()

    hd.Draw("colz")
    hd.SetTitle("Data")
    canv.Print(file_name+"(")

    hmc.Draw("colz")
    hd.SetTitle("MC")
    canv.Print(file_name)

    hrat = hd.Clone()
    hrat.Divide(hmc)
    hrat.SetTitle("Data/MC")
    hrat.Draw("colztext")
    r.gStyle.SetPaintTextFormat("0.3f")
    hrat.GetZaxis().SetRangeUser(0.7, 1.3)
    canv.Print(file_name+")")

def main():

    HTbins = ["200_275","275_325","325_375","375_475","475_575","575_675","675_775","775_875","875_975","975_1075","1075"]
    ht_scheme = ["incl", "excl"][1]
    n_j = ["le3j", "ge4j", "ge2j"][:2]
    n_b = ["eq0b", "eq1b", "eq2b", "eq3b", "ge0b", "ge1b"][:-2]
    selec = ["OneMuon", "DiMuon"][0]
    my_path = "Root_Files_22Oct_noMHTMET_v0"
    my_title = "MHT_MET_vs_alphaT_"

    if ht_scheme == "incl":
        my_iter = [HTbins]
        ht_string = "inclHT"
    elif ht_scheme == "excl":
        my_iter = HTbins

    out_root_file = r.TFile.Open("out/MHT_MET_vs_alphaT.root", "RECREATE")

    for htbin in my_iter:
        if ht_scheme == "excl":
            ht_string = htbin
        for nb, nj in product(n_b, n_j):
            h_data = grab_plots(f_path = "/Users/chrislucas/SUSY/AnalysisCode/rootfiles/%s/Muon_Data.root" % my_path,
                                sele = selec, h_title = my_title, njet = nj,
                                btag = nb, ht_bins = htbin)
            
            mc_hists = dict.fromkeys( mc_samples() )
            mc_total = None
            for samp in mc_hists:
                this_hist = grab_plots(f_path = "/Users/chrislucas/SUSY/AnalysisCode/rootfiles/%s/Muon_%s.root" % (my_path, samp),
                                            sele = selec, h_title = my_title+jet_string(nj),
                                            btag = nb, ht_bins = htbin)


                if not mc_total:
                    mc_total = this_hist.Clone()
                else:
                    mc_total.Add(this_hist)

            mhtmet_vs_at_analysis(hd = h_data, hmc = mc_total,
                                    file_name = "out/mht_met_vs_alphaT_mu_%s_%s_%s.pdf" % (nb, nj, ht_string))

            h_data.SetName("mht_met_vs_alphaT_data_mu_%s_%s_%s" % (nb, nj, ht_string))
            mc_total.SetName("mht_met_vs_alphaT_mc_mu_%s_%s_%s" % (nb, nj, ht_string))
            out_root_file.cd()
            h_data.Write()
            mc_total.Write()

    out_root_file.Close()

            # h_rats = {}

            # # loop through all the alphaT bins
            # for i in range(1, 7):
            #     projx_mc = mc_total.ProjectionX("_pfy", i, i)

            #     projx_data = h_data.ProjectionX("_pfydata", i, i)

            #     ratio = draw_ratio(projx_data, projx_mc, o_name = "out/mhtovmet_aTge%s_%s_%s.pdf" % (alphat_thresh(i), nj, nb))
            #     h_rats[i] = ratio

            # c1 = r.TCanvas()
            # for i in range(2, 7):
            #     h_tmp = h_rats[i].Clone()
            #     h_tmp.Divide(h_rats[1])
            #     h_tmp.Draw()
            #     h_tmp.SetTitle("Ratio(%s) / Ratio(%s)" % (alphat_thresh(i), alphat_thresh(1)))
            #     c1.Print("mhtovmet_ratioCompare_aTge%s_%s_%s.pdf" % (alphat_thresh(i), nj, nb))

if __name__ == '__main__':
    main()