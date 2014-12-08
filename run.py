from collections import OrderedDict
from uncertainties import ufloat
import ROOT as r
import numpy as np
import math
import bisect

################################################################################
# CONFIG #######################################################################
################################################################################

#input_file = "./input/Muon_Trigger_noSecondJetThresh.root"
input_file = "./input/Muon_Data_Adam.root"

offset = 1.e-6

diff = False

ht_bins = [200,275,325,375,475,575,675,775,875,975,1075]
#yaxis_thresholds = [160.]*len(ht_bins)
yaxis_thresholds = [0.65,0.60,]+[0.55]*(len(ht_bins)-2)

#samples = OrderedDict([("phot","Photon")])

samples = OrderedDict([("mu","OneMuon"),])#("mm","DiMuon"),])
njet_bins = OrderedDict([("le3j",["2","3"]),("ge4j",["4","5"]),])#("ge2j",["all"]),])
bjet_bins = ["ge0b"]#,"ge1b","eq0b"]#,"eq1b","eq2b","eq3b","ge4b",]#"ge2b",]
cats = OrderedDict([
        #(("le3j","eq0b"),("2--3","0")),
        # (("le3j","eq1b"),("2--3","1")),
        #(("le3j","ge1b"),("2--3","$\geq 1$")),
        # (("le3j","eq2b"),("2--3","2")),
        #(("ge4j","eq0b"),("$\geq 4$","0")),
        # (("ge4j","eq1b"),("$\geq 4$","1")),
        #(("ge4j","ge1b"),("$\geq 4$","$\geq 1$")),
        # ((("ge4j","eq2b"),("$\geq 4$","2")),
        # ((("ge4j","eq3b"),("$\geq 4$","3")),
        # ((("ge4j","ge4b"),("$\geq 4$","$\geq 4$")),
        ])

#samples = OrderedDict([("mu","OneMuon")])
#njet_bins = OrderedDict([("le3j",["2","3"]),("ge4j",["4","5"])])
#bjet_bins = ["ge0b"]
#cats = OrderedDict([(("le3j","ge0b"),("2--3","$\geq 0$")),
#                    (("ge4j","ge0b"),("$\geq 4$","$\geq 0$")),
#                    ])

ht_binning = [200,275,325,375,475,575,675,775,875,975,1075]
#ht_binning = [ x*1. for x in range(0,1125,25) ]
yaxis_binning = [ x/100. for x in range(30,71,1) ]
#yaxis_binning = [ x*1. for x in range(0,210,10) ]

print "ht_binning: ",ht_binning
print "yaxis_binning: ",yaxis_binning

################################################################################
# FUNCTIONS ####################################################################
################################################################################

def rebin( input, xbinning, ybinning, title="" ) : 
    xbinning = [ x*1. for x in xbinning ]
    ybinning = [ x*1. for x in ybinning ]
    output = r.TH2D(input.GetTitle()+"_rebinned" if title == "" else title,
                    title, 
                    len(xbinning)-1, 
                    np.array(xbinning), 
                    len(ybinning)-1, 
                    np.array(ybinning) )
    xaxis = input.GetXaxis()
    yaxis = input.GetYaxis()
    for j in range(yaxis.GetNbins()+2) : 
        for i in range(xaxis.GetNbins()+2) : 
            output.Fill( xaxis.GetBinCenter(i),
                         yaxis.GetBinCenter(j),
                         input.GetBinContent(i,j) )
    return output

def cumu( input, ht_low=None, ht_high=None, at_low=None, at_high=None, title="" ) : 

    output = r.TH2D( input.GetTitle()+"_cumulative" if title == "" else title,
                     title,
                     input.GetNbinsX(),
                     np.array(input.GetXaxis().GetXbins()),
                     input.GetNbinsY(),
                     np.array(input.GetYaxis().GetXbins()) )

    xaxis = output.GetXaxis()
    yaxis = output.GetYaxis()

#    for j in range(yaxis.GetNbins()+2) : 
#        for i in range(xaxis.GetNbins()+2) : 
#            output.Fill( xaxis.GetBinCenter(i),
#                         yaxis.GetBinCenter(j),
#                         input.Integral(i,xaxis.GetNbins()+1,
#                                        j,yaxis.GetNbins()+1) )

    # If not set, use full range including overflow bins
    ht_l = xaxis.FindBin(ht_low+offset) if ht_low is not None else 0
    ht_h = xaxis.FindBin(ht_high+offset) if ht_high is not None else xaxis.GetNbins()+1
    at_l = yaxis.FindBin(at_low+offset) if at_low is not None else 0
    at_h = yaxis.FindBin(at_high+offset) if at_high is not None else yaxis.GetNbins()+1

    # If set, add overflow bins if required
    if ht_l == 1 : ht_l = 0
    if ht_h == xaxis.GetNbins() : ht_h = xaxis.GetNbins()+1
    if at_l == 1 : at_l = 0
    if at_h == yaxis.GetNbins() : at_h = yaxis.GetNbins()+1

    for iat in range( at_l, at_h+1 ) :
        for iht in range( ht_l, ht_h+1 ) :
            output.Fill( xaxis.GetBinCenter(iht),
                         yaxis.GetBinCenter(iat),
                         input.Integral(iht,ht_h,iat,at_h) )

    return output

def make_plot( numer_diff, denom_diff, 
               numer_cumu, denom_cumu, 
               iht_bin, ht_bin, ht_upper,
               title="" ) :

    r.gROOT.SetBatch(r.kTRUE)
    c1 = r.TCanvas()
    c1.cd()
    pad1 = r.TPad("pad1","",0,0,1,1)
    pad2 = r.TPad("pad2","",0,0,1,1)
    pad2.SetFillStyle(4000)
    pad1.Draw()
    pad1.cd()
    #pad1.SetLogy()
    #c1.SetLogy()
    r.gStyle.SetOptStat(0)

    numer_diff_proj = numer_diff.ProjectionY( "numer_diff_proj"+title, 
                                              numer_diff.GetXaxis().FindBin(ht_bin+offset), 
                                              ( numer_diff.GetXaxis().FindBin(ht_upper+offset) \
                                                    if ht_upper is not None else -1 ) )

    denom_diff_proj = denom_diff.ProjectionY( "denom_diff_proj"+title, 
                                              denom_diff.GetXaxis().FindBin(ht_bin+offset), 
                                              ( denom_diff.GetXaxis().FindBin(ht_upper+offset) \
                                                    if ht_upper is not None else -1 ) )

    numer_cumu_proj = numer_cumu.ProjectionY( "numer_cumu_proj"+title, 
                                              numer_cumu.GetXaxis().FindBin(ht_bin+offset), 
                                              ( numer_cumu.GetXaxis().FindBin(ht_upper+offset) \
                                                    if ht_upper is not None else -1 ) )

    denom_cumu_proj = denom_cumu.ProjectionY( "denom_cumu_proj"+title, 
                                              denom_cumu.GetXaxis().FindBin(ht_bin+offset), 
                                              ( denom_cumu.GetXaxis().FindBin(ht_upper+offset) \
                                                    if ht_upper is not None else -1 ) )
    
    if False : 

        denom_diff_proj.Draw("hist")
        denom_diff_proj.SetLineWidth(2)
        denom_diff_proj.SetLineColor(r.kBlack)

        denom_diff_proj.SetTitle("Eff"+title)
        denom_diff_proj.GetXaxis().SetTitle("#alpha_{T}")
        denom_diff_proj.GetYaxis().SetTitle("Events")
        denom_diff_proj.GetYaxis().SetRangeUser(1.,denom_diff_proj.GetMaximum()*10.)

        numer_diff_proj.Draw("histsame")
        numer_diff_proj.SetLineWidth(2)
        numer_diff_proj.SetLineStyle(2)
        numer_diff_proj.SetLineColor(r.kGreen)

    else :

        denom_cumu_proj.Draw("hist")
        denom_cumu_proj.SetLineWidth(2)
        denom_cumu_proj.SetLineColor(r.kBlack)

        denom_cumu_proj.SetTitle("Eff"+title)
        denom_cumu_proj.GetXaxis().SetTitle("#alpha_{T}")
        denom_cumu_proj.GetYaxis().SetTitle("Events")
        denom_cumu_proj.GetYaxis().SetRangeUser(numer_cumu_proj.GetMinimum()*0.9,
                                                denom_cumu_proj.GetMaximum()*1.1)

        numer_cumu_proj.Draw("histsame")
        numer_cumu_proj.SetLineWidth(2)
        numer_cumu_proj.SetLineColor(r.kGreen)

    pad1.Update()
    pad1.Modified()
    c1.cd()

    # Compute the pad range with suitable margins
    ymin = 0.
    ymax = 1.1
    dy = (ymax-ymin)/0.8 # 10 per cent margins top and bottom
    xmin = denom_cumu_proj.GetXaxis().GetXmin()
    xmax = denom_cumu_proj.GetXaxis().GetXmax()
    dx = (xmax-xmin)/0.8 # 10 per cent margins left and right
    pad2.Range( xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy )
    pad2.Draw()
    pad2.cd()

    graph_diff = r.TGraphAsymmErrors( numer_diff_proj, denom_diff_proj, "cp" )
    graph_cumu = r.TGraphAsymmErrors( numer_cumu_proj, denom_cumu_proj, "cp" )

    vals_diff = sorted([ (at,eff,eh,el) for at,eff,eh,el in zip( [ graph_diff.GetX()[x] for x in range(graph_diff.GetN()) ],
                                                                 [ graph_diff.GetY()[x] for x in range(graph_diff.GetN()) ],
                                                                 [ graph_diff.GetEYhigh()[x] for x in range(graph_diff.GetN()) ],
                                                                 [ graph_diff.GetEYlow()[x] for x in range(graph_diff.GetN()) ] ) ], key = lambda val : val[0] )
    index_diff = bisect.bisect_right( [ x[0] for x in vals_diff ], yaxis_thresholds[iht_bin]+offset )
    eff_diff = (vals_diff[index_diff][1:]) if index_diff < len(vals_diff) else (np.nan,np.nan,np.nan)

    vals_cumu = sorted([ (at,eff,eh,el) for at,eff,eh,el in zip( [ graph_cumu.GetX()[x] for x in range(graph_cumu.GetN()) ],
                                                                 [ graph_cumu.GetY()[x] for x in range(graph_cumu.GetN()) ],
                                                                 [ graph_cumu.GetEYhigh()[x] for x in range(graph_cumu.GetN()) ],
                                                                 [ graph_cumu.GetEYlow()[x] for x in range(graph_cumu.GetN()) ] ) ], key = lambda val : val[0] )
    index_cumu = bisect.bisect_right( [ x[0] for x in vals_cumu ], yaxis_thresholds[iht_bin]+offset )
    eff_cumu = (vals_cumu[index_cumu][1:]) if index_cumu < len(vals_cumu) else (np.nan,np.nan,np.nan)
    
    graph_diff.Draw("psame")
    graph_diff.SetLineWidth(2)
    graph_diff.SetLineColor(r.kRed)
    graph_diff.SetMarkerColor(r.kRed)
    graph_cumu.Draw("psame")
    graph_cumu.SetLineWidth(2)
    graph_cumu.SetLineColor(r.kBlue)
    graph_cumu.SetMarkerColor(r.kBlue)

    pad2.Update()
    pad2.Modified()

    # draw axis on the right side of the pad
    axis = r.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
    #axis.SetLabelColor(r.kRed)
    #axis.SetTitleColor(r.kRed)
    axis.Draw("same")
    axis.SetTitle("Efficiency")
    c1.Modified()
    c1.Update()
    c1.SaveAs("plots/eff"+title+".pdf")
    r.gROOT.SetBatch(r.kFALSE)
    r.SetOwnership( c1, False ) 
    #raw_input("")
    #quit()

    return eff_diff,eff_cumu

def run() :
    
    file = r.TFile(input_file) 

    count = file.Get("count_total/count")
    print "Number of events processed: {:.0f}".format( count.GetEntries() )

    # Extract numerator and denominator histograms
    numer = OrderedDict()
    denom = OrderedDict()
    for ratio in ["numer","denom"] :
        for sample in samples.items() :
            for bjet_bin in bjet_bins :
                for iht_bin,ht_bin in enumerate(ht_bins) :
                    ht_upper = ht_bins[iht_bin+1] if iht_bin+1 < len(ht_bins) else None

                    for njet_bin in njet_bins.items() :

                        entries = 0
                        key = (sample[0],njet_bin[0],bjet_bin,ht_bin)
                        for njet in njet_bin[1] :
                            
                            ht_upper_str = "_{:.0f}".format( ht_upper ) if ht_upper is not None else ""
                            name = "{:s}_{:s}_{:.0f}{:s}/{:s}_{:s}".format( sample[1], 
                                                                            bjet_bin, 
                                                                            ht_bin,
                                                                            ht_upper_str,
                                                                            ratio,
                                                                            njet )
                            histo = file.Get(name)
                            entries += histo.GetEntries()
                            if ratio == "numer" :
                                if key in numer.keys() : numer[key].Add(histo)
                                else : numer[key] = r.TH2D(histo)
                                # Inclusive histogram
                                incl = (sample[0],"ge2j","ge0b","0")
                                if incl in numer.keys() : numer[incl].Add(histo)
                                else : numer[incl] = r.TH2D(histo)
                            if ratio == "denom" :
                                if key in denom.keys() : denom[key].Add(histo)
                                else : denom[key] = r.TH2D(histo)
                                # Inclusive histogram
                                incl = (sample[0],"ge2j","ge0b","0")
                                if incl in denom.keys() : denom[incl].Add(histo)
                                else : denom[incl] = r.TH2D(histo)

                        print "histo={:s} sample={:s} bjet={:s} ht={:.0f} njet={:s}, entries={:.0f}".format( ratio,
                                                                                                             sample[0], 
                                                                                                             bjet_bin, 
                                                                                                             ht_bin, 
                                                                                                             njet_bin[0],
                                                                                                             entries )

    # Rebin histograms, produce cumulative histograms, calculate efficiencies
    effs_diff = OrderedDict()
    effs_cumu = OrderedDict()
    for sample in samples.items() :
        for bjet_bin in bjet_bins :
            for iht_bin,ht_bin in enumerate(ht_bins) :
                ht_upper = ht_bins[iht_bin+1] if iht_bin+1 < len(ht_bins) else None
                for njet_bin in njet_bins.items() :

                    key = (sample[0],njet_bin[0],bjet_bin,ht_bin)

                    ht_upper_str = "_{:.0f}".format( ht_upper ) if ht_upper is not None else ""
                    title = "_{:s}_{:s}_{:s}_{:.0f}{:s}".format( sample[0], 
                                                                 njet_bin[0],
                                                                 bjet_bin, 
                                                                 ht_bin,
                                                                 ht_upper_str,
                                                                 )
                    
                    if key in numer.keys() and key in denom.keys() :
                        numer_diff = rebin( numer[key], ht_binning, yaxis_binning, "numer_diff"+title )
                        denom_diff = rebin( denom[key], ht_binning, yaxis_binning, "denom_diff"+title )
                        numer_diff.Sumw2()
                        denom_diff.Sumw2()
                        numer_cumu = cumu( numer_diff, title="numer_cumu"+title, ht_low=ht_bin, ht_high=ht_upper )
                        denom_cumu = cumu( denom_diff, title="denom_cumu"+title, ht_low=ht_bin, ht_high=ht_upper )
                        numer_cumu.Sumw2()
                        denom_cumu.Sumw2()
                        ratio_diff = numer_diff.Clone()
                        ratio_diff.SetName("ratio_diff"+title)
                        ratio_diff.SetTitle("ratio_diff"+title)
                        ratio_diff.Divide(numer_diff,denom_diff)
                        ratio_cumu = numer_cumu.Clone()
                        ratio_cumu.SetName("ratio_cumu"+title)
                        ratio_cumu.SetTitle("ratio_cumu"+title)
                        ratio_cumu.Divide(numer_cumu,denom_cumu)

                        #print 
                        #print [ "{:.3f}".format(ratio_cumu.GetYaxis().GetBinLowEdge(x)) for x in range(ratio_cumu.GetYaxis().GetNbins()) ] 
                        #print [ "{:.3f}".format(ratio_cumu.GetYaxis().GetBinCenter(x+1)) for x in range(ratio_cumu.GetYaxis().GetNbins()) ] 
                        #print [ "{:.3f}".format(ratio_cumu.GetBinContent(iht_bin+1,x+1)) for x in range(ratio_cumu.GetYaxis().GetNbins()) ] 
                        
                        eff_diff,eff_cumu = make_plot( numer_diff, denom_diff,
                                                       numer_cumu, denom_cumu,
                                                       iht_bin, ht_bin, ht_upper, 
                                                       title ) 
                        
                        effs_diff[key] = eff_diff
                        effs_cumu[key] = eff_cumu

                        print "Efficiency: sample={:s} njet={:s} bjet={:s} ht={:.0f}{:s} eff={:5.3f} + {:5.3f} - {:5.3f}".format( sample[0], 
                                                                                                                                  njet_bin[0],
                                                                                                                                  bjet_bin,
                                                                                                                                  ht_bin,
                                                                                                                                  ht_upper_str,
                                                                                                                                  effs_cumu[key][0],
                                                                                                                                  effs_cumu[key][1],
                                                                                                                                  effs_cumu[key][2],
                                                                                                                                  )

    #quit()

    effs_cumu[("Yossof","le3j","ge0b",200)] = (0.818,0.004,0.004)
    effs_cumu[("Yossof","le3j","ge0b",275)] = (0.952,0.004,0.004)
    effs_cumu[("Yossof","le3j","ge0b",325)] = (0.979,0.003,0.003)
    effs_cumu[("Yossof","le3j","ge0b",375)] = (0.992,0.002,0.002)
    effs_cumu[("Yossof","le3j","ge0b",475)] = (0.998,0.003,0.003)

    effs_cumu[("Yossof","ge4j","ge0b",200)] = (0.789,0.004,0.004)
    effs_cumu[("Yossof","ge4j","ge0b",275)] = (0.900,0.013,0.013)
    effs_cumu[("Yossof","ge4j","ge0b",325)] = (0.956,0.010,0.010)
    effs_cumu[("Yossof","ge4j","ge0b",375)] = (0.987,0.007,0.007)
    effs_cumu[("Yossof","ge4j","ge0b",475)] = (0.996,0.007,0.007)

    for ht_bin in [575,675,775,875,975,1075] :
        effs_cumu[("Yossof","le3j","ge0b",ht_bin)] = (1.000,0.000,0.000)
        effs_cumu[("Yossof","ge4j","ge0b",ht_bin)] = (1.000,0.000,0.000)

    print effs_cumu.keys()

    with open('tables/output.tex','w') as f:

        f.write( "\\documentclass[]{article}")
        f.write( "\\usepackage{lscape}")
        f.write( "\\begin{document}")
        f.write( "\\clearpage")
        f.write( "\\begin{landscape}")
        f.write( "\\begin{center}")
        f.write( "\\begin{table}[h!]")
        f.write( "\\caption{Signal trigger efficiencies.}" )
        f.write( "\\centering")
        f.write( "\\scriptsize")
        f.write( "\\begin{tabular}{lll"+"l"*len(ht_binning)+"}")
        f.write( "\\hline")
        f.write( "\\hline")
        f.write( "$N_{\\textrm{jet}}$ & $N_{\\textrm{b}}$ & Sample & \\multicolumn{"+str(len(ht_binning))+"}{c}{$H_{\\textrm{T}}$ (GeV)} \\\\")
        f.write( "\\cline{4-14}")
        f.write( " & & & "+" & ".join( [str(x) for x in ht_binning] )+"\\\\")

        f.write( "\\hline")
        f.write( "2--3 & $\geq 0$ & Yossof &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("Yossof","le3j","ge0b",x)][0],
                                                                                                     effs_cumu[("Yossof","le3j","ge0b",x)][1],
                                                                                                     effs_cumu[("Yossof","le3j","ge0b",x)][2] ) \
                                                           for x in ht_binning ])+"\\\\")
        f.write( "$\geq 4$ & $\geq 0$ & Yossof &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("Yossof","ge4j","ge0b",x)][0],
                                                                                                         effs_cumu[("Yossof","ge4j","ge0b",x)][1],
                                                                                                         effs_cumu[("Yossof","ge4j","ge0b",x)][2] ) \
                                                               for x in ht_binning ])+"\\\\")

        f.write( "\\hline")
        f.write( "2--3 & $\geq 0$ & $\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mu","le3j","ge0b",x)][0],
                                                                                                            effs_cumu[("mu","le3j","ge0b",x)][1],
                                                                                                            effs_cumu[("mu","le3j","ge0b",x)][2] ) \
                                                                  for x in ht_binning ])+"\\\\")
        f.write( "$\geq 4$ & $\geq 0$ & $\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mu","ge4j","ge0b",x)][0],
                                                                                                                effs_cumu[("mu","ge4j","ge0b",x)][1],
                                                                                                                effs_cumu[("mu","ge4j","ge0b",x)][2] ) \
                                                                      for x in ht_binning ])+"\\\\")

        if "mm" in samples.keys() :
            f.write( "2--3 & $\geq 0$ & $\\mu\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mm","le3j","ge0b",x)][0],
                                                                                                                    effs_cumu[("mm","le3j","ge0b",x)][1],
                                                                                                                    effs_cumu[("mm","le3j","ge0b",x)][2] ) \
                                                                          for x in ht_binning ])+"\\\\")
            f.write( "$\geq 4$ & $\geq 0$ & $\\mu\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mm","ge4j","ge0b",x)][0],
                                                                                                                        effs_cumu[("mm","ge4j","ge0b",x)][1],
                                                                                                                        effs_cumu[("mm","ge4j","ge0b",x)][2] ) \
                                                                              for x in ht_binning ])+"\\\\")

        f.write( "\\hline")
        for cat in cats.items() :
            header = cat[1][0] + " & " + cat[1][1] + " & " 
#            val,errh,errl = None,None,None
#            if diff :
#                tmp1 = effs_cumu[("mu","le3j","ge0b",x)][0]
#                tmp1h = effs_cumu[("mu","le3j","ge0b",x)][1]
#                tmp1l = effs_cumu[("mu","le3j","ge0b",x)][2]
#                tmp2 = effs_cumu[("mu",cat[0][0],cat[0][1],x)][0]
#                tmp2h = effs_cumu[("mu",cat[0][0],cat[0][1],x)][1]
#                tmp2l = effs_cumu[("mu",cat[0][0],cat[0][1],x)][2]
#                val = tmp2 - tmp1
#                errh = math.sqrt( abs( tmp2h*tmp2h - tmp1h*tmp1h ) )
#                errl = math.sqrt( abs( tmp2l*tmp2l - tmp1l*tmp1l ) )
#            else :
#                val = effs_cumu[("mu",cat[0][0],cat[0][1],x)][0]
#                errh = effs_cumu[("mu",cat[0][0],cat[0][1],x)][1]
#                errl = effs_cumu[("mu",cat[0][0],cat[0][1],x)][2]
#            f.write( header,"$\\mu$ + jets &"," & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( val, errh, errl ) for x in ht_binning ]),"\\\\"
            f.write( header+"$\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mu",cat[0][0],cat[0][1],x)][0],
                                                                                                     effs_cumu[("mu",cat[0][0],cat[0][1],x)][1],
                                                                                                     effs_cumu[("mu",cat[0][0],cat[0][1],x)][2] 
                                                                                                     ) for x in ht_binning ])+"\\\\")

        if "mm" in samples.keys() :
            for cat in cats.items() :
                header = cat[1][0] + " & " + cat[1][1] + " & " 
#                val,errh,errl = None,None,None
#                if diff :
#                    tmp1 = effs_cumu[("mm","le3j","ge0b",x)][0]
#                    tmp1h = effs_cumu[("mm","le3j","ge0b",x)][1]
#                    tmp1l = effs_cumu[("mm","le3j","ge0b",x)][2]
#                    tmp2 = effs_cumu[("mm",cat[0][0],cat[0][1],x)][0]
#                    err2h = effs_cumu[("mm",cat[0][0],cat[0][1],x)][1]
#                    err2l = effs_cumu[("mm",cat[0][0],cat[0][1],x)][2]
#                    val = tmp2 - tmp1
#                    errh = math.sqrt( abs( tmp2h*tmp2h - tmp1h*tmp1h ) )
#                    errl = math.sqrt( abs( tmp2l*tmp2l - tmp1l*tmp1l ) )
#                else :
#                    val = effs_cumu[("mm",cat[0][0],cat[0][1],x)][0]
#                    errh = effs_cumu[("mm",cat[0][0],cat[0][1],x)][1]
#                    errl = effs_cumu[("mm",cat[0][0],cat[0][1],x)][2]
#                f.write( header,"$\\mu\\mu$ + jets &"," & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( val, errh, errl ) for x in ht_binning ]),"\\\\"
                f.write( header+"$\\mu\\mu$ + jets &"+" & ".join(["${:4.2f}^{{+{:4.2f}}}_{{-{:4.2f}}}$".format( effs_cumu[("mm",cat[0][0],cat[0][1],x)][0],
                                                                                                             effs_cumu[("mm",cat[0][0],cat[0][1],x)][1],
                                                                                                             effs_cumu[("mm",cat[0][0],cat[0][1],x)][2] 
                                                                                                             ) for x in ht_binning ])+"\\\\")

        f.write( "\\hline")
        f.write( "\\hline")
        f.write( "\\end{tabular}")
        f.write( "\\end{table}")
        f.write( "\\end{center}")
        f.write( "\\end{landscape}")
        f.write( "\\end{document}")

#    quit()

    print 
    key = ("mu","ge2j","ge0b","0")
    print incl

    numer_orig = numer[key]
    denom_orig = denom[key]
    
    print "numer ",numer_orig.Print()
    print "denom ",denom_orig.Print()

    numer_diff = rebin( numer_orig, ht_binning, yaxis_binning, "numer_mu_ge2j_ge0b_HT0_rebinned" )
    denom_diff = rebin( denom_orig, ht_binning, yaxis_binning, "denom_mu_ge2j_ge0b_HT0_rebinned" )
    
    print "numer_diff ",numer_diff.Print()
    print "denom_diff ",denom_diff.Print()

    numer_cumu = cumu( numer_diff, title="numer_mu_ge2j_ge0b_HT0_cumulative" )
    denom_cumu = cumu( denom_diff, title="denom_mu_ge2j_ge0b_HT0_cumulative" )
    
    print "numer_cumu ",numer_cumu.Print()
    print "denom_cumu ",denom_cumu.Print()

    eff_diff = numer_diff.Clone()
    eff_diff.SetName("eff_mu_ge2j_ge0b_HT0_diff")
    eff_diff.SetTitle("eff_mu_ge2j_ge0b_HT0_diff")
    eff_diff.Divide(numer_diff,denom_diff)

    eff_cumu = numer_cumu.Clone()
    eff_cumu.SetName("eff_mu_ge2j_ge0b_HT0_cumu")
    eff_cumu.SetTitle("eff_mu_ge2j_ge0b_HT0_cumu")
    eff_cumu.Divide(numer_cumu,denom_cumu)

    print "eff_diff ",eff_diff.Print()
    print "eff_cumu ",eff_cumu.Print()

#    numer.Draw()
#    raw_input("")
#    denom.Draw()
#    raw_input("")

#    numer_diff.Draw()
#    raw_input("")
#    denom_diff.Draw()
#    raw_input("")

#    numer_cumu.Draw()
#    raw_input("")
#    denom_cumu.Draw()
#    raw_input("")

#    eff_diff.Draw()
#    raw_input("")
#    eff_cumu.Draw()
#    raw_input("")

################################################################################
# EXECUTE ######################################################################
################################################################################

if __name__ == '__main__': 
    run()
