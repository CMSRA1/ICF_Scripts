settings = {}
settings["bins"] = ["200","275","325","375","475","575","675","775","875","975","1075"]

Btag_Efficiencies = {'Had_Z0':{},'Muon_Z0':{},'DiMuon':{},'Had_Z2':{},'Muon_Z2':{},'Photon':{},'DiLepton_Z0':{},'DiLepton_Z2':{}}
bins = tuple(settings["bins"])

dict_entries = ('Btag_Efficiency','Mistag_Efficiency','Ctag_Efficiency','Btag_Error','Mistag_Error','Ctag_Error')

for key in Btag_Efficiencies:
  # print key
  # print Btag_Efficiencies[key]
  Btag_Efficiencies[key] = dict.fromkeys(bins)
  print Btag_Efficiencies[key]
  for a in bins: 
    Btag_Efficiencies[key][a] = dict.fromkeys(dict_entries,0)
  print Btag_Efficiencies["Had_Z2"]