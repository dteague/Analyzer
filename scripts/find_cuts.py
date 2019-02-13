#!/usr/bin/env python
import argparse
import os
import json
from pprint import pprint



parser = argparse.ArgumentParser(description='Scrap configs to get cuts.')
parser.add_argument('--dir', dest='filespace', 
                   default="PartDet",  help='File space of cut files')

args = parser.parse_args()

configdir = args.filespace
if not os.path.isdir(configdir):
    print "Directory doesn't exist!"
    exit(1)

cutloc = "%s/Cuts.in" % (configdir)
if not os.path.exists(cutloc):
    print "No config files in directory!"
    exit(1)

filemap = {"NRecoMuon1":["Muon_info.json","Muon1"], "NRecoMuon2":["Muon_info.json","Muon2"],
           "NRecoElectron1":["Electron_info.json","Elec1"], "NRecoElectron2":["Electron_info.json","Elec2"],
           "NRecoJet1":["Jet_info.json","Jet1"], "NRecoJet2":["Jet_info.json","Jet2"],
           "NRecoBJet":["Jet_info.json","BJet"], "NRecoWJet":["FatJet_info.json","WJet"],
           "METCut": ["Run_info.json", "Run"], "NRecoTriggers1": ["Run_info.json", "Run"], 
           "NLeptonPair":["LeptonPair_info.json", "GenericPair"]
           }


    
cutfile = open(cutloc)
for line in cutfile:
    line = line.strip()
    slash_loc = line.find('//')
    if slash_loc != -1:
        line = line[0:slash_loc].strip()
    cuts = line.split()
    if len(cuts) != 3 or (cuts[1] == '0' and cuts[2] == '-1'): continue

    print "%s from %s to %s" %(cuts[0], cuts[1], cuts[2])
    print "-"*20
    if cuts[0] not in filemap:
        print
        continue
    partfile = "%s%s" % (configdir, filemap[cuts[0]][0])
    
    with open(partfile, 'r') as f2:
        clean_file=''
        for f2line in f2:
            f2line = f2line.strip()
            slash_loc = f2line.find('//')
            if slash_loc != -1:
                f2line = f2line[0:slash_loc].strip()
            clean_file += f2line
        data = json.loads(clean_file)
        f2.close()
        subdata = data[filemap[cuts[0]][1]]
    for item in subdata:
        if type(subdata[item]) == type(True) and subdata[item]:
            subitem = item[7:] + "Cut"
            extra = "true"
            if subitem in subdata:
                extra = subdata[subitem]

            print "%s:  " % item,  extra
        if item == "EtaCut" or item == "PtCut":
            print "%s:  " % item,  subdata[item]
    print "-"*20
    print
