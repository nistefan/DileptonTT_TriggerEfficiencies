#!/bin/env python
__author__ = 'Nicole Stefanov'
__date__ = '2020-06-09'
         
#from array import array
#import collections
#import json
#from collections import defaultdict
#from collections import Counter
import sys, os
import math
#from math import sin, cos, pi, sqrt
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
#import numpy
import numpy as np
#import matplotlib.pyplot as plt          
#from numpy import matrix
#from numpy import linalg
#from numpy import array
#import matplotlib.cm as cm
#from numpy import linalg as LA
#from sympy import *
#LinAlgError = numpy.linalg.linalg.LinAlgError
from operator import itemgetter
#import heapq
import time
#import copy
import ROOT as r
from ROOT import TFile, TTree, TH1F, TObject, Math
from ROOT import gROOT, AddressOf
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed, kGreen
r.gROOT.SetBatch(1)
r.gROOT.LoadMacro('vecUtilsM4D2.h'+'++') #'__PWD__/vecUtilsM4D__NUMBER__.h'+'++') #./vecUtilsM4D8.h'+'++')                                  
lv = r.Math.LorentzVector(r.Math.PtEtaPhiM4D('float'))
#llv = r.Math.LorentzVector(r.Math.PtEtaPhiM4D('float'))                                                                        
tlv = r.TLorentzVector
DeltaPhi = r.Math.VectorUtil.DeltaPhi
InvariantMass=r.Math.VectorUtil.InvariantMass
DeltaR = r.Math.VectorUtil.DeltaR
Angle =r.Math.VectorUtil.Angle
#r.gROOT.LoadMacro('./dileptonSolver/DileptonAnalyticalSolver.cc+')
#solver = r.llsolver.DileptonAnalyticalSolver()
import subprocess

#gROOT.ProcessLine(
#"struct MyStruct {\
#   Float_t     fMyInt1;\
#   Float_t     fMyInt2;\
#   Float_t     fMyInt3;\
#   Float_t     fMyInt4;\
#   Float_t     fMyInt5;\
#   Float_t     fMyInt6;\
#   Float_t     fMyInt7;\
#   Float_t     fMyInt8;\
#};" );



verbose = True


def createEfficiencyHistogram(): #histoTag):
     nbins_subleadPt_2D = 3;
     edges_subleadPt_2D = [20.0, 50.0, 90.0, 200.0]
#     // === 09-13-18 (bin optimization) ===

     nbins_pt_2D = 4;
     edges_pt_2D = [20.0, 50.0, 80.0, 120.0, 200.0]

 #    // === 09-11-18 (matching AN2016_392) ===
     nbins_pt = 6;
     edges_pt = [20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0]
#     hnamePre, htitlePre = 'h_ratioMETPHIoverMETNeut_', ' h_ratioMETPHIoverMETNeut'
#     titleX, titleY = 'RatioMETMediator_over_METNeut', 'Events'
#     h_ratioMETPHIoverMETNeut_x = r.TH1F('h_ratioMETPHIoverMETNeut_x', 'RatioMETPHIoverMETNeut_x;'+titleX+';'+titleY, 100, -20, 20)
#     h_ratioMETPHIoverMETNeut_y = r.TH1F('h_ratioMETPHIoverMETNeut_y', 'RatioMETPHIoverMETNeut_y;'+titleX+';'+titleY, 100, -20, 20)
#     h_ratioMETPHIoverMETNeutcomplete = r.TH1F('h_ratioMETPHIoverMETNeutcomplete', 'RatioMETPHIoverMETNeutcomplete;'+titleX+';'+titleY, nbins_pt, edges_pt)
     hnamePre, htitlePre ='h_lepton_pt_', 'h_lepton_pt_'
     titleX, titleY = 'Pt', 'Events'
     h_lepton_pt_MET_F = r.TH1F('h_lepton_pt_MET_F', 'lepton_pt_MET_F;'+titleX+';'+titleY, nbins_pt, edges_pt)
     h_lepton_pt_MET_E = r.TH1F('h_lepton_pt_MET_E', 'lepton_pt_MET_E;'+titleX+';'+titleY, nbins_pt, edges_pt) 
     h_lepton_pt_MET_D = r.TH1F('h_lepton_pt_MET_D', 'lepton_pt_MET_D;'+titleX+';'+titleY, nbins_pt, edges_pt)
     h_lepton_pt_MET_C = r.TH1F('h_lepton_pt_MET_C', 'lepton_pt_MET_C;'+titleX+';'+titleY, nbins_pt, edges_pt)
     h_lepton_pt_MET_B = r.TH1F('h_lepton_pt_MET_B', 'lepton_pt_MET_B;'+titleX+';'+titleY, nbins_pt, edges_pt)
     h_lepton_pt_TT = r.TH1F('h_lepton_pt_TT', 'lepton_pt_TT;'+titleX+';'+titleY, nbins_pt, edges_pt)

def createStartingCodeForEfficiencyHistogram(h_sethistogramName, setBinningIntervals, setTitleX, setTitleY):
#     edges_subleadPt_2D = np.array([20.0, 50.0, 90.0, 200.0])
      #     // === 09-13-18 (bin optimization) ===
#     edges_pt_2D = np.array([20.0, 50.0, 80.0, 120.0, 200.0])
#     edges_pt = np.array([20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0])
     #hnamePre, htitlePre ='h_', 'h_'
# example
#     titleX, titleY = '1st lepton Pt', 'Events'
     titleX, titleY = setTitleX, setTitleY
#     example
#     h_leadingLepton_pt_MET_D = r.TH1F('h_lepton_pt_MET_D', 'lepton_pt_MET_D;'+titleX+';'+titleY, len(edges_pt)-1, edges_pt) #20., 200.) #nbins_pt, edges_pt)
     h_sethistogramName = r.TH1F(''+str(h_sethistogramName)+'', ''+str(h_sethistogramName)[2:]+';'+titleX+';'+titleY, len(setBinningIntervals)-1, setBinningIntervals)

#     titleX, titleY = '2nd lep Pt', 'Events'
#     h_subleadingLepton_pt_MET_D = r.TH1F('h_subleadingLepton_pt_MET_D', '2nd_lepton_pt_MET_D;'+titleX+';'+titleY, 6, 20., 200.)# nbins_pt, edges_pt)                     
     return titleX, titleY, h_sethistogramName



 
def writeEfficiencyHistogram():
        c = r.TCanvas('DM_ratio_to_met','solver comparison DM_ratio_to_met')
        c.Divide(2)
        c.cd(1)
        h_DMx_ratio_to_Metx.Draw()
        c.cd(2)
        h_DMy_ratio_to_Mety.Draw()
        for ext in ['png'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
        for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
        h_DMx_ratio_to_Metx.Write()
        h_DMy_ratio_to_Mety.Write()


def stringNameInVariableOut(stringName, variableValue):
     stringName = variableValue
     return stringName


## Vorsicht aender das fuer 2017/8; eta 2.5 anpassen!! Nur jetzt zum Testen.
def applyMuonCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_relIso):

    if (lepton_isMuon == 1) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (lepton_relIso < 0.25):
        return True
    else: return False

#def applyMuonCuts(muon_CutBasedID_Medium, lepton, 


def applyElectronCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_scEta):
    if (lepton_isMuon == 0) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (abs(lepton_scEta) < 1.444 or abs(lepton_scEta) > 1.566):
        return True
    else: return False


def applyJetCuts(jet_pt, jet_eta, jet_puid):
     if (jet_pt > 20.) and abs(jet_eta) < 2.4 and jet_puid == 7:
          return True
     else: return False


#def applyPileupReweighting():
 #    dataPileupFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/pileup/pileupDATA__Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1__max120_n120__MinBias69200.root"
     
## Getting SF from root file with ROOT. Just has been used for consistency check. Not anymore.
def getLeptonSF(SF_FileHisto, sfForWhichxAxisValue, sfForWhichyAxisValue):
     xBin = SF_FileHisto.GetXaxis().FindBin(sfForWhichxAxisValue)
     yBin = SF_FileHisto.GetYaxis().FindBin(sfForWhichyAxisValue)

     # Taking care of underflow (xBin = 0) bins, setting to first bin ...
     if xBin == 0:
          xBin = 1
     # ... and overflow (xBin > SF_FileHisto.GetNbinsX()) bins, setting to max. possible bin of that histogram
     elif xBin > SF_FileHisto.GetNbinsX():
          xBin = SF_FileHisto.GetNbinsX()

     # Taking care of underflow (yBin = 0) bins, setting to first bin ...
     if yBin == 0:
          yBin = 1
     # ... and overflow (yBin > SF_FileHisto.GetNbinsY()) bins, setting to max. possible bin of that histogram
     elif yBin > SF_FileHisto.GetNbinsY():
          yBin = SF_FileHisto.GetNbinsY()
     SFvalue = SF_FileHisto.GetBinContent(xBin, yBin)
#     print(xBin)
#     print(yBin)
     if not SFvalue:
          SFvalue = 1.
     return SFvalue

###### convert TH2D or TH2F histogram into python xD list, origninally made for SF handling without ROOT, but this works for general case, too - of course.
def create_sfArray(inputSFarrayFileRootPath, whichSFhistogram):
     inputSFarrayFileRootPath_OpenIt = r.TFile.Open(inputSFarrayFileRootPath, "read")
     inputSFarrayFileRootPath_GetHisto = inputSFarrayFileRootPath_OpenIt.Get(whichSFhistogram)
     gatheredSFweights = []
     for item_sfRangeX in range(1, inputSFarrayFileRootPath_GetHisto.GetNbinsX()+1):
          getOneYBinRange = []
          for item_sfRangeY in range(1, inputSFarrayFileRootPath_GetHisto.GetNbinsY()+1):
               getOneYBin=[]
               getOneYBin = [inputSFarrayFileRootPath_GetHisto.GetYaxis().GetBinLowEdge(item_sfRangeY), inputSFarrayFileRootPath_GetHisto.GetBinContent(item_sfRangeX, item_sfRangeY)]
               getOneYBinRange.append(getOneYBin)
               print(str(inputSFarrayFileRootPath_GetHisto.GetXaxis().GetBinUpEdge(item_sfRangeX))+"|"+str(inputSFarrayFileRootPath_GetHisto.GetYaxis().GetBinLowEdge(item_sfRangeY))+" =sf= "+str(inputSFarrayFileRootPath_GetHisto.GetBinContent(item_sfRangeX, item_sfRangeY)))
          gatheredSFweights.append([inputSFarrayFileRootPath_GetHisto.GetXaxis().GetBinLowEdge(item_sfRangeX), getOneYBinRange])
          
     print(gatheredSFweights)
     inputSFarrayFileRootPath_OpenIt.Close()
     return gatheredSFweights


def leptonSF(sfArray, sfForWhichxAxisValue, sfForWhichyAxisValue):
     sf = 1.
     count_item_X = 0
     count_item_Y = 0
     for item_X in range(len(sfArray)):
          if sfArray[item_X][0] < sfForWhichxAxisValue:
               count_item_X = item_X
          else:
               break
     for item_Y in range(len(sfArray[count_item_X][1])):
          if sfArray[count_item_X][1][item_Y][0] < sfForWhichyAxisValue:
               count_item_Y = item_Y
          else:
               break
     sf = sfArray[count_item_X][1][count_item_Y][1]
     return sf


def splitListintoSubListsOfPieceSize(ListToSplit, PieceSize):
    for item_splitListintoSubListsOfPieceSize in range(0, len(ListToSplit), PieceSize):
        yield ListToSplit[item_splitListintoSubListsOfPieceSize:item_splitListintoSubListsOfPieceSize + PieceSize]


def errorSqrtSumOfsquaredWeights(someArrayWithNumbers):
     error = 0.
     for item_someArrayWithNumbers in someArrayWithNumbers:
          error += (item_someArrayWithNumbers**2)
     return np.sqrt(error)
     

RootfilesToRunOver =["/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/met_run2017D.root",
#                     "/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/ttbarsignalplustau_fromDilepton.root"
#/nfs/dust/cms/user/stefanon/analyses2017/ckoraka_TriggerSF/Converted_Ntuples_For_TiggerSF_Script/met_run2017D.root",
#"/nfs/dust/cms/group/topcmsdesy/ntuple13tev/2017_09_05_TAG_V037/ttbarsignalplustau.root",
]

whichSubDirectoriesToCreate = ["DoubleEl_OR__X__allMET",
"DoubleEl_OR__X__allMET_METHIGH",
"DoubleEl_OR__X__allMET_METLOW",
"DoubleEl_OR__X__allMET_NJETSHIGH",
"DoubleEl_OR__X__allMET_NJETSLOW",
"DoubleEl_OR__X__allMET_NPVHIGH",
"DoubleEl_OR__X__allMET_NPVLOW",
"DoubleMu_OR__X__allMET",
"DoubleMu_OR__X__allMET_METHIGH",
"DoubleMu_OR__X__allMET_METLOW",
"DoubleMu_OR__X__allMET_NJETSHIGH",
"DoubleMu_OR__X__allMET_NJETSLOW",
"DoubleMu_OR__X__allMET_NPVHIGH",
"DoubleMu_OR__X__allMET_NPVLOW",
"EMu_OR__X__allMET",
"EMu_OR__X__allMET_METHIGH",
"EMu_OR__X__allMET_METLOW",
"EMu_OR__X__allMET_NJETSHIGH",
"EMu_OR__X__allMET_NJETSLOW",
"EMu_OR__X__allMET_NPVHIGH",
"EMu_OR__X__allMET_NPVLOW",
"HLT_DoubleEl_OR",
"HLT_DoubleMu_OR",
"HLT_EMu_OR",
"DoubleEl__X__allMET", 
"DoubleMu__X__allMET",
"EMu__X__allMET",
"HLT_DoubleEl",
"HLT_DoubleMu",
"HLT_EMu",
"corr2D"]


Trigger_2017_se = [
                "HLT_Ele27_WPTight_Gsf_v_",
                "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_",
                "HLT_Ele35_WPTight_Gsf_v_",
                "HLT_Ele38_WPTight_Gsf_v_",
                "HLT_Ele40_WPTight_Gsf_v_",
                "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_",
                "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_",
]


Trigger_2017_smu = [
                # smu Triggers
                "HLT_IsoMu27_v_",
                "HLT_IsoMu24_eta2p1_v_",
]



Trigger_2017_mumu = [
# mumu Triggers
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_",
]


Trigger_2017_ee = [
# ee Triggers
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_",
                #    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_",
                "HLT_DoubleEle33_CaloIdL_MW_v_",
                "HLT_DoubleEle25_CaloIdL_MW_v_",
                "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_",
]


Trigger_2017_emu = [
# emu Triggers   
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_",
]

Trigger_2017_met = [
# met Triggers
                "HLT_PFMET120_PFMHT120_IDTight_v_",
                "HLT_PFMET250_HBHECleaned_v_",
                "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_",
                "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_",
                "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_"
]


Trigger_2017 = []
Trigger_2017 += Trigger_2017_se
Trigger_2017 += Trigger_2017_smu
Trigger_2017 += Trigger_2017_mumu
Trigger_2017 += Trigger_2017_ee
Trigger_2017 += Trigger_2017_emu
Trigger_2017 += Trigger_2017_met

define_counts = []
dictTrigger_2017 = {}
for item_Trigger_2017 in Trigger_2017:
     globals()[item_Trigger_2017] = False
#     exec("%s = %r" % (item_Trigger_2017, False))
     define_counts.append("count_"+str(item_Trigger_2017))
#     dictTrigger_2017[item_Trigger_2017] = False
     
for item_define_counts in define_counts:
     globals()[item_define_counts] = 0
#     exec("%s = %d" % (item_define_counts, 0))

nan = "_TestSeptemberRecoMllbB1"

### ==========================================
###          For pileup reweighting
### ==========================================
pileupReweighting = False
if pileupReweighting:
     dataPileupFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/pileup/pileupDATA__Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1__max120_n120__MinBias69200.root"
     inputdataPileupFile =  r.TFile.Open(dataPileupFile, "read")
     inputdataPileupHisto = inputdataPileupFile.Get("pileup")
     integralPileupHisto = inputdataPileupHisto.Integral()
#     print("integralPileupHisto")
 #    print(integralPileupHisto)
     #for item_binrange in range(0, 20):
      #    print("inputdataPileupHisto.GetXaxis().FindBin("+str(item_binrange))
       #   print(inputdataPileupHisto.GetXaxis().FindBin(item_binrange))
        #  print(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(item_binrange)))
  #   binContentArrayOfdataPileupHisto = []
     weightedBinContentArrayOfdataPileupHisto = []
    
     for itemBinContent in range(0, inputdataPileupHisto.GetNbinsX()+1):
          weightedBinContentArrayOfdataPileupHisto.append(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(itemBinContent))/integralPileupHisto)
    #      binContentArrayOfdataPileupHisto.append(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(itemBinContent)))
#     print("weightedBinContentArrayOfdataPileupHisto")
 #    print(weightedBinContentArrayOfdataPileupHisto)
  #   print(len(weightedBinContentArrayOfdataPileupHisto))
   #  print("binContentArrayOfdataPileupHisto")
    # print(binContentArrayOfdataPileupHisto)
     inputdataPileupFile.Close()
     MCPileupFile = "/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/ttbarsignalplustau_fromDilepton.root"
     inputMCPileupFile = r.TFile.Open(MCPileupFile, "read")
     inputMCPilepuFileHisto = inputMCPileupFile.Get("PileupMCTemplateMaker/MC_TrueNIntBX0")
     integralMCPilepuFileHisto = inputMCPilepuFileHisto.Integral()
#     print("integralMCPilepuFileHisto")
#     print(integralMCPilepuFileHisto)
#     binContentArrayOfMCPileupHisto = []
     weightedBinContentArrayOfMCPileupHisto = []
     #for item_binrange in range(0, 20):
      #    print("inputMCPilepuFileHisto.GetXaxis().FindBin("+str(item_binrange))
       #   print(inputMCPilepuFileHisto.GetXaxis().FindBin(item_binrange))
        #  print(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(item_binrange)))
     for itemBinContent in range(0, inputMCPilepuFileHisto.GetNbinsX()+1):
          weightedBinContentArrayOfMCPileupHisto.append(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(itemBinContent))/integralMCPilepuFileHisto)
 #         binContentArrayOfMCPileupHisto.append(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(itemBinContent)))
#     print("weightedBinContentArrayOfMCPileupHisto")
#     print(weightedBinContentArrayOfMCPileupHisto)
#     print(len(weightedBinContentArrayOfMCPileupHisto))
#     print("binContentArrayOfMCPileupHisto")
#     print(binContentArrayOfMCPileupHisto)
     inputMCPileupFile.Close()

     if len(weightedBinContentArrayOfMCPileupHisto) != len(weightedBinContentArrayOfdataPileupHisto): 
          sys.exit("Something went wrong for the pileup reweighting step. \n len(binContentArrayOfMCPileupHisto) != len(binContentArrayOfdataPileupHisto)")

     pileupWeights = []
     for item_run in range(len(weightedBinContentArrayOfMCPileupHisto)):
          if weightedBinContentArrayOfMCPileupHisto[item_run] != 0.:
               pileupWeights.append(weightedBinContentArrayOfdataPileupHisto[item_run]/weightedBinContentArrayOfMCPileupHisto[item_run])
          else: pileupWeights.append(1.)

     # Nope, problem with divding by 0.:
     #pileupWeights = [item_weightedBinContentArrayOfdataPileupHisto/item_weightedBinContentArrayOfMCPileupHisto for item_weightedBinContentArrayOfdataPileupHisto, item_weightedBinContentArrayOfMCPileupHisto in zip(weightedBinContentArrayOfdataPileupHisto, weightedBinContentArrayOfMCPileupHisto)]

     print("pileupWeights")
     print(pileupWeights)
     #print(pileupWeights[15])
     #print(pileupWeights[18])


ElectronIDsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/2017_ElectronMVA90.root"
#inputElectronIDsfFile = r.TFile.Open(ElectronIDsfFile, "read")
#inputElectronIDsfFileHisto = inputElectronIDsfFile.Get("EGamma_SF2D")
inputElectronIDhistogram = "EGamma_SF2D"
#gatheredElectronIDsfWeights = []
#getOneXBinRange = []

#for item_sfRangeX in range(1, inputElectronIDsfFileHisto.GetNbinsX()+1):
 #    getOneYBinRange = []
  #   print("Range of Xbin ("+str(inputElectronIDsfFileHisto.GetXaxis().GetBinLowEdge(item_sfRangeX))+" to "+str(inputElectronIDsfFileHisto.GetXaxis().GetBinUpEdge(item_sfRangeX)))
   #  for item_sfRangeY in range(1, inputElectronIDsfFileHisto.GetNbinsY()+1):
    #      getOneYBin=[]
     #     print("Range of Ybin ("+str(inputElectronIDsfFileHisto.GetYaxis().GetBinLowEdge(item_sfRangeY))+" to "+str(inputElectronIDsfFileHisto.GetYaxis().GetBinUpEdge(item_sfRangeY)))
      #    print("For bins ("+str(item_sfRangeX)+"|"+str(item_sfRangeY)+") gives SF = "+str(inputElectronIDsfFileHisto.GetBinContent(item_sfRangeX, item_sfRangeY)))
       #   getOneYBin = [inputElectronIDsfFileHisto.GetYaxis().GetBinLowEdge(item_sfRangeY), inputElectronIDsfFileHisto.GetBinContent(item_sfRangeX, item_sfRangeY)]
        #  getOneYBinRange.append(getOneYBin)
          
    # gatheredElectronIDsfWeights.append([inputElectronIDsfFileHisto.GetXaxis().GetBinLowEdge(item_sfRangeX), getOneYBinRange])
     
#inputElectronIDsfFile.Close()
#print(gatheredElectronIDsfWeights)

#sys.exit("D")
gatheredElectronIDsfWeights = create_sfArray(ElectronIDsfFile, inputElectronIDhistogram)

###### Example Test
#lep_Eta = -1.957742
#lep_Pt = 44.2188
#exampleSf = leptonSF(gatheredElectronIDsfWeights, lep_Eta, lep_Pt)
#print("exampleSf")
#print(exampleSf)

ElectronRECOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"
gatheredElectronRECOsfWeights = create_sfArray(ElectronRECOsfFile, inputElectronIDhistogram)


###### Example Test1
#exampleSf1 = leptonSF(gatheredElectronRECOsfWeights, lep_Eta, lep_Pt)
#print("exampleSf1")
#print(exampleSf1)

MuonIDsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ID_syst.root"
inputMuonIDhistogram = "NUM_MediumID_DEN_genTracks_pt_abseta"
gatheredMuonIDsfWeights = create_sfArray(MuonIDsfFile, inputMuonIDhistogram)



###### Example Test2
#exampleMuonSf = leptonSF(gatheredMuonIDsfWeights, lep_Pt, abs(lep_Eta))
#print("exampleMuonSf")
#print(exampleMuonSf)

MuonISOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ISO_syst.root"
inputMuonIsohistogram = "NUM_LooseRelIso_DEN_MediumID_pt_abseta"
gatheredMuonIsosfWeights = create_sfArray(MuonISOsfFile, inputMuonIsohistogram)

###### Example Test3
#exampleMuonSf1 = leptonSF(gatheredMuonIsosfWeights, lep_Pt, abs(lep_Eta))
#print("exampleMuonSf1")
#print(exampleMuonSf1)







#test = (44.2188,-0.957742,2.36075,-0.010981)
#test1 = (3.,-2.97742,2.36075,-0.010981) 

#def leptonSF(gatheredElectronIDsfWeights, sfForWhichxAxisValue, sfForWhichyAxisValue):
#sf = 1.
#count_item_X = 0
#count_item_Y = 0
#for item_X in range(len(gatheredElectronIDsfWeights)):
       #   print("Next event")
       #   checkIfTrue = gatheredElectronIDsfWeights[item_X][0] < test1[1]
        #  print(checkIfTrue)
         # if not checkIfTrue: break
#          elif checkIfTrue: count_item_X +=1
#          if gatheredElectronIDsfWeights[item_X][0] < test1[1]:
          #     print(gatheredElectronIDsfWeights[item_X][0])
           #    print(test[1])
#               count_item_X = item_X
            #   print(gatheredElectronIDsfWeights[count_item_X])
             #  print("count_item_X")
              # print(count_item_X)
               #print("item_X")
               #print(item_X)
#          else:       
#               break
#print("item_X")
#print(item_X)
#print("count_item_X = "+str(count_item_X))




#for item_Y in range(len(gatheredElectronIDsfWeights[count_item_X][1])):
 #                   print("---------------------Next event------------------")
  #                  print(item_Y)
                    #print(test[0])
                    #print(gatheredElectronIDsfWeights[count_item_X][1][count_item_Y][0])
#                    checkIfTrue1 = gatheredElectronIDsfWeights[count_item_X][1][item_Y][0] < test1[0]
 #                   print(checkIfTrue1)
  #                  if not checkIfTrue1: break
   #                 whichDouble = gatheredElectronIDsfWeights[count_item_X][1][item_Y][0]
    #                whichToCompare = test[0]
     #               print(whichDouble)
      #              print(whichToCompare)
       #             if whichDouble <= whichToCompare: print("Is true")
        #            else: break
 #                   if gatheredElectronIDsfWeights[count_item_X][1][item_Y][0] < test1[0]:
  #                       count_item_Y = item_Y
   #                      print(count_item_Y)
    #                     print(gatheredElectronIDsfWeights[count_item_X][1][count_item_Y])
     #                    print(gatheredElectronIDsfWeights[count_item_X][1][count_item_Y][0])
      #                   print(test[0])                         
#                    else: break
#print("count_item_Y = "+str(count_item_Y))
#sf = gatheredElectronIDsfWeights[count_item_X][1][count_item_Y][1]
#print("sf = "+str(sf))



#sys.exit("Here we stop test run")
### ========================================== 
###       For Lepton SF application
### ========================================== 
LeptonSFapplication = False
if LeptonSFapplication:
     ElectronIDsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/2017_ElectronMVA90.root"
     inputElectronIDsfFile = r.TFile.Open(ElectronIDsfFile, "read")
     inputElectronIDsfFileHisto = inputElectronIDsfFile.Get("EGamma_SF2D")
     
     ElectronRECOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"
     inputElectronRECOsfFile = r.TFile.Open(ElectronRECOsfFile, "read")
     inputElectronRECOsfFileHisto = inputElectronRECOsfFile.Get("EGamma_SF2D")
     
     MuonIDsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ID_syst.root"
     inputMuonIDsfFile = r.TFile.Open(MuonIDsfFile, "read")
     inputMuonIDsfFileHisto = inputMuonIDsfFile.Get("NUM_MediumID_DEN_genTracks_pt_abseta")
     

     MuonISOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ISO_syst.root"
     inputMuonISOsfFile = r.TFile.Open(MuonISOsfFile, "read")
     inputMuonISOsfFileHisto = inputMuonISOsfFile.Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta")
     #test = (44.2188,-0.957742,2.36075,-0.010981)
     #test = (501.,-2.57742,2.36075,-0.010981)
#test = (44.2188,-0.957742,2.36075,-0.010981)

#print(test[0])
#index_eta = 0
#index_pteta = 0


     #getSFvalue = getLeptonSF(inputElectronIDsfFileHisto, test[0], test[1]) 
     #print(getSFvalue)

     #inputElectronIDsfFile.Close()

     #sys.exit("Stopping test run.")


edges_subleadPt_2D = np.array([20.0, 50.0, 90.0, 200.0])
#     // === 09-13-18 (bin optimization) ===                                                                                                                         

nbins_pt_2D = 4;
edges_pt_2D = np.array([20.0, 50.0, 80.0, 120.0, 200.0])


nbins_pt = 6;
edges_pt = np.array([20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0])


hnamePre, htitlePre ='h_', 'h_'
titleX, titleY = '1st lepton Pt', 'Events'
h_leadingLepton_pt_MET_D = r.TH1F('h_lepton_pt_MET_D', 'lepton_pt_MET_D;'+titleX+';'+titleY, len(edges_pt)-1, edges_pt) #20., 200.) #nbins_pt, edges_pt)



titleX, titleY = '2nd lep Pt', 'Events'
h_subleadingLepton_pt_MET_D = r.TH1F('h_subleadingLepton_pt_MET_D', '2nd_lepton_pt_MET_D;'+titleX+';'+titleY, 6, 20., 200.)# nbins_pt, edges_pt)




### Change parameters as needed
years_considered = ["2016", "2017", "2018"]

RunPart_considered = ["A", "B", "C", "D", "E", "F", "G", "H"]  

#strength = ["lead", "subleading"]

strength = ["1st", "2nd"]

lepton = ["Electron", "Muon"]

oneLeptonFeature = ["Pt", "Eta", "relIso"]

elseFeature = ["dileptonMass", "Met", "PrimaryVertex_count", "BTag_count", "Jet_count"]

whichDileptonTriggerGroup = ["HLT_DoubleMuon", "HLT_DoubleElectron", "HLT_ElectronMuon", 
                             "HLT_DoubleMuon_WITH_SLtrigger", "HLT_DoubleElectron_WITH_SLtrigger", "HLT_ElectronMuon_WITH_SLtrigger", "HLT_allMET"]

whichTriggerGroup = ["HLT_DoubleMuon__X__HLT_allMET", "HLT_DoubleElectron__X__HLT_allMET", "HLT_ElectronMuon__X__HLT_allMET",
                     "HLT_DoubleMuon_WITH_SLtrigger__X__HLT_allMET", "HLT_DoubleElectron_WITH_SLtrigger__X__HLT_allMET", "HLT_ElectronMuon_WITH_SLtrigger__X__HLT_allMET"]

whichTriggerGroup += whichDileptonTriggerGroup


#which2DHistogramsCreatingForFiring = []


#for item_whichTriggerGroup in whichTriggerGroup:
 #    for item_Trigger2017 in Trigger2017:
  #        which2DHistogramsCreatingForFiring.append(item_whichTriggerGroup+"_vs_"+item_Trigger2017)

whichTriggerPassed = ["HLT_allMET", "HLT_allMET_SingleElectronCase", "HLT_allMET_SingleMuonCase", "HLT_allMET_DoubleMuonCase", "HLT_allMET_DoubleElectronCase",
                      "HLT_allMET_ElectronMuonCase", "HLT_DoubleMuon", "HLT_DoubleElectron", "HLT_ElectronMuon",
                      "HLT_DoubleMuon_AND_HLT_allMET", "HLT_DoubleElectron_AND_HLT_allMET", "HLT_ElectronMuon_AND_HLT_allMET"]

whichLeptonTriggerPassed = ["seTrigger", "smuTrigger", "eeTrigger", "emuTrigger", "mumuTrigger"]

which1DHistogramsToCreate = []

which1DHistogramsToCreateSingleElectronTrigger = []
which1DHistogramsToCreateSingleMuonTrigger = []
which1DHistogramsToCreateDoubleElectronTrigger = []
which1DHistogramsToCreateDoubleMuonTrigger = []
which1DHistogramsToCreateElectronMuonTrigger = []

which2DFeatures = ["1st_Electron_Pt_vs_Eta", "2nd_Electron_Pt_vs_Eta", "1st_Electron_Eta_vs_2nd_Electron_Eta", "1st_Electron_Pt_vs_2nd_Electron_Pt",
                             "1st_Muon_Pt_vs_Eta", "2nd_Muon_Pt_vs_Eta", "1st_Muon_Eta_vs_2nd_Muon_Eta", "1st_Muon_Pt_vs_2nd_Muon_Pt", 
                             "1st_Electron_Pt_vs_1st_Muon_Pt", "1st_Electron_Eta_vs_1st_Muon_Eta"
                  ]

which2DHistogramsToCreate = []
which2DHistogramsToCreateSingleElectronTrigger = []
which2DHistogramsToCreateSingleMuonTrigger = []
which2DHistogramsToCreateDoubleElectronTrigger = []
which2DHistogramsToCreateDoubleMuonTrigger = []
which2DHistogramsToCreateElectronMuonTrigger = []


strength_lepton_oneLeptonFeature = []

addingFeature = False

['1st_Electron_Pt', '1st_Muon_Pt', '2nd_Electron_Pt', '2nd_Muon_Pt', '1st_Electron_Eta', '1st_Muon_Eta', '2nd_Electron_Eta', '2nd_Muon_Eta',]
for itemwhichTriggerPassed in whichTriggerPassed:
     for item_elseFeature in elseFeature:
          which1DHistogramsToCreate.append(itemwhichTriggerPassed+"_"+item_elseFeature)
          print(itemwhichTriggerPassed+"_"+item_elseFeature)
     for item_oneLeptonFeature in oneLeptonFeature:
          for item_strength in strength:
               for item_lepton in lepton:
#                    if "relIso" not in itemoneLeptonFeature:
 #                        print(itemstrength+"_"+item_lepton+"_"+item_oneLeptonFeature)
  #                       strength_lepton_oneLeptonFeature.append(itemstrength+"_"+item_lepton+"_"+item_oneLeptonFeature)
#                    print(itemwhichTriggerPassed+"_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
                    which1DHistogramsToCreate.append(itemwhichTriggerPassed+"_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
#                    which1DHistogramsToCreateSingleElectronTrigger.append(itemwhichTriggerPassed+"_seTrigger_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
 #                   which1DHistogramsToCreateSingleMuonTrigger.append(itemwhichTriggerPassed+"_smuTrigger_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
  #                  which1DHistogramsToCreateDoubleElectronTrigger.append(itemwhichTriggerPassed+"_eeTrigger_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
   #                 which1DHistogramsToCreateDoubleMuonTrigger.append(itemwhichTriggerPassed+"_mumuTrigger_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
    #                which1DHistogramsToCreateElectronMuonTrigger.append(itemwhichTriggerPassed+"_emuTrigger_"+item_strength+"_"+item_lepton+"_"+item_oneLeptonFeature)
#                    print(itemstrength+"_"+item_lepton+"_"+item_oneLeptonFeature)

for itemwhichTriggerPassed in whichTriggerPassed:
     for item_which2DFeatures in which2DFeatures:
          which2DHistogramsToCreate.append(itemwhichTriggerPassed+"_"+item_which2DFeatures)
     #     which2DHistogramsToCreateSingleElectronTrigger.append(itemwhichTriggerPassed+"_seTrigger_"+item_which2DFeatures)
      #    which2DHistogramsToCreateSingleMuonTrigger.append(itemwhichTriggerPassed+"_smuTrigger_"+item_which2DFeatures)
       #   which2DHistogramsToCreateDoubleElectronTrigger.append(itemwhichTriggerPassed+"_eeTrigger_"+item_which2DFeatures)
        #  which2DHistogramsToCreateDoubleMuonTrigger.append(itemwhichTriggerPassed+"_mumuTrigger_"+item_which2DFeatures)
         # which2DHistogramsToCreateElectronMuonTrigger.append(itemwhichTriggerPassed+"_emuTrigger_"+item_which2DFeatures)



 #    for itemstrength in strength:
  #       for item_oneLeptonFeature in oneLeptonFeature:
   #           if "Pt" in item_oneLeptonFeature or "Eta" in item_oneLeptonFeature:
    #               for item_lepton in lepton:
     #                   for item_oneLeptonFeature2 in oneLeptonFeature:
       #                      if "Pt" in item_oneLeptonFeature2 or "Eta" in item_oneLeptonFeature2:
      #                            print(itemstrength+"_"+item_lepton+"_"+item_oneLeptonFeature+"_vs_"+item_oneLeptonFeature2)
                   
#     for item1stloop in strength_lepton_oneLeptonFeature:
 #         for item_2ndloop in strength_lepton_oneLeptonFeature:
  #             if item1stloop != item_2ndloop:
   #                 if not ("1st_Electron" in item1stloop and "2nd_Muon" in item_2ndloop):
    #                     print(item1stloop+"_vs_"+item_2ndloop)
     #                    which2DHistogramsToCreate.append(item1stloop+"_vs_"+item_2ndloop)

#channels = ["ee", "emu", "mumu"]
#whatType = ["data", "MC"]
#strength = ["1st", "2nd"]
#whatCombination = ["pt_eta_iso", "njet_nvertex_met"]

#whatArraysToCreate = []

#for item_channels in channels:
 #    for item_whatType in whatType:
  #        for item_whatCombination in whatCombination:
   #            whatArraysToCreate.append(item_channels+"_Channel_"+item_whatType+
#for item in Trigger_2017:
 #    globals()[item_Trigger_2017] = False



ee_Channel_data_1st_Lepton_pt_eta_iso = []
ee_Channel_data_2nd_Lepton_pt_eta_iso = []
ee_Channel_data_njet_nvertex_met = []

emu_Channel_data_1st_Lepton_pt_eta_iso = []
emu_Channel_data_2nd_Lepton_pt_eta_iso = []
emu_Channel_data_njet_nvertex_met = []

mumu_Channel_data_1st_Lepton_pt_eta_iso = []
mumu_Channel_data_2nd_Lepton_pt_eta_iso = []
mumu_Channel_data_njet_nvertex_met = []

ee_Channel_MC_1st_Lepton_pt_eta_iso = []
ee_Channel_MC_2nd_Lepton_pt_eta_iso = []
ee_Channel_MC_njet_nvertex_met = []

emu_Channel_MC_1st_Lepton_pt_eta_iso = []
emu_Channel_MC_2nd_Lepton_pt_eta_iso = []
emu_Channel_MC_njet_nvertex_met = []

mumu_Channel_MC_1st_Lepton_pt_eta_iso = []
mumu_Channel_MC_2nd_Lepton_pt_eta_iso = []
mumu_Channel_MC_njet_nvertex_met = []


for inputFileName in RootfilesToRunOver:            
            inputFile = r.TFile.Open(inputFileName, "read")
            #inputTree = inputFile.Get("ttHTreeMaker/worldTree")
            inputTree=inputFile.Get("writeNTuple/NTuple")
            events = inputTree.GetEntries()
            print(events)
            str2="."
            bis = inputFileName.find(str2)
            textname=inputFileName[:bis]
            print(textname)

            if "/" in inputFileName:
                textnamePre = inputFileName.split("/")
                textname = textnamePre[-1].replace(".root", "")
            print(textname)
            isMC = True
            if "run" in textname or "Run" in textname: 
                isMC = False
                for itemRunPart_considered in RunPart_considered:
    #                 print(itemRunPart_considered)
                     if itemRunPart_considered in textname: 
                          isRun = itemRunPart_considered
                          break
            year = ""
            for itemyears_considered in years_considered:
     #            print(itemyears_considered)
                 if itemyears_considered in textname:
                      year = itemyears_considered
                      break
            print("isRun is "+isRun)
            print("isMC is "+str(isMC))
            print("year is "+str(year))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~  Read in correct triggers which should be used for the considered data/MC   ~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Trigger_met = []
            Trigger_SingleElectron = []
            Trigger_SingleMuon = []
            Trigger_DoubleMuon = []
            Trigger_DoubleElectron = []
            Trigger_ElectronMuon = []
            Trigger_met = eval("Trigger_"+year+"_met")
            Trigger_SingleElectron = eval("Trigger_"+year+"_se")
            Trigger_SingleMuon = eval("Trigger_"+year+"_smu")
            Trigger_DoubleMuon = eval("Trigger_"+year+"_mumu")
            Trigger_DoubleElectron = eval("Trigger_"+year+"_ee")
            Trigger_ElectronMuon = eval("Trigger_"+year+"_emu")
            ### ~~~~~~~~~~~~~~~ Run dependent trigger handling for 2017 data ~~~~~~~~~~~~~~~~~~~~~~
            if year == "2017":
                 if isRun == "B": 
                      Trigger_DoubleMuon = []
                      Trigger_DoubleMuon = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_"]
                 elif isRun == "C" or isRun == "D":
                      Trigger_DoubleMuon = []
                      Trigger_DoubleMuon = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_"]
                 elif isRun == "E" or isRun == "F":
                      Trigger_DoubleMuon = []
                      Trigger_DoubleMuon = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_"]
                      Trigger_SingleMuon = []
                      Trigger_SingleMuon = ["HLT_IsoMu27_v_"]
            ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            TriggerToBeConsidered = []
            TriggerToBeConsidered += Trigger_SingleElectron
            TriggerToBeConsidered += Trigger_SingleMuon
            TriggerToBeConsidered += Trigger_DoubleMuon
            TriggerToBeConsidered += Trigger_DoubleElectron
            TriggerToBeConsidered += Trigger_ElectronMuon
            TriggerToBeConsidered += Trigger_met
            

            which2DHistogramsCreatingForFiring = []
            for item_whichDileptonTriggerGroup in whichDileptonTriggerGroup:
                 for item_TriggerToBeConsidered in TriggerToBeConsidered:
                      which2DHistogramsCreatingForFiring.append(item_whichDileptonTriggerGroup+"_vs_"+item_TriggerToBeConsidered)

            channel = ""
            #sys.exit("Here we stop test run")
            directory = textname +nan+"/"
             #   textname +="_"+nStart+"_"+nStopMinus1
#            nan = "_TestSeptemberRecoMllbB"#"_fixedWTops_RecoNTuplesDeltaRCut___DIRECTORYSUFFIX__"                                                  
#    nan ="_fixedWTops_NTuplesRecoLevelIfSonn"                                                                                              
            #directoryPlot2D = directory + "EventPlot2D/"
            directory = textname +nan+"/"
            if not os.path.exists(directory):
                os.makedirs(directory)
            for item_whichSubDirectoriesToCreate in whichSubDirectoriesToCreate:
                if not os.path.exists(directory +item_whichSubDirectoriesToCreate):
                    os.makedirs(directory +item_whichSubDirectoriesToCreate)
            outputFileName = textname + nan + ".root"
            rootdirectory =directory + outputFileName
            newFile = r.TFile.Open(rootdirectory, "RECREATE")
            #newTree = r.TTree( 'NTuple', 'Just A Tree without Clone' )

            nEntries = inputTree.GetEntries()
            countevents =0

            count_HLT_DoubleEl_ = 0
            count_HLT_DoubleMuon_ = 0
            count_HLT_EMu_ = 0
            count_HLT_SingleMuon_ = 0
            count_HLT_SingleElectron_ = 0
            count_HLT_METTriggers_ = 0

            
            for i in range(100): #nEntries): #10): #int(nStart), int(nStopMinus1)):
                #    for i in range(nEntries): #events): #xrange(nEntries): #in range(events):
                #for i in range(nEntries):
                #Call range of events
                #for i in range(50):
                #if i == 0:
                inputTree.GetEntry(i)
                #print "---------------------------------"
                #print("Event:", i)
                #print "---------------------------------"
                countevents +=1
                if(countevents%1000 == 0): print("Event:", i)
                it = inputTree

#                passSLCuts_e = False
 #               passSLCuts_mu = False
              #  passDLCuts_mumu = False
               # passDLCuts_ee = False
                #passDLCuts_emu = False

#                for item in dictTrigger_2017.keys():
 #                    print(item+": "+str(dictTrigger_2017.get(item)))

  #              for item in dictTrigger_2017.keys():
   #                  dictTrigger_2017[item] = False
               # for item_Trigger_2017 in Trigger_2017:
                #     globals()[item_Trigger_2017] = False
#                     exec("%s = %r" % (item_Trigger_2017, False))
                 #    print(str(item_Trigger_2017)+": "+str(eval(item_Trigger_2017)))
#                eeChannel_1st_Lepton_pt_eta_iso=[]
 #               eeChannel_2nd_Lepton_pt_eta_iso=[]
  #              eeChannel_njet_nvertex_met = []
#
 #               emuChannel_1st_Lepton_pt_eta_iso=[]
  #              emuChannel_2nd_Lepton_pt_eta_iso=[]
   #             emuChannel_njet_nvertex_met = []
#
 #               mumuChannel_1st_Lepton_pt_eta_iso=[]
  #              mumuChannel_2nd_Lepton_pt_eta_iso=[]
   #             mumuChannel_njet_nvertex_met = []

                for item_TriggerToBeConsidered in TriggerToBeConsidered:
                     globals()[item_TriggerToBeConsidered] = False
                     print(str(item_TriggerToBeConsidered)+": "+str(eval(item_TriggerToBeConsidered)))

                HLT_DoubleEl_fired = False
                HLT_DoubleMuon_fired = False
                HLT_EMu_fired = False
                HLT_SingleMuon_fired = False
                HLT_SingleElectron_fired = False
                HLT_METTriggers_fired = False
#                seTriggerfired_OR = False
 #               emuTriggerfired_OR = False
  #              eeTriggerfired_OR = False
   #             mumuTriggerfired_OR = False
    #            smuTriggerfired_OR = False



################### 1. Check all Trigger results of interest
################### ===========================================================

                for item_TriggerToBeConsidered in TriggerToBeConsidered:
                   print("inputTree."+item_TriggerToBeConsidered[:-1])
                   print(eval("inputTree."+item_TriggerToBeConsidered[:-1]))
                   if eval("inputTree."+item_TriggerToBeConsidered[:-1]) == 1:
                        globals()[item_TriggerToBeConsidered] = True
                        globals()["count_"+str(item_TriggerToBeConsidered)] += 1

################## 2. Read in results from 1. and check category trigger result
################## ============================================================

                #### MET trigger
                for item_Trigger_met in Trigger_met:
                     print(item_Trigger_met)
                     print(eval(item_Trigger_met))
                     if eval(item_Trigger_met) == True:
                           count_HLT_METTriggers_ += 1
                           HLT_METTriggers_fired = True
                           break

             #################   Don't proceed checking this event if no MET trigger fired.
                if not HLT_METTriggers_fired: 
                     contine

                #### Double Electron trigger
                for item_Trigger_DoubleElectron in Trigger_DoubleElectron:
                     print(item_Trigger_DoubleElectron)
                     print(eval(item_Trigger_DoubleElectron))
                     if eval(item_Trigger_DoubleElectron) == True:
                          count_HLT_DoubleEl_ += 1
                          HLT_DoubleEl_fired = True
                          break

                #### Double Muon trigger
                for item_Trigger_DoubleMuon in Trigger_DoubleMuon:
                     print(item_Trigger_DoubleMuon)
                     print(eval(item_Trigger_DoubleMuon))
                     if eval(item_Trigger_DoubleMuon) == True:
                            count_HLT_DoubleMuon_ += 1
                            HLT_DoubleMuon_fired = True
                            break

                #### Electron Muon trigger
                for item_Trigger_ElectronMuon in Trigger_ElectronMuon:
                     print(item_Trigger_ElectronMuon)
                     print(eval(item_Trigger_ElectronMuon))
                     if eval(item_Trigger_ElectronMuon) == True:
                          count_HLT_EMu_ += 1
                          HLT_EMu_fired = True
                          break

                #### Single Electron trigger
                for item_Trigger_SingleElectron in Trigger_SingleElectron:
                     print(item_Trigger_SingleElectron)
                     print(eval(item_Trigger_SingleElectron))
                     if eval(item_Trigger_SingleElectron) == True:
                            count_HLT_SingleElectron_ += 1
                            HLT_SingleElectron_fired = True
                            break

                #### Single Muon trigger
                for item_Trigger_SingleMuon in Trigger_SingleMuon:
                     print(item_Trigger_SingleMuon)
                     print(eval(item_Trigger_SingleMuon))
                     if eval(item_Trigger_SingleMuon) == True:
                          count_HLT_SingleMuon_ += 1
                          HLT_SingleMuon_fired = True
                          break


#######################################  Get Lepton Cut Results ##############################################################                
                nMuons = 0
                nElectrons = 0
                #lepton_pt_array = inputTree.eve.lepton_pt_
                #lepton_isMuon_array = inputTree.eve.lepton_isMuon_
                #lepton_eta_array = inputTree.eve.lepton_eta_

                ## ----------------------------- electrons/ muon cuts start
                
                get_indicesAndPt_of_passed_lepton = []
                get_indicesAndPt_of_passed_muons = []
                get_indicesAndPt_of_passed_electrons = []
                
                
                #print(inputTree.muon_PFIso[0])
                #print(len(inputTree.muon_PFIso))
#                sys.exit("Stopping test")

#                for item_lepton in range(len(inputTree.muon_PFIso)):
                
#                print(inputTree.leptons.size())
#                for item_lepton in inputTree.leptons_fCoordinates_fPt:
 #                    print(item_lepton)
 #               sys.exit("Stopping test")
   #             for i in inputTree.leptons_:
  #                  print(item_lepton)
    #                oneLepton = lv()
     #               oneLepton = inputTree.leptons_[item_lepton]
      #              print(oneLepton.Pt())
                numberOfJets = int(inputTree.jets.size())
                print(numberOfJets)
                numberOfPrimaryVertices = int(inputTree.NPV_all)
                print(numberOfPrimaryVertices)
                met = inputTree.met_Puppi.Pt()
                
                for item_lepton in range(inputTree.leptons.size()):
#################
######################### If Muon, let's check if inputTree.passes the chosen selections
######################### ==========================================================
                    print("inputTree.muon_PFIso["+str(item_lepton)+"]")
                    print(inputTree.muon_PFIso[item_lepton])
                    print("Is inputTree.muon_PFIso[item_lepton] >= 76?")
                    print(ord(inputTree.muon_PFIso[item_lepton]) >= 76)
                    if(inputTree.muon_CutBasedID_Medium[item_lepton] == "Y" and inputTree.leptons[item_lepton].Pt() > 15. and (abs(inputTree.leptons[item_lepton].Eta()) < 2.4) and ord(inputTree.muon_PFIso[item_lepton]) >= 76):
                              nMuons += 1
                              get_indicesAndPt_of_passed_muons.append((inputTree.leptons[item_lepton].Pt(), item_lepton))
                              get_indicesAndPt_of_passed_lepton.append((inputTree.leptons[item_lepton].Pt(), item_lepton, "M"))
                              print("Passed Muon with "+str(inputTree.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(inputTree.muon_CutBasedID_Medium[item_lepton]))
                              
######################### If Electron, let's check if it passes the chosen selections
######################### ============================================================

                    elif(inputTree.eleID_MVA_Iso_90[item_lepton] == 1 and inputTree.leptons[item_lepton].Pt() > 15. and (abs(inputTree.leptons[item_lepton].Eta()) < 2.4)):
                         ### removing ECAL transistion region 1.4442 < abs(eta) < 1.566
                         if (abs(inputTree.leptons[item_lepton].Eta()) < 1.442) or (abs(inputTree.leptons[item_lepton].Eta()) > 1.566):
                              nElectrons += 1
                              get_indicesAndPt_of_passed_lepton.append((inputTree.leptons[item_lepton].Pt(), item_lepton, "E"))
                              get_indicesAndPt_of_passed_electrons.append((inputTree.leptons[item_lepton].Pt(), item_lepton))
                              print("Passed Electron with "+str(inputTree.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(inputTree.eleID_MVA_Iso_90[item_lepton]))

                ## ----------------------------- electrons / muon cuts end
                
################### Get leading and subleading lepton index if available
################### =======================================================
                sortingLeptons = sorted(get_indicesAndPt_of_passed_lepton, reverse=True)
                # Don't check this event further on if event failed Dilepton cuts
                if len(sortingLeptons) < 2:
                     continue
                elif len(sortingLeptons) > 1: 
                    leadingLeptonIndex = sortingLeptons[0][1]
                    leadLeptonisMuon = sortingLeptons[0][2]
                    print("leadingLeptonIndex = "+str(leadingLeptonIndex))
                    if len(sortingLeptons) > 1:
                         subleadingLeptonIndex = sortingLeptons[1][1]
                         subleadingLeptonisMuon = sortingLeptons[1][2]
                         print("subleadingLeptonIndex = "+str(subleadingLeptonIndex))
                    print(sortingLeptons)
                    
#                    leadingLepton = tlv()
                    leadingLepton = lv()
                    if len(sortingLeptons) > 0:
                         leadingLepton = inputTree.leptons[leadingLeptonIndex]
#                         leadingLepton.SetPtEtaPhiE(inputTree.eve.lepton_pt_[leadingLeptonIndex], inputTree.eve.lepton_eta_[leadingLeptonIndex], inputTree.eve.lepton_phi_[leadingLeptonIndex], inputTree.eve.lepton_e_[leadingLeptonIndex])

                    subleadingLepton = lv()
#                    subleadingLepton = tlv()
                    if len(sortingLeptons) > 1:
                         subleadingLepton = inputTree.leptons[subleadingLeptonIndex]
#                         subleadingLepton.SetPtEtaPhiE(inputTree.eve.lepton_pt_[subleadingLeptonIndex], inputTree.eve.lepton_eta_[subleadingLeptonIndex], inputTree.eve.lepton_phi_[subleadingLeptonIndex], inputTree.eve.lepton_e_[subleadingLeptonIndex])
                    #dilepton = tlv()

                    dilepton = lv()
                    if len(sortingLeptons) > 1:
                         dilepton = leadingLepton + subleadingLepton

                    ### Get leading and subleading muon index if available
                    sortingMuons = sorted(get_indicesAndPt_of_passed_muons, reverse=True)
                    if len(sortingMuons) > 0:
                         leadingMuonIndex = sortingMuons[0][1]
                         print("leadingMuonIndex = "+str(leadingMuonIndex))
                         if len(sortingMuons) > 1:
                              subleadingMuonIndex = sortingMuons[1][1]
                              print("subleadingMuonIndex = "+str(subleadingMuonIndex))
                    ### Get leading and subleading electron index if available
                    sortingElectrons = sorted(get_indicesAndPt_of_passed_electrons, reverse=True)
                    if len(sortingElectrons) > 0:
                         leadingElectronIndex = sortingElectrons[0][1]
                         print("leadingElectronIndex = "+str(leadingElectronIndex))
                         if len(sortingElectrons) > 1:
                              subleadingElectronIndex = sortingElectrons[1][1]
                              print("subleadingElectronIndex = "+str(subleadingElectronIndex))

                    ## ----------------------------- categorisation into channels
                    
                    if leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15. and dilepton.M() > 20.:
                         # if this passes we have an OS muon pair --> mumu Channel
                         if (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -169):
                              MuonTriggersfired = False
                              if HLT_SingleMuon_fired or HLT_DoubleMuon_fired:
                                   MuonTriggersfired = True
                              if not isMC:
                                   
                                   mumu_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   mumu_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   mumu_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, MuonTriggersfired, HLT_SingleMuon_fired])
                              else:
                                   mumu_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   mumu_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   mumu_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, MuonTriggersfired, HLT_SingleMuon_fired])

                         # if this passes we have an OS electron pair --> ee Channel
                         elif (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -121):
                              ElectronTriggersfired = False
                              if HLT_SingleElectron_fired or HLT_DoubleEl_fired:
                                   ElectronTriggersfired = True 
                              if not isMC:
                                   ee_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   ee_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   ee_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, ElectronTriggersfired, HLT_SingleElectron_fired])
                              else:
                                   ee_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   ee_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   ee_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, ElectronTriggersfired, HLT_SingleElectron_fired])
                         # if this passes we have an OS electron & muon --> emu Channel
                         elif (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -143):
                              ElectronMuonTriggersfired = False
                              SingleTriggersfired = False
                              if HLT_SingleElectron_fired or HLT_SingleMuon_fired or HLT_EMu_fired:
                                   ElectronMuonTriggersfired = True
                              if HLT_SingleElectron_fired or HLT_SingleMuon_fired:
                                   SingleTriggersfired = True
                              if not isMC:
                                   emu_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   emu_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   emu_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, ElectronMuonTriggersfired, SingleTriggersfired])
                              else:
                                   emu_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                   emu_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                   emu_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, ElectronMuonTriggersfired, SingleTriggersfired])
                    ####### DL muon
#                   elif nMuons == 2 and nElectrons == 0 
#                               and leadLeptonisMuon and inputTree.eve.lepton_pt_[leadingLeptonIndex] >= 25. 
#                               and subleadingLeptonisMuon and inputTree.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
                    #elif nMuons == 2 and nElectrons == 0 and leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
#                    if nMuons > 1 and leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
 #                        checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)          
 #                        # Nope, THIS code will decide inputTree.itself AND does all channels in one go and does not waste resources. Run 3 is calling! 
                         #runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "mumu"))
  #                       if checkDileptonMassPass and (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -169): #and runOnCorrectSamples:
  #                            passDLCuts_mumu = True
   #                           mumuChannel_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
    #                          mumuChannel_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
     #                         mumuChannel_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met]) 

                    ####### DL ee
      #              elif nMuons == 0 and nElectrons == 2 and leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
#inputTree.eve.lepton_pt_[leadingLeptonIndex] >= 25. and inputTree.eve.lepton_pt_[subleadingLeptonIndex] >= 15.: 
#                   elif nMuons == 0 and nElectrons == 2
 #                              and not leadLeptonisMuon and inputTree.eve.lepton_pt_[leadingLeptonIndex] >= 25.
  #                             and not subleadingLeptonisMuon and inputTree.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
       #                  checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)
        #                 runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "ee"))
         #                if checkDileptonMassPass and (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex]  == -121) and runOnCorrectSamples:
          #                    passDLCuts_ee = True
           #                   eeChannel_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
            #                  eeChannel_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
             #                 eeChannel_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met])

                    ####### DL emu
              #      elif nMuons == 1 and nElectrons ==1:
               ##      if (leadLeptonisMuon == "E" and subleadingLeptonisMuon == "M") or (leadLeptonisMuon == "M" and subleadingLeptonisMuon == "E"):
                #          if leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
#                #          if inputTree.eve.lepton_pt_[leadingLeptonIndex] >= 25. and inputTree.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
                 #              checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)
                  #             runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "emu"))
                   #            if checkDileptonMassPass and (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex]  == -143) and runOnCorrectSamples:
                    #                passDLCuts_emu = True
                     #               emuChannel_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                      #              emuChannel_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                       #             emuChannel_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met])
                ## ----------------------------- category cut end
                #get_indicesAndPt_of_passed_jets = []
                #for item_jet in range(len(inputTree.eve.jet_pt_)):
                 #    jet_pt = inputTree.eve.jet_pt_[item_jet]
                  #   jet_eta = inputTree.eve.jet_eta_[item_jet]
                   #  jet_phi = inputTree.eve.jet_phi_[item_jet]
                    # jet_puid = inputTree.eve.jet_puid_[item_jet]
                     #if applyJetCuts(jet_pt, jet_eta, jet_puid):
                      #    get_indicesAndPt_of_passed_jets.append((jet_pt, item_jet))


#                dilepton = tlv()
 #               dilepton = leadingLepton + subleadingLepton
  #              if dilepton.M() > 20. and dilepton.M() :
   #                 dileptonMassCutpassed = True
    #            else: dileptonMassCutpassed = False



                        
################## Check if at least one trigger of considered trigger group fired; 
##################                                  if yes, accept and go out of loop for ...
################## ============================================================================
            
            
#                if passDLCuts_ee and HLT_DoubleEl_fired:
 #                h_leadingLepton_pt_MET_D.Fill(leadingLepton.Pt())
  #               h_subleadingLepton_pt_MET_D.Fill(subleadingLepton.Pt())
            
#            sortingLeptons = sorted(, reverse=True)
#
 #           c = r.TCanvas('Lepton_pt', 'lepton pt')
  #          c.Divide(2)
   #         c.cd(1)
    #        h_leadingLepton_pt_MET_D.Draw()
     #       c.cd(2)
#            h_subleadingLepton_pt_MET_D.Draw()
 #           for ext in ['png'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
  #          for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
   #         h_leadingLepton_pt_MET_D.Write()
    #        h_subleadingLepton_pt_MET_D.Write()

            newFile.Close()
            inputFile.Close()
            if LeptonSFapplication:
                 inputElectronIDsfFile.Close()
                 inputElectronRECOsfFile.Close()
                 inputMuonIDsfFile.Close()
                 inputMuonISOsfFile.Close()

            if verbose:
                print("count_HLT_DoubleEl_")
                print(count_HLT_DoubleEl_)
                print("HLT_Ele27_WPTight_Gsf_v_")
                print(HLT_Ele27_WPTight_Gsf_v_)
                print("count_HLT_Ele27_WPTight_Gsf_v_")
                print(count_HLT_Ele27_WPTight_Gsf_v_)
                print("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_")
                print(HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_)
                print("count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_")
                print(count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_)
                print("HLT_Ele35_WPTight_Gsf_v_")
                print(HLT_Ele35_WPTight_Gsf_v_)
                print("count_HLT_Ele35_WPTight_Gsf_v_")
                print(count_HLT_Ele35_WPTight_Gsf_v_)
                print("HLT_Ele38_WPTight_Gsf_v_")
                print(HLT_Ele38_WPTight_Gsf_v_)
                print("count_HLT_Ele38_WPTight_Gsf_v_")
                print(count_HLT_Ele38_WPTight_Gsf_v_)
                print("HLT_Ele40_WPTight_Gsf_v_")
                print(HLT_Ele40_WPTight_Gsf_v_)
                print("count_HLT_Ele40_WPTight_Gsf_v_")
                print(count_HLT_Ele40_WPTight_Gsf_v_)
                print("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_")
                print(HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_)
                print("count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_")
                print(count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_)
                print("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_")
                print(HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_)
                print("count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_")
                print(count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_)
                print("HLT_IsoMu27_v_")
                print(HLT_IsoMu27_v_)
                print("count_HLT_IsoMu27_v_")
                print(count_HLT_IsoMu27_v_)
                print("HLT_IsoMu24_eta2p1_v_")
                print(HLT_IsoMu24_eta2p1_v_)
                print("count_HLT_IsoMu24_eta2p1_v_")
                print(count_HLT_IsoMu24_eta2p1_v_)
                print("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_")
                print(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_)
                print("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_")
                print(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_)
                print("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_")
                print(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_)
                print("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_")
                print(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_)
                print("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_")
                print(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_)
                print("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_")
                print(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_)
                print("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                print(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                print("count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                print(count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                print("HLT_DoubleEle33_CaloIdL_MW_v_")
                print(HLT_DoubleEle33_CaloIdL_MW_v_)
                print("count_HLT_DoubleEle33_CaloIdL_MW_v_")
                print(count_HLT_DoubleEle33_CaloIdL_MW_v_)
                print("HLT_DoubleEle25_CaloIdL_MW_v_")
                print(HLT_DoubleEle25_CaloIdL_MW_v_)
                print("count_HLT_DoubleEle25_CaloIdL_MW_v_")
                print(count_HLT_DoubleEle25_CaloIdL_MW_v_)
                print("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_")
                print(HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_)
                print("count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_")
                print(count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_)
                print("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                print(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                print("count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                print(count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                print("HLT_PFMET120_PFMHT120_IDTight_v_")
                print(HLT_PFMET120_PFMHT120_IDTight_v_)
                print("count_HLT_PFMET120_PFMHT120_IDTight_v_")
                print(count_HLT_PFMET120_PFMHT120_IDTight_v_)
                print("HLT_PFMET250_HBHECleaned_v_")
                print(HLT_PFMET250_HBHECleaned_v_)
                print("count_HLT_PFMET250_HBHECleaned_v_")
                print(count_HLT_PFMET250_HBHECleaned_v_)
                print("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_")
                print(HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_)
                print("count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_")
                print(count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_)
                print("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_")
                print(HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_)
                print("count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_")
                print(count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_)
                print("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_")
                print(HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_)
                print("count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_")
                print(count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_)
                print("count_HLT_METTriggers_")
                print(count_HLT_METTriggers_)
                print("count_HLT_DoubleEl_")
                print(count_HLT_DoubleEl_)

print("ee_Channel_data_1st_Lepton_pt_eta_iso")
print(ee_Channel_data_1st_Lepton_pt_eta_iso)




ee_Channel_data_1st_Lepton_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0])]
#[indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=itemgetter(1))]
print(ee_Channel_data_1st_Lepton_pt_Order)
ee_Channel_data_1st_Lepton_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][1])]
#[indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=itemgetter(2))]
print(ee_Channel_data_1st_Lepton_eta_Order)
ee_Channel_data_1st_Lepton_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][2])]
#[indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=itemgetter(3))]
print(ee_Channel_data_1st_Lepton_iso_Order)


#==============================================================
#==============================================================
#### DETERMINING SYSTEMATIC uncertainties using the same array:
#==============================================================
#==============================================================

#### by subset selections on NUMBER OF JETS:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff number of jets < 3:
ee_Channel_data_1st_Lepton_pt_Order_smaller3Njets = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][0] <3]
### ordering of events according to 1st lepton pt iff number of jets >= 3:
ee_Channel_data_1st_Lepton_pt_Order_3AndMoreNjets = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][0] >= 3]


#### by subset selections on NUMBER OF VERTICES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff number of vertices < 30: 
ee_Channel_data_1st_Lepton_pt_Order_smaller30Nvertices = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][1] <30]
### ordering of events according to 1st lepton pt iff number of vertices >= 30:
ee_Channel_data_1st_Lepton_pt_Order_30AndMoreNvertices = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][1] >= 30]


#### by subset selections on MET:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff MET < 80 GeV:
ee_Channel_data_1st_Lepton_pt_Order_METsmaller80 = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][2] <80.]
### ordering of events according to 1st lepton pt iff MET >= 80 GeV:
ee_Channel_data_1st_Lepton_pt_Order_METlarger80 = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][2] >= 80.]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("ee_Channel_data_2nd_Lepton_pt_eta_iso")
print(ee_Channel_data_2nd_Lepton_pt_eta_iso)
ee_Channel_data_2nd_Lepton_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][0])]
#[indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=itemgetter(1))]                                                               
print(ee_Channel_data_2nd_Lepton_pt_Order)
ee_Channel_data_2nd_Lepton_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][1])]
print(ee_Channel_data_2nd_Lepton_eta_Order)
ee_Channel_data_2nd_Lepton_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][2])]
print(ee_Channel_data_2nd_Lepton_iso_Order)

print("ee_Channel_data_njet_nvertex_met")
print(ee_Channel_data_njet_nvertex_met)
ee_Channel_data_njet_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][0])]
print(ee_Channel_data_njet_Order)
ee_Channel_data_nvertex_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][1])]
print(ee_Channel_data_nvertex_Order)
ee_Channel_data_met_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][2])]
print(ee_Channel_data_met_Order)



print("emu_Channel_data_1st_Lepton_pt_eta_iso")
print(emu_Channel_data_1st_Lepton_pt_eta_iso)
print("emu_Channel_data_2nd_Lepton_pt_eta_iso")
print(emu_Channel_data_2nd_Lepton_pt_eta_iso)
print("emu_Channel_data_njet_nvertex_met")
print(emu_Channel_data_njet_nvertex_met)

print("mumu_Channel_data_1st_Lepton_pt_eta_iso")
print(mumu_Channel_data_1st_Lepton_pt_eta_iso)
print("mumu_Channel_data_2nd_Lepton_pt_eta_iso")
print(mumu_Channel_data_2nd_Lepton_pt_eta_iso)
print("mumu_Channel_data_njet_nvertex_met")
print(mumu_Channel_data_njet_nvertex_met)


## Having seen how the wind blows, let's write the remaining ones in shorter way ;)
#==================================================================
#==================================================================
#### DETERMINING remaining SYSTEMATIC uncertainties in ee Channel:
#==================================================================
#================================================================== 

whatToDivide = []

channel = ["ee", "emu", "mumu"] 
whichType = ["data", "MC"]

whichLepton = ["1st_Lepton", "2nd_Lepton"]
 
#whatToOrder

for item_channel in channel:
     for item_whichType in whichType:
          item_channel+"_Channel_"+item_whichType+
          
whatToDivide = [#"ee_Channel_data_1st_Lepton_pt_Order",
                "ee_Channel_data_1st_Lepton_eta_Order",
                "ee_Channel_data_1st_Lepton_iso_Order",
                "ee_Channel_data_2nd_Lepton_pt_Order",
                "ee_Channel_data_2nd_Lepton_eta_Order"
                "ee_Channel_data_2nd_Lepton_iso_Order",
                ]

inWhichSubsetsToDivide = ["smaller3Njets", "3AndMoreNjets", "smaller30Nvertices", "30AndMoreNvertices", "METsmaller80", "METlarger80"]


whichCommandToExecute =[]

for item_whatToDivide in whatToDivide:
     for item_inWhichSubsetsToDivide in inWhichSubsetsToDivide:
          if "_pt_" in item_whatToDivide: item_elem = "0"
          elif "_eta_" in item_whatToDivide: item_elem = "1"
          elif "_iso_" in item_whatToDivide: item_elem = "2"
          if not item_elem: sys.exit("Something went wrong in assigning Systematic uncertainties in whatToDivide loop. Please check.")
          if "Njets" in item_inWhichSubsetsToDivide: 
               item_indices_ordered = "0"
               if "smaller3" in item_inWhichSubsetsToDivide:
                    item_conditional = "<3]"
               elif "3AndMore" in item_inWhichSubsetsToDivide:
                    item_conditional = ">= 3]"
          elif "Nvertices" in item_inWhichSubsetsToDivide: 
               item_indices_ordered = "1"
               if "smaller30" in item_inWhichSubsetsToDivide:
                    item_conditional = "<30]"
               elif "30AndMore" in item_inWhichSubsetsToDivide:
                    item_conditional = ">= 30]"
          elif "MET" in item_inWhichSubsetsToDivide: 
               item_indices_ordered = "2"
               if "smaller80" in item_inWhichSubsetsToDivide:
                    item_conditional = "<80.]"
               elif "larger80" in item_inWhichSubsetsToDivide:
                    item_conditional = ">= 80.]"
          if not item_indices_ordered: sys.exit("Something went wrong in assigning Systematic uncertainties in inWhichSubsetsToDivide loop. Please check.")
          if not item_conditional: sys.exit("Something went wrong in assigning Systematic uncertainties for item_conditional. Please check.")
          whichCommandToExecute.append(item_whatToDivide+"_"+item_inWhichSubsetsToDivide+" = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1]["+item_elem+"]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]]["+item_indices_ordered+"] "+item_conditional)
          

print(whichCommandToExecute)
          
#for item_TriggerToBeConsidered in TriggerToBeConsidered:
 #                  print("inputTree."+item_TriggerToBeConsidered[:-1])
  #                 print(eval("inputTree."+item_TriggerToBeConsidered[:-1]))
   #                if eval("inputTree."+item_TriggerToBeConsidered[:-1]) == 1:
    #                    globals()[item_TriggerToBeConsidered] = True
     #                   globals()["count_"+str(item_TriggerToBeConsidered)] += 1





