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
#from operator import itemgetter
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

####  https://gist.github.com/rmanzoni/d31db4d0c692fc5a9580
####  No weighted SF TEfficiency for frequentists.
####  https://root-forum.cern.ch/t/tefficiency-with-weighted-events/18067/2
####  https://root.cern.ch/doc/master/classTEfficiency.html on 2020-06-24

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




## Vorsicht aender das fuer 2017/8; eta 2.5 anpassen!! Nur jetzt zum Testen.
def applyMuonCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_relIso):

    if (lepton_isMuon == 1) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (lepton_relIso < 0.25):
        return True
    else: return False

#def applyMuonCuts(muon_CutBasedID_Medium, lepton, 


def applyElectronCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_scEta):
    if (lepton_isMuon == 0) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (abs(lepton_scEta) < 1.4442 or abs(lepton_scEta) > 1.566):
        return True
    else: return False


def applyJetCuts(jet_pt, jet_eta, jet_puid):
     if (jet_pt > 20.) and abs(jet_eta) < 2.4 and jet_puid == 7:
          return True
     else: return False


#def applyPileupReweighting():
 #    dataPileupFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/pileup/pileupDATA__Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1__max120_n120__MinBias69200.root"
     

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

###### convert TH2D or TH2F histogram into python xD list, origninally made for SFs, but this works for general case, too.
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
               print(str(inputSFarrayFileRootPath_GetHisto.GetXaxis().GetBinUpEdge(item_sfRangeX))+"|"+str(inputSFarrayFileRootPath_GetHisto.GetYaxis().GetBinLowEdge(item_sfRangeY))+" =sf= "+str(item_sfRangeY))
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
#def getLeptonSF(SF_FileHisto, lepton_pt, lepton_eta):
#     xBin = SF_FileHisto.GetXaxis().FindBin(lepton_eta)
#     yBin = SF_FileHisto.GetYaxis().FindBin(lepton_pt)
#     if xBin == 0:
#          xBin = 1
#     elif xBin > SF_FileHisto.GetNbinsX():
#          xBin = SF_FileHisto.GetNbinsX()
#     if yBin == 0:
#          yBin = 1
#     elif yBin > SF_FileHisto.GetNbinsY():
#          yBin = SF_FileHisto.GetNbinsY()
#     SFvalue = SF_FileHisto.GetBinContent(xBin, yBin)
#     print(xBin)
#     print(yBin)
#     if not SFvalue:
#          SFvalue = 1.
#     return SFvalue

RootfilesToRunOver =["/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/met_run2017D.root"
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

nan = "_TestSeptemberRecoMllbB"

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

#     print("pileupWeights")
 #    print(pileupWeights)
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

###### Example
lep_Eta = -1.957742
lep_Pt = 44.2188
exampleSf = leptonSF(gatheredElectronIDsfWeights, lep_Eta, lep_Pt)
print("exampleSf")
print(exampleSf)

ElectronRECOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"
gatheredElectronRECOsfWeights = create_sfArray(ElectronRECOsfFile, inputElectronIDhistogram)


exampleSf1 = leptonSF(gatheredElectronRECOsfWeights, lep_Eta, lep_Pt)
print("exampleSf1")
print(exampleSf1)

MuonIDsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ID_syst.root"
inputMuonIDhistogram = "NUM_MediumID_DEN_genTracks_pt_abseta"
gatheredMuonIDsfWeights = create_sfArray(MuonIDsfFile, inputMuonIDhistogram)
exampleMuonSf = leptonSF(gatheredMuonIDsfWeights, lep_Pt, abs(lep_Eta))
print("exampleMuonSf")
print(exampleMuonSf)

MuonISOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/muon/94X/RunBCDEF_SF_ISO_syst.root"
inputMuonIsohistogram = "NUM_LooseRelIso_DEN_MediumID_pt_abseta"
gatheredMuonIsosfWeights = create_sfArray(MuonISOsfFile, inputMuonIsohistogram)
exampleMuonSf1 = leptonSF(gatheredMuonIsosfWeights, lep_Pt, abs(lep_Eta))
print("exampleMuonSf1")
print(exampleMuonSf1)







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


#print("ddddddddddddddddddddddddddddddddd\n")

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


years_considered = ["2016", "2017", "2018"]

RunPart_considered = ["A", "B", "C", "D", "E", "F", "G", "H"]  

strength = ["lead", "subleading"]

lepton = ["Electron", "Muon"]

oneLeptonFeature = ["Pt", "Eta", "relIso"]

elseFeature = ["dileptonMass", "Met", "numberPrimaryVertex", "numberBTags", "numberOfJets"]

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
                for item_RunPart_considered in RunPart_considered:
    #                 print(item_RunPart_considered)
                     if item_RunPart_considered in textname: 
                          isRun = item_RunPart_considered
                          break
#                if "A" in textname: isRun = "A"
 #               elif "B" in textname: isRun = "B"
  #              elif "C" in textname: isRun = "C" 
   #             elif "D" in textname: isRun = "D"
    #            elif "E" in textname: isRun = "E"
     #           elif "F" in textname: isRun = "F"
      #          elif "G" in textname: isRun = "G"
       #         elif "H" in textname: isRun = "H"
            year = ""
            for item_years_considered in years_considered:
     #            print(item_years_considered)
                 if item_years_considered in textname:
                      year = item_years_considered
                      break
#            if "2016" in textname:
 #                year = "2016"
  #          elif "2017" in textname:
   #              year = "2017"
    #        elif "2018" in textname:
     #            year = "2018"
            print("isRun is "+isRun)
            print("isMC is "+str(isMC))
            print("year is "+str(year))
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
            nEntries = inputTree.GetEntries()
            countevents =0

            count_HLT_Ele27_WPTight_Gsf_v_ = 0
            count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_ = 0
            count_HLT_Ele35_WPTight_Gsf_v_ = 0
            count_HLT_Ele38_WPTight_Gsf_v_ = 0
            count_HLT_Ele40_WPTight_Gsf_v_ = 0
            count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = 0
            count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_ = 0

            # smu Triggers                                                                                                                                                                                                     
            count_HLT_IsoMu27_v_ = 0
            count_HLT_IsoMu24_eta2p1_v_ = 0

            # mumu Triggers                                                                                                                                                                                                    
            count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = 0
            count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = 0
            count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = 0

            # ee Triggers                                                                                                                                                                                                      
            count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0
            #    count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0    
            count_HLT_DoubleEle33_CaloIdL_MW_v_ = 0
            count_HLT_DoubleEle25_CaloIdL_MW_v_ = 0
            count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_ = 0

            # emu Triggers                                                                                                                                                                                                     
            count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0
            count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0
            count_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0
            count_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0

            # met Triggers                                                                                                                                                                                                     
            count_HLT_PFMET120_PFMHT120_IDTight_v_ = 0
            count_HLT_PFMET250_HBHECleaned_v_ = 0
            count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_ = 0
            count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_ = 0
            count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_ = 0

            count_HLT_DoubleEl_ = 0
            count_HLT_DoubleMuon_ = 0
            count_HLT_EMu_ = 0
            count_HLT_SingleMuon_ = 0
            count_HLT_SingleElectron_ = 0
            count_HLT_METTriggers_ = 0

            
            for i in range(10):  #nEntries): #10): #int(nStart), int(nStopMinus1)):
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

                passSLCuts_e = False
                passSLCuts_mu = False
                passDLCuts_mumu = False
                passDLCuts_ee = False
                passDLCuts_emu = False

                HLT_Ele27_WPTight_Gsf_v_ = False
                HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_ = False
                HLT_Ele35_WPTight_Gsf_v_ = False
                HLT_Ele38_WPTight_Gsf_v_ = False
                HLT_Ele40_WPTight_Gsf_v_ = False
                HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = False
                HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_ = False
    
                # smu Triggers                                                                                                                            
                HLT_IsoMu27_v_ = False
                HLT_IsoMu24_eta2p1_v_ = False
                
                # mumu Triggers                                                                                                                           
                HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = False
                HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = False
                HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = False
                
                # ee Triggers                                                                                                                             
                HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = False
                #    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = False                                                                             
                HLT_DoubleEle33_CaloIdL_MW_v_ = False
                HLT_DoubleEle25_CaloIdL_MW_v_ = False
                HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_ = False
                
                # emu Triggers                                                                                                                            
                HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = False
                HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = False
                HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = False
                HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = False
                
                # met Triggers                                                                                                                            
                HLT_PFMET120_PFMHT120_IDTight_v_ = False
                HLT_PFMET250_HBHECleaned_v_ = False
                HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_ = False
                HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_ = False
                HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_ = False

                HLT_DoubleEl_ = False
                HLT_DoubleMuon_ = False
                HLT_EMu_ = False
                HLT_SingleMuon_ = False
                HLT_SingleElectron_ = False
                HLT_METTriggers_ = False


                #######################################  Get Lepton Cut Results ##############################################################





                
                nMuons = 0
                nElectrons = 0
                #lepton_pt_array = it.eve.lepton_pt_
                #lepton_isMuon_array = it.eve.lepton_isMuon_
                #lepton_eta_array = it.eve.lepton_eta_

                ## ----------------------------- electrons/ muon cuts start
                
                get_indicesAndPt_of_passed_lepton = []
                get_indicesAndPt_of_passed_muons = []
                get_indicesAndPt_of_passed_electrons = []
                #print(it.muon_PFIso[0])
                #print(len(it.muon_PFIso))
#                sys.exit("Stopping test")

#                for item_lepton in range(len(it.muon_PFIso)):
                
#                print(it.leptons.size())
#                for item_lepton in it.leptons_fCoordinates_fPt:
 #                    print(item_lepton)
 #               sys.exit("Stopping test")
   #             for i in it.leptons_:
  #                  print(item_lepton)
    #                oneLepton = lv()
     #               oneLepton =it.leptons_[item_lepton]
      #              print(oneLepton.Pt())

                for item_lepton in range(it.leptons.size()):

######################### If Muon, let's check if it passes the chosen selections
######################### ==========================================================
                    print("it.muon_PFIso["+str(item_lepton)+"]")
                    print(it.muon_PFIso[item_lepton])
                    print("Is it.muon_PFIso[item_lepton] >= 76?")
                    print(ord(it.muon_PFIso[item_lepton]) >= 76)
                    if(it.muon_CutBasedID_Medium[item_lepton] == "Y" and it.leptons[item_lepton].Pt() > 15. and (abs(it.leptons[item_lepton].Eta()) < 2.4) and ord(it.muon_PFIso[item_lepton]) >= 76):
                              nMuons += 1
                              get_indicesAndPt_of_passed_muons.append((it.leptons[item_lepton].Pt(), item_lepton))
                              get_indicesAndPt_of_passed_lepton.append((it.leptons[item_lepton].Pt(), item_lepton, "M"))
                              print("Passed Muon with "+str(it.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(it.muon_CutBasedID_Medium[item_lepton]))

######################### If Electron, let's check if it passes the chosen selections
######################### ============================================================

                    elif(it.eleID_MVA_Iso_90[item_lepton] == 1 and it.leptons[item_lepton].Pt() > 15. and (abs(it.leptons[item_lepton].Eta()) < 2.4)):
                              nElectrons += 1
                              get_indicesAndPt_of_passed_lepton.append((it.leptons[item_lepton].Pt(), item_lepton, "E"))
                              get_indicesAndPt_of_passed_electrons.append((it.leptons[item_lepton].Pt(), item_lepton))
                              print("Passed Electron with "+str(it.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(it.eleID_MVA_Iso_90[item_lepton]))

                ## ----------------------------- electrons / muon cuts end
                
################### Get leading and subleading lepton index if available
################### =======================================================
                sortingLeptons = sorted(get_indicesAndPt_of_passed_lepton, reverse=True)
                if len(sortingLeptons) > 0: 
                    leadingLeptonIndex = sortingLeptons[0][1]
                    leadLeptonisMuon = sortingLeptons[0][2]
                    print("leadingLeptonIndex = "+str(leadingLeptonIndex))
                    if len(sortingLeptons) > 1:
                         subleadingLeptonIndex = sortingLeptons[1][1]
                         subleadingLeptonisMuon = sortingLeptons[1][2]
                         print("subleadingLeptonIndex = "+str(subleadingLeptonIndex))
#                else: continue
                    print(sortingLeptons)

#                    leadingLepton = tlv()
                    leadingLepton = lv()
                    if len(sortingLeptons) > 0:
                         leadingLepton = it.leptons[leadingLeptonIndex]
#                         leadingLepton.SetPtEtaPhiE(it.eve.lepton_pt_[leadingLeptonIndex], it.eve.lepton_eta_[leadingLeptonIndex], it.eve.lepton_phi_[leadingLeptonIndex], it.eve.lepton_e_[leadingLeptonIndex])

                    subleadingLepton = lv()
#                    subleadingLepton = tlv()
                    if len(sortingLeptons) > 1:
                         subleadingLepton = it.leptons[subleadingLeptonIndex]
#                         subleadingLepton.SetPtEtaPhiE(it.eve.lepton_pt_[subleadingLeptonIndex], it.eve.lepton_eta_[subleadingLeptonIndex], it.eve.lepton_phi_[subleadingLeptonIndex], it.eve.lepton_e_[subleadingLeptonIndex])
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

                    ## ----------------------------- category cuts start
                
                    ####### SL muon
                    if nMuons == 1 and nElectrons == 0 and leadingLepton.Pt() >= 26.  and it.lepPfIso[leadingLeptonIndex] < 0.15:
#                   if nMuons == 1 and nElectrons == 0 and leadLeptonisMuon and it.eve.lepton_pt_[leadingLeptonIndex] >= 26.  and it.eve.lepton_relIso_[leadingLeptonIndex] < 0.15:
                         passSLCuts_mu = True

                    ####### SL electron
                    elif nElectrons == 1 and nMuons == 0 and leadingLepton.Pt() >= 30.:
#                   elif nElectrons == 1 and nMuons == 0 and not leadLeptonisMuon and it.eve.lepton_pt_[leadingLeptonIndex] >= 30.:
                         passSLCuts_e = True

                    ####### DL muon
#                   elif nMuons == 2 and nElectrons == 0 
#                               and leadLeptonisMuon and it.eve.lepton_pt_[leadingLeptonIndex] >= 25. 
#                               and subleadingLeptonisMuon and it.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
                    elif nMuons == 2 and nElectrons == 0 and leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
                         checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)          
                         runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "mumu"))
                         if checkDileptonMassPass and (it.lepPdgId[leadingLeptonIndex] * it.lepPdgId[subleadingLeptonIndex] == -169) and runOnCorrectSamples:
                              passDLCuts_mumu = True

                    ####### DL ee
                    elif nMuons == 0 and nElectrons == 2 and leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
#it.eve.lepton_pt_[leadingLeptonIndex] >= 25. and it.eve.lepton_pt_[subleadingLeptonIndex] >= 15.: 
#                   elif nMuons == 0 and nElectrons == 2
 #                              and not leadLeptonisMuon and it.eve.lepton_pt_[leadingLeptonIndex] >= 25.
  #                             and not subleadingLeptonisMuon and it.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
                         checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)
                         runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "ee"))
                         if checkDileptonMassPass and (it.lepPdgId[leadingLeptonIndex] * it.lepPdgId[subleadingLeptonIndex]  == -121) and runOnCorrectSamples:
                              passDLCuts_ee = True
                    ####### DL emu
                    elif nMuons == 1 and nElectrons ==1:
                     if (leadLeptonisMuon == "E" and subleadingLeptonisMuon == "M") or (leadLeptonisMuon == "M" and subleadingLeptonisMuon == "E"):
                          if leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 15.:
#                          if it.eve.lepton_pt_[leadingLeptonIndex] >= 25. and it.eve.lepton_pt_[subleadingLeptonIndex] >= 15.:
                               checkDileptonMassPass = (dilepton.M() > 20.) #and (dilepton.M() < 76. or dilepton.M() > 106.)
                               runOnCorrectSamples = isMC or (not isMC and (channel == "" or channel == "emu"))
                               if checkDileptonMassPass and (it.lepPdgId[leadingLeptonIndex] * it.lepPdgId[subleadingLeptonIndex]  == -143) and runOnCorrectSamples:
                                    passDLCuts_emu = True
                ## ----------------------------- category cut end
                #get_indicesAndPt_of_passed_jets = []
                #for item_jet in range(len(it.eve.jet_pt_)):
                 #    jet_pt = it.eve.jet_pt_[item_jet]
                  #   jet_eta = it.eve.jet_eta_[item_jet]
                   #  jet_phi = it.eve.jet_phi_[item_jet]
                    # jet_puid = it.eve.jet_puid_[item_jet]
                     #if applyJetCuts(jet_pt, jet_eta, jet_puid):
                      #    get_indicesAndPt_of_passed_jets.append((jet_pt, item_jet))


#                dilepton = tlv()
 #               dilepton = leadingLepton + subleadingLepton
  #              if dilepton.M() > 20. and dilepton.M() :
   #                 dileptonMassCutpassed = True
    #            else: dileptonMassCutpassed = False



################### Check all Trigger results of interest 
################### ======================================================= 

                #### Single Electron triggers

                if(it.HLT_Ele27_WPTight_Gsf_v == 1): 
                         count_HLT_Ele27_WPTight_Gsf_v_ += 1 
                         HLT_Ele27_WPTight_Gsf_v_ = True 

                if(it.HLT_Ele32_WPTight_Gsf_L1DoubleEG_v == 1): 
                         count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_ += 1 
                         HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_ = True 

                if(it.HLT_Ele35_WPTight_Gsf_v == 1): 
                         count_HLT_Ele35_WPTight_Gsf_v_ += 1 
                         HLT_Ele35_WPTight_Gsf_v_ = True 

                if(it.HLT_Ele38_WPTight_Gsf_v == 1): 
                         count_HLT_Ele38_WPTight_Gsf_v_ += 1 
                         HLT_Ele38_WPTight_Gsf_v_ = True 

                if(it.HLT_Ele40_WPTight_Gsf_v == 1): 
                         count_HLT_Ele40_WPTight_Gsf_v_ += 1 
                         HLT_Ele40_WPTight_Gsf_v_ = True 

                if(it.HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v == 1): 
                         count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ += 1 
                         HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = True 

                if(it.HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v == 1): 
                         count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_ += 1 
                         HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_ = True 

                ######## Single Muon triggers
                if(it.HLT_IsoMu27_v == 1): 
                         count_HLT_IsoMu27_v_ += 1 
                         HLT_IsoMu27_v_ = True 

                if(it.HLT_IsoMu24_eta2p1_v == 1): 
                         count_HLT_IsoMu24_eta2p1_v_ += 1 
                         HLT_IsoMu24_eta2p1_v_ = True 

                ######## Double Muon Triggers
                if(it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v == 1): 
                         count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ += 1 
                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = True 

                if(it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v == 1): 
                         count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ += 1 
                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = True 

                if(it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v == 1): 
                         count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ += 1 
                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = True 

                ########## Double Electron Triggers
                if(it.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1): 
                         count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ += 1 
                         HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = True 

                if(it.HLT_DoubleEle33_CaloIdL_MW_v == 1): 
                         count_HLT_DoubleEle33_CaloIdL_MW_v_ += 1 
                         HLT_DoubleEle33_CaloIdL_MW_v_ = True 

                if(it.HLT_DoubleEle25_CaloIdL_MW_v == 1): 
                         count_HLT_DoubleEle25_CaloIdL_MW_v_ += 1 
                         HLT_DoubleEle25_CaloIdL_MW_v_ = True 

                if(it.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v == 1): 
                         count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_ += 1 
                         HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_ = True 

                ########## EMu Triggers
                if(it.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1): 
                         count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ += 1 
                         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = True 

                if(it.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v == 1): 
                         count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ += 1 
                         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = True 

                if(it.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v == 1): 
                         count_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ += 1 
                         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = True 

                if(it.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v == 1): 
                         count_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ += 1 
                         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = True 


                ########### MET Triggers
                if(it.HLT_PFMET120_PFMHT120_IDTight_v == 1): 
                         count_HLT_PFMET120_PFMHT120_IDTight_v_ += 1 
                         HLT_PFMET120_PFMHT120_IDTight_v_ = True 

                if(it.HLT_PFMET250_HBHECleaned_v == 1): 
                         count_HLT_PFMET250_HBHECleaned_v_ += 1 
                         HLT_PFMET250_HBHECleaned_v_ = True 

                if(it.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v == 1): 
                         count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_ += 1 
                         HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_ = True 

                if(it.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v == 1): 
                         count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_ += 1 
                         HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_ = True 

                if(it.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v == 1): 
                         count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_ += 1 
                         HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_ = True 

                 ##### Check if at least one Double Electron trigger fired
                if(it.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1 or it.HLT_DoubleEle33_CaloIdL_MW_v == 1 or it.HLT_DoubleEle25_CaloIdL_MW_v == 1 or it.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v == 1):
                         count_HLT_DoubleEl_ += 1 
                         HLT_DoubleEl_ = True 
                        

                ##### Check if at least one Double Muon trigger fired
                if(it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v == 1 or it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v == 1 or it.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v == 1):
                         count_HLT_DoubleMuon_ += 1 
                         HLT_DoubleMuon_ = True 

                ##### Check if at least one EMu trigger fired
                if(it.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1 or it.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v == 1 or 
                   it.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v == 1 or it.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v == 1):
                         count_HLT_EMu_ += 1 
                         HLT_EMu_ = True 
                ##### Check if at least one Single Muon trigger fired
                if(it.HLT_IsoMu27_v == 1 or it.HLT_IsoMu24_eta2p1_v == 1):
                         count_HLT_SingleMuon_ += 1 
                         HLT_SingleMuon_ = True 

                ##### Check if at least one Single Electron trigger fired
                if(it.HLT_Ele27_WPTight_Gsf_v == 1 or it.HLT_Ele32_WPTight_Gsf_L1DoubleEG_v == 1 or it.HLT_Ele35_WPTight_Gsf_v == 1
                   or it.HLT_Ele38_WPTight_Gsf_v == 1 or it.HLT_Ele40_WPTight_Gsf_v == 1 or it.HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v == 1 
                   or it.HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v == 1):
                #if(it.HLT_SingleElectron_ == 1): 
                         count_HLT_SingleElectron_ += 1 
                         HLT_SingleElectron_ = True 

                ##### Check if at least one MET trigger fired
                if(it.HLT_PFMET120_PFMHT120_IDTight_v == 1 or it.HLT_PFMET250_HBHECleaned_v == 1 or it.HLT_PFHT500_PFMET100_PFMHT100_IDTight_v == 1
                   or it.HLT_PFHT700_PFMET85_PFMHT85_IDTight_v == 1 or it.HLT_PFHT800_PFMET75_PFMHT75_IDTight_v == 1):
                         count_HLT_METTriggers_ += 1 
                         HLT_METTriggers_ = True 
            

                         
                if passDLCuts_ee and HLT_DoubleEl_:
                 h_leadingLepton_pt_MET_D.Fill(leadingLepton.Pt())
                 h_subleadingLepton_pt_MET_D.Fill(subleadingLepton.Pt())
            c = r.TCanvas('Lepton_pt', 'lepton pt')
            c.Divide(2)
            c.cd(1)
            h_leadingLepton_pt_MET_D.Draw()
            c.cd(2)
            h_subleadingLepton_pt_MET_D.Draw()
            for ext in ['png'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
            for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+nan+textname+'.'+ext)
            h_leadingLepton_pt_MET_D.Write()
            h_subleadingLepton_pt_MET_D.Write()

            inputFile.Close()
            if LeptonSFapplication:
                 inputElectronIDsfFile.Close()
                 inputElectronRECOsfFile.Close()
                 inputMuonIDsfFile.Close()
                 inputMuonISOsfFile.Close()

            if verbose:
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
        
