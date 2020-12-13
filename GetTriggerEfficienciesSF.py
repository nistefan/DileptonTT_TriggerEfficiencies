#!/bin/env python
__author__ = 'Nicole Stefanov'
__date__ = '2020-06-09'
         
#from array import array
#0;115;0cimport collections
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
import scipy.stats

import matplotlib.pyplot as plt          
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
from ROOT import TFile, TTree, TAxis, TH1F, TObject, Math
from ROOT import gROOT, AddressOf
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed, kGreen
from ROOT import TCanvas, TGraphErrors
from ROOT import gROOT

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

r.gStyle.SetOptStat(000000) 
verbose = False

# ============== Binning used ===================

#met_binning = [0., 20., 40., 60., 80., 100., 125., 150., 175., 200.]
met_binning = [20., 40., 60., 80., 100., 125., 150., 175., 200.]
#electron>_leading_Pt_binning = [20., 50., 80., 120., 200.]

#electron_subleading_Pt_binning = [20., 50., 90., 200.]
#lepton_Pt_binning = [ 20., 30., 40., 60., 80., 100., 200.]
#lepton_Pt_binning = [ 30., 40., 60., 80., 100., 200.]

lepton_Pt_binning = [20., 40., 60., 80., 100., 125., 150., 175., 200.]

#lepton_Pt_binning = [50., 80., 120., 200.]
#lepton_Pt_binning = [20., 30., 40., 60., 80., 100., 200.]

#lepton_2D_etaByEta_binning = [0.6, 1.4, 2.4]

lepton_2D_etaByEta_binning = [0.6, 1.4, 2.4] 
#[0.3, 0.6, 1.2, 1.7, 2.4]

#leadingLepton_eta_2D = [0, 0.4, 0.9, 1.5, 2.4]

#subleadingLepton_eta_2D = [0, 0.4, 0.9, 2.4]

#electron_eta_binning = [-2.4, -2.1, -1.566, -1.4442, -1.0, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 1.0, 1.4442, 1.566, 2.1, 2.4]

electron_eta_binning = [-2.1, -1.566, -1.4442, -1.0, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 1.0, 1.4442, 1.566, 2.1, 2.4] 

muon_eta_binning = [-2.1, -1.8, -1.5, -1.2, -0.9, -0.5, -0.2, 0.0, 0.2, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4]
#muon_eta_binning = [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.5, -0.2, 0.0, 0.2, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4]

#nJets_binning =  [1., 2., 3., 4., 5., 6., 7., 8.]
#nJets_binning =  [0., 1., 2., 3., 4., 5., 6., 7., 8.]
nJets_binning =  [1, 2, 3, 4, 5, 6, 7, 8, 9]



#electron_eta_binning = [-2.4, -2.1, -1.566, -1.4442, -1.0, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 1.0, 1.566, 2.1, 2.4]

#muon_eta_binning = [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.5, -0.2, 0.0, 0.2, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4]

#nJets_binning =  [0., 1., 2., 3., 4., 5., 6., 7., 8.]

nVertex_binning  = [float(i) for i in range(5,65,5)]

#nVertex_alternative_binning = [float(i) for i in range(0,26)]
nameRootoutputDirectory = "MuonMedium_ElectronMva90Iso_CoarseBinning40"
nameRootoutputCsVName = nameRootoutputDirectory+"/"+nameRootoutputDirectory+".txt"
nameRootoutput = "MuonMedium_ElectronMva90Iso_CoarseBinning40"

directory = nameRootoutputDirectory+"/"

if not os.path.exists(directory):
     os.makedirs(directory)


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
#def applyMuonCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_relIso):

 #   if (lepton_isMuon == 1) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (lepton_relIso < 0.25):
  #      return True
   # else: return False


#def applyElectronCuts(lepton_isMuon, lepton_pt, lepton_eta, lepton_isTight, lepton_scEta):
 #   if (lepton_isMuon == 0) and (lepton_pt > 15.) and (abs(lepton_eta) < 2.4) and (lepton_isTight == 1) and (abs(lepton_scEta) < 1.444 or abs(lepton_scEta) > 1.566):
  #      return True
   # else: return False


#def applyJetCuts(jet_pt, jet_eta, jet_puid):
 #    if (jet_pt > 20.) and abs(jet_eta) < 2.4 and jet_puid == 7:
  #        return True
   #  else: return False


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
#               print(str(inputSFarrayFileRootPath_GetHisto.GetXaxis().GetBinUpEdge(item_sfRangeX))+"|"+str(inputSFarrayFileRootPath_GetHisto.GetYaxis().GetBinLowEdge(item_sfRangeY))+" =sf= "+str(inputSFarrayFileRootPath_GetHisto.GetBinContent(item_sfRangeX, item_sfRangeY)))
          gatheredSFweights.append([inputSFarrayFileRootPath_GetHisto.GetXaxis().GetBinLowEdge(item_sfRangeX), getOneYBinRange])
          
#     print(gatheredSFweights)
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

def getSplit(l, n):
    stepsize = int(round(len(l)/float(n)))
    m =  []
    g = 0
   # print("stepsize")
    #print(stepsize)
    for iqw in range(stepsize, len(l)-stepsize, stepsize):
       # print(i)
       # print(i-stepsize)
        m.append(l[iqw-stepsize:iqw])
        #print(l[i-stepsize:i])                         
        g = iqw
    if len(l)%stepsize < 0.5*stepsize:
        m.append(l[g:])
    else:
        m.append(l[g:g+stepsize])
        m.append(l[g+stepsize:])
    #print(m)                                
    return m

def splitIntoBins(anyArray):
     numberOfBins = 1
     if len(anyArray) > 60 :
          numberOfBins = 6
     elif len(anyArray) > 10:
          numberOfBins = 3
     splitanyArray = np.array_split(anyArray, numberOfBins)
     return splitanyArray


def splitIntoBinsDependent(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value):
     outputArray = [np.array([], dtype=np.float32)]*len(alreadyExistingDataSplitIndexArray)
     
     for item_alreadyOrderedMCArray in alreadyOrderedMCArray:
          bin = 0
          for item_alreadyExistingDataSplitIndexArray in alreadyExistingDataSplitIndexArray:
               if alreadyExistingDataArray[item_alreadyExistingDataSplitIndexArray[-1]][which_value] > alreadyExistingMCArray[item_alreadyOrderedMCArray][which_value]:
           #         print("datapoint")
            #        print(alreadyExistingDataArray[item_alreadyExistingDataSplitIndexArray[-1]][which_value])
             #       print("Belongs to this bin")
              #      print(bin)
                    outputArray[bin] = np.append(outputArray[bin], item_alreadyOrderedMCArray)
                    break
               elif bin == len(outputArray)-1:
               #     print(bin == len(outputArray)-1)
                #    print("datapoint")
                 #   print(alreadyExistingDataArray[item_alreadyExistingDataSplitIndexArray[-1]][which_value])
                  #  print("Belongs to this bin")
                   # print(bin)
                    outputArray[bin] = np.append(outputArray[bin], item_alreadyOrderedMCArray)
                    break
               bin += 1
    # print(outputArray)
     return outputArray




def splitIntoBinsDependentPt(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value = 0):
     return splitIntoBinsDependent(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value)


def splitIntoBinsDependentEta(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value = 1):
     return splitIntoBinsDependent(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value)


def splitIntoBinsDependentIso(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value = 2):
     return splitIntoBinsDependent(alreadyOrderedMCArray, alreadyExistingDataSplitIndexArray, alreadyExistingDataArray, alreadyExistingMCArray, which_value)

def getoutValues(splitIndexArray, dataOrMCValueArray, which_value):
     getoutValuesArray = []
     for i_splitIndexArray in splitIndexArray:
          getoutValuesArray.append(dataOrMCValueArray[i_splitIndexArray][which_value])
     return getoutValuesArray


def getoutIndexArrayForSplit(splitIndexArray, dataOrMCValueArray, which_value, valuesAccordingToSplit):
     getoutIndexToSplitArray = []
     getoutValuesArray = []
     for i_splitIndexArray in splitIndexArray:
          getoutValuesArray.append(dataOrMCValueArray[i_splitIndexArray][which_value])
     for i_valuesAccordingToSplit in valuesAccordingToSplit:
          try:
               oneIndexFound = next(x[0] for x in enumerate(getoutValuesArray) if x[1] > i_valuesAccordingToSplit)
       #        print("i_valuesAccordingToSplit")
        #       print(i_valuesAccordingToSplit)
         #      print("oneIndexFound")
          #     print(oneIndexFound)
          except:
               oneIndexFound = -1
          if oneIndexFound != -1:
                    getoutIndexToSplitArray.append(oneIndexFound)
     
     
     #f len(getoutIndexToSplitArray) > 0:
#     print(getoutValuesArray)
 #    print(getoutIndexToSplitArray)
     return getoutIndexToSplitArray

#          [next(x[0] for x in enumerate(eval(what)) if x[1] > 20.), next(x[0] for x in enumerate(eval(what)) if x[1] > 50.), next(x[0] for x in enumerate(eval(what)) if x[1] > 80.), next(x[0] for x in enumerate(eval(what)) if x[1] > 120.), next(x[0] for x in enumerate(eval(what)) if x[1] > 200.)]
def getoutIndexArrayForSplitABS(splitIndexArray, dataOrMCValueArray, which_value, valuesAccordingToSplit):
     getoutIndexToSplitArray = []
     getoutValuesArray = []
     for i_splitIndexArray in splitIndexArray:
          getoutValuesArray.append(dataOrMCValueArray[i_splitIndexArray][which_value])
     for i_valuesAccordingToSplit in valuesAccordingToSplit:
          try:
               oneIndexFound = next(x[0] for x in enumerate(getoutValuesArray) if abs(x[1]) > i_valuesAccordingToSplit)
          except:
               oneIndexFound = -1
          if oneIndexFound != -1:
               getoutIndexToSplitArray.append(oneIndexFound)
     #f len(getoutIndexToSplitArray) > 0:                                                                                                                                                                                                                                                                                                           
     return getoutIndexToSplitArray


def getoutIndexArrayFor1DSplit(arrayToBeSplit, array_binning, whichIndex):

     getoutSplitIndex = []

     for item_array_binning in array_binning:
          #print(item_array_binning)
          nn = [indices_ordered[0] for indices_ordered in sorted(enumerate(arrayToBeSplit), key=lambda elem: elem[1][whichIndex]) 
                if arrayToBeSplit[indices_ordered[0]][whichIndex] < item_array_binning and not any(indices_ordered[0] in sublist for sublist in getoutSplitIndex)]
          #print(nn)
          getoutSplitIndex.append(nn)

     nn = [indices_ordered[0] for indices_ordered in sorted(enumerate(arrayToBeSplit), key=lambda elem: elem[1][whichIndex])
                if arrayToBeSplit[indices_ordered[0]][whichIndex] >= array_binning[-1] and not any(indices_ordered[0] in sublist for sublist in getoutSplitIndex)]
     getoutSplitIndex.append(nn)
   #  print(getoutSplitIndex)
     return getoutSplitIndex



def getoutIndexArrayFor2DSplitEtaABS(lepton1_array, lepton2_array, lepton1_binning, lepton2_binning):
     if len(lepton1_array) != len(lepton2_array):
          sys.exit("leading lepton array and subleading lepton array must be of equal length, but they are not.\n Please check what went wrong.")
     
     getoutSplitIndex2D = []

     for item_lepton1_binning in lepton1_binning:
       #   print("item_lepton1_binning: "+str(item_lepton1_binning))
          for item_lepton2_binning in lepton2_binning:
        #     print("item_lepton2_binning: "+str(item_lepton2_binning))
             nn = [indices_ordered[0] for indices_ordered in sorted(enumerate(lepton1_array), key=lambda elem: elem[1][1])
                   if (abs(lepton1_array[indices_ordered[0]][1]) < item_lepton1_binning and abs(lepton2_array[indices_ordered[0]][1]) < item_lepton2_binning 
                       and not any(indices_ordered[0] in sublist for sublist in getoutSplitIndex2D))]
         #    print(nn)
             getoutSplitIndex2D.append(nn)


     #print(getoutSplitIndex2D)
     return getoutSplitIndex2D

def getoutIndexArrayFor2DSplitEtaABSupperSystematics(lepton1_array, lepton2_array, lepton1_binning, lepton2_binning, njet_vertex_metArray, whichIndexInNjetVertexMetArray, cutValueSystematics):
     if len(lepton1_array) != len(lepton2_array):
          sys.exit("leading lepton array and subleading lepton array must be of equal length, but they are not.\n Please check what went wrong.")
     getoutSplitIndex2D = []
     for item_lepton1_binning in lepton1_binning:
          #print("item_lepton1_binning: "+str(item_lepton1_binning))
          for item_lepton2_binning in lepton2_binning:
          #   print("item_lepton2_binning: "+str(item_lepton2_binning))
             nn = [indices_ordered[0] for indices_ordered in sorted(enumerate(lepton1_array), key=lambda elem: elem[1][1])
                   if (abs(lepton1_array[indices_ordered[0]][1]) < item_lepton1_binning and abs(lepton2_array[indices_ordered[0]][1]) < item_lepton2_binning
                       and njet_vertex_metArray[indices_ordered[0]][whichIndexInNjetVertexMetArray] > cutValueSystematics
                       and not any(indices_ordered[0] in sublist for sublist in getoutSplitIndex2D))]
            # print(nn)
             getoutSplitIndex2D.append(nn)
     #print(getoutSplitIndex2D)                                                                                                                                                              
     return getoutSplitIndex2D

def getoutIndexArrayFor2DSplitEtaABlowerSystematics(lepton1_array, lepton2_array, lepton1_binning, lepton2_binning, njet_vertex_metArray, whichIndexInNjetVertexMetArray, cutValueSystematics):
     if len(lepton1_array) != len(lepton2_array):
          sys.exit("leading lepton array and subleading lepton array must be of equal length, but they are not.\n Please check what went wrong.")
     getoutSplitIndex2D = []
     for item_lepton1_binning in lepton1_binning:
          #print("item_lepton1_binning: "+str(item_lepton1_binning))                                                                                                                                                                                                           
          for item_lepton2_binning in lepton2_binning:
          #   print("item_lepton2_binning: "+str(item_lepton2_binning))                                                                                                                                                                                                        
             nn = [indices_ordered[0] for indices_ordered in sorted(enumerate(lepton1_array), key=lambda elem: elem[1][1])
                   if (abs(lepton1_array[indices_ordered[0]][1]) < item_lepton1_binning and abs(lepton2_array[indices_ordered[0]][1]) < item_lepton2_binning
                       and njet_vertex_metArray[indices_ordered[0]][whichIndexInNjetVertexMetArray] <= cutValueSystematics
                       and not any(indices_ordered[0] in sublist for sublist in getoutSplitIndex2D))]
             # print(nn)  
             getoutSplitIndex2D.append(nn)
     #print(getoutSplitIndex2D)
     return getoutSplitIndex2D




def getSplitOfIndexArray(toBesplitIndexArray, indexAccordingToSplit):
     getoutSplitArray = []
     #print("toBesplitIndexArray")
     #print(toBesplitIndexArray)
    # print("indexAccordingToSplit")
    # print(indexAccordingToSplit)
     for i in range(0, len(indexAccordingToSplit)):
        
          if i == 0:
#             print("Starting")
 #            print([toBesplitIndexArray[:indexAccordingToSplit[i]]])
  #           print(toBesplitIndexArray[indexAccordingToSplit[i]:indexAccordingToSplit[i+1]])
             if indexAccordingToSplit[i] != 0:
                  getoutSplitArray.append(toBesplitIndexArray[:indexAccordingToSplit[i]])
             getoutSplitArray.append(toBesplitIndexArray[indexAccordingToSplit[i]:indexAccordingToSplit[i+1]])
          elif i == len(indexAccordingToSplit)-1:
   #          print("Ending")
    #         print([toBesplitIndexArray[indexAccordingToSplit[i]:]])
             getoutSplitArray.append(toBesplitIndexArray[indexAccordingToSplit[i]:])
          else:   
             getoutSplitArray.append(toBesplitIndexArray[indexAccordingToSplit[i]:indexAccordingToSplit[i+1]])
     #        print(toBesplitIndexArray[indexAccordingToSplit[i]:indexAccordingToSplit[i+1]])
     #print("getoutSplitArray")
     #print(getoutSplitArray)
     return getoutSplitArray



def validateRange(splitIndexArray,  dataOrMCValueArray, which_value):
     getValues = []
     for item_splitIndexArray in splitIndexArray:
          one_split_piece_of_dataOrMCValueArray = []
          for value_in_item_splitIndexArray in item_splitIndexArray:
               #print(value_in_item_splitIndexArray)
               one_split_piece_of_dataOrMCValueArray.append(dataOrMCValueArray[int(value_in_item_splitIndexArray)][which_value])
          getValues.append(one_split_piece_of_dataOrMCValueArray)
     return getValues

def binArrayGetOutTheValuesFromSplitIndices(splitIndexArray, dataOrMCValueArray, which_value):
     getBinnedArrayFromDataOrMCValueArray = []
     for item_splitIndexArray in splitIndexArray:
#          print(item_splitIndexArray)
#          item_splitIndexArray = item_splitIndexArray.astype(int)
          one_split_piece_of_dataOrMCValueArray = []
          for value_in_item_splitIndexArray in item_splitIndexArray:
#              print(dataOrMCValueArray[value_in_item_splitIndexArray][which_value])
#              print("value_in_item_splitIndexArray")
#              print(value_in_item_splitIndexArray)
              one_split_piece_of_dataOrMCValueArray.append(dataOrMCValueArray[value_in_item_splitIndexArray][which_value])
          getBinnedArrayFromDataOrMCValueArray.append(one_split_piece_of_dataOrMCValueArray)
#     print(getBinnedArrayFromDataOrMCValueArray)
     return getBinnedArrayFromDataOrMCValueArray



def errorSqrtSumOfsquaredWeights(someArrayWithNumbers):
     error = 0.
     for item_someArrayWithNumbers in someArrayWithNumbers:
          error += (item_someArrayWithNumbers**2)
     return np.sqrt(error)
     
def binsErrorSqrtSumOfSquaredWeights(someArrayOfBinArrays):
     errorBinArray = []
     for item_someArrayOfBinArrays in someArrayOfBinArrays:
          error = 0.
          error += np.sqrt(sum([iqq*iqq for iqq in item_someArrayOfBinArrays]))
          errorBinArray.append(error)
     return errorBinArray

def binsErrorSqrtLen(someArrayOfBinArrays):
     errorBinArray = []
     for item_someArrayOfBinArrays in someArrayOfBinArrays:
          error = 0.
          error += np.sqrt(len(item_someArrayOfBinArrays))
          errorBinArray.append(error)
     return errorBinArray

def errorSqrtLen(someArrayWithNumbers):
     error = 0.
     for item_someArrayWithNumbers in someArrayWithNumbers:
          error += len(someArrayWithNumbers)
     return np.sqrt(error)

def getScaleFactorAndUncertainty( efficiencyMC, uncertaintyMC, efficiencyData, uncertaintyData):
     if efficiencyMC > 0.:
          scaleFactor = float(efficiencyData)/efficiencyMC

          uncertaintyScaleFactor = np.sqrt( (uncertaintyMC/float(efficiencyMC))**2 + (uncertaintyData/float(efficiencyData))**2) * scaleFactor
     else:
          scaleFactor = uncertaintyScaleFactor = 0.

     return scaleFactor, uncertaintyScaleFactor


def getScaleFactorAndUncertaintyFromTupel( efficiencyAndUncertaintyTupelMC, efficiencyAndUncertaintyTupelData):
     if efficiencyAndUncertaintyTupelMC[0] > 0.:
          scaleFactor = float(efficiencyAndUncertaintyTupelData[0])/efficiencyAndUncertaintyTupelMC[0]
          if efficiencyAndUncertaintyTupelData[0] > 0.:
               uncertaintyScaleFactor = np.sqrt( (efficiencyAndUncertaintyTupelMC[1]/float(efficiencyAndUncertaintyTupelMC[0]))**2 + (efficiencyAndUncertaintyTupelData[1]/float(efficiencyAndUncertaintyTupelData[0]))**2) * scaleFactor
          else: 
               uncertaintyScaleFactor = 0.

     else:
          scaleFactor = uncertaintyScaleFactor = 0.

     return scaleFactor, uncertaintyScaleFactor


def getScaleFactorAndUncertaintyPerBin(efficiencyAndUncertaintyTupelMC, efficiencyAndUncertaintyTupelData):
     getScaleFactorAndUncertaintyForEachBin_Array = []
     if len(efficiencyAndUncertaintyTupelData) != len(efficiencyAndUncertaintyTupelMC):
      #    print(efficiencyAndUncertaintyTupelMC)
       #   print(efficiencyAndUncertaintyTupelData)
          sys.exit("ERROR! len(efficiencyData) != len(efficiencyMC):in getScaleFactorAndUncertaintyPerBin( efficiencyMC, uncertaintyMC, efficiencyData, uncertaintyData). \n Something went heavily wrong. Please check.")
     for one_bin in range(0, len(efficiencyAndUncertaintyTupelData)):
          getScaleFactorAndUncertaintyForEachBin_Array.append(getScaleFactorAndUncertaintyFromTupel(efficiencyAndUncertaintyTupelMC[one_bin], efficiencyAndUncertaintyTupelData[one_bin]))
          
     return getScaleFactorAndUncertaintyForEachBin_Array


def getScaleFactorSystematicsPerBin(scaleFactorNominal, arrayOfScaleFactorsUsedForSystemtics):
     numberOfBinsInNominal = len(scaleFactorNominal)
     getSystematicDiffsPerBin = []
     for i in arrayOfScaleFactorsUsedForSystemtics:
          if len(i) != numberOfBinsInNominal:
               sys.exit("ERROR! len(i) != numberOfBinsInNominal in getScaleFactorSystematics(scaleFactorNominal, arrayOfScaleFactorsUsedForSystemtics). \n Something went heavily wrong. Please check.")
     for oneBinIndex in range(0, numberOfBinsInNominal):
          scaleFactorNominalInOneBinIndex = scaleFactorNominal[oneBinIndex][0]
          getAllDiffs = []
          for oneArray in arrayOfScaleFactorsUsedForSystemtics:
               if oneArray[oneBinIndex][0] != 0. and scaleFactorNominalInOneBinIndex != 0.:
                    diff = abs(oneArray[oneBinIndex][0] - scaleFactorNominalInOneBinIndex)
               else:
                    diff = 0.
               getAllDiffs.append(diff)
          getSystematicDiffsPerBin.append(max(getAllDiffs))

     return getSystematicDiffsPerBin

     
def getScaleFactorSystematics(scaleFactorNominal, arrayOfScaleFactorsUsedForSystemtics):
     numberOfBinsInNominal = len(scaleFactorNominal)
     getSystematicDiffsPerBin = []
     for i in arrayOfScaleFactorsUsedForSystemtics:
          if len(i) != numberOfBinsInNominal:
               sys.exit("ERROR! len(i) != numberOfBinsInNominal in getScaleFactorSystematics(scaleFactorNominal, arrayOfScaleFactorsUsedForSystemtics). \n Something went heavily wrong. Please check.")
     
     scaleFactorNominalInOneBinIndex = scaleFactorNominal[0]
     getAllDiffs = []
     for oneArray in arrayOfScaleFactorsUsedForSystemtics:
          if oneArray[0] != 0. and scaleFactorNominalInOneBinIndex != 0.:
                    diff = abs(oneArray[0] - scaleFactorNominalInOneBinIndex)
          else:
                    diff = 0.
     getAllDiffs.append(diff)

     return max(getAllDiffs)

def getScaleFactorSystematicAndStatisticalUncertaintyPerBin(scaleFactorNominal, arrayScaleFactorSystematicUncertainty, alpha):
     if len(arrayScaleFactorSystematicUncertainty) != len(scaleFactorNominal):
          sys.exit("ERROR! len(arrayScaleFactorSystematicUncertainty) != len(scaleFactorNominal) in getScaleFactorSystematicAndStatisticalUncertaintyPerBin(scaleFactorNominal, arrayScaleFactorSystematicUncertainty, alphaUncertainty). \n Something went heavily wrong. Please check.")
     alphaUncertainty = abs(1.-alpha)
     scaleFactorNominalWithStatisticalAndSystematicalUncertainties = []
     for i in range(0, len(scaleFactorNominal)):
          combinedUncertainty = (scaleFactorNominal[i][1])**2 + (arrayScaleFactorSystematicUncertainty[i])**2 + alphaUncertainty**2
          scaleFactorNominalWithStatisticalAndSystematicalUncertainties.insert(i, (scaleFactorNominal[i][0], np.sqrt(combinedUncertainty)))
     return scaleFactorNominalWithStatisticalAndSystematicalUncertainties

def getScaleFactorSystematicAndStatisticalUncertainty(scaleFactorNominal, scaleFactorSystematicUncertainty, alpha):
     alphaUncertainty = abs(1.-alpha)
     combinedUncertainty = (scaleFactorNominal[1])**2 + (scaleFactorSystematicUncertainty)**2 + alphaUncertainty**2
     scaleFactorNominalWithStatisticalAndSystematicalUncertainties = (scaleFactorNominal[0], np.sqrt(combinedUncertainty))
     return scaleFactorNominalWithStatisticalAndSystematicalUncertainties

          
# Taken from https://gist.github.com/DavidWalz/8538435
def clopper_pearson(k,n,alpha=0.32):
    """                                                                                                          
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval                                         
    alpha confidence intervals for a binomial distribution of k expected successes on n trials                   
    Clopper Pearson intervals are a conservative estimate.                                                       
    """
    lo = scipy.stats.beta.ppf(alpha/2, k, n-k+1)
    hi = scipy.stats.beta.ppf(1 - alpha/2, k+1, n-k)
    return lo, hi
#------------------------------------------------------------

def getEffAndErrorByClopperPearson(counts, shots):
     if shots > 0:
#          print("counts "+str(counts))
 #         print("shots "+str(shots))
          eff = float(counts)/shots
          effDown, effUp = clopper_pearson(int(counts), int(shots))
     
          efferrDown = eff - effDown 
          efferrUp = effUp -eff
          
          effErr = max(efferrDown, efferrUp)
          
     else: eff = effErr = 0 
     return eff, effErr


def binArrayGetOutEffValuesFromSplitIndices(splitIndexArray, dataOrMCValueArray, dataOrMCNjetVertexMetArray, which_value):
     getBinnedArrayFromDataOrMCValueArray = []
#     getBinnedArrayFromDataOrMCNjetVertexMetArray = []
     get_HLT_DoubleMuon_fired_Weights_perBin = []
     getAllCounts_HLT_DoubleMuon_fired_perBin = []
     getAdditional_HLT_SingleMuon_fired_Weights_perBin = []
     getAdditionalCounts_HLT_SingleMuon_fired_perBin = []
     getlepton_weight_perBin = []
     getAllShots_perBin = []
     getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin = []
     getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = []
     getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin = []
     getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin = []
     c_all = 0.
     d_all = 0.
     d_all_OR = 0.
     e_all = 0.
     d_all_weights = 0.
     e_all_weights = 0.
     all_weights = 0.
     for item_splitIndexArray in splitIndexArray:
          one_split_piece_all_lepton_weight = []
          one_split_piece_of_dataOrMCValueArray = []
          #one_split_piece_of_dataOrMCNjetVertexMetArray = []
         # one_split_piece_HLT_DoubleMuon_fired_getlepton_weight = []
         # one_split_piece_HLT_SingleMuon_fired_getlepton_weight = []
       #   one_split_piece_HLT_DoubleMuon_fired = []
        #  one_split_piece_HLT_SingleMuon_fired = []
         # one_split_piece_item_splitIndexArray = []
          c = len(item_splitIndexArray)
          c_all += c
          d = 0.
          e = 0.
          e_weights = 0.
          d_weights = 0.
          lepton_weights = 0.
          for value_in_item_splitIndexArray in item_splitIndexArray:
           #   one_split_piece_of_dataOrMCNjetVertexMetArray.append(dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray])
              one_split_piece_of_dataOrMCValueArray.append(dataOrMCValueArray[value_in_item_splitIndexArray][which_value])
              if dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][3] == True:
                   d += 1.
                   d_weights += dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1]
                   #one_split_piece_HLT_DoubleMuon_fired_getlepton_weight.append(dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1])
              elif dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][3] == False and dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][4] == True:
                   e +=1.
                   e_weights += dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1]
                   #one_split_piece_HLT_SingleMuon_fired_getlepton_weight.append(dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1])
              one_split_piece_all_lepton_weight.append(dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1])
              lepton_weights += dataOrMCNjetVertexMetArray[value_in_item_splitIndexArray][-1]
          countHLT_DoubleMuon_fired_perBin = d
          countHLT_DoubleMuonfired_OR_HLT_SingleMuon_fired_perBin = d + e
          d_all += d
          d_all_OR += countHLT_DoubleMuonfired_OR_HLT_SingleMuon_fired_perBin
          e_all += e
          sum_HLT_DoubleMuon_fired_lepton_weight_perBin = d_weights
          sum_HLT_DoubleMuon_fired_OR_HLT_SingleMuon_fired_lepton_weight_perBin = d_weights + e_weights
          get_HLT_DoubleMuon_fired_Weights_perBin.append(sum_HLT_DoubleMuon_fired_lepton_weight_perBin)
          d_all_weights += sum_HLT_DoubleMuon_fired_lepton_weight_perBin
          getAllCounts_HLT_DoubleMuon_fired_perBin.append(countHLT_DoubleMuon_fired_perBin)
          getAdditional_HLT_SingleMuon_fired_Weights_perBin.append(e_weights)
          e_all_weights += e_weights
          getAdditionalCounts_HLT_SingleMuon_fired_perBin.append(e)
          sum_Of_All_lepton_weights_in_one_split_piece = lepton_weights
          getlepton_weight_perBin.append(sum_Of_All_lepton_weights_in_one_split_piece)
          all_weights += sum_Of_All_lepton_weights_in_one_split_piece
          getAllShots_perBin.append(c)
          getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin.append(getEffAndErrorByClopperPearson(countHLT_DoubleMuon_fired_perBin, c))
          getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin.append(getEffAndErrorByClopperPearson(countHLT_DoubleMuonfired_OR_HLT_SingleMuon_fired_perBin, c))
          getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin.append(getEffAndErrorByClopperPearson(sum_HLT_DoubleMuon_fired_lepton_weight_perBin, sum_Of_All_lepton_weights_in_one_split_piece))
          getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin.append(getEffAndErrorByClopperPearson(sum_HLT_DoubleMuon_fired_OR_HLT_SingleMuon_fired_lepton_weight_perBin, sum_Of_All_lepton_weights_in_one_split_piece))
          #getBinnedArrayFromDataOrMCNjetVertexMetArray.append(one_split_piece_of_dataOrMCNjetVertexMetArray)
          getBinnedArrayFromDataOrMCValueArray.append(one_split_piece_of_dataOrMCValueArray)
#     print(getBinnedArrayFromDataOrMCValueArray)
#     print("sum(getAllCounts_HLT_DoubleMuon_fired_perBin)")
#     print(sum(getAllCounts_HLT_DoubleMuon_fired_perBin))
#     print("c_all")
#     print(c_all)
     getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL = getEffAndErrorByClopperPearson(d_all, c_all)
#     print("d_all_OR: "+str(d_all_OR))
#     print("d_all_weights: "+str(d_all_weights))
#     print("(d_all_weights+e_all_weights) : "+str((d_all_weights+e_all_weights)))
#     print("all_weights: "+str(all_weights))
     getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL = getEffAndErrorByClopperPearson(d_all_OR, c_all)
     getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted = getEffAndErrorByClopperPearson(d_all_weights, all_weights)
     getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted = getEffAndErrorByClopperPearson((d_all_weights+e_all_weights), all_weights)
#     print("sum(getAllCounts_HLT_DoubleMuon_fired_perBin)")
#     print(sum(getAllCounts_HLT_DoubleMuon_fired_perBin))
#     print("sum(getAllShots_perBin)")
#     print(sum(getAllShots_perBin))
     
     return getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL, getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted,  getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL, getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted,  getAllShots_perBin, getAllCounts_HLT_DoubleMuon_fired_perBin, getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin, getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin, getAdditionalCounts_HLT_SingleMuon_fired_perBin, getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin, getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin, getlepton_weight_perBin, get_HLT_DoubleMuon_fired_Weights_perBin, getAdditional_HLT_SingleMuon_fired_Weights_perBin, getBinnedArrayFromDataOrMCValueArray #, getBinnedArrayFromDataOrMCNjetVertexMetArray



def tripleGraphPlot(directroy, binning_array, plotName, graph1, graph2, graph3):

     getDiffIntervals = np.array([float(abs(binning_array[i+1] - binning_array[i])/2) for i in range(0, len(binning_array)-1)])

     print(getDiffIntervals)
     print(len(getDiffIntervals))
     c = r.TCanvas(str(plotName), ' ')
     n = len(getDiffIntervals)+2
     x  = np.array(binning_array)

     ex = getDiffIntervals #np.array( [  0.05,  0.1, 0.07, 0.07, 0.04, 0.05, 0.06, 0.07, 0.08, 0.05 ] )                                                                                                                                                  
     y = np.array([i[0] for i in graph1])
     ey = np.array([i[1] for i in graph1])
#     print(y)
#     print(ey)
#     print(x)
#     print(ex)

     y1 = np.array([i[0] for i in graph2])
     ey1 = np.array([i[1] for i in graph2])

     z1 = np.array([i[0] for i in graph3])
     ez1 = np.array([i[1] for i in graph3])
     #y  = np.array(  [     1,  2.9,  5.6,  7.4,  9.0,  9.6,  8.7,  6.3,  4.5,    1 ] )
     #ey = np.array(  [  0.8,  0.7,  0.6,  0.5,  0.4,  0.4,  0.5,  0.6,  0.7,  0.8  ] )
     gr = TGraphErrors( n, x, y, ex, ey )
     gr.SetTitle( 'TGraphErrors Example' )
     #gr.SetMarkerColor( 4 )
     #gr.SetMarkerStyle( 21 )
     axis = gr.GetYaxis()
     #axis.SetLimits(0., 1.2)
     #gr.GetHistogram().SetMinimum(0.)
     #gr.GetHistogram().SetMaximium(1.2)
     gr.Draw('AP' )
     gr.GetYaxis().SetRangeUser(0., 1.2)
#     gr.GetXaxis().SetRangeUser(-2.7, 2.7)
     #gr.Draw('AP' )
     #c.Update()
     gr1 = TGraphErrors( n, x, y1, ex, ey1 )
     gr1.SetLineColor(2)
     gr1.SetMarkerColor( 2 )
     gr1.Draw('Psame' )
     #c.Update()                                                                                                                                                                                                                                         
     gr2 = TGraphErrors( n, x, z1, ex, ez1 )
     gr2.SetLineColor(4)
     gr2.SetMarkerColor( 4)
     gr2.Draw('Psame' )
     for ext in ['png'] : c.SaveAs(directory+str(plotName)+'.'+ext)
     for ext in ['pdf'] : c.SaveAs(directory+str(plotName)+'.'+ext)


RootfilesToRunOver =[
     "/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017B.root",
     "/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017C.root",
     "/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017D.root",
     "/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017E.root",
     "/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017F.root"
#"/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/met_run2017D.root",
#                     "/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/ttbarsignalplustau_fromDilepton.root"
#/nfs/dust/cms/user/stefanon/analyses2017/ckoraka_TriggerSF/Converted_Ntuples_For_TiggerSF_Script/met_run2017D.root",
#"/nfs/dust/cms/group/topcmsdesy/ntuple13tev/2017_09_05_TAG_V037/ttbarsignalplustau.root",
]

MC_SampleName = ["/pnfs/desy.de/cms/tier2/store/user/nstefano/22_10_2020_V07/ttbarsignalplustau_fromDilepton.root"
]

RootfilesToRunOver += MC_SampleName


#    ______                                                       
#   /\__  _\       __                                             
#   \/_/\ \/ _ __ /\_\     __      __      __   _ __              
#      \ \ \/\`'__\/\ \  /'_ `\  /'_ `\  /'__`\/\`'__\            
#       \ \ \ \ \/ \ \ \/\ \L\ \/\ \L\ \/\  __/\ \ \/             
#        \ \_\ \_\  \ \_\ \____ \ \____ \ \____\\ \_\             
#         \/_/\/_/   \/_/\/___L\ \/___L\ \/____/ \/_/             
#                          /\____/ /\____/                        
#                          \_/__/  \_/__/                         
#      ____________________________________________________       
#       ____________________________________________________      

# made with https://patorjk.com/software/taag/#p=display&f=Larry%203D&t=Trigger
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
pileupReweighting = True
if pileupReweighting:
     dataPileupFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/pileup/pileupDATA__Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1__max120_n120__MinBias69200.root"
     inputdataPileupFile =  r.TFile.Open(dataPileupFile, "read")
     inputdataPileupHisto = inputdataPileupFile.Get("pileup")
     integralPileupHisto = inputdataPileupHisto.Integral()
#     print("integralPileupHisto")
#     print(integralPileupHisto)
#     for item_binrange in range(0, 20):
 #         print("inputdataPileupHisto.GetXaxis().FindBin("+str(item_binrange))
  #        print(inputdataPileupHisto.GetXaxis().FindBin(item_binrange))
   #       print(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(item_binrange)))
    # binContentArrayOfdataPileupHisto = []
     weightedBinContentArrayOfdataPileupHisto = []
    
     for itemBinContent in range(0, inputdataPileupHisto.GetNbinsX()+1):
          weightedBinContentArrayOfdataPileupHisto.append(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(itemBinContent))/integralPileupHisto)
     #     binContentArrayOfdataPileupHisto.append(inputdataPileupHisto.GetBinContent(inputdataPileupHisto.GetXaxis().FindBin(itemBinContent)))
#     print("weightedBinContentArrayOfdataPileupHisto")
#     print(weightedBinContentArrayOfdataPileupHisto)
#     print(len(weightedBinContentArrayOfdataPileupHisto))
#     print("binContentArrayOfdataPileupHisto")
#     print(binContentArrayOfdataPileupHisto)
     inputdataPileupFile.Close()
     MCPileupFile = "/pnfs/desy.de/cms/tier2/store/user/nstefano/02_05_2020_V06/ttbarsignalplustau_fromDilepton.root"
     inputMCPileupFile = r.TFile.Open(MCPileupFile, "read")
     inputMCPilepuFileHisto = inputMCPileupFile.Get("PileupMCTemplateMaker/MC_TrueNIntBX0")
     integralMCPilepuFileHisto = inputMCPilepuFileHisto.Integral()
#     print("integralMCPilepuFileHisto")
#     print(integralMCPilepuFileHisto)
#     binContentArrayOfMCPileupHisto = []
     weightedBinContentArrayOfMCPileupHisto = []
#     for item_binrange in range(0, 20):
#          print("inputMCPilepuFileHisto.GetXaxis().FindBin("+str(item_binrange))
#          print(inputMCPilepuFileHisto.GetXaxis().FindBin(item_binrange))
#          print(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(item_binrange)))
     for itemBinContent in range(0, inputMCPilepuFileHisto.GetNbinsX()+1):
          weightedBinContentArrayOfMCPileupHisto.append(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(itemBinContent))/integralMCPilepuFileHisto)
#          binContentArrayOfMCPileupHisto.append(inputMCPilepuFileHisto.GetBinContent(inputMCPilepuFileHisto.GetXaxis().FindBin(itemBinContent)))
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

#     print(pileupWeights[15])
#     print(pileupWeights[18])

#sys.exit()
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
#lep_Eta = 1.9238435389
#lep_Pt = 114.65415955
#exampleSf = leptonSF(gatheredElectronIDsfWeights, lep_Eta, lep_Pt)
#print("exampleSf")
#print(exampleSf)

#lep1_Eta = 1.4279511513
#lep1_Pt = 111.87143707
#exampleSfsecond = leptonSF(gatheredElectronIDsfWeights, lep1_Eta, lep1_Pt)
#print("exampleSfsecond")
#print(exampleSfsecond)


ElectronRECOsfFile = "/nfs/dust/cms/user/stefanon/analyses2017/fixes/CMSSW_9_4_17/src/TopAnalysis/Configuration/analysis/common/data/electron/94X/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"
gatheredElectronRECOsfWeights = create_sfArray(ElectronRECOsfFile, inputElectronIDhistogram)


###### Example Test1
#exampleSfReco = leptonSF(gatheredElectronRECOsfWeights, lep_Eta, lep_Pt)
#print("exampleSfReco")
#print(exampleSfReco)

#exampleSfRecoSecond = leptonSF(gatheredElectronRECOsfWeights, lep1_Eta, lep1_Pt)
#print("exampleSfRecoSecond")
#print(exampleSfRecoSecond)

#print("First lepton weight")
#print(exampleSf*exampleSfReco)
#print("Second lepton weight")
#print(exampleSfsecond*exampleSfRecoSecond)

#print("Final lepton weight")
#print(exampleSf*exampleSfReco*exampleSfsecond*exampleSfRecoSecond)

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



#sys.exit()



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

#['1st_Electron_Pt', '1st_Muon_Pt', '2nd_Electron_Pt', '2nd_Muon_Pt', '1st_Electron_Eta', '1st_Muon_Eta', '2nd_Electron_Eta', '2nd_Muon_Eta',]

for itemwhichTriggerPassed in whichTriggerPassed:
     for item_elseFeature in elseFeature:
          which1DHistogramsToCreate.append(itemwhichTriggerPassed+"_"+item_elseFeature)
          #print(itemwhichTriggerPassed+"_"+item_elseFeature)
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

#
#  _                    _   __       _____                        _   ______               _       
# | |                  | | /_ |  _  |  __ \                      | | |  ____|             | |      
# | |     _____   _____| |  | | (_) | |__) |___  ___ ___  _ __ __| | | |____   _____ _ __ | |_ ___ 
# | |    / _ \ \ / / _ \ |  | |     |  _  // _ \/ __/ _ \| '__/ _` | |  __\ \ / / _ \ '_ \| __/ __|
# | |___|  __/\ V /  __/ |  | |  _  | | \ \  __/ (_| (_) | | | (_| | | |___\ V /  __/ | | | |_\__ \
# |______\___| \_/ \___|_|  |_| (_) |_|  \_\___|\___\___/|_|  \__,_| |______\_/ \___|_| |_|\__|___/
#                                                                                                  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                  
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


ee_Channel_MC_trigger_correlation = []
emu_Channel_MC_trigger_correlation = []
mumu_Channel_MC_trigger_correlation = []

ee_Channel_data_trigger_correlation = []
emu_Channel_data_trigger_correlation = []
mumu_Channel_data_trigger_correlation = []



year = ""

for inputFileName in RootfilesToRunOver:            
            inputFile = r.TFile.Open(inputFileName, "read")
            #inputTree = inputFile.Get("ttHTreeMaker/worldTree")
            inputTree=inputFile.Get("writeNTuple/NTuple")
            events = inputTree.GetEntries()
            printNumberOfFiles =str(inputFileName)+" has got "+str(events)+" to be evaluated"
            print(printNumberOfFiles)
            with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write(printNumberOfFiles)
                 nameRootoutput.write("\n")
            #print(events)
            str2="."
            bis = inputFileName.find(str2)
            textname=inputFileName[:bis]
            #print(textname)

            if "/" in inputFileName:
                textnamePre = inputFileName.split("/")
                textname = textnamePre[-1].replace(".root", "")
            #print(textname)
            isMC = True
            if "run" in textname or "Run" in textname: 
                isMC = False
                for itemRunPart_considered in RunPart_considered:
    #                 print(itemRunPart_considered)
                     if itemRunPart_considered in textname: 
                          isRun = itemRunPart_considered
                          break
            
            for itemyears_considered in years_considered:
     #            print(itemyears_considered)
                 if itemyears_considered in textname:
                      year = itemyears_considered
                      break
            #print("isRun is "+isRun)
            #print("isMC is "+str(isMC))
            #print("year is "+str(year))
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
            

            ee_OnlyTriggerDileptonWithSingleLepton_fired = 0
            emu_OnlyTriggerDileptonWithSingleLepton_fired = 0
            mumu_OnlyTriggerDileptonWithSingleLepton_fired = 0

            ee_MET_AND_TriggerDileptonWithSingleLepton_fired = 0
            emu_MET_AND_TriggerDileptonWithSingleLepton_fired = 0
            mumu_MET_AND_TriggerDileptonWithSingleLepton_fired = 0

            ee_Only_MET_Trigger_fired = 0
            emu_Only_MET_Trigger_fired = 0
            mumu_Only_MET_Trigger_fired = 0

            ee_allEvents = 0
            emu_allEvents = 0
            mumu_allEvents = 0

            which2DHistogramsCreatingForFiring = []
            for item_whichDileptonTriggerGroup in whichDileptonTriggerGroup:
                 for item_TriggerToBeConsidered in TriggerToBeConsidered:
                      which2DHistogramsCreatingForFiring.append(item_whichDileptonTriggerGroup+"_vs_"+item_TriggerToBeConsidered)

            channel = ""

            nEntries = inputTree.GetEntries()
            countevents =0

            count_HLT_DoubleEl_ = 0
            count_HLT_DoubleMuon_ = 0
            count_HLT_EMu_ = 0
            count_HLT_SingleMuon_ = 0
            count_HLT_SingleElectron_ = 0
            count_HLT_METTriggers_ = 0

            if not isMC:
                 allEntries = nEntries
            elif isMC:
                 allEntries = nEntries #10000
            for i in range(allEntries): #nEntries): #10): #int(nStart), int(nStopMinus1)):
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
                if(countevents%100000 == 0): print("Event:", i)


                for item_TriggerToBeConsidered in TriggerToBeConsidered:
                     globals()[item_TriggerToBeConsidered] = False
                     #print(str(item_TriggerToBeConsidered)+": "+str(eval(item_TriggerToBeConsidered)))

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
                   #print("inputTree."+item_TriggerToBeConsidered[:-1])
                   #print(eval("inputTree."+item_TriggerToBeConsidered[:-1]))
                   if eval("inputTree."+item_TriggerToBeConsidered[:-1]) == 1:
                        globals()[item_TriggerToBeConsidered] = True
                        globals()["count_"+str(item_TriggerToBeConsidered)] += 1

################## 2. Read in results from 1. and check category trigger result
################## ============================================================

                #### MET trigger
                for item_Trigger_met in Trigger_met:
                     #print(item_Trigger_met)
                     #print(eval(item_Trigger_met))
                     if eval(item_Trigger_met) == True:
                           count_HLT_METTriggers_ += 1
                           HLT_METTriggers_fired = True
                           break


                #### Double Electron trigger
                for item_Trigger_DoubleElectron in Trigger_DoubleElectron:
                     #print(item_Trigger_DoubleElectron)
                     #print(eval(item_Trigger_DoubleElectron))
                     if eval(item_Trigger_DoubleElectron) == True:
                          count_HLT_DoubleEl_ += 1
                          HLT_DoubleEl_fired = True
                          break

                #### Double Muon trigger
                for item_Trigger_DoubleMuon in Trigger_DoubleMuon:
                     #print(item_Trigger_DoubleMuon)
                     #print(eval(item_Trigger_DoubleMuon))
                     if eval(item_Trigger_DoubleMuon) == True:
                            count_HLT_DoubleMuon_ += 1
                            HLT_DoubleMuon_fired = True
                            break

                #### Electron Muon trigger
                for item_Trigger_ElectronMuon in Trigger_ElectronMuon:
                     #print(item_Trigger_ElectronMuon)
                     #print(eval(item_Trigger_ElectronMuon))
                     if eval(item_Trigger_ElectronMuon) == True:
                          count_HLT_EMu_ += 1
                          HLT_EMu_fired = True
                          break

                #### Single Electron trigger
                for item_Trigger_SingleElectron in Trigger_SingleElectron:
                     #print(item_Trigger_SingleElectron)
                     #print(eval(item_Trigger_SingleElectron))
                     if eval(item_Trigger_SingleElectron) == True:
                            count_HLT_SingleElectron_ += 1
                            HLT_SingleElectron_fired = True
                            break

                #### Single Muon trigger
                for item_Trigger_SingleMuon in Trigger_SingleMuon:
                     #print(item_Trigger_SingleMuon)
                     #print(eval(item_Trigger_SingleMuon))
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
                #print(numberOfJets)
                numberOfPrimaryVertices = int(inputTree.NPV_all)
                #print(numberOfPrimaryVertices)
                pileupSF = 1.
                theException = []
                try:
                     if pileupReweighting:
                          pileupSF = pileupWeights[numberOfPrimaryVertices]
                except:
                     theException.append(i)
                     pileupSF = 1.
                met = inputTree.met_Puppi.Pt()
                
                for item_lepton in range(inputTree.leptons.size()):
#################
######################### If Muon, let's check if inputTree.passes the chosen selections
######################### ==========================================================
                    #print("inputTree.muon_PFIso["+str(item_lepton)+"]")
                    #print(inputTree.muon_PFIso[item_lepton])
                    #print("Is inputTree.muon_PFIso[item_lepton] >= 76?")
                    #print(ord(inputTree.muon_PFIso[item_lepton]) >= 76)
                    if(inputTree.muon_CutBasedID_Medium[item_lepton] == "Y" and inputTree.leptons[item_lepton].Pt() > 20. and (abs(inputTree.leptons[item_lepton].Eta()) < 2.4) and ord(inputTree.muon_PFIso[item_lepton]) >= 76):
                              nMuons += 1
                              get_indicesAndPt_of_passed_muons.append((inputTree.leptons[item_lepton].Pt(), item_lepton))
                              get_indicesAndPt_of_passed_lepton.append((inputTree.leptons[item_lepton].Pt(), item_lepton, "M"))
                              #print("Passed Muon with "+str(inputTree.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(inputTree.muon_CutBasedID_Medium[item_lepton]))
                              
######################### If Electron, let's check if it passes the chosen selections
######################### ============================================================

                    elif(inputTree.eleID_MVA_Iso_90[item_lepton] == 1 and inputTree.leptons[item_lepton].Pt() > 20. and (abs(inputTree.leptons[item_lepton].Eta()) < 2.4)):
                         ### removing ECAL transistion region 1.4442 < abs(eta) < 1.566
                         if (abs(inputTree.leptons[item_lepton].Eta()) < 1.442) or (abs(inputTree.leptons[item_lepton].Eta()) > 1.566):
                              nElectrons += 1
                              get_indicesAndPt_of_passed_lepton.append((inputTree.leptons[item_lepton].Pt(), item_lepton, "E"))
                              get_indicesAndPt_of_passed_electrons.append((inputTree.leptons[item_lepton].Pt(), item_lepton))
                              #print("Passed Electron with "+str(inputTree.leptons[item_lepton].Pt())+" and "+str(item_lepton)+" and "+str(inputTree.eleID_MVA_Iso_90[item_lepton]))

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
                    #print("leadingLeptonIndex = "+str(leadingLeptonIndex))
                    if len(sortingLeptons) > 1:
                         subleadingLeptonIndex = sortingLeptons[1][1]
                         subleadingLeptonisMuon = sortingLeptons[1][2]
                     #    print("subleadingLeptonIndex = "+str(subleadingLeptonIndex))
                    #print(sortingLeptons)
                    
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
                         #print("leadingMuonIndex = "+str(leadingMuonIndex))
                         if len(sortingMuons) > 1:
                              subleadingMuonIndex = sortingMuons[1][1]
                          #    print("subleadingMuonIndex = "+str(subleadingMuonIndex))
                    ### Get leading and subleading electron index if available
                    sortingElectrons = sorted(get_indicesAndPt_of_passed_electrons, reverse=True)
                    if len(sortingElectrons) > 0:
                         leadingElectronIndex = sortingElectrons[0][1]
                         #print("leadingElectronIndex = "+str(leadingElectronIndex))
                         if len(sortingElectrons) > 1:
                              subleadingElectronIndex = sortingElectrons[1][1]
                          #    print("subleadingElectronIndex = "+str(subleadingElectronIndex))

                    ## ----------------------------- categorisation into channels
                    lepton_weight = 1.
                    if leadingLepton.Pt() >= 25. and subleadingLepton.Pt() >= 20. and dilepton.M() > 20.:
                         # if this passes we have an OS muon pair --> mumu Channel
                         if (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -169) and met > 40.:
                              #if not HLT_METTriggers_fired:
                              mumu_allEvents += 1
                              if HLT_DoubleMuon_fired or HLT_SingleMuon_fired:
                                        mumu_OnlyTriggerDileptonWithSingleLepton_fired += 1
                              
                              if HLT_METTriggers_fired:
                                   mumu_Only_MET_Trigger_fired += 1
                                   if HLT_DoubleMuon_fired or HLT_SingleMuon_fired:
                                        mumu_MET_AND_TriggerDileptonWithSingleLepton_fired += 1
                                   #else:
                                    #    mumu_Only_MET_Trigger_fired += 1

                                   lepton_weight = leptonSF(gatheredMuonIDsfWeights, leadingLepton.Pt(), abs(leadingLepton.Eta()))*leptonSF(gatheredMuonIsosfWeights, leadingLepton.Pt(), abs(leadingLepton.Eta()))*leptonSF(gatheredMuonIDsfWeights, subleadingLepton.Pt(),  abs(subleadingLepton.Eta()))*leptonSF(gatheredMuonIsosfWeights, subleadingLepton.Pt(), abs(subleadingLepton.Eta()))
                                   if numberOfJets == 0 or numberOfJets == 1:
                                             print(numberOfJets)
                                             print(numberOfPrimaryVertices)
                                             print(met)
                                             print(lepton_weight)

                                   if not isMC:                                   
                                        mumu_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        mumu_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        mumu_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_DoubleMuon_fired, HLT_SingleMuon_fired, lepton_weight])
                                   else:
                                        mumu_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        mumu_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        mumu_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_DoubleMuon_fired, HLT_SingleMuon_fired, pileupSF*lepton_weight])

                         # if this passes we have an OS electron pair --> ee Channel
                         elif (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -121) and met > 40:
                              #if not HLT_METTriggers_fired:
                              ee_allEvents +=1
                              if HLT_DoubleEl_fired or HLT_SingleElectron_fired:
                                   ee_OnlyTriggerDileptonWithSingleLepton_fired += 1
                              if HLT_METTriggers_fired:
                                   ee_Only_MET_Trigger_fired += 1
                                   if HLT_DoubleEl_fired or HLT_SingleElectron_fired:
                                        ee_MET_AND_TriggerDileptonWithSingleLepton_fired += 1
#                                   else:
 #                                       ee_Only_MET_Trigger_fired += 1

                                   lepton_weight = leptonSF(gatheredElectronIDsfWeights, leadingLepton.Eta(), leadingLepton.Pt())*leptonSF(gatheredElectronRECOsfWeights, leadingLepton.Eta(), leadingLepton.Pt())*leptonSF(gatheredElectronIDsfWeights, subleadingLepton.Eta(), subleadingLepton.Pt())*leptonSF(gatheredElectronRECOsfWeights, subleadingLepton.Eta(), subleadingLepton.Pt())
                                   if numberOfJets == 0 or numberOfJets == 1:
                                             print(numberOfJets)
                                             print(numberOfPrimaryVertices)
                                             print(met)
                                             print(lepton_weight)

                                   if not isMC:
                                        ee_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        ee_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        ee_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_DoubleEl_fired, HLT_SingleElectron_fired, lepton_weight])
                                        #                                   ee_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, ElectronTriggersfired, HLT_SingleElectron_fired, lepton_weight])
                                        #if not HLT_DoubleEl_fired:
                                        #    if HLT_SingleElectron_fired:
                                        #        print(ee_Channel_data_1st_Lepton_pt_eta_iso)
                                        #       print(ee_Channel_data_2nd_Lepton_pt_eta_iso)
                                       #      print(ee_Channel_data_njet_nvertex_met)
                                        #     print("FFFFFFfOOOOOOOOOOOOOOOOuuuuunnnnnnd!")
                                         #    sys.exit()
                                   else:
                                        ee_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        ee_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        ee_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_DoubleEl_fired, HLT_SingleElectron_fired, pileupSF*lepton_weight])
                                        #if not HLT_DoubleEl_fired:
                                        #    
                                        # if HLT_SingleElectron_fired:
                                        #  print(ee_Channel_MC_1st_Lepton_pt_eta_iso)
                                       # print(ee_Channel_MC_2nd_Lepton_pt_eta_iso)
                                        #print(ee_Channel_MC_njet_nvertex_met)
                                        #print("FFFFFFfOOOOOOOOOOOOOOOOuuuuunnnnnnd!")
                                        #sys.exit()

                         # if this passes we have an OS electron & muon --> emu Channel
                         elif (inputTree.lepPdgId[leadingLeptonIndex] * inputTree.lepPdgId[subleadingLeptonIndex] == -143):
                              #ElectronMuonTriggersfired = False
                              SingleTriggersfired = False
                              
                              #if HLT_SingleElectron_fired or HLT_SingleMuon_fired or HLT_EMu_fired:
                               #    ElectronMuonTriggersfired = True

                              if HLT_SingleElectron_fired or HLT_SingleMuon_fired:
                                   SingleTriggersfired = True
                                   
                              emu_allEvents += 1
#                              if not HLT_METTriggers_fired:
                              if HLT_EMu_fired or SingleTriggersfired:
                                        emu_OnlyTriggerDileptonWithSingleLepton_fired += 1

                              if HLT_METTriggers_fired:
                                   emu_Only_MET_Trigger_fired += 1

                                   if HLT_EMu_fired or SingleTriggersfired:
                                        emu_MET_AND_TriggerDileptonWithSingleLepton_fired += 1
#                                   else:
 #                                       emu_Only_MET_Trigger_fired += 1
                                        
                                   muon_weight = 1.
                                   electron_weight = 1.

                                   if abs(inputTree.lepPdgId[leadingLeptonIndex]) == 11:
                                        electron_weight = leptonSF(gatheredElectronIDsfWeights, leadingLepton.Eta(), leadingLepton.Pt())*leptonSF(gatheredElectronRECOsfWeights, leadingLepton.Eta(), leadingLepton.Pt())
                                        muon_weight = leptonSF(gatheredMuonIDsfWeights, subleadingLepton.Pt(),  abs(subleadingLepton.Eta()))*leptonSF(gatheredMuonIsosfWeights, subleadingLepton.Pt(), abs(subleadingLepton.Eta()))

                                   elif abs(inputTree.lepPdgId[subleadingLeptonIndex]) == 11:
                                        electron_weight = leptonSF(gatheredElectronIDsfWeights, subleadingLepton.Eta(), subleadingLepton.Pt())*leptonSF(gatheredElectronRECOsfWeights, subleadingLepton.Eta(), subleadingLepton.Pt())
                                        muon_weight = leptonSF(gatheredMuonIDsfWeights, leadingLepton.Pt(), abs(leadingLepton.Eta()))*leptonSF(gatheredMuonIsosfWeights, leadingLepton.Pt(), abs(leadingLepton.Eta()))
 
                                   if numberOfJets == 0 or numberOfJets == 1:
                                             print(numberOfJets)
                                             print(numberOfPrimaryVertices)
                                             print(met)
                                             print(lepton_weight)

                                   if not isMC:
                                        emu_Channel_data_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        emu_Channel_data_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        emu_Channel_data_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_EMu_fired, SingleTriggersfired, muon_weight*electron_weight])
                                   else:
                                        emu_Channel_MC_1st_Lepton_pt_eta_iso.append([leadingLepton.Pt(), leadingLepton.Eta(), inputTree.lepPfIso[leadingLeptonIndex]])
                                        emu_Channel_MC_2nd_Lepton_pt_eta_iso.append([subleadingLepton.Pt(), subleadingLepton.Eta(), inputTree.lepPfIso[subleadingLeptonIndex]])
                                        emu_Channel_MC_njet_nvertex_met.append([numberOfJets, numberOfPrimaryVertices, met, HLT_EMu_fired, SingleTriggersfired, pileupSF*muon_weight*electron_weight])

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

           
            inputFile.Close()
            if LeptonSFapplication:
                 inputElectronIDsfFile.Close()
                 inputElectronRECOsfFile.Close()
                 inputMuonIDsfFile.Close()
                 inputMuonISOsfFile.Close()
            if isMC:
                 ee_Channel_MC_trigger_correlation.append([ee_allEvents, ee_OnlyTriggerDileptonWithSingleLepton_fired, ee_MET_AND_TriggerDileptonWithSingleLepton_fired, ee_Only_MET_Trigger_fired])
                 emu_Channel_MC_trigger_correlation.append([emu_allEvents, emu_OnlyTriggerDileptonWithSingleLepton_fired, emu_MET_AND_TriggerDileptonWithSingleLepton_fired, emu_Only_MET_Trigger_fired])
                 mumu_Channel_MC_trigger_correlation.append([mumu_allEvents, mumu_OnlyTriggerDileptonWithSingleLepton_fired, mumu_MET_AND_TriggerDileptonWithSingleLepton_fired, mumu_Only_MET_Trigger_fired])
            else:
                 ee_Channel_data_trigger_correlation.append([ee_allEvents, ee_OnlyTriggerDileptonWithSingleLepton_fired, ee_MET_AND_TriggerDileptonWithSingleLepton_fired, ee_Only_MET_Trigger_fired])
                 emu_Channel_data_trigger_correlation.append([emu_allEvents, emu_OnlyTriggerDileptonWithSingleLepton_fired, emu_MET_AND_TriggerDileptonWithSingleLepton_fired, emu_Only_MET_Trigger_fired])
                 mumu_Channel_data_trigger_correlation.append([mumu_allEvents, mumu_OnlyTriggerDileptonWithSingleLepton_fired, mumu_MET_AND_TriggerDileptonWithSingleLepton_fired, mumu_Only_MET_Trigger_fired])

            if verbose:
                with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                                     nameRootoutput.write("count_HLT_DoubleEl_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_DoubleEl_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele27_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele27_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele27_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele27_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele35_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele35_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele35_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele35_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele38_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele38_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele38_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele38_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele40_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele40_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele40_WPTight_Gsf_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele40_WPTight_Gsf_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_IsoMu27_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_IsoMu27_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_IsoMu27_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_IsoMu27_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_IsoMu24_eta2p1_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_IsoMu24_eta2p1_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_IsoMu24_eta2p1_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_IsoMu24_eta2p1_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_DoubleEle33_CaloIdL_MW_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_DoubleEle33_CaloIdL_MW_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_DoubleEle33_CaloIdL_MW_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_DoubleEle33_CaloIdL_MW_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_DoubleEle25_CaloIdL_MW_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_DoubleEle25_CaloIdL_MW_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_DoubleEle25_CaloIdL_MW_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_DoubleEle25_CaloIdL_MW_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_PFMET120_PFMHT120_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_PFMET120_PFMHT120_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_PFMET120_PFMHT120_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_PFMET120_PFMHT120_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_PFMET250_HBHECleaned_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_PFMET250_HBHECleaned_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_PFMET250_HBHECleaned_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_PFMET250_HBHECleaned_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_METTriggers_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_METTriggers_)
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write("count_HLT_DoubleEl_")
                                     nameRootoutput.write("\n")
                                     nameRootoutput.write(count_HLT_DoubleEl_)

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

#print("\n\nee_Channel_data_1st_Lepton_pt_eta_iso")
#print(ee_Channel_data_1st_Lepton_pt_eta_iso)
#print("\n\nee_Channel_MC_1st_Lepton_pt_eta_iso")
#print(ee_Channel_MC_1st_Lepton_pt_eta_iso)


#print("\n\nee_Channel_data_njet_nvertex_met")
#print(ee_Channel_data_njet_nvertex_met)
#print("\n\nee_Channel_MC_njet_nvertex_met")
#print(ee_Channel_MC_njet_nvertex_met)

#sys.exit()



#ee_Channel_data_1st_Lepton_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0])]
#print(ee_Channel_data_1st_Lepton_pt_Order)
#ee_Channel_data_1st_Lepton_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][1])]
#print(ee_Channel_data_1st_Lepton_eta_Order)
#ee_Channel_data_1st_Lepton_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][2])]
#print(ee_Channel_data_1st_Lepton_iso_Order)


#==============================================================
#==============================================================
#### DETERMINING SYSTEMATIC uncertainties using the same array:
#==============================================================
#==============================================================

#### by subset selections on NUMBER OF JETS:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff number of jets < 3:
#ee_Channel_data_1st_Lepton_pt_Order_smaller3Njets = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][0] <3]
### ordering of events according to 1st lepton pt iff number of jets >= 3:
#ee_Channel_data_1st_Lepton_pt_Order_3AndMoreNjets = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][0] >= 3]


#### by subset selections on NUMBER OF VERTICES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff number of vertices < 30: 
#ee_Channel_data_1st_Lepton_pt_Order_smaller30Nvertices = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][1] <30]
### ordering of events according to 1st lepton pt iff number of vertices >= 30:
#ee_Channel_data_1st_Lepton_pt_Order_30AndMoreNvertices = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][1] >= 30]


#### by subset selections on MET:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ordering of events according to 1st lepton pt iff MET < 80 GeV:
#ee_Channel_data_1st_Lepton_pt_Order_METsmaller80 = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][2] <80.]
### ordering of events according to 1st lepton pt iff MET >= 80 GeV:
#ee_Channel_data_1st_Lepton_pt_Order_METlarger80 = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_1st_Lepton_pt_eta_iso), key=lambda elem: elem[1][0]) if ee_Channel_data_njet_nvertex_met[indices_ordered[0]][2] >= 80.]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#print("\n\nee_Channel_data_2nd_Lepton_pt_eta_iso")
#print(ee_Channel_data_2nd_Lepton_pt_eta_iso)

#print("\n\nee_Channel_MC_2nd_Lepton_pt_eta_iso")
#print(ee_Channel_MC_2nd_Lepton_pt_eta_iso)

#ee_Channel_data_2nd_Lepton_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][0])]
#[indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=itemgetter(1))]                                                               
#print("ee_Channel_data_2nd_Lepton_pt_Order")
#print(ee_Channel_data_2nd_Lepton_pt_Order)

#split_ee_Channel_data_2nd_Lepton_pt_Order = splitListintoSubListsOfPieceSize(ee_Channel_data_2nd_Lepton_pt_Order, 2)
#for item_split_ee_Channel_data_2nd_Lepton_pt_Order in split_ee_Channel_data_2nd_Lepton_pt_Order:
#     print(item_split_ee_Channel_data_2nd_Lepton_pt_Order)
#ee_Channel_data_2nd_Lepton_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][1])]
#print(ee_Channel_data_2nd_Lepton_eta_Order)
#ee_Channel_data_2nd_Lepton_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_2nd_Lepton_pt_eta_iso), key=lambda elem: elem[1][2])]
#print(ee_Channel_data_2nd_Lepton_iso_Order)

#print("\n\nee_Channel_data_njet_nvertex_met")
#print(ee_Channel_data_njet_nvertex_met)

#print("\n\nee_Channel_MC_njet_nvertex_met")
#print(ee_Channel_MC_njet_nvertex_met)



#ee_Channel_data_njet_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][0])]
#print(ee_Channel_data_njet_Order)
#ee_Channel_data_nvertex_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][1])]
#print(ee_Channel_data_nvertex_Order)
#ee_Channel_data_met_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate(ee_Channel_data_njet_nvertex_met), key=lambda elem: elem[1][2])]
#print("ee_Channel_data_met_Order")
#print(ee_Channel_data_met_Order)

verbose = False
if verbose:
     print("\n\nemu_Channel_data_1st_Lepton_pt_eta_iso")
     print(emu_Channel_data_1st_Lepton_pt_eta_iso)
     print("\n\nemu_Channel_data_2nd_Lepton_pt_eta_iso")
     print(emu_Channel_data_2nd_Lepton_pt_eta_iso)
     print("\n\nemu_Channel_data_njet_nvertex_met")
     print(emu_Channel_data_njet_nvertex_met)

     print("\n\nmumu_Channel_data_1st_Lepton_pt_eta_iso")
     print(mumu_Channel_data_1st_Lepton_pt_eta_iso)
     print("\n\nmumu_Channel_data_2nd_Lepton_pt_eta_iso")
     print(mumu_Channel_data_2nd_Lepton_pt_eta_iso)
     print("\n\nmumu_Channel_data_njet_nvertex_met")
     print(mumu_Channel_data_njet_nvertex_met)

     print("\n\nemu_Channel_MC_1st_Lepton_pt_eta_iso")
     print(emu_Channel_MC_1st_Lepton_pt_eta_iso)
     print("\n\nemu_Channel_MC_2nd_Lepton_pt_eta_iso")
     print(emu_Channel_MC_2nd_Lepton_pt_eta_iso)
     print("\n\nemu_Channel_MC_njet_nvertex_met")
     print(emu_Channel_MC_njet_nvertex_met)

     print("\n\nmumu_Channel_MC_1st_Lepton_pt_eta_iso")
     print(mumu_Channel_MC_1st_Lepton_pt_eta_iso)
     print("\n\nmumu_Channel_MC_2nd_Lepton_pt_eta_iso")
     print(mumu_Channel_MC_2nd_Lepton_pt_eta_iso)
     print("\n\nmumu_Channel_MC_njet_nvertex_met")
     print(mumu_Channel_MC_njet_nvertex_met)
## Having seen how the wind blows, let's write the remaining ones in shorter way ;)
#==================================================================
#==================================================================
#
#  _                    _   ___                             _                 ______               _       
# | |                  | | |__ \   _      /\               | |               |  ____|             | |      
# | |     _____   _____| |    ) | (_)    /  \   _ __   __ _| |_   _ _______  | |____   _____ _ __ | |_ ___ 
# | |    / _ \ \ / / _ \ |   / /        / /\ \ | '_ \ / _` | | | | |_  / _ \ |  __\ \ / / _ \ '_ \| __/ __|
# | |___|  __/\ V /  __/ |  / /_   _   / ____ \| | | | (_| | | |_| |/ /  __/ | |___\ V /  __/ | | | |_\__ \
# |______\___| \_/ \___|_| |____| (_) /_/    \_\_| |_|\__,_|_|\__, /___\___| |______\_/ \___|_| |_|\__|___/
#                                                              __/ |                                       
#                                                             |___/
#  _____       _                      _              ______  _      _      _                 _             
# |  __ \     | |                    (_)            |  ____|/ _|/ _(_)    (_)               (_)            
# | |  | | ___| |_ ___ _ __ _ __ ___  _ _ __   ___  | |__  | |_| |_ _  ___ _  ___ _ __   ___ _  ___  ___   
# | |  | |/ _ \ __/ _ \ '__| '_ ` _ \| | '_ \ / _ \ |  __| |  _|  _| |/ __| |/ _ \ '_ \ / __| |/ _ \/ __|  
# | |__| |  __/ ||  __/ |  | | | | | | | | | |  __/ | |____| | | | | | (__| |  __/ | | | (__| |  __/\__ \  
# |_____/ \___|\__\___|_|  |_| |_| |_|_|_| |_|\___| |______|_| |_| |_|\___|_|\___|_| |_|\___|_|\___||___/  
#                                                                                                          
# ------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------


whichCommandToExecute =[]

whatToDivide = []
channel = ["ee", "emu", "mumu"]
whichType = ["data", "MC"]
whichLepton = ["1st_Lepton", "2nd_Lepton"]
whichleptonFeature =["pt", "eta", "iso"]
#inWhichSubsetsToDivide = ["smaller3Njets", "3AndMoreNjets", "smaller30Nvertices", "30AndMoreNvertices", "METsmaller80", "METlarger80"]                                                               
#whichOtherFeature = ["njet", "nvertex", "met"]                   
#inWhichSubsetsToDivide = [["smaller3Njets", "<3"], ["3AndMoreNjets", ">=3"], ["smaller30Nvertices", "< 30"], ["30AndMoreNvertices", ">= 30"], ["METsmaller80", "< 80."], ["METlarger80", ">= 80."]]

### Let's analyse events having survived ...

for item_channel in channel:
    for item_whichType in whichType:
     #   print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
     #   print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
     #   print(len(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")))
        
             
        
# Ordered according to number of jets
#        print(str(item_channel)+"_Channel_"+item_whichType+"_njet_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][0])]")
  #      exec(str(item_channel)+"_Channel_"+item_whichType+"_njet_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][0])]")
 #       print(eval(str(item_channel)+"_Channel_"+item_whichType+"_njet_Order"))
  #      print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_njet_Order")))
        prefix = "Njets_"
        which_binning = nJets_binning
        which_Order ="_njet_Order"
#        what ="getoutValues("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0)"
#        print(what)                                                                                                                                                                                                                                                                                                                                                               
#        print(eval(what))
#        what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0, which_binning)"
#        print(what1)                                                                                                                                                                                                                                                                                                                                                              
#        print(eval(what1))                                                                                                                                                                                                                                                                                                                                                        
#        what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+"eval(what1))"   
#        print(what2)                                                                                                                                                                                                                                                                                                                                                               
#        print(eval(what2))  

        
#        what1 = "getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0, "+str(which_binning)+")"
 #       what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+what1+")"
        what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, which_binning, 0)"
 #       print(what2)
  #      print(eval(what2))

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0)")
        
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
        
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))


#------------------------------------------------------------------------ Systematics for varying N Jets ------------------------        


        prefix = "Eta_2D_syst_njet_larger4_"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABSupperSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0, 4)"
#        print(what2)
#        print(eval(what2))
#        for item_what2 in eval(what2):
#             print(item_what2)
#             for subitem_item_what2 in item_what2:
#                  print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])
        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))         
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

 
        prefix = "Eta_2D_syst_njet_smallerEqual4_"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABlowerSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0, 4)"
#        print(what2)
#        print(eval(what2))
#        for item_what2 in eval(what2):
#             print(item_what2)
#             for subitem_item_what2 in item_what2:
#                  print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])
        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))   
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")


#--------------------------------------------------------------------------- Ordered according to number of vertices

     #   exec(str(item_channel)+"_Channel_"+item_whichType+"_nvertex_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][1])]")

        prefix = "NVertex_"
        which_binning = nVertex_binning
        which_Order ="_nvertex_Order"
#        if item_channel == "emu":
  #           print(str(item_channel)+"_Channel_"+item_whichType+"_nvertex_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][1])]")
   #          print(eval(str(item_channel)+"_Channel_"+item_whichType+"_nvertex_Order"))
    #         print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_nvertex_Order")))
       #      what ="getoutValues("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)"
        #     print(what)
         #    print(eval(what))
      #       what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1, which_binning)"
     #        print(what1)
    #         print(eval(what1))
   #          what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+"eval(what1))"
  #           print(what2)
 #            print(eval(what2))
#

  #      what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1, "+str(which_binning)+")"
 #       what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+what1+")"
        what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, which_binning, 1)"
#        print(what2)                                                                                                                                                                                                                                                                                                                                                                 
 #       print(eval(what2))
#arrayToBeSplit, array_binning, whichIndex)
        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
        

        test_Nvertex_Splitting = False
        if test_Nvertex_Splitting:
             print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
             print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
             print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))

                     
             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))
             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
             print(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
             
             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))

        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

        
        



#------------------------------------------------------------------------ Systematics for varying number of vertices 

        prefix = "Eta_2D_syst_nvertex_larger30_"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABSupperSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1, 30)"

#        print(what2)
#        print(eval(what2))
#        for item_what2 in eval(what2):
#             print(item_what2)
#             for subitem_item_what2 in item_what2:
#                  print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])
        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))         
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")
        
 
        prefix = "Eta_2D_syst_nvertex_smallerEqual30"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABlowerSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1, 30)"
 #       print(what2)
 #       print(eval(what2))
 #       for item_what2 in eval(what2):
 #            print(item_what2)
 #            for subitem_item_what2 in item_what2:
 #                 print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))  
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

#----------------------------------------------------------------------  Ordered according to MET
#        print(str(item_channel)+"_Channel_"+item_whichType+"_met_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][2])]")
        exec(str(item_channel)+"_Channel_"+item_whichType+"_met_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met), key=lambda elem: elem[1][2])]")
#        print(eval(str(item_channel)+"_Channel_"+item_whichType+"_met_Order"))
#        print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_met_Order")))
 
        prefix = "MET_"
        which_binning = met_binning
        which_Order ="_met_Order"
        testingMET = False
        if testingMET:
             what ="getoutValues("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2)"
             print(what)
             print(eval(what))
             what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2, which_binning)"
             print(what1)
             print(eval(what1))
             what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+"eval(what1))"
             print(what2)
             print(eval(what2))


             what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+which_Order+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2, "+str(which_binning)+")"
             what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+which_Order+","+what1+")"
             
             exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                  prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2)")

             
             print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
             print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))

             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
              
             print(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType))

             print(prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType)
             print(eval(prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType))
             print("NNNNNNNNNNNNNEEEEEEEEEEEEEEWWWWWWWWWWWWWWWWWWWWW")

        what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, which_binning, 2)"
        if testingMET: 
             print(what2)
             print(eval(what2))
#             print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
 #            print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2)")
        

         #print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
         #print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
         #print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
         #print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
         #print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
         #print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

        if testingMET:

              print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
              print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
              
              print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
              
              print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))         

              print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
         
              print(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
              
              print(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType))

              print(prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType)
              print(eval(prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType))


  
#---------------------------------------------------------------------  Systematics for varying number MET

        prefix = "Eta_2D_syst_met_larger100_"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABSupperSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2, 100)"
#        print(what2)
#        print(eval(what2))
#        for item_what2 in eval(what2):
#             print(item_what2)
#             for subitem_item_what2 in item_what2:
#                  print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))         
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

 
        prefix = "Eta_2D_syst_met_smallerEqual100"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABlowerSystematics("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+", "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2, 100)"
#        print(what2)
#        print(eval(what2))
#        for item_what2 in eval(what2):
#             print(item_what2)
#             for subitem_item_what2 in item_what2:
#                  print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")
#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))  
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")

  
      # Ordered according to abs(eta) 2D:                                                                                                                                                                                                                                                                                                                                               
        prefix = "Eta_2D_"
        item_eta_binning = lepton_2D_etaByEta_binning
#        print("------------------------------------------------####22222DDDDDD")
        what2 = "getoutIndexArrayFor2DSplitEtaABS("+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso, "+ str(item_eta_binning)+", "+str(item_eta_binning)+")"
 #       print(what2)
 #       print(eval(what2))
 #       for item_what2 in eval(what2):
 #            print(item_what2)
 #            for subitem_item_what2 in item_what2:
 #                 print([eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]"), eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso["+str(subitem_item_what2)+"]")])

        exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
             prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")

#        print(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"))
#        print(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_1st_Lepton_pt_eta_iso"))
#        print(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso")
#        print(eval(item_channel+"_Channel_"+item_whichType+"_2nd_Lepton_pt_eta_iso"))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType))
#        print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
#        print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType))
        with open(nameRootoutputCsVName, 'a') as nameRootoutput:
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)
             nameRootoutput.write("\n")
             nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType)))
             nameRootoutput.write("\n")
             nameRootoutput.write("\n")
             nameRootoutput.write("------------------------------------------------------------------------------------\n")



        for item_whichLepton in whichLepton:
# -----------------------------------------------------------------   Ordered according to pt  
            #print(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][0])]")
  #          exec(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][0])]")
            #print(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order"))
            #print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order")))

            
#            print("------------------------------------------------####1111y")
            
#            what ="getoutValues("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order,"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 0)"
#            print(what)
#            print(eval(what))
                        
#            what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order,"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 0, lepton_Pt_binning)"
#            print(what1)
#            print(eval(what1))
#            what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order,"+"eval(what1))"            
            
#            print(what2)
#            print(eval(what2))
 #           sys.exit()
#            eff = eval(what2)
            prefix = "Pt_"
       #     print(binArrayGetOutTheValuesFromSplitIndices(eval(what2), eval(item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso"), 0))
            #getoutIndexArrayForSplit(splitIndexArray, dataOrMCValueArray, which_value, valuesAccordingToSplit)
        #    print(binArrayGetOutEffValuesFromSplitIndices(eval(what2), eval(item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso"), eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"), 0))
#            what1 = "getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order,"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 0, "+str(lepton_Pt_binning)+")"
 #           what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_Order,"+what1+")" 
            
            what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, which_binning, 0)"
#            print(what2)
 #           print(eval(what2))

            exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
                  prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+ 
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+ 
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
                  prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 0)")
            

#            (Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_ee_Channel_data_1st_Lepton,
#             Pt_getAllShots_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getAllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getAdditionalCounts_HLT_SingleMuon_fired_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getlepton_weight_perBin_ee_Channel_data_1st_Lepton,
#             Pt_get_HLT_DoubleMuon_fired_Weights_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getAdditional_HLT_SingleMuon_fired_Weights_perBin_ee_Channel_data_1st_Lepton,
#             Pt_getBinnedArrayFromDataOrMCValueArray_ee_Channel_data_1st_Lepton) = binArrayGetOutEffValuesFromSplitIndices(getSplitOfIndexArray(ee_Channel_data_1st_Lepton_pt_Order, getoutIndexArrayForSplit(ee_Channel_data_1st_Lepton_pt_Order, ee_Channel_data_1st_Lepton_pt_eta_iso, 0, [20., 30., 40., 60., 80., 100., 200.])), ee_Channel_data_1st_Lepton_pt_eta_iso, ee_Channel_data_njet_nvertex_met, 0)
            with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("------------------------------------------------------------------------------------\n")


#            print("Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
 #           print(eval("Pt_getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton))
  #          print("Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
   #         print(eval("Pt_getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton))

                 
            # Ordered according to eta
#            print(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][1])]")
            #exec(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][1])]")
#            print(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order"))
#            print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order")))
  #          print("splitIntoBins("+str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order")
   #         print(splitIntoBins(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_eta_Order")))
            which_Order = "_eta_Order"
            item_eta_binning = []
            if item_channel == 'ee':
                 item_eta_binning = electron_eta_binning
            else: 
                 item_eta_binning = muon_eta_binning
#            print("------------------------------------------------####1111y")

#            what ="getoutValues("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+which_Order+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 1)"
#            print(what)
#            print(eval(what))

#            what1 ="getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+which_Order+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 1, item_eta_binning)"
#            print(what1)
#            print(eval(what1))
#            what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+which_Order+","+"eval(what1))"

#            print(what2)
#            print(eval(what2))

            prefix = "Eta_"

            what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, item_eta_binning, 1)"
#            print(what2)
 #           print(eval(what2))

#            print(binArrayGetOutEffValuesFromSplitIndices(eval(what2), eval(item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso"), eval(item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met"), 1))
 #           what1 = "getoutIndexArrayForSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+which_Order+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, 1, "+str(item_eta_binning)+")"
  #          what2 = "getSplitOfIndexArray("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+which_Order+","+what1+")"
            exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
                  prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
            prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
                 prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 1)")

            
#            print(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
 #           print(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton))
  #          print(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
   #         print(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton))
            with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval(prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("------------------------------------------------------------------------------------\n")


            # Ordered according to iso
            # prefix = "Iso_"
            
#            what2 = "getoutIndexArrayFor1DSplit("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, which_binning, 2)"
#            print(what2)                                                                                                                                                                                                                                                                                                                                                                  
 #           print(eval(what2))      
#            exec("("+prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
 #                 prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
  #          prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
   #         prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
    #        prefix+"getAllShots_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
     #       prefix+"getAllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
      #      prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
       #     prefix+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
 #           prefix+"getAdditionalCounts_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
  #          prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
   #         prefix+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
    #        prefix+"getlepton_weight_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
     ##       prefix+"get_HLT_DoubleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
       #     prefix+"getAdditional_HLT_SingleMuon_fired_Weights_perBin_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+","+
        #         prefix+"getBinnedArrayFromDataOrMCValueArray_"+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+") = binArrayGetOutEffValuesFromSplitIndices("+what2+","+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso, "+item_channel+"_Channel_"+item_whichType+"_njet_nvertex_met, 2)")


            
    #        print(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][2])]")
  #          exec(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_iso_Order = [indices_ordered[0] for indices_ordered in sorted(enumerate("+item_channel+"_Channel_"+item_whichType+"_"+item_whichLepton+"_pt_eta_iso), key=lambda elem: elem[1][2])]")
     #       print(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_iso_Order"))
      #      print(len(eval(str(item_channel)+"_Channel_"+item_whichType+"_"+item_whichLepton+"_iso_Order")))
           


#
#  _                    _   ____             _____           _        ______         _                 
# | |                  | | |___ \   _       / ____|         | |      |  ____|       | |                
# | |     _____   _____| |   __) | (_)     | (___   ___ __ _| | ___  | |__ __ _  ___| |_ ___  _ __ ___ 
# | |    / _ \ \ / / _ \ |  |__ <           \___ \ / __/ _` | |/ _ \ |  __/ _` |/ __| __/ _ \| '__/ __|
# | |___|  __/\ V /  __/ |  ___) |  _       ____) | (_| (_| | |  __/ | | | (_| | (__| || (_) | |  \__ \
# |______\___| \_/ \___|_| |____/  (_)     |_____/ \___\__,_|_|\___| |_|  \__,_|\___|\__\___/|_|  |___/
#
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                                                                                  


#availablePlots = ["MET_", "Pt_", "Eta_", "Njets_", "NVertex_", "Eta_2D_", "Eta_2D_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_njet_smallerEqual4_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_nvertex_larger30_", "Eta_2D_syst_nvertex_smallerEqual30", "Eta_2D_syst_met_larger100_", "Eta_2D_syst_met_smallerEqual100" ]

#availableEffPlots = ["MET_", "Pt_", "Eta_", "Njets_", "NVertex_", "Eta_2D_"]

availablePlotsForSystematicUncertaintyDetermination = ["Eta_2D_", "Eta_2D_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_njet_smallerEqual4_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_nvertex_larger30_", "Eta_2D_syst_nvertex_smallerEqual30", "Eta_2D_syst_met_larger100_", "Eta_2D_syst_met_smallerEqual100" ]

availableLeptonLeadSubleadPlots = ["Pt_", "Eta_"]

availablePlots = ["MET_", "Njets_", "NVertex_", "Eta_2D_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_njet_smallerEqual4_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_nvertex_larger30_", "Eta_2D_syst_nvertex_smallerEqual30", "Eta_2D_syst_met_larger100_", "Eta_2D_syst_met_smallerEqual100"]


channel = ["ee", "emu", "mumu"]
ifweightedOrNot = ["weighted_"]
ifweightedBin = ["", "weighted_"]

whichLepton = ["1st_Lepton", "2nd_Lepton"]


## leading or subleading lepton plots

for prefix in  availableLeptonLeadSubleadPlots:
     for item_channel in channel:
          for item_whichLepton in whichLepton:
               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_data_"+item_whichLepton+")")

               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton))


               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")



               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton))


               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")

               
               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_weighted_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
               
               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton))


               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton))


               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton+" = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_MC_"+item_whichLepton+", "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_data_"+item_whichLepton+")")
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+item_channel+"_Channel_"+item_whichLepton)))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


## All lepton plots or lepton 2D 

for prefix in availablePlots:
  for item_channel in channel:
#     for item_Lepton in whichLeptonConfigured:
          for ifweighted in ifweightedOrNot:

               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel_data)") 

               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel"))
               
               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               
               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_data)")
               
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print(str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_MC")
               print(eval(str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_MC"))
               print(str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_data")
               print(eval(str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_data"))

               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel_data)")
               
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(ifweighted)+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               
               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyFromTupel("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel_data)")
               
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_data)")

               
               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
               print(str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_MC")
               print(eval(str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_MC"))
               print(str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_data")
               print(eval(str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel_data"))

               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel_data)")


               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")

               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel")
               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel_data)")


               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")

               print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")

               exec("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel = getScaleFactorAndUncertaintyPerBin("+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel_MC, "+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel_data)")
               

               print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"))

               with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")


#for item_channel in channel:
#     for item_Lepton in whichLepton:

#          Total_Histogramm_Eff_OnlyDoubleMuon =
#          Total_Histogramm_Eff_OnlyDoubleMuon_weighted = 
#          Total_Histogramm_Eff_DoubleMuonOrSingleLeptonTrigger = 
#          Total_Histogramm_Eff_DoubleMuonOrSingleLeptonTrigger_weighted = 
          

#          Total_Histogramm_ScaleFactor_OnlyDoubleMuon =
#          Total_Histogramm_ScaleFactor_OnlyDoubleMuon_weighted =
#          Total_Histogramm_ScaleFactor_DoubleMuonOrSingleLeptonTrigger =
#          Total_Histogramm_ScaleFactor_DoubleMuonOrSingleLeptonTrigger_weighted =
# propOfUncertainty(y, u1, x1, u2, x2):
#          getScaleFactorAndUncertainty
#def getScaleFactorAndUncertainty(uncertaintyMC, efficiencyMC, uncertaintyData, efficiencyData):



#
#  _                    _   _  _        _____           _                       _   _        _    _                     _        _       _   _           
# | |                  | | | || |  _   / ____|         | |                     | | (_)      | |  | |                   | |      (_)     | | (_)          
# | |     _____   _____| | | || |_(_) | (___  _   _ ___| |_ ___ _ __ ___   __ _| |_ _  ___  | |  | |_ __   ___ ___ _ __| |_ __ _ _ _ __ | |_ _  ___  ___ 
# | |    / _ \ \ / / _ \ | |__   _|    \___ \| | | / __| __/ _ \ '_ ` _ \ / _` | __| |/ __| | |  | | '_ \ / __/ _ \ '__| __/ _` | | '_ \| __| |/ _ \/ __|
# | |___|  __/\ V /  __/ |    | |  _   ____) | |_| \__ \ ||  __/ | | | | | (_| | |_| | (__  | |__| | | | | (_|  __/ |  | || (_| | | | | | |_| |  __/\__ \
# |______\___| \_/ \___|_|    |_| (_) |_____/ \__, |___/\__\___|_| |_| |_|\__,_|\__|_|\___|  \____/|_| |_|\___\___|_|   \__\__,_|_|_| |_|\__|_|\___||___/
#                                              __/ |                                                                                                     
#                                             |___/                                                                                                      
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::..:::..


print("ee_Channel_MC_trigger_correlation")
print(ee_Channel_MC_trigger_correlation)
print("emu_Channel_MC_trigger_correlation")
print(emu_Channel_MC_trigger_correlation)
print("mumu_Channel_MC_trigger_correlation")
print(mumu_Channel_MC_trigger_correlation)

print("ee_Channel_data_trigger_correlation")
print(ee_Channel_data_trigger_correlation)
print("emu_Channel_data_trigger_correlation")
print(emu_Channel_data_trigger_correlation)
print("mumu_Channel_data_trigger_correlation")
print(mumu_Channel_data_trigger_correlation)

arrayCorrelations = [ee_Channel_MC_trigger_correlation, emu_Channel_MC_trigger_correlation, mumu_Channel_MC_trigger_correlation]
#ee_Channel_MC_trigger_correlation.append([events, ee_OnlyTriggerDileptonWithSingleLepton_fired, ee_MET_AND_TriggerDileptonWithSingleLepton_fired, ee_Only_MET_Trigger_fired])

for correlation in arrayCorrelations:
     a=float(correlation[0][0])
     b=float(correlation[0][1])
     c=float(correlation[0][2])
     d=float(correlation[0][3])
     
     e=float(b)/a
     f=float(c)/a
     g=float(d)/a

     h = e*g
     j = h/f
     if correlation == arrayCorrelations[0]:
          ee_alpha = j
     elif correlation == arrayCorrelations[1]:
          emu_alpha = j
     elif correlation == arrayCorrelations[2]:
          mumu_alpha = j
     print(j)


with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                  nameRootoutput.write("ee_alpha")
                  nameRootoutput.write("\n")
                  nameRootoutput.write(str(ee_alpha))
                  nameRootoutput.write("\n")
                  nameRootoutput.write("\n")
                  nameRootoutput.write("emu_alpha")
                  nameRootoutput.write("\n")
                  nameRootoutput.write(str(emu_alpha))
                  nameRootoutput.write("\n")
                  nameRootoutput.write("\n")
                  nameRootoutput.write("mumualpha")
                  nameRootoutput.write("\n")
                  nameRootoutput.write(str(mumu_alpha))
                  nameRootoutput.write("\n")
                  nameRootoutput.write("\n")

print("ee_alpha")
print(ee_alpha)

print("emu_alpha")
print(emu_alpha) 

print("mumu_alpha")
print(mumu_alpha) 


availablePlotsForSystematicUncertaintyDetermination = ["Eta_2D_syst_njet_larger4_", "Eta_2D_syst_njet_smallerEqual4_", "Eta_2D_syst_njet_larger4_", "Eta_2D_syst_nvertex_larger30_", "Eta_2D_syst_nvertex_smallerEqual30", "Eta_2D_syst_met_larger100_", "Eta_2D_syst_met_smallerEqual100" ]

availablePlotsNominal = "Eta_2D_"


for item_channel in channel:
       #eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel")
       #eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel")
       arrayOfTotalScalFactorsForSystematics = "["
       arrayOfScaleFactorsUsedForSystemtics = "["
       for prefix in availablePlotsForSystematicUncertaintyDetermination:
            arrayOfTotalScalFactorsForSystematics += "ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel, "
            arrayOfScaleFactorsUsedForSystemtics += "ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel, "
       arrayOfTotalScalFactorsForSystematics += "]"
       arrayOfScaleFactorsUsedForSystemtics += "]"
       
       print(arrayOfScaleFactorsUsedForSystemtics)
       prefix = availablePlotsNominal
       if item_channel == "ee":
            ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = eval("getScaleFactorSystematicsPerBin("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel, "+arrayOfScaleFactorsUsedForSystemtics+")")
            ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL = eval("getScaleFactorSystematics("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel, "+arrayOfTotalScalFactorsForSystematics+")")
       elif item_channel == "emu":
            emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = eval("getScaleFactorSystematicsPerBin("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel, "+arrayOfScaleFactorsUsedForSystemtics+")")
            emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL = eval("getScaleFactorSystematics("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel, "+arrayOfTotalScalFactorsForSystematics+")")
       elif item_channel == "mumu":
            mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = eval("getScaleFactorSystematicsPerBin("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel, "+arrayOfScaleFactorsUsedForSystemtics+")")
            mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL = eval("getScaleFactorSystematics("+"ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel, "+arrayOfTotalScalFactorsForSystematics+")")
prefix = availablePlotsNominal
for item_channel in channel:
     print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel")
     print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel") )
     print("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel")
     print(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"))

print("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
print("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
print("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)            

print("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
print(ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL)
print("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
print(emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL)
print("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
print(mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL)

with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")

item_channel = "ee"
ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = getScaleFactorSystematicAndStatisticalUncertaintyPerBin(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"), ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin, ee_alpha)
ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired = getScaleFactorSystematicAndStatisticalUncertainty(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel"), ee_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL, ee_alpha)
item_channel = "emu"
emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = getScaleFactorSystematicAndStatisticalUncertaintyPerBin(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"), emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin, emu_alpha)
emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired = getScaleFactorSystematicAndStatisticalUncertainty(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel"), emu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL, emu_alpha)

item_channel = "mumu"
mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin = getScaleFactorSystematicAndStatisticalUncertaintyPerBin(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"), mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin, mumu_alpha)
mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired = getScaleFactorSystematicAndStatisticalUncertainty(eval("ScaleFactorUncertainty_"+str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL_"+str(item_channel)+"_Channel"), mumu_getScaleFactorSytematicUncertaintiesHLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_TOTAL, mumu_alpha)


print("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
print("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
print(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired)       
       
print("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
print("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
print(emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired)

print("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
print(mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
print("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
print(mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired)

with open(nameRootoutputCsVName, 'a') as nameRootoutput:
                 nameRootoutput.write("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")
                 nameRootoutput.write("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")
                 nameRootoutput.write("\n")
                 nameRootoutput.write(str(eval("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired")))
                 nameRootoutput.write("\n")
                 nameRootoutput.write("\n")



#  _                    _   _____       _____  _       _       
# | |                  | | | ____|  _  |  __ \| |     | |      
# | |     _____   _____| | | |__   (_) | |__) | | ___ | |_ ___ 
# | |    / _ \ \ / / _ \ | |___ \      |  ___/| |/ _ \| __/ __|
# | |___|  __/\ V /  __/ |  ___) |  _  | |    | | (_) | |_\__ \
# |______\___| \_/ \___|_| |____/  (_) |_|    |_|\___/ \__|___/
#...............................................................                                                              

outputFileName = nameRootoutputDirectory+ ".root"
rootdirectory =directory + outputFileName
newFile = r.TFile.Open(rootdirectory, "RECREATE")

complete_lepton_2D_etaByEta_binning = [0.]

complete_lepton_2D_etaByEta_binning += lepton_2D_etaByEta_binning

binning_2dEta = np.array(complete_lepton_2D_etaByEta_binning)


#binning_2dEta = np.array([0, 0.3, 0.6, 1.2, 1.7, 2.4])
hnamePre, htitlePre = 'h_ee_lepton2DEta', 'lepton2DEta'
titleX, titleY = '#eta_{leadinglepton}', '#eta_{subleadinglepton}'
h_ee_lepton2DEta = r.TH2F('h_ee_lepton2DEta', htitlePre+';'+titleX+';'+titleY,  len(binning_2dEta)-1, binning_2dEta, len(binning_2dEta)-1, binning_2dEta)
#for i in range(0, len(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)):
hnamePre, htitlePre = 'h_mumu_lepton2DEta', 'lepton2DEta'
h_mumu_lepton2DEta = r.TH2F('h_mumu_lepton2DEta', htitlePre+';'+titleX+';'+titleY,  len(binning_2dEta)-1, binning_2dEta, len(binning_2dEta)-1, binning_2dEta)
hnamePre, htitlePre = 'h_emu_lepton2DEta', 'lepton2DEta'
h_emu_lepton2DEta = r.TH2F('h_emu_lepton2DEta', htitlePre+';'+titleX+';'+titleY,  len(binning_2dEta)-1, binning_2dEta, len(binning_2dEta)-1, binning_2dEta)

## Interface for ROOT users 
c = 0
#print("ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
#print(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
for i in range(1, len(binning_2dEta)):
     
     for j in range(1, len(binning_2dEta)):
#          print(str(i)+','+str(j))
#          print(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
#          print(ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          h_ee_lepton2DEta.SetBinContent(i, j, ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
          h_ee_lepton2DEta.SetBinError(i, j, ee_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          c+=1
c = r.TCanvas('ee_2D_Eta_lepton', 'leading and subleading lepton eta')
h_ee_lepton2DEta.Draw('colz texte')
for ext in ['png'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
h_ee_lepton2DEta.Write()


c = 0
#print("mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
#print(mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
for i in range(1, len(binning_2dEta)):
     for j in range(1, len(binning_2dEta)):
#          print(str(i)+','+str(j))
#          print(mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
#          print(mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          h_mumu_lepton2DEta.SetBinContent(i, j, mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
          h_mumu_lepton2DEta.SetBinError(i, j, mumu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          c+=1
c = r.TCanvas('mumu_2D_Eta_lepton', 'leading and subleading lepton eta')
h_mumu_lepton2DEta.Draw('colz texte')
for ext in ['png'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
h_mumu_lepton2DEta.Write() 

c = 0
#print("emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin")
#print(emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin)
for i in range(1, len(binning_2dEta)):
     for j in range(1, len(binning_2dEta)):
#          print(str(i)+','+str(j))
#          print(emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
#          print(emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          h_emu_lepton2DEta.SetBinContent(i, j, emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][0])
          h_emu_lepton2DEta.SetBinError(i, j, emu_ScaleFactorPlusSystematicAndStatisticalUncertainties_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin[c][1])
          c+=1
c = r.TCanvas('emu_2D_Eta_lepton', 'leading and subleading lepton eta')
h_emu_lepton2DEta.Draw('colz texte')
for ext in ['png'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
for ext in ['pdf'] : c.SaveAs(directory + c.GetName()+'_'+nameRootoutputDirectory+'.'+ext)
h_emu_lepton2DEta.Write()


#newTree = r.TTree( 'NTuple', 'Just A Tree without Clone' )          
newFile.Close()



#tripleGraphPlot(binning_array, plotName, graph1, graph2, graph3)

#Eta_getEfficNjets_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel_dataiency_AllCounts_HLT_DoubleMuon_fired_perBin_mumu_Channel_data_1st_Lepton
#Eta_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_mumu_Channel_MC_1st_Lepton
#ScaleFactorUncertainty_Eta_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_mumu_Channel_1st_Lepton

#Njets_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel_data
#Njets_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel_MC


#ScaleFactorUncertainty_Njets_getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_ee_Channel


met_binning = [0., 20., 40., 60., 80., 100., 125., 150., 175., 200.]
nJets_binning =  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
nVertex_binning  = [float(i) for i in range(0,65,5)]

available1DPlots = ["MET_", "Njets_", "NVertex_"]
lepton_Pt_binning = [20., 50., 80., 120., 200.] 
muon_eta_binning = [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.5, -0.2, 0.0, 0.2, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4]
electron_eta_binning =  [-2.4, -2.1, -1.566, -1.4442, -1.0, -0.6, -0.3, -0.1, 0.1, 0.3, 0.6, 1.0, 1.4442, 1.566, 2.1, 2.4]

for item_channel in channel:
       for prefix in available1DPlots:
            if "MET" in prefix:
                 binning_array = met_binning
            elif "Njets" in prefix:
                 binning_array = nJets_binning
            elif "NVertex" in prefix:
                 binning_array = nVertex_binning
            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data"), eval(whichEffs+"_MC"), eval(whichScale))
            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data"), eval(whichEffs+"_MC"), eval(whichScale))
            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data"), eval(whichEffs+"_MC"), eval(whichScale))
            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data"), eval(whichEffs+"_MC"), eval(whichScale))

available1DPlots = ["Pt_", "Eta_"]


#["MET_", "Njets_", "NVertex_",
 
for item_channel in channel:
       for prefix in available1DPlots:
            if "Pt" in prefix:
                 binning_array = lepton_Pt_binning
            elif "Eta" in prefix:
                 if "ee" in item_channel:
                      binning_array = electron_eta_binning
                 elif "mu" in item_channel:
                      binning_array = muon_eta_binning
            
            
            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
          #  if "ee" in item_channel:
#                 exec(whichEffs+"_data_1st_Lepton = "+whichEffs+"_data_1st_Lepton.insert(3, (0.,0.))")
 #                print(eval(whichEffs+"_data_1st_Lepton"))
  #               exec(whichEffs+"_MC_1st_Lepton = "+whichEffs+"_MC_1st_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_MC_1st_Lepton"))
    #             exec(whichScale+"_1st_Lepton = "+whichScale+"_1st_Lepton.insert(3, (0., 0.))")
     #            print(eval(whichScale+"_1st_Lepton"))
                 
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_1st_Lepton"), eval(whichEffs+"_MC_1st_Lepton"), eval(whichScale+"_1st_Lepton"))
            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
           # if "ee" in item_channel:
      #           exec(whichEffs+"_data_1st_Lepton = "+whichEffs+"_data_1st_Lepton.insert(3, (0.,0.))")
       #          print(eval(whichEffs+"_data_1st_Lepton"))
        #         exec(whichEffs+"_MC_1st_Lepton = "+whichEffs+"_MC_1st_Lepton.insert(3, (0.,0.))")
         #        print(eval(whichEffs+"_MC_1st_Lepton"))
          #       exec(whichScale+"_1st_Lepton = "+whichScale+"_1st_Lepton.insert(3, (0., 0.))")
           #      print(eval(whichScale+"_1st_Lepton"))

            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_1st_Lepton"), eval(whichEffs+"_MC_1st_Lepton"), eval(whichScale+"_1st_Lepton")) 
            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            #if "ee" in item_channel:
             #    exec(whichEffs+"_data_1st_Lepton = "+whichEffs+"_data_1st_Lepton.insert(3, (0.,0.))")
              ##   print(eval(whichEffs+"_data_1st_Lepton"))
                # exec(whichEffs+"_MC_1st_Lepton = "+whichEffs+"_MC_1st_Lepton.insert(3, (0.,0.))")
                # print(eval(whichEffs+"_MC_1st_Lepton"))
        #         exec(whichScale+"_1st_Lepton = "+whichScale+"_1st_Lepton.insert(3, (0., 0.))")
         #        print(eval(whichScale+"_1st_Lepton"))

            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_1st_Lepton"), eval(whichEffs+"_MC_1st_Lepton"), eval(whichScale+"_1st_Lepton")) 
            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
 #           if "ee" in item_channel:
  #               exec(whichEffs+"_data_1st_Lepton = "+whichEffs+"_data_1st_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_data_1st_Lepton"))
    ##             exec(whichEffs+"_MC_1st_Lepton = "+whichEffs+"_MC_1st_Lepton.insert(3, (0.,0.))")
      #           print(eval(whichEffs+"_MC_1st_Lepton"))
       #          exec(whichScale+"_1st_Lepton = "+whichScale+"_1st_Lepton.insert(3, (0., 0.))")
        #         print(eval(whichScale+"_1st_Lepton"))

            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_1st_Lepton"), eval(whichEffs+"_MC_1st_Lepton"), eval(whichScale+"_1st_Lepton"))


            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
 #           if "ee" in item_channel:
  #               exec(whichEffs+"_data_2nd_Lepton = "+whichEffs+"_data_2nd_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_data_2nd_Lepton"))
    #             exec(whichEffs+"_MC_2nd_Lepton = "+whichEffs+"_MC_2nd_Lepton.insert(3, (0.,0.))")
     #            print(eval(whichEffs+"_MC_2nd_Lepton"))
      #           exec(whichScale+"_2nd_Lepton = "+whichScale+"_2nd_Lepton.insert(3, (0., 0.))")
       #          print(eval(whichScale+"_2nd_Lepton"))

            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_2nd_Lepton"), eval(whichEffs+"_MC_2nd_Lepton"), eval(whichScale+"_2nd_Lepton"))
            whichEffs = str(prefix)+"getEfficiency_AllCounts_HLT_DoubleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
 #           if "ee" in item_channel:
  #               exec(whichEffs+"_data_2nd_Lepton = "+whichEffs+"_data_2nd_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_data_2nd_Lepton"))
    #             exec(whichEffs+"_MC_2nd_Lepton = "+whichEffs+"_MC_2nd_Lepton.insert(3, (0.,0.))")
     #            print(eval(whichEffs+"_MC_2nd_Lepton"))
      #           exec(whichScale+"_2nd_Lepton = "+whichScale+"_2nd_Lepton.insert(3, (0., 0.))")
       #          print(eval(whichScale+"_2nd_Lepton"))

            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_2nd_Lepton"), eval(whichEffs+"_MC_2nd_Lepton"), eval(whichScale+"_2nd_Lepton"))
            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
 #           if "ee" in item_channel:
  #               exec(whichEffs+"_data_2nd_Lepton = "+whichEffs+"_data_2nd_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_data_2nd_Lepton"))
    #             exec(whichEffs+"_MC_2nd_Lepton = "+whichEffs+"_MC_2nd_Lepton.insert(3, (0.,0.))")
     #            print(eval(whichEffs+"_MC_2nd_Lepton"))
       #          exec(whichScale+"_2nd_Lepton = "+whichScale+"_2nd_Lepton.insert(3, (0., 0.))")
      #           print(eval(whichScale+"_2nd_Lepton"))


            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_2nd_Lepton"), eval(whichEffs+"_MC_2nd_Lepton"), eval(whichScale+"_2nd_Lepton"))
 #           if "ee" in item_channel:
  #               exec(whichEffs+"_data_2nd_Lepton = "+whichEffs+"_data_2nd_Lepton.insert(3, (0.,0.))")
   #              print(eval(whichEffs+"_data_2nd_Lepton"))
    #             exec(whichEffs+"_MC_2nd_Lepton = "+whichEffs+"_MC_2nd_Lepton.insert(3, (0.,0.))")
     #            print(eval(whichEffs+"_MC_2nd_Lepton"))
      #           exec(whichScale+"_2nd_Lepton = "+whichScale+"_2nd_Lepton.insert(3, (0., 0.))")
       #          print(eval(whichScale+"_2nd_Lepton"))

            whichEffs = str(prefix)+"getEfficiency_HLT_DoubleMuonPlusAdditional_HLT_SingleMuon_fired_weighted_perBin_"+str(item_channel)+"_Channel"
            whichScale = "ScaleFactorUncertainty_"+whichEffs
            tripleGraphPlot(directory, binning_array, str(whichEffs), eval(whichEffs+"_data_2nd_Lepton"), eval(whichEffs+"_MC_2nd_Lepton"), eval(whichScale+"_2nd_Lepton"))

sys.exit()
