#include <iostream>
#include "rootFitHeaders.h"
#include "commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include "PsetCollection.h"
#include "CMS_lumi_raaCent.C"
#include "tdrstyle.C"
#include "SONGKYO.h"

using namespace std;
using namespace RooFit;
void tempPlot(
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-2.4, float yHigh=2.4,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool fixParameters=0  )
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  

  using namespace RooFit;
  gStyle->SetEndErrorSize(0);
 
  TString SignalCB = "Double";

  float massLow = 8; 
  float massHigh = 14;

  float yBinLow = -2.4; 
  float yBinHigh = 2.4;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = 60;
  int   nYBin  = 60;
  //int   nMassBin  = (massHigh-massLow)*10;
  TFile* f1;
  if      ( collId == kPPDATA) f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/skimmedFiles/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20164251755_3964bbec2f15f2cf9baa0676644690f40cee27c4.root");
  else if ( collId == kAADATA) f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20164272229_95c28a5bdf107c32b9e54843b8c85939ffe1aa23.root");
  else if ( collId == kPADATA) f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/pPb5TeV/skimmedFiles/all_2.root");
  //else if ( collId == kPADATA) f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/skimmedFiles/yskimPA_OpSign_20173281728_unIdentified.root");
  else if ( collId == kAADATAPeri) f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/skimmedFiles/yskimPbPb_PeripheralPD_Trig-L1DoubleMu0Peripheral_OpSign_EP-OppositeHF_20164272252_95c28a5bdf107c32b9e54843b8c85939ffe1aa23.root");
  else if ( collId == kPPMCUps1S) f1 = new TFile("skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20163251233_2b58ba03c4751c9d10cb9d60303271ddd6e1ba3a.root");
  else if ( collId == kAAMCUps1S) f1 = new TFile("skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20163251233_2b58ba03c4751c9d10cb9d60303271ddd6e1ba3a.root");
 
  if(collId == kAADATAPeri) collId =2; 
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && y>-2.4 && y<2.4",ptLow, ptHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f)", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) )
    kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d) && ( abs(abs(dphiEp2/3.141592)-0.5)>%.3f && abs(abs(dphiEp2/3.141592)-0.5)<%.3f )",cLow, cHigh, dphiEp2Low, dphiEp2High);
  
  
  TTree* tree = (TTree*) f1->Get("mm");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  //RooWorkspace *ws = new RooWorkspace(Form("workspace_%s",kineLabel.Data()));
  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("y")->setRange(yBinLow, yBinHigh);
  ws->var("mass")->Print();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"),MarkerSize(.8));

  RooPlot* myPlot2 = ws->var("y")->frame(nYBin); // bins
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataHist"));

  myPlot->SetFillStyle(4000);
  myPlot->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot->GetYaxis()->SetTitleOffset(1.77);
  myPlot->GetYaxis()->CenterTitle();
  myPlot->GetYaxis()->SetTitleSize(0.05);
  myPlot->GetYaxis()->SetLabelSize(0.045) ;
  myPlot->GetXaxis()->SetRangeUser(8,14);
  myPlot->GetXaxis()->SetTitleSize(0);
  myPlot->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot->GetXaxis()->SetTitleOffset(1.22) ;
  myPlot->GetXaxis()->SetLabelOffset(0.018) ;
  myPlot->GetXaxis()->SetLabelSize(0.045) ;
  myPlot->GetXaxis()->SetTitleSize(0.057) ;
  myPlot->GetXaxis()->CenterTitle();

  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.77);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.05);
  myPlot2->GetYaxis()->SetLabelSize(0.045) ;
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->GetXaxis()->SetTitle("y_{Lab}");
  myPlot2->GetXaxis()->SetTitleOffset(1.22) ;
  myPlot2->GetXaxis()->SetLabelOffset(0.018) ;
  myPlot2->GetXaxis()->SetLabelSize(0.045) ;
  myPlot2->GetXaxis()->SetTitleSize(0.057) ;
  myPlot2->GetXaxis()->CenterTitle();

  /*float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.056;
  float text_size = 19;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 2.5 && ptHigh==5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S) 
  {
      drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
      drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
      drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
    drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  } 
 */ 
//  drawText(Form("Signal Function : %s CB", SignalCB.Data() ), 0.55,0.54,1,14);

  TCanvas* c1 =  new TCanvas("canvas1","My plots1",600,600);
  TCanvas* c2 =  new TCanvas("canvas2","My plots2",600,600);
  c1->cd();
  myPlot->Draw();
  c2->cd();
  myPlot2->Draw();

  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  CMS_lumi_raaCent(c1, 3 ,33);
  CMS_lumi_raaCent(c2, 3 ,33);

  c1->Update();
  c2->Update();
  c1->Draw();
  c2->Draw();

} 
 
