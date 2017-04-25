#include <iostream>
#include <sstream>
#include "TROOT.h"
#include "commonUtility.h"
#include "SONGKYO.h"
#include "cutsAndBin.h"
#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "CMS_lumi_raaCent.C"
#include "tdrstyle_raaCent.C"
#include "SONGKYO.h"

using namespace std;

void plot_kinematic_test()
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 3; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  TFile *f1 = new TFile("skimmedFiles/all_2.root","read");
  //TFile *f1 = new TFile("skimmedFiles/yskimPA_OpSign_2017329147_unIdentified.root","read");
  //TFile *f1 = new TFile("skimmedFiles/yskimPA_OpSign_20173311049_unIdentified_tight.root","read");
  TTree *tr = (TTree*) f1->Get("mm");
  TTree *hftr = (TTree*) f1->Get("hf");
  DiMuon dm;
  Int_t Ntracks;
  Float_t SumET_HFplusEta4,SumET_HFminusEta4,SumET_HFplus,SumET_HFminus,SumET_HF;
  tr->SetBranchAddress("mm",&dm.run);
  hftr->SetBranchAddress("SumET_HFplusEta4",&SumET_HFplusEta4);
  hftr->SetBranchAddress("SumET_HFminusEta4",&SumET_HFminusEta4);
  hftr->SetBranchAddress("SumET_HFplus",&SumET_HFplus);
  hftr->SetBranchAddress("SumET_HFminus",&SumET_HFminus);
  hftr->SetBranchAddress("SumET_HF",&SumET_HF);
  hftr->SetBranchAddress("Ntracks",&Ntracks);


  TH1D *hPtD = new TH1D("hPtD",";p_{T}^{#mu#mu} (GeV/c);Counts",220,0,130);
  TH1D *hPtS = new TH1D("hPtS",";p_{T}^{#mu} (GeV/c);",80,0,30);
  TH1D *hrapDL = new TH1D("hrapDL",";y_{lab}",50,-2.4,2.4);
  TH1D *hphiDL = new TH1D("hphiDL",";#phi^{#mu#mu} (rad);",50,-3.15,3.15);
  TH1D *hphiSL = new TH1D("hphiSL",";#phi^{#mu} (rad);",50,-3.15,3.15);
  TH1D *hetaDL = new TH1D("hetaDL",";#eta^{#mu#mu};",50,-2.4,2.4);
  TH1D *hetaSL = new TH1D("hetaSL",";#eta^{#mu};",50,-2.4,2.4);
  TH1D *hrapSL = new TH1D("hrapSL","",50,-2.5,2.5);
  TH2D *hPtRap = new TH2D("hPtRap",";y_{lab};p_{T} (GeV/c)",100,-2.4,2.4,100,0,9);
  TH2D *hPtMass = new TH2D("hPtMass",";m_{#mu^{+}#mu^{-}} (GeV/c^{2});p_{T}^{#mu#mu} (GeV/c)",100,8,14,140,0,130);
  hPtRap->GetXaxis()->SetTitleSize(0.07);
  hPtRap->GetYaxis()->SetTitleSize(0.07);
  hPtRap->GetXaxis()->SetTitleOffset(0.65);
  hPtRap->GetYaxis()->SetTitleOffset(0.65);

  TH2D *hETPlMin = new TH2D("hETPlMin",";SumET_HFplusEta4;SumET_HFminusEta4",200,0,100,200,0,100);
  TH2D *hETtrack = new TH2D("hETtrack",";SumET_HF;Ntracks",200,0,400,200,0,400);
  hETPlMin->GetXaxis()->SetTitleSize(0.07);
  hETPlMin->GetYaxis()->SetTitleSize(0.07);
  hETPlMin->GetXaxis()->SetTitleOffset(0.65);
  hETPlMin->GetYaxis()->SetTitleOffset(0.65);
  hETtrack->GetXaxis()->SetTitleSize(0.07);
  hETtrack->GetYaxis()->SetTitleSize(0.07);
  hETtrack->GetXaxis()->SetTitleOffset(0.65);
  hETtrack->GetYaxis()->SetTitleOffset(0.65);
  

  hftr->Draw("SumET_HFplusEta4:SumET_HFminusEta4>>hETPlMin","","colz");
  hftr->Draw("SumET_HF:Ntracks>>hETtrack","","colz");

  for(int i=0;i<tr->GetEntries();i++)
  {
    tr->GetEntry(i);
    if(dm.mass>=8 && dm.mass<=14 && dm.pt1>=4 && dm.pt2>=4 && dm.pt>=0 && TMath::Abs(dm.y)<2.4)
    {
      hPtD->Fill(dm.pt);
      hPtS->Fill(dm.pt1);
      hPtS->Fill(dm.pt2);
      hrapDL->Fill(dm.y);
      hphiDL->Fill(dm.phi);
      hphiSL->Fill(dm.phi1);
      hphiSL->Fill(dm.phi2);
      hetaDL->Fill(dm.eta);
      hetaSL->Fill(dm.eta1);
      hetaSL->Fill(dm.eta2);
      hPtRap->Fill(dm.y,dm.pt);
      hPtMass->Fill(dm.mass,dm.pt);
    }
  }

  

  SetHistStyle(hPtD,8,0); 
  SetHistStyle(hPtS,8,0); 
  SetHistStyle(hphiDL,8,0); 
  SetHistStyle(hphiSL,8,0); 
  SetHistStyle(hrapDL,8,0); 
  SetHistStyle(hetaDL,8,0); 
  SetHistStyle(hetaSL,8,0); 
  hPtD->GetXaxis()->SetTitleFont(132);
  hPtD->GetYaxis()->SetTitleFont(132);
  hPtS->GetXaxis()->SetTitleFont(132);
  hPtS->GetYaxis()->SetTitleFont(132);
  hrapDL->GetXaxis()->SetTitleFont(132);
  hrapDL->GetYaxis()->SetTitleFont(132);
  hphiDL->GetXaxis()->SetTitleFont(132);
  hphiDL->GetYaxis()->SetTitleFont(132);
  hphiSL->GetXaxis()->SetTitleFont(132);
  hphiSL->GetYaxis()->SetTitleFont(132);
  hetaDL->GetXaxis()->SetTitleFont(132);
  hetaDL->GetYaxis()->SetTitleFont(132);
  hetaSL->GetXaxis()->SetTitleFont(132);
  hetaSL->GetYaxis()->SetTitleFont(132);

  hPtRap->GetXaxis()->CenterTitle();
  hPtRap->GetYaxis()->CenterTitle();
  hPtMass->GetXaxis()->CenterTitle();
  hPtMass->GetYaxis()->CenterTitle();

  //hPtD->Scale(1./hPtD->GetEntries());
  hPtS->Scale(1./hPtS->GetEntries());
  hrapDL->Scale(1./hrapDL->GetEntries());
  hphiDL->Scale(1./hphiDL->GetEntries());
  hphiSL->Scale(1./hphiSL->GetEntries());
  hetaDL->Scale(1./hetaDL->GetEntries());
  hetaSL->Scale(1./hetaSL->GetEntries());

  hrapDL->GetYaxis()->SetRangeUser(0,0.05);
  hphiDL->GetYaxis()->SetRangeUser(0,0.03);
  hphiSL->GetYaxis()->SetRangeUser(0,0.03);
 
   
  hPtD->GetXaxis()->SetLabelSize(0.045);
  hPtD->GetYaxis()->SetLabelSize(0.045);
  hPtD->GetXaxis()->SetLabelOffset(0.011);
  hPtD->GetYaxis()->SetLabelOffset(0.011);
  hPtD-> GetXaxis()->SetTitleSize(0.05);
  hPtD-> GetYaxis()->SetTitleSize(0.05);
  hPtD->GetXaxis()->SetTitleOffset(1.11);
  hPtD->GetYaxis()->SetTitleOffset(1.38);

  hPtMass->GetXaxis()->SetLabelSize(0.045);
  hPtMass->GetYaxis()->SetLabelSize(0.045);
  hPtMass->GetXaxis()->SetLabelOffset(0.011);
  hPtMass->GetYaxis()->SetLabelOffset(0.011);
  hPtMass->GetXaxis()->SetTitleSize(0.05);
  hPtMass->GetYaxis()->SetTitleSize(0.05);
  hPtMass->GetXaxis()->SetTitleOffset(1.11);
  hPtMass->GetYaxis()->SetTitleOffset(1.38);

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  gPad->SetTopMargin(0.06);
  c1->SetTicks(1,1);
  gPad->SetLogy();
  hPtD->Draw("pe");
  double text_x = 0.5;
  double text_y = 0.75;
  double text_y_diff = 0.05;
  double text_size = 18;
  drawText("8<m_{#mu^{+}#mu^{-}}<14 GeV/c^{2}",text_x,text_y,1,text_size);
  drawText("p_{T}^{#mu}>4 GeV/c",text_x,text_y-text_y_diff,1,text_size);
  drawText("|y_{lab}^{#mu#mu}|<2.4",text_x,text_y-2*text_y_diff,1,text_size);
  //  c1->cd(7);
//  gPad->SetBottomMargin(0.15);
//  hETPlMin->Draw("colz");
  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  gPad->SetTopMargin(0.06);
  gPad->SetRightMargin(0.165);
  hPtMass->Draw("colz");
  //drawText("60<m_{#mu^{+}#mu^{-}}<120 GeV/c^{2}",text_x,text_y,1,text_size);
  drawText("8<m_{#mu^{+}#mu^{-}}<14 GeV/c^{2}",text_x,text_y,1,text_size);
  drawText("p_{T}^{#mu}>4 GeV/c",text_x,text_y-text_y_diff,1,text_size);
  drawText("|y_{lab}^{#mu#mu}|<2.4",text_x,text_y-2*text_y_diff,1,text_size);
  jumSun(pdgMass.Y1S,0,pdgMass.Y1S,130,2,3);
  jumSun(pdgMass.Y2S,0,pdgMass.Y2S,130,2,3);
  jumSun(pdgMass.Y3S,0,pdgMass.Y3S,130,2,3);
  drawText("m_{#Upsilon(1S)}",0.268,0.893,2,18);
  drawText("m_{#Upsilon(2S)}",0.346,0.893,2,18);
  drawText("m_{#Upsilon(3S)}",0.438,0.893,2,18);

  CMS_lumi_raaCent( c1, iPeriod, iPos );
  CMS_lumi_raaCent( c2, iPeriod, iPos );
  /*gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.145);
  gPad->SetRightMargin(0.165);
  gPad->SetTopMargin(0.05);
  */
  c1->Update();
  c2->Update();
  c1->SaveAs("DiMuPt_extend.pdf");
  c2->SaveAs("DiMu_PtvsMass.pdf");
}

