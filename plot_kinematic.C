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

using namespace std;

void plot_kinematic()
{
  TFile *f1 = new TFile("skimmedFiles/all_2.root","read");
  //TFile *f1 = new TFile("skimmedFiles/yskimPA_OpSign_2017329147_unIdentified.root","read");
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


  TH1D *hPtD = new TH1D("hPtD",";p_{T}^{#mu#mu} (GeV/c);",80,0,30);
  TH1D *hPtS = new TH1D("hPtS",";p_{T}^{#mu} (GeV/c);",80,0,30);
  TH1D *hrapDL = new TH1D("hrapDL",";y_{lab}",50,-2.4,2.4);
  TH1D *hphiDL = new TH1D("hphiDL",";#phi^{#mu#mu} (rad);",50,-3.15,3.15);
  TH1D *hphiSL = new TH1D("hphiSL",";#phi^{#mu} (rad);",50,-3.15,3.15);
  TH1D *hetaDL = new TH1D("hetaDL",";#eta^{#mu#mu};",50,-2.4,2.4);
  TH1D *hetaSL = new TH1D("hetaSL",";#eta^{#mu};",50,-2.4,2.4);
  TH1D *hrapSL = new TH1D("hrapSL","",50,-2.5,2.5);
  TH2D *hPtRap = new TH2D("hPtRap",";y_{lab};p_{T} (GeV/c)",100,-2.4,2.4,100,0,9);
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
    if(dm.mass>=8 && dm.mass<=14 && dm.pt1>=4 && dm.pt2>=4 && dm.pt<=30 && dm.pt>=0 && TMath::Abs(dm.y)<2.4)
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

  hPtD->Scale(1./hPtD->GetEntries());
  hPtS->Scale(1./hPtS->GetEntries());
  hrapDL->Scale(1./hrapDL->GetEntries());
  hphiDL->Scale(1./hphiDL->GetEntries());
  hphiSL->Scale(1./hphiSL->GetEntries());
  hetaDL->Scale(1./hetaDL->GetEntries());
  hetaSL->Scale(1./hetaSL->GetEntries());

  hrapDL->GetYaxis()->SetRangeUser(0,0.05);
  hphiDL->GetYaxis()->SetRangeUser(0,0.03);
  hphiSL->GetYaxis()->SetRangeUser(0,0.03);
 
  cout << "entry : " << hPtD->GetEntries() << endl; 

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",1300,700);
  c1->Divide(4,2);
  c1->cd(1);
  gPad->SetLogy();
  gPad->SetBottomMargin(0.15);
  hPtD->Draw("pe");
  c1->cd(2);
  gPad->SetLogy();
  gPad->SetBottomMargin(0.15);
  hPtS->Draw("pe");
  c1->cd(3);
  gPad->SetBottomMargin(0.15);
  hrapDL->Draw("pe");
  c1->cd(4);
  gPad->SetBottomMargin(0.15);
  hphiDL->Draw("pe");
  c1->cd(5);
  gPad->SetBottomMargin(0.15);
  hphiSL->Draw("pe");
  c1->cd(6);
  gPad->SetBottomMargin(0.15);
  hPtRap->Draw("colz");
//  c1->cd(7);
//  gPad->SetBottomMargin(0.15);
//  hETPlMin->Draw("colz");
  c1->cd(7);
  gPad->SetBottomMargin(0.15);
  hETtrack->Draw("colz");

  c1->SaveAs("pPb5TeV_kinematicPlot.pdf");
}

