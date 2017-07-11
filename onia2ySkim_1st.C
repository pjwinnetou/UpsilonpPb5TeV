#include <ctime>

#include "cutsAndBin.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include <TLorentzVector.h>
#include "TriggerManipulation.h" 
#include "commonUtility.h"
static const long MAXTREESIZE = 10000000000;

TString getDayAndTime();

void onia2ySkim_1st( int nevt = -1,
		 int fileID = kPADATA,
		 int trigId=kL1DoubleMuOpen2016, 
		 int epSelection = kEPOppositeHF, 
		 bool saveTracks=false, 
		 TString skimVersion="unIdentified", 
		 bool DiMuSign = false
		 ) 
{
  using namespace std;
  
  bool isMC = false; 
  if ( (fileID == kPPMC) || (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) || (fileID == kAAMC) || (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )
    isMC = true;

  TChain *mytree = new TChain("myTree");

  TString fname;  TString fname1;  TString fname2;  TString fname3;  TString fname4;   TString fname5;
  
  if (fileID == kPPDATA) {
    fname = "/home/samba.old/UpsilonAnalysis/tempfiles/Upsilon_pPb/data/RD2013_pp.root";
    mytree->Add(fname.Data());
  }  
  else if (fileID == kPADATA) {
    fname = "/home/samba.old/UpsilonAnalysis/tempfiles/Upsilon_pPb/data/RD2013_pa_1st_run_merged.root"; // in lxplus
    mytree->Add(fname.Data());
  }

  cout << endl << "*==*==*==*==*==*==*==* INPUT FILE *==*==*==*==*==*==*==*==*" << endl;

  TFile *f1 = new TFile(Form("%s",fname.Data()),"read");
  if (f1->IsZombie()) { cout << "*** KYO : CANNOT open the root file!! Macro terminated ***" << endl; return;} 
  cout <<"* file ::" << fname << endl;
  
  
  // Same or Opposite sign event
  TString fdimusign;
  if(!DiMuSign) fdimusign = "OpSign";
  else if(DiMuSign) fdimusign = "SSign";

  cout << endl;
  cout << "*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*" << endl;
  cout << " Sign of selecting dimuons : " << fdimusign.Data() << endl;
  cout << endl;

  

  // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
  TFile* newfile;
  if (fileID == kPPDATA) {
    newfile = new TFile(Form("skimmedFiles/yskimPP_L1DoubleMu0PD_%s_%s_%s.root",fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  if (fileID == kPADATA) {
    newfile = new TFile(Form("skimmedFiles/yskimPA1st_%s_%s_%s.root",fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
   
   
  // import the tree to the RooDataSet
  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Ntracks;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_QQ_trig[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[200];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Ntracks;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_QQ_mupl_isHighPurity[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_isHighPurity[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_isHighPurity;   //!
  TBranch        *b_Reco_QQ_mumi_isHighPurity;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_isHighPurity", Reco_QQ_mupl_isHighPurity, &b_Reco_QQ_mupl_isHighPurity);
  mytree->SetBranchAddress("Reco_QQ_mumi_isHighPurity", Reco_QQ_mumi_isHighPurity, &b_Reco_QQ_mumi_isHighPurity);


  
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  //  mytree->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  //  mytree->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_QQ_mupl_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_normChi2_global[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
  TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  mytree->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_QQ_mupl_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nMuValHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
  TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_QQ_mupl_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_StationsMatched[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
  TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  mytree->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_QQ_mupl_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dxy[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dxy;   //!
  TBranch        *b_Reco_QQ_mumi_dxy;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  mytree->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_QQ_mupl_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dz[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[200];   //[Reco_mu_size]
  Float_t         SumET_HFplusEta4;   //
  Float_t         SumET_HFminusEta4;   //
  Float_t         SumET_HFplus;   //
  Float_t         SumET_HFminus;   //
  Float_t         SumET_HF;   //
  TBranch        *b_Reco_QQ_mupl_dz;   //!
  TBranch        *b_Reco_QQ_mumi_dz;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  TBranch        *b_SumET_HFplusEta4;   //!
  TBranch        *b_SumET_HFminusEta4;   //!
  TBranch        *b_SumET_HFplus;   //!
  TBranch        *b_SumET_HFminus;   //!
  TBranch        *b_SumET_HF;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  mytree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  mytree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  mytree->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
  mytree->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
  mytree->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
  mytree->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  Int_t           Reco_QQ_mupl_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_QQ_mupl_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_mu_TMOneStaTight[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
  TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  mytree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_QQ_mupl_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nPixWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_QQ_mupl_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  Int_t           Reco_QQ_mumi_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  Float_t         Reco_QQ_mupl_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  Float_t         Reco_QQ_mumi_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  dimuon tree 
  ////////////////////////////////////////////////////////////////////////
  DiMuon dm;
  TTree *mmTree = new TTree("mm","diMuonPairs");
  TTree *hfTree = new TTree("hf","hf");
  mmTree->SetMaxTreeSize(MAXTREESIZE);
  mmTree->Branch("mm",&dm.run,branchString.Data());
  hfTree->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
  hfTree->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
  hfTree->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
  hfTree->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
  hfTree->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
  hfTree->Branch("Ntracks",&Ntracks,"Ntracks/I");
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","pt weight", 0, 10000,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *evtWeight);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  Event selection tree 
  ////////////////////////////////////////////////////////////////////////
  TH1D* hEvent = new TH1D("hFilter","",20 ,0.5, 21.5);

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;
  
  
  int HLTPASS = 0;
  int DIMUTRIGPASS = 0;
  int DIMUIDPASS = 0;
  int DIMUPASS_all = 0;
  
  int DIMUPASS_CENT1 = 0;
  int DIMUPASS_CENT2 = 0;
  int DIMUPASS_CENT3 = 0;
  int DIMUPASS_CENT4 = 0;
  int DIMUPASS_CENT5 = 0;
  int DIMUPASS_CENT6 = 0;
  int DIMUPASS_CENT7 = 0;
  int DIMUPASS_CENT8 = 0;
  int DIMUPASS_CENT9 = 0;

  double ch_pluseta4,ch_minuseta4,ch_plus,ch_minus;
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;
    
    ///////////////////////////
    ///// Call the values /////
    ///////////////////////////
    mytree->GetEntry(iev);
        
    hEvent->GetXaxis()->SetBinLabel(1,"Events total");     hEvent->Fill(1);  

    if ( !((HLTriggers&1)==1))  
      continue; // trigger selection

    HLTPASS++;

    dm.clear();
    dm.run = runNb;
    dm.lumi = LS ;
    dm.event = eventNb ;
    dm.vz = zVtx;
    dm.cBin = Centrality ;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      hEvent->GetXaxis()->SetBinLabel(4,"Di-muons Total");      hEvent->Fill(4);
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      
      if ( !((Reco_QQ_trig[irqq]&1)==1) )
        continue; // trigger selection
      hEvent->GetXaxis()->SetBinLabel(7,"Di-muons trig");      hEvent->Fill(7);

      DIMUTRIGPASS++;


      if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta() ) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta()) )   )
        continue;      
      hEvent->GetXaxis()->SetBinLabel(5,"Di-muons Accep");      hEvent->Fill(5);

      bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mupl_nPixWMea[irqq] > 1) &&
          (Reco_QQ_mupl_dxy[irqq]<0.3) &&
          (Reco_QQ_mupl_dz[irqq]<30.)	 //&&  (Reco_QQ_mupl_isHighPurity[irqq]==true) 
          ) ; 

      bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mumi_nPixWMea[irqq] > 1) &&
          (Reco_QQ_mumi_dxy[irqq]<0.3) &&
          (Reco_QQ_mumi_dz[irqq]<30.)  //&& (Reco_QQ_mumi_isHighPurity[irqq]==true)
          ) ; 

      bool muplHighPtCut = ( (Reco_QQ_mupl_nMuValHits[irqq]>0) &&
          (Reco_QQ_mupl_StationsMatched[irqq]>1) &&
          (Reco_QQ_mupl_ptErr_global[irqq]/mupl_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mupl_dxy[irqq]<0.2) &&
          (Reco_QQ_mupl_dz[irqq] <0.5) &&
          (Reco_QQ_mupl_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) 
          );
      bool mumiHighPtCut = ( (Reco_QQ_mumi_nMuValHits[irqq]>0) &&
          (Reco_QQ_mumi_StationsMatched[irqq]>1) &&
          (Reco_QQ_mumi_ptErr_global[irqq]/mumi_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mumi_dxy[irqq]<0.2) &&
          (Reco_QQ_mumi_dz[irqq] <0.5) &&
          (Reco_QQ_mumi_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) 
          );


      if ( !(muplSoft && mumiSoft) ) 
        continue;   
      hEvent->GetXaxis()->SetBinLabel(8,"Di-muons mu ID");      hEvent->Fill(8);

      DIMUIDPASS++;

      // vertex probablitiy cut:
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;
      hEvent->GetXaxis()->SetBinLabel(6,"Di-muons Vtx prob.");      hEvent->Fill(6);

      if( !DiMuSign )
      {
        if ( Reco_QQ_sign[irqq] != 0 ) continue;
      }
      else if(DiMuSign)
      {
        if( Reco_QQ_sign[irqq] == 0) continue;
      }

      hEvent->GetXaxis()->SetBinLabel(9,"Di-muoons charge sign");      hEvent->Fill(9);


      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();
      dm.phi    = JP_Reco->Phi();
      dm.y      = -JP_Reco->Rapidity();
      dm.eta      = -JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = -mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = -mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();

      dm.oniaIndex = irqq;
      mmTree->Fill();

      massVar->setVal( (double)dm.mass ) ;
      ptVar->setVal(   (double)dm.pt   ) ;
      yVar->setVal(    (double)dm.y    ) ;
      pt1Var->setVal(  (double)dm.pt1  ) ;
      dphiEp2Var->setVal(   (double)dm.dphiEp2  ) ;
      eta1Var->setVal( (double)dm.eta1 ) ;
      pt2Var->setVal(  (double)dm.pt2  ) ;
      eta2Var->setVal( (double)dm.eta2 ) ;
      ep2Var->setVal( (double)dm.ep2 ) ;
      cBinVar->setVal( (double)dm.cBin ) ;
      evtWeight->setVal( (double)dm.weight ) ;
      dataSet->add( *argSet);
    } // end of dimuon loop
     
      ch_pluseta4 = SumET_HFplusEta4; 
      
      ch_minuseta4 = SumET_HFminusEta4; 
      ch_plus = SumET_HFplus; 
      ch_minus = SumET_HFminus;

      SumET_HFplusEta4 = ch_minuseta4; 
      SumET_HFminusEta4 = ch_pluseta4; 
      SumET_HFplus = ch_minus; 
      SumET_HFminus = ch_plus; 
      hfTree->Fill();
  } //end of event loop
  
  
  cout << endl;
  cout << endl;
  cout << "******************" << endl;
  cout << "**** SUMMARY *****" << endl;
  cout << endl;

  cout << "1. # of events passing HLT : " << HLTPASS << endl;
  cout << "2. # of dimuons passing trigger matching after pass 1 : " << DIMUTRIGPASS << endl;
  cout << "3. # of dimuons passing muon ID cut after pass 1-2 & acceptance cut : " << DIMUIDPASS << endl;
  cout << "Cent All : # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for all central events (0-100%) " << DIMUPASS_all << endl;
  cout << "Cent 1. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (0-5%)   : " << DIMUPASS_CENT1 << endl;
  cout << "Cent 2. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (5-10%)  : " << DIMUPASS_CENT2 << endl;
  cout << "Cent 3. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (10-20%) : " << DIMUPASS_CENT3 << endl;
  cout << "Cent 4. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (20-30%) : " << DIMUPASS_CENT4 << endl;
  cout << "Cent 5. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (30-40%) : " << DIMUPASS_CENT5 << endl;
  cout << "Cent 6. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (40-50%) : " << DIMUPASS_CENT6 << endl;
  cout << "Cent 7. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (50-60%) : " << DIMUPASS_CENT7 << endl;
  cout << "Cent 8. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (60-70%) : " << DIMUPASS_CENT8 << endl;
  cout << "Cent 9. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (70-80%) : " << DIMUPASS_CENT9 << endl;
  cout << "******************" << endl;
  cout << "******************" << endl;
  cout << endl;
  cout << endl;

  dataSet->Write();
  mmTree->Write();  // Don't need to call Write() for trees
  hEvent->Write();
  hfTree->Write();
  newfile->Close();
} 



TString getDayAndTime() 
{ 
  time_t currentTime;
  struct tm *localTime;
  
  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time
  
  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%d%d%d%d%d",Year,Month,Day,Hour,Min);
}
