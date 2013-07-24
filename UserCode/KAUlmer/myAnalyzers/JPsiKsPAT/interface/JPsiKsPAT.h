// -*- C++ -*-
//
// Package:    JPsiKsPAT
// Class:      JPsiKsPAT
// 
/**\class JPsiKsPAT JPsiKsPAT.cc myAnalyzers/JPsiKsPAT/src/JPsiKsPAT.cc

 Description: <one line class summary>
Make rootTuple for b->s mu mu reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer
//         Created:  Mon Apr 21 09:53:19 MDT 2008
// $Id: JPsiKsPAT.h,v 1.19 2011/05/08 17:29:49 kaulmer Exp $
//
//

#ifndef _JPsiKsPAT_h
#define _JPsiKsPAT_h

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/interface/JPsiKsPAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class JPsiKsPAT : public edm::EDAnalyzer {
public:
  explicit JPsiKsPAT(const edm::ParameterSet&);
  ~JPsiKsPAT();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool const HasGoodME11(reco::Muon const& muon, double const dxdzCut) const;
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;
  
  // ----------member data ---------------------------
  std::string hlTriggerResults_;
  std::string v0Producer;
  std::string v0Type;
  std::string vtxSample;
  std::string genParticles_;
  std::string v0Collection_;
  std::string muonType;
  std::string muonTypeForPAT;
  bool doMC_;
  bool onlyCount_;
  TTree*      tree_;
  TTree*      treeMC_;
  TTree*      treeMC2_;
  TTree*      treeMC3_;
  TTree*      treeMC4_;
  TTree*      treeMC5_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;

  unsigned int        l1_mu3, l1_2mu3, l1_muOpen, l1_mu0;
  unsigned int        hlt_mu3, hlt_mu5, hlt_mu7, hlt_mu9, hlt_2mu0, hlt_2mu0_quark,  hlt_2mu3, hlt_2mu0L2, hlt_2mu3JPsi, hlt_BJPsiMuMu;
  unsigned int        hlt_mu0trk0, hlt_mu3trk0, hlt_mu0trkmu0, hlt_mu3trkmu0, hlt_mu0trkmu0OST, hlt_mu3trkmu0OST, hlt_mu0trkmu0OST_tight;
  unsigned int        hlt_L1muOpen, hlt_L12muOpen, hlt_L12muOpenTight;
  unsigned int        nB;
  float               priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  std::vector<float>       *priRfVtxX, *priRfVtxY, *priRfVtxZ, *priRfVtxXE, *priRfVtxYE, *priRfVtxZE, *priRfVtxCL;
  std::vector<int>         *priRfNTrkDif;
  std::vector<float>       *bMass, *bVtxCL, *bPx, *bPy, *bPz;
  std::vector<double>      *bPxE, *bPyE, *bPzE;
  std::vector<float>       *bctau, *bctau2D, *bctauBS, *bctauMPV, *bctauMPV2D, *bctauMPVBS, *bctauRf, *bctauMPVRf;
  std::vector<float>       *bctauE, *bctau2DE, *bctauBSE, *bctauMPVE, *bctauMPV2DE, *bctauMPVBSE, *bctauRfE, *bctauMPVRfE;
  std::vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  std::vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;
  std::vector<float>       *bResMass, *bResVtxCL, *bResPx, *bResPy, *bResPz;
  std::vector<float>       *bResDecayVtxX, *bResDecayVtxY, *bResDecayVtxZ;
  std::vector<float>       *bResDecayVtxXE, *bResDecayVtxYE, *bResDecayVtxZE;
  std::vector<float>       *VMass, *VVtxCL, *VPx, *VPy, *VPz;
  std::vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  std::vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  std::vector<float>       *JMass, *JVtxCL, *JPx, *JPy, *JPz;
  std::vector<float>       *JDecayVtxX, *JDecayVtxY, *JDecayVtxZ;
  std::vector<float>       *JDecayVtxXE, *JDecayVtxYE, *JDecayVtxZE;
  std::vector<int>         *JmuOL;
  std::vector<float>       *mumPx, *mumPy, *mumPz, *mumD0, *mumD0E, *mumC2;
  std::vector<int>         *mumCat, *mumME1, *mumAngT, *mumNHits, *mumNPHits, *mumTrigL1Open1mu, *mumTrigL1Open2mu, *mumTrigL1Open2muTight, *mumTrig2mu0, *mumTrig2mu3, *mumTrig2mu0_quark;
  std::vector<int>         *mumTrig2mu0L2, *mumTrigmu0trk0, *mumTrigmu3trk0, *mumTrigmu0trkmu0, *mumTrigmu3trkmu0;
  std::vector<float>       *mupPx, *mupPy, *mupPz, *mupD0, *mupD0E, *mupC2;
  std::vector<int>         *mupCat, *mupME1, *mupAngT, *mupNHits, *mupNPHits, *mupTrigL1Open1mu, *mupTrigL1Open2mu, *mupTrigL1Open2muTight, *mupTrig2mu0, *mupTrig2mu3, *mupTrig2mu0_quark;
  std::vector<int>         *mupTrig2mu0L2, *mupTrigmu0trk0, *mupTrigmu3trk0, *mupTrigmu0trkmu0, *mupTrigmu3trkmu0;
  std::vector<float>       *VTrkpMass, *VTrkpPx, *VTrkpPy, *VTrkpPz, *VTrkpD0, *VTrkpD0E;
  std::vector<float>       *VTrkmMass, *VTrkmPx, *VTrkmPy, *VTrkmPz, *VTrkmD0, *VTrkmD0E;
  std::vector<float>       *bResTrkPx, *bResTrkPy, *bResTrkPz, *bResTrkD0, *bResTrkD0E;
  std::vector<int>         *bResTrkChg;
  int                 genKsPsi, genKsPsi2, genKstarpPsi, genLambdaPsi, feedup, feeddown;
  float               genBPt, genBeta, genBy, genJPsiPt, genJPsieta;
  float               genMupEta, genMupPt, genMupP, genMupPhi, genMumEta, genMumPt, genMumP, genMumPhi;
  int                 muAcc, muTrig, weight;
  int                 genHLT_2mu0, genHLT_2mu0_quark, genHLT_2mu3, genHLT_2mu0L2, genHLT_mu0trk0, genHLT_mu3trk0, genHLT_mu0trkmu0;
  int                 genHLT_mu3trkmu0, genHLT_mu0trkmu0OST, genHLT_mu0trkmu0OST_tight;
  int                 genHLT_mu3trkmu0OST, genHLT_L1muOpen, genHLT_L12muOpen, genHLT_L12muOpenTight;

  std::vector<float> *bMass2, *bMass3, *bMass4, *bMass5;
  std::vector<float> *bPx2, *bPx3, *bPx4, *bPx5;
  std::vector<float> *bPy2, *bPy3, *bPy4, *bPy5;
  std::vector<float> *bPz2, *bPz3, *bPz4, *bPz5;
  std::vector<float> *bDx2, *bDx3, *bDx4, *bDx5;
  std::vector<float> *bDy2, *bDy3, *bDy4, *bDy5;
  std::vector<float> *bDz2, *bDz3, *bDz4, *bDz5;


  float       truePriVtxX, truePriVtxY, truePriVtxZ;
  float       trueBPx, trueBPy, trueBPz;
  float       trueBDecayVtxX, trueBDecayVtxY, trueBDecayVtxZ;
  float       trueBResPx, trueBResPy, trueBResPz;
  float       trueBResDecayVtxX, trueBResDecayVtxY, trueBResDecayVtxZ;
  float       trueVPx, trueVPy, trueVPz;
  float       trueVDecayVtxX, trueVDecayVtxY, trueVDecayVtxZ;
  float       trueJPx, trueJPy, trueJPz;
  float       trueJDecayVtxX, trueJDecayVtxY, trueJDecayVtxZ;
  float       trueMumPx, trueMumPy, trueMumPz, trueMumD0;
  float       trueMupPx, trueMupPy, trueMupPz, trueMupD0;
  float       genJPt, genJP;
  float       trueVTrkpPx, trueVTrkpPy, trueVTrkpPz, trueVTrkpD0;
  float       trueVTrkmPx, trueVTrkmPy, trueVTrkmPz, trueVTrkmD0;
  float       trueBResTrkPx, trueBResTrkPy, trueBResTrkPz, trueBResTrkD0;
  int         trueBResTrkChg;
  std::vector<int> truthMatch, truthKs, truthPsi;
  std::vector<int> truthMatchPAT, truthKsPAT, truthPsiPAT;
  int         prompt;

};

#endif
