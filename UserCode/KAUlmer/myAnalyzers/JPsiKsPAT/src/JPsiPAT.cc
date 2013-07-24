// -*- C++ -*-
//
// Package:    JPsiPAT
// Class:      JPsiPAT
// 
/**\class JPsiPAT JPsiPAT.cc myAnalyzers/JPsiPAT/src/JPsiPAT.cc

 Description: <one line class summary>
Make rootTuple for b->s JPsi(mu+mu-) analyses

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer
//         Created:  Wed May  7 13:15:04 MDT 2008
// $Id: JPsiPAT.cc,v 1.3 2011/05/08 17:30:06 kaulmer Exp $
//
//


// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/interface/JPsiPAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/CLHEP/interface/Migration.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <utility>

//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsiPAT::JPsiPAT(const edm::ParameterSet& iConfig)
  :
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
  vtxSample( iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices")) ), 
  genParticles_( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  doMC_( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  onlyCount_( iConfig.getUntrackedParameter<bool>("onlyCount",false) ),

  tree_(0), treeMC_(0), l1_mu3(0), l1_2mu3(0), l1_muOpen(0), l1_mu0(0),
  hlt_mu3(0), hlt_mu5(0), hlt_mu7(0), hlt_mu9(0), hlt_2mu0(0), hlt_2mu3(0), hlt_2mu3JPsi(0), hltBJPsiMuMu(0), 
  hlt_L1muOpen(0), hlt_L12muOpen(0), nB(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priRfVtxX(0), priRfVtxY(0), priRfVtxZ(0), priRfVtxXE(0), priRfVtxYE(0), priRfVtxZE(0), priRfVtxCL(0), 
  priRfNTrkDif(0),
  bMass(0), bVtxCL(0), bPx(0), bPy(0), bPz(0), bPxE(0), bPyE(0), bPzE(0), 
  bctau(0), bctau2D(0), bctauBS(0), bctauMPV(0), bctauMPVBS(0), bctauRf(0), bctauMPVRf(0),
  bctauE(0), bctau2DE(0), bctauBSE(0), bctauMPVE(0), bctauMPVBSE(0), bctauRfE(0), bctauMPVRfE(0),
  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bResMass(0), bResVtxCL(0), bResPx(0), bResPy(0), bResPz(0),
  bResDecayVtxX(0), bResDecayVtxY(0), bResDecayVtxZ(0), bResDecayVtxXE(0), bResDecayVtxYE(0), bResDecayVtxZE(0),
  VMass(0), VVtxCL(0), VPx(0), VPy(0), VPz(0),
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0),
  VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  JMass(0), JVtxCL(0), JPx(0), JPy(0), JPz(0),
  JDecayVtxX(0), JDecayVtxY(0), JDecayVtxZ(0), JDecayVtxXE(0), JDecayVtxYE(0), JDecayVtxZE(0), JmuOL(0),
  mumPx(0), mumPy(0), mumPz(0), mumD0(0), mumD0E(0), mumC2(0), mumCat(0), mumME1(0), mumAngT(0), mumNHits(0), mumNPHits(0),
  mumTrigL1Open(0), mumTrig2mu0(0), mumTrig2mu3(0), mumTrigL1OpenD(0),
  mumTrig2mu0D(0), mumTrig2mu3D(0),
  mupPx(0), mupPy(0), mupPz(0), mupD0(0), mupD0E(0), mupC2(0), mupCat(0), mupME1(0), mupAngT(0), mupNHits(0), mupNPHits(0),
  mupTrigL1Open(0), mupTrig2mu0(0), mupTrig2mu3(0), mupTrigL1OpenD(0),
  mupTrig2mu0D(0), mupTrig2mu3D(0),
  VTrkpMass(0), VTrkpPx(0), VTrkpPy(0), VTrkpPz(0), 
  VTrkpD0(0), VTrkpD0E(0), 
  VTrkmMass(0), VTrkmPx(0), VTrkmPy(0), VTrkmPz(0), 
  VTrkmD0(0), VTrkmD0E(0), 
  bResTrkPx(0), bResTrkPy(0), bResTrkPz(0), 
  bResTrkD0(0), bResTrkD0E(0),bResTrkChg(0), 
  genKsPsi(0), genKstarpPsi(0), genLambdaPsi(0), feedup(0), feeddown(0),

  bMass2(0), bMass3(0), bMass4(0), bMass5(0), bPx2(0), bPx3(0), bPx4(0), bPx5(0),   bPy2(0), bPy3(0), bPy4(0), bPy5(0), 
  bPz2(0), bPz3(0), bPz4(0), bPz5(0), 
  bDx2(0), bDx3(0), bDx4(0), bDx5(0), bDy2(0), bDy3(0), bDy4(0), bDy5(0), bDz2(0), bDz3(0), bDz4(0), bDz5(0),

  truePriVtxX(0), truePriVtxY(0), truePriVtxZ(0), trueBPx(0), trueBPy(0), trueBPz(0), trueBDecayVtxX(0), trueBDecayVtxY(0), trueBDecayVtxZ(0),
  trueBResPx(0), trueBResPy(0), trueBResPz(0), trueBResDecayVtxX(0), trueBResDecayVtxY(0), trueBResDecayVtxZ(0),
  trueVPx(0), trueVPy(0), trueVPz(0), trueVDecayVtxX(0), trueVDecayVtxY(0), trueVDecayVtxZ(0),
  trueJPx(0), trueJPy(0), trueJPz(0), trueJDecayVtxX(0), trueJDecayVtxY(0), trueJDecayVtxZ(0),
  trueMumPx(0), trueMumPy(0), trueMumPz(0), trueMumD0(0), trueMupPx(0), trueMupPy(0), trueMupPz(0), trueMupD0(0),
  trueVTrkpPx(0), trueVTrkpPy(0), trueVTrkpPz(0), trueVTrkpD0(0),
  trueVTrkmPx(0), trueVTrkmPy(0), trueVTrkmPz(0), trueVTrkmD0(0),
  trueBResTrkPx(0), trueBResTrkPy(0), trueBResTrkPz(0), trueBResTrkD0(0), trueBResTrkChg(0),
  truthMatch(0), truthKs(0), truthPsi(0), truthMatchPAT(0), truthKsPAT(0), truthPsiPAT(0), prompt(0)

{
   //now do what ever initialization is needed
}


JPsiPAT::~JPsiPAT()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void


JPsiPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using std::vector;
   using namespace edm;
   using namespace reco;
   using namespace std;

   // get event content information

   ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

   //get genParticles  
   Handle<GenParticleCollection> genParticles;
   if (doMC_) {
     iEvent.getByLabel(genParticles_, genParticles);
   }

   // first get HLT results
   edm::Handle<edm::TriggerResults> hltresults;
   try {
     iEvent.getByLabel(hlTriggerResults_,hltresults);
   }
   catch ( ... ) {
     cout << "Couldn't get handle on HLT Trigger!" << endl;
   }
   if (!hltresults.isValid()) {
     cout << "No Trigger Results!" << endl;
   }
   else {
     int ntrigs=hltresults->size();
     if (ntrigs==0){
       cout << "No trigger name given in TriggerResults of the input " << endl;
     } 
     
     // get hold of trigger names - based on TriggerResults object!
     //edm::TriggerNames triggerNames_;
     //triggerNames_.init(*hltresults); 
     
     const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
     
     for (int itrig=0; itrig< ntrigs; itrig++) {
       TString trigName = triggerNames_.triggerName(itrig);
       int hltflag = (*hltresults)[itrig].accept();
       //int wasrun  = (*hltresults)[itrig].wasrun();
       //cout << "Trigger " <<  trigName << " was passed = " <<  hltflag << endl;
       if (trigName=="HLT_DoubleMu3_BJPsi") hltBJPsiMuMu = hltflag;
       if (trigName=="HLT_Mu3") hlt_mu3 = hltflag;
       if (trigName=="HLT_Mu5") hlt_mu5 = hltflag;
       if (trigName=="HLT_Mu7") hlt_mu7 = hltflag;
       if (trigName=="HLT_Mu9") hlt_mu9 = hltflag;  
       if (trigName=="HLT_DoubleMu0") hlt_2mu0 = hltflag;
       if (trigName=="HLT_DoubleMu3") hlt_2mu3 = hltflag;
       if (trigName=="HLT_DoubleMu3_JPsi") hlt_2mu3JPsi = hltflag;
       if (trigName=="HLT_L1MuOpen") hlt_L1muOpen = hltflag;
       if (trigName=="HLT_L1DoubleMuOpen") hlt_L12muOpen = hltflag;
     }
   }


   // get L1 trigger info
   
   edm::ESHandle<L1GtTriggerMenu> menuRcd;
   iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
   const L1GtTriggerMenu* menu = menuRcd.product();
   
   edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
   iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
   const DecisionWord dWord = gtRecord->decisionWord();  
   
   if ( menu->gtAlgorithmResult( "L1_SingleMu3", dWord) )  l1_mu3 = 1;
   if ( menu->gtAlgorithmResult( "L1_DoubleMu3", dWord) )  l1_2mu3 = 1;
   if ( menu->gtAlgorithmResult( "L1_SingleMuOpen", dWord) )  l1_muOpen = 1;
   if ( menu->gtAlgorithmResult( "L1_SingleMu0", dWord) )  l1_mu0 = 1;


   // get primary vertex
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(vtxSample, recVtxs);
   unsigned int nVtxTrks = 0;
   //reco::VertexCollection::const_iterator bestVtx = recVtxs->begin();
   reco::Vertex bestVtx = *(recVtxs->begin());

   //get primary with beamspot constraint
   Handle<reco::VertexCollection> recVtxsBS;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", recVtxsBS);
   
   nVtxTrks = 0;
   //reco::VertexCollection::const_iterator bestVtxBS = recVtxsBS->begin();
   reco::Vertex bestVtxBS = *(recVtxsBS->begin());
   
   //get beamspot
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
   else cout << "No beam spot available from EventSetup" << endl;
   
   //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // try reconstruction without fitting modules
   //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   Handle<vector<VertexCompositeCandidate> > theV0Handle;
   iEvent.getByLabel("generalV0Candidates", "Kshort", theV0Handle);

   Handle< vector<pat::GenericParticle> > thePATTrackHandle;
   iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);

   Handle< vector<pat::Muon> > thePATMuonHandle;
   iEvent.getByLabel("cleanPatMuons", thePATMuonHandle);
   

   cout << "event has " << theV0Handle->size() << "Ks and " << thePATMuonHandle->size() << "muons." << endl; 

   float pi = 3.14159265;

   if ( thePATMuonHandle->size()>=2 && !onlyCount_) {

       for ( std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin();
	     iMuonP != thePATMuonHandle->end(); ++iMuonP ) {
	 //check for mu+ first
	 if (iMuonP->charge() == 1) {
	   const pat::Muon *patMuonP = &(*iMuonP);
	   TrackRef glbTrackP = iMuonP->track();

	   if ( glbTrackP.isNull() ) {
	     cout << "continue due to no track ref" << endl;
	     continue;
	   }

	   //next check for mu-
	   for ( std::vector<pat::Muon>::const_iterator iMuonM = thePATMuonHandle->begin();
		 iMuonM != thePATMuonHandle->end(); ++iMuonM ) {
	     if (iMuonM->charge() == -1) {
	       const pat::Muon *patMuonM = &(*iMuonM);
	       //cout << "for muM: isGlobalMuon = " << iMuonM->isGlobalMuon() << "  isTrackerMuon = " << iMuonM->isTrackerMuon() << "  isStandAloneMuon = " << iMuonM->isStandAloneMuon() << "  isCaloMuon = " << iMuonM->isCaloMuon() << endl;
	       TrackRef glbTrackM = iMuonM->track();
	       if ( glbTrackM.isNull() ) {
		 cout << "continue from no track ref" << endl;
		 continue;
	       }

	       cout << "have 2 good oppositely charged muons. " << endl;
	       
	       TransientTrack muon1TT(glbTrackP, &(*bFieldHandle) );
               TransientTrack muon2TT(glbTrackM, &(*bFieldHandle) );





               //Creating a KinematicParticleFactory
	       KinematicParticleFactoryFromTransientTrack pFactory;
	       
	       //The mass of a muon and the insignificant mass sigma 
	       //to avoid singularities in the covariance matrix.
	       ParticleMass muon_mass = 0.10565837; //pdg mass
	       float muon_sigma = muon_mass*1.e-6;

	       
	       //initial chi2 and ndf before kinematic fits.
	       float chi = 0.;
	       float ndf = 0.;
	       vector<RefCountedKinematicParticle> muonParticles;
	       muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	       muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));



	       KinematicParticleVertexFitter fitter;   

	       RefCountedKinematicTree psiVertexFitTree;
	       psiVertexFitTree = fitter.fit(muonParticles); 
	       if (!psiVertexFitTree->isValid()) {
		 std::cout << "caught an exception in the psi vertex fit" << std::endl;
		 continue; 
	       }
	       psiVertexFitTree->movePointerToTheTop();
	       
	       RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	       RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

	       if ( psi_vFit_vertex_noMC->chiSquared() < 0 ) cout << "negative chisq from psi fit" << endl;

               RefCountedKinematicTree vertexFitTree = psiVertexFitTree;

	       vertexFitTree->movePointerToTheFirstChild();
	       RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
	       cout << "mass of mu1 = " << mu1CandMC->currentState().mass() << " and charge = " << mu1CandMC->currentState().particleCharge() << endl;
	       vertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
	       cout << "mass of mu2 = " << mu2CandMC->currentState().mass() << " and charge = " << mu2CandMC->currentState().particleCharge() << endl;

	       // get mu+ and mu- momenta from final B fit
	       KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMupKP;
	       KinematicParameters psiMumKP;
	       
	       if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	       if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	       if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
	       if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 
	       
	       //get cateogories for muon candidates
	       if (iMuonP->isGlobalMuon()) {
	         if (iMuonP->isTrackerMuon()) mupCategory = 1;
                 else mupCategory = 8;
	       }
	       else if (iMuonP->isTrackerMuon()) {
	         if (muon::isGoodMuon(*iMuonP, muon::TrackerMuonArbitrated)) {
	           if (muon::isGoodMuon(*iMuonP, muon::TMOneStationLoose)) mupCategory = 4;
		   else mupCategory = 3;
		 } else mupCategory = 2;
	       }
	       else if (iMuonP->isStandAloneMuon()) mupCategory = 5;
	       else if (iMuonP->isCaloMuon()) mupCategory = 6;
	       else mupCategory = 7;
	       
	       if (iMuonM->isGlobalMuon()) {
	         if (iMuonM->isTrackerMuon()) mumCategory = 1;
	         else mumCategory = 8;        
	       } 
	       else if (iMuonM->isTrackerMuon()) {
	         if (muon::isGoodMuon(*iMuonM, muon::TrackerMuonArbitrated)) {
	           if (muon::isGoodMuon(*iMuonM, muon::TMOneStationLoose)) mumCategory = 4;
		   else mumCategory = 3;
		 } else mumCategory = 2;
	       }
	       else if (iMuonM->isStandAloneMuon()) mumCategory = 5;
	       else if (iMuonM->isCaloMuon()) mumCategory = 6;
	       else mumCategory = 7;

               const reco::Muon *recoMuonM = patMuonM;
	       const reco::Muon *recoMuonP = patMuonP;

               mumME1Clean = HasGoodME11(*recoMuonM,2.);
               mupME1Clean = HasGoodME11(*recoMuonP,2.);

	       
	       cout << "filling new candidate" << endl;
	       
	       // fill candidate variables now
	       
	       
	       JMass->push_back( psi_vFit_noMC->currentState().mass() ); 
	       JDecayVtxX->push_back( psi_vFit_vertex_noMC->position().x() );
	       JDecayVtxY->push_back( psi_vFit_vertex_noMC->position().y() );
	       JDecayVtxZ->push_back( psi_vFit_vertex_noMC->position().z() );
	       JDecayVtxXE->push_back( sqrt(psi_vFit_vertex_noMC->error().cxx()) );
	       JDecayVtxYE->push_back( sqrt(psi_vFit_vertex_noMC->error().cyy()) );
	       JDecayVtxZE->push_back( sqrt(psi_vFit_vertex_noMC->error().czz()) );
	       JVtxCL->push_back( ChiSquaredProbability((double)(psi_vFit_vertex_noMC->chiSquared()),(double)(psi_vFit_vertex_noMC->degreesOfFreedom())) );
	       JPx->push_back( psiMumKP.momentum().x() + psiMupKP.momentum().x() );
	       JPy->push_back( psiMumKP.momentum().y() + psiMupKP.momentum().y() );
	       JPz->push_back( psiMumKP.momentum().z() + psiMupKP.momentum().z() );
	       JmuOL->push_back( muon::overlap( *recoMuonM, *recoMuonP, 1.0, 1.0 ) );

	       mumPx->push_back(psiMumKP.momentum().x());
	       mumPy->push_back(psiMumKP.momentum().y());
	       mumPz->push_back(psiMumKP.momentum().z());
	       mumD0->push_back( glbTrackM->d0() );
	       mumD0E->push_back(glbTrackM->d0Error() );
	       mumC2->push_back( glbTrackM->normalizedChi2() );
	       mumCat->push_back( mumCategory );
               mumME1->push_back(mumME1Clean);
               mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMLastStationAngTight) );
	       mumNHits->push_back( glbTrackM->numberOfValidHits() );
               mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
	       mupPx->push_back(psiMupKP.momentum().x());
	       mupPy->push_back(psiMupKP.momentum().y());
	       mupPz->push_back(psiMupKP.momentum().z());
	       mupD0->push_back( glbTrackP->d0() );
	       mupD0E->push_back(glbTrackP->d0Error() );
	       mupC2->push_back( glbTrackP->normalizedChi2() );
	       mupCat->push_back( mupCategory );
               mupME1->push_back(mupME1Clean);
               mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMLastStationAngTight) );
	       mupNHits->push_back( glbTrackP->numberOfValidHits() );
               mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );	 

               mupTrigL1Open->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty());
	       mupTrig2mu0->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty());
	       mupTrig2mu3->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty());
               mumTrigL1Open->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty());
	       mumTrig2mu0->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty());
	       mumTrig2mu3->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty());



               /////////////////////////////////////////////////////////////////////////////////////////
               // check trigger objects
               // compare to trigger match from David
               //////////////////////////////////////////////////
   
               edm::InputTag tag2("hltTriggerSummaryAOD::HLT");
               edm::Handle<trigger::TriggerEvent> trgEvent;
               iEvent.getByLabel(tag2,trgEvent);

               const trigger::size_type numFilterObjects(trgEvent->sizeFilters());
               //cout << "Used Processname: " << trgEvent->usedProcessName() << endl;
               //cout << "Number of TriggerFilterObjects: " << numFilterObjects << endl;
               //cout << "The TriggerFilterObjects: #, tag" << endl;
               //for ( trigger::size_type i = 0; i != numFilterObjects; ++i ) 
               //  cout << i << " " << trgEvent->filterTag(i).encode() << endl;

               bool passed2mu0 = false;
               bool passed2mu3 = false;
               for ( trigger::size_type i = 0; i != numFilterObjects; ++i  ) {
                 string module = trgEvent->filterTag(i).encode();
                 if (TString(module) == "hltDiMuonL3PreFiltered0::HLT" )
                   passed2mu0 = true;
                 if (TString(module) == "hltDiMuonL3PreFiltered::HLT" )
                   passed2mu3 = true;   
               }

               // Get the L1 candidates
               vector<trigger::TriggerObject> L1cands;
               L1cands.clear();
               const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());
               for(uint j=0; j < TOC.size(); j++) {
                 if ( abs(TOC[j].id()) == 13 ) {
                   L1cands.push_back(TOC[j]);
                   //cout << "L1 id: " << TOC[j].id() << " and eta: " << TOC[j].eta() << " and phi = " << TOC[j].phi() << endl;
                 }
               }
    
               // Get the HLT candidates
               vector<trigger::TriggerObject> hlt2mu0Cands;
               hlt2mu0Cands.clear();
               vector<trigger::TriggerObject> hlt2mu3Cands;
               hlt2mu3Cands.clear();
               InputTag tag2mu0 = InputTag("hltDiMuonL3PreFiltered0::HLT");
               InputTag tag2mu3 = InputTag("hltDiMuonL3PreFiltered::HLT");
               if ( passed2mu0 ) {
                 const trigger::Keys& KEYS2mu0(trgEvent->filterKeys(trgEvent->filterIndex(tag2mu0)));
                 for(uint j=0; j < KEYS2mu0.size(); j++) {
                   trigger::size_type hltf=KEYS2mu0[j];
                   if ( abs(TOC[hltf].id()) == 13 ) {
                     hlt2mu0Cands.push_back(TOC[hltf]);
                     //cout << "HLT n: " << j << " id: " << TOC[hltf].id() << " and eta: " << TOC[hltf].eta() << " and phi = " << TOC[hltf].phi() << 
                     //  " from tag = " << "hltDiMuonL3PreFiltered0::HLT" << endl;
                   }
                 }     
               }
               if ( passed2mu3 ) {
                 const trigger::Keys& KEYS2mu3(trgEvent->filterKeys(trgEvent->filterIndex(tag2mu3)));
                 for(uint j=0; j < KEYS2mu3.size(); j++) {
                   trigger::size_type hltf=KEYS2mu3[j];
                   if ( abs(TOC[hltf].id()) == 13 ) {
                     hlt2mu3Cands.push_back(TOC[hltf]);
                     //cout << "HLT n: " << j << " id: " << TOC[hltf].id() << " and eta: " << TOC[hltf].eta() << " and phi = " << TOC[hltf].phi() << 
                     //  " from tag = " << "hltDiMuonL3PreFiltered::HLT" << endl;
                   }
                 }
               }   
  
  
               // Then the objects are in L1cands, hlt2mu0Cands and hlt2mu3Cands

               bool mupTrigMatchL1open = false;
               bool mupTrigMatch2mu0 = false;
               bool mupTrigMatch2mu3 = false;
               bool mumTrigMatchL1open = false;
               bool mumTrigMatch2mu0 = false;
               bool mumTrigMatch2mu3 = false;

               //cout << "size of L1 cands = " << L1cands.size() << endl;
               for ( size_t candL = 0; candL < L1cands.size(); candL++ ) {
                 double eta= L1cands[candL].eta();
                 double phi= L1cands[candL].phi();
                 //cout << "L1 eta,phi = " << eta << "," << phi << endl;
                 if ( deltaR(eta, phi, psiMumKP.momentum().eta() ,psiMumKP.momentum().phi() ) < 0.3 ) mumTrigMatchL1open = true;
                 if ( deltaR(eta, phi, psiMupKP.momentum().eta() ,psiMupKP.momentum().phi() ) < 0.3 ) mupTrigMatchL1open = true; 
               }
   
               //cout << "mum eta,phi = " << psiMumKP.momentum().eta() << "," << psiMumKP.momentum().phi() << endl;
               //cout << "mup eta,phi = " << psiMupKP.momentum().eta() << "," << psiMupKP.momentum().phi() << endl;
               //cout << "mum L1open trig match David = " << mumTrigMatchL1open << endl;
               //cout << "mup L1open trig match David = " << mupTrigMatchL1open << endl;

               mumTrigL1OpenD->push_back( mumTrigMatchL1open );
               mupTrigL1OpenD->push_back( mupTrigMatchL1open );


               //cout << "size of hlt 2mu0 cands = " << hlt2mu0Cands.size() << endl;
               for ( size_t candNum = 0; candNum < hlt2mu0Cands.size(); candNum++ ) {
                 double eta= hlt2mu0Cands[candNum].eta();
                 double phi= hlt2mu0Cands[candNum].phi();   
                 //cout << "HLT 2mu0 eta,phi = " << eta << "," << phi << endl;
                 if ( deltaR(eta, phi, psiMumKP.momentum().eta() ,psiMumKP.momentum().phi() ) < 0.05 ) mumTrigMatch2mu0 = true;
                 if ( deltaR(eta, phi, psiMupKP.momentum().eta() ,psiMupKP.momentum().phi() ) < 0.05 ) mupTrigMatch2mu0 = true; 
               }

               //cout << "mum 2mu0 trig match David = " << mumTrigMatch2mu0 << endl;
               //cout << "mup 2mu0 trig match David = " << mupTrigMatch2mu0 << endl;
               mumTrig2mu0D->push_back( mumTrigMatch2mu0 );
               mupTrig2mu0D->push_back( mupTrigMatch2mu0 );

               //cout << "size of hlt 2mu3 cands = " << hlt2mu3Cands.size() << endl;
               for ( size_t candNum = 0; candNum < hlt2mu3Cands.size(); candNum++ ) {
                 double eta= hlt2mu3Cands[candNum].eta();
                 double phi= hlt2mu3Cands[candNum].phi();   
                 //cout << "HLT 2mu3 eta,phi = " << eta << "," << phi << endl;
                 if ( deltaR(eta, phi, psiMumKP.momentum().eta() ,psiMumKP.momentum().phi() ) < 0.05 ) mumTrigMatch2mu3 = true;
                 if ( deltaR(eta, phi, psiMupKP.momentum().eta() ,psiMupKP.momentum().phi() ) < 0.05 ) mupTrigMatch2mu3 = true; 
               }
 
               //cout << "mum 2mu3 trig match David = " << mumTrigMatch2mu3 << endl;
               //cout << "mup 2mu3 trig match David = " << mupTrigMatch2mu3 << endl;
               mumTrig2mu3D->push_back( mumTrigMatch2mu3 );
               mupTrig2mu3D->push_back( mupTrigMatch2mu3 );

	       nB++;
	       
	       
	       
	
	       /////////////////////////////////////////////////
	

	       muonParticles.clear();


	     
	   }
	 }
       } 
     } 
   } // if V0Handle > 0 and muHandle > 1

   //////////////////////////////////////////////////////
   //////// get truth information from genParticles only for events with a B candidate
   
   priVtxX = bestVtx.x();
   priVtxY = bestVtx.y();
   priVtxZ = bestVtx.z();
   priVtxXE = bestVtx.xError();
   priVtxYE = bestVtx.yError();
   priVtxZE = bestVtx.zError();
   priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
   if (nB > 0 && doMC_) {

     genKsPsi = -1; genKstarpPsi = -1; genLambdaPsi = -1; prompt = 1; feedup = -1; feeddown = -1;

     cout << "Size of genParticle collection is " << genParticles->size() << endl;
     
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       // check if any of our signals were generated
       
       const Candidate & BCand = (*genParticles)[ k ];

       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
	 // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found B0";
	 int ipsi(-1), iks(-1);
	 bool wrong = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   //cout << "B0 daughter " << i << " has id = " << genDau->pdgId() << endl;
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1),  ipi1(-1), ipi2(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     if ( genGDau->pdgId()==13 && genDau->pdgId()==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && genDau->pdgId()==443 ) imu2 = j;
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
	     for ( uint m = 0; m < genGDau->numberOfDaughters(); m++){
	       const Candidate * genGGDau = genGDau->daughter(m);
	       cout << " ===" << genGGDau->pdgId();
	       if ( genGGDau->pdgId()==211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi1 = m;
	       if ( genGGDau->pdgId()==-211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi2 = m;
	       if ( genGDau->pdgId()==310 && abs(genGGDau->pdgId())!=211 && genGGDau->pdgId()!=22 ) wrong = true; 
	     }
           }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=311 && genDau->pdgId()!=22 ) wrong = true;
           if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi1!=-1 && ipi2!=-1&&!wrong) iks = i;
	 }
	 
	 if (ipsi!=-1&&iks!=-1&&!wrong)
	   cout << " found genKsPsi";
	 
	 cout << endl;
	 
	 if (ipsi!=-1&&iks!=-1) {
	   
	   genKsPsi = 1;
	   //write out info from daughters
	   const Candidate * genpsi =  BCand.daughter(ipsi);
	   const Candidate * genks =  BCand.daughter(iks)->daughter(0);
	   
	   const Candidate *BTemp = &BCand;
	   while ( (BTemp->vx() == genpsi->vx()) && (BTemp->vy() == genpsi->vy()) ) {
	     if (BTemp->numberOfMothers() > 0) {
	       BTemp = BTemp->mother(0);
	     } else {
	       cout << "Can't find mother of oscillating B!" << endl;
	       BTemp = 0;
	     }
	   }
	   
	   truePriVtxX = BTemp->vx();
	   truePriVtxY = BTemp->vy();
	   truePriVtxZ = BTemp->vz();
	   
	   trueBPx = BCand.px();
	   trueBPy = BCand.py();
	   trueBPz = BCand.pz();
	   
	   trueJPx = genpsi->px();
	   trueJPy = genpsi->py();
	   trueJPz = genpsi->pz();
	   trueBDecayVtxX = genpsi->vx();
	   trueBDecayVtxY = genpsi->vy();
	   trueBDecayVtxZ = genpsi->vz();

	   cout << "true B decay vertex = (" << genpsi->vx() << "," << genpsi->vy() << "," << genpsi->vz() << ")" << endl;

	   //cout << "getting info for mu+ and mu- with imup = " << imu1 << " and imum = " << imu2 << endl;

	   for (uint j=0; j<genpsi->numberOfDaughters(); j++) {
	     cout << "daughter " << j << " from genpsi has id = " << genpsi->daughter(j)->pdgId() << endl;
	     if (genpsi->daughter(j)->pdgId()==13) { //13 is a mu-
	       trueMumPx = genpsi->daughter(j)->px();
	       trueMumPy = genpsi->daughter(j)->py();
	       trueMumPz = genpsi->daughter(j)->pz();
	       trueJDecayVtxX = genpsi->daughter(j)->vx();
	       trueJDecayVtxY = genpsi->daughter(j)->vy();
	       trueJDecayVtxZ = genpsi->daughter(j)->vz();
	     }
	     if (genpsi->daughter(j)->pdgId()==-13) { //-13 is a mu+
	       trueMupPx = genpsi->daughter(j)->px();
	       trueMupPy = genpsi->daughter(j)->py();
	       trueMupPz = genpsi->daughter(j)->pz();
	     }
	   }

	   trueVPx = genks->px();
	   trueVPy = genks->py();
	   trueVPz = genks->pz();
	   
	   for (uint j=0; j<genks->numberOfDaughters(); j++) {
	     if ( genks->daughter(j)->charge()>0 ) {
	       trueVTrkpPx = genks->daughter(j)->px();
	       trueVTrkpPy = genks->daughter(j)->py();
	       trueVTrkpPz = genks->daughter(j)->pz();
	       trueVDecayVtxX = genks->daughter(j)->vx();
	       trueVDecayVtxY = genks->daughter(j)->vy();
	       trueVDecayVtxZ = genks->daughter(j)->vz();
	     }
	     if ( genks->daughter(j)->charge()<0 ) {
	       trueVTrkmPx = genks->daughter(j)->px();
	       trueVTrkmPy = genks->daughter(j)->py();
	       trueVTrkmPz = genks->daughter(j)->pz();
	     }
	   }
	
	   cout << "moving to truth match check from reco objects by hand" << endl;
	   
	   /////////////////////////////////////////////////////////////////////////
	   // determine MC truth

	   // calculate true eta and phi for all tracks
	   float trueMupPhi = atan(trueMupPy/trueMupPx);
	   if ( trueMupPx < 0 && trueMupPy < 0 ) trueMupPhi -= pi;
	   if ( trueMupPx < 0 && trueMupPy > 0 ) trueMupPhi += pi;
	   float trueMupP = sqrt( trueMupPx*trueMupPx +  trueMupPy*trueMupPy +  trueMupPz*trueMupPz );
	   float trueMupEta = 0.5*log( (trueMupP + trueMupPz)/(trueMupP - trueMupPz) );
	   
	   float trueMumPhi = atan(trueMumPy/trueMumPx);
	   if ( trueMumPx < 0 && trueMumPy < 0 ) trueMumPhi -= pi;
	   if ( trueMumPx < 0 && trueMumPy > 0 ) trueMumPhi += pi;
	   float trueMumP = sqrt( trueMumPx*trueMumPx +  trueMumPy*trueMumPy +  trueMumPz*trueMumPz );
	   float trueMumEta = 0.5*(log( (trueMumP + trueMumPz)/(trueMumP - trueMumPz) ) );
	   
	   float truePipPhi = atan(trueVTrkpPy/trueVTrkpPx);
	   if ( trueVTrkpPx < 0 && trueVTrkpPy < 0 ) truePipPhi -= pi;
	   if ( trueVTrkpPx < 0 && trueVTrkpPy > 0 ) truePipPhi += pi;
	   float truePipP = sqrt( trueVTrkpPx*trueVTrkpPx +  trueVTrkpPy*trueVTrkpPy +  trueVTrkpPz*trueVTrkpPz );
	   float truePipEta = 0.5*log( (truePipP + trueVTrkpPz)/(truePipP - trueVTrkpPz) );
	   
	   float truePimPhi = atan(trueVTrkmPy/trueVTrkmPx);
	   if ( trueVTrkmPx < 0 && trueVTrkmPy < 0 ) truePimPhi -= pi;
	   if ( trueVTrkmPx < 0 && trueVTrkmPy > 0 ) truePimPhi += pi;
	   float truePimP = sqrt( trueVTrkmPx*trueVTrkmPx +  trueVTrkmPy*trueVTrkmPy +  trueVTrkmPz*trueVTrkmPz );
	   float truePimEta = 0.5*log( (truePimP + trueVTrkmPz)/(truePimP - trueVTrkmPz) );
	   
	   //cout << "=======================" << endl;
	   //cout << "For true B muP eta, phi = " << trueMupEta << "," << trueMupPhi << endl;
	   //cout << "For true B muM eta, phi = " << trueMumEta << "," << trueMumPhi << endl;
	   //cout << "For true B piP eta, phi = " << truePipEta << "," << truePipPhi << endl;
	   //cout << "For true B piM eta, phi = " << truePimEta << "," << truePimPhi << endl;	 
	   
	   float RcutMu = 0.02;
	   float RcutPi = 0.10;
	   float RcutVtx = 10.;
	   
	   truthMatch.clear(); truthKs.clear(); truthPsi.clear();
	   
	   for (uint i = 0; i<mupPx->size(); i++) {
	     //loop to check all B candidates found

	     bool istrueMup = false;
	     bool istrueMum = false;
	     bool istruePip = false;
	     bool istruePim = false;
	     bool istrueKs = false;
	     bool istruePsi = false;
	     bool istrueB = false;
	     
	     // calculate eta and phi for all tracks in B candidate
	     float mupPhi = atan(mupPy->at(i)/mupPx->at(i));
	     if ( mupPx->at(i) < 0 && mupPy->at(i) < 0 ) mupPhi -= pi;
	     if ( mupPx->at(i) < 0 && mupPy->at(i) > 0 ) mupPhi += pi;
	     float mupP = sqrt( mupPx->at(i)*mupPx->at(i) +  mupPy->at(i)*mupPy->at(i) +  mupPz->at(i)*mupPz->at(i) );
	     float mupEta = 0.5*log( (mupP + mupPz->at(i))/(mupP - mupPz->at(i)) );
	     
	     float mumPhi = atan(mumPy->at(i)/mumPx->at(i));
	     if ( mumPx->at(i) < 0 && mumPy->at(i) < 0 ) mumPhi -= pi;
	     if ( mumPx->at(i) < 0 && mumPy->at(i) > 0 ) mumPhi += pi;
	     float mumP = sqrt( mumPx->at(i)*mumPx->at(i) +  mumPy->at(i)*mumPy->at(i) +  mumPz->at(i)*mumPz->at(i) );
	     float mumEta = 0.5*log( (mumP + mumPz->at(i))/(mumP - mumPz->at(i)) );	 
	     
	     float pipPhi = atan(VTrkpPy->at(i)/VTrkpPx->at(i));
	     if ( VTrkpPx->at(i) < 0 && VTrkpPy->at(i) < 0 ) pipPhi -= pi;
	     if ( VTrkpPx->at(i) < 0 && VTrkpPy->at(i) > 0 ) pipPhi += pi;
	     float pipP = sqrt( VTrkpPx->at(i)*VTrkpPx->at(i) +  VTrkpPy->at(i)*VTrkpPy->at(i) +  VTrkpPz->at(i)*VTrkpPz->at(i) );
	     float pipEta = 0.5*log( (pipP + VTrkpPz->at(i))/(pipP - VTrkpPz->at(i)) );
	     
	     float pimPhi = atan(VTrkmPy->at(i)/VTrkmPx->at(i));
	     if ( VTrkmPx->at(i) < 0 && VTrkmPy->at(i) < 0 ) pimPhi -= pi;
	     if ( VTrkmPx->at(i) < 0 && VTrkmPy->at(i) > 0 ) pimPhi += pi;
	     float pimP = sqrt( VTrkmPx->at(i)*VTrkmPx->at(i) +  VTrkmPy->at(i)*VTrkmPy->at(i) +  VTrkmPz->at(i)*VTrkmPz->at(i) );
	     float pimEta = 0.5*log( (pimP + VTrkmPz->at(i))/(pimP - VTrkmPz->at(i)) );
	     
	     //cout << "For reco B muP eta, phi = " << mupEta << "," << mupPhi << endl;
	     //cout << "For reco B muM eta, phi = " << mumEta << "," << mumPhi << endl;
	     //cout << "For reco B piP eta, phi = " << pipEta << "," << pipPhi << endl;
	     //cout << "For reco B piM eta, phi = " << pimEta << "," << pimPhi << endl;	 
	     
	     float deltaRmup = sqrt( (mupEta-trueMupEta)*(mupEta-trueMupEta) +  (mupPhi-trueMupPhi)*(mupPhi-trueMupPhi) );
	     if ( deltaRmup < RcutMu ) istrueMup = true;
	     
	     float deltaRmum = sqrt( (mumEta-trueMumEta)*(mumEta-trueMumEta) +  (mumPhi-trueMumPhi)*(mumPhi-trueMumPhi) ) ;
	     if ( deltaRmum < RcutMu ) istrueMum = true;
	     
	     float deltaRpip = sqrt( (pipEta-truePipEta)*(pipEta-truePipEta) +  (pipPhi-truePipPhi)*(pipPhi-truePipPhi) );
	     if ( deltaRpip < RcutPi ) istruePip = true;
	     
	     float deltaRpim = sqrt( (pimEta-truePimEta)*(pimEta-truePimEta) +  (pimPhi-truePimPhi)*(pimPhi-truePimPhi) );
	     if ( deltaRpim < RcutPi ) istruePim = true;
	     
	     //cout << "deltaR for mup = " << deltaRmup << ", mum = " << deltaRmum << ", deltaRpip = " << deltaRpip << ", deltaRpim = " << deltaRpim << endl;

	     //check Ks vertex position truth match
	     float deltaRksvtx = sqrt( (trueVDecayVtxX - VDecayVtxX->at(i))*
				       (trueVDecayVtxX - VDecayVtxX->at(i)) +
				       (trueVDecayVtxY - VDecayVtxY->at(i))*
				       (trueVDecayVtxY - VDecayVtxY->at(i)) +
				       (trueVDecayVtxZ - VDecayVtxZ->at(i))*
				       (trueVDecayVtxZ - VDecayVtxZ->at(i)) );	     

	     if ( istrueMup & istrueMum ) istruePsi = true;
	     if ( istruePip & istruePim && (deltaRksvtx<RcutVtx) ) istrueKs = true;
	     if ( istruePsi & istrueKs ) istrueB = true;

	     if (istruePsi) {
	       //cout << "true Psi from reco from cand " << i << endl;
	       truthPsi.push_back(1);
	     } else truthPsi.push_back(-1);
	     if (istrueKs) {
	       //cout << "true Ks from reco from cand " << i << endl;
	       truthKs.push_back(1);
	     } else truthKs.push_back(-1);
	     if (istrueB) {
	       //cout << "true B from reco from cand " << i << endl;
	       truthMatch.push_back(1);
	     } else truthMatch.push_back(-1);
	     
	   }
	   
	 }
	 
	 ///////////////////////////////////////////////////////////////////////
	   
       } // closes if (BCand == B0ID)
       
       //check for B->JPsiK*+(kspi) decay   
       if ( abs(BCand.pdgId())==521 ) {
	 cout << "found B+";
	 int ipsi(-1), ikstp(-1);
	 bool wrong = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1),  ik0(-1), ipi(-1), iks(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     if ( genGDau->pdgId()==13 && abs(genDau->pdgId())==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && abs(genDau->pdgId())==443 ) imu2 = j;
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 )
	       wrong = true;
	     if (abs(genGDau->pdgId())==311 && abs(genDau->pdgId())==323 ) {
	       //K*+ decays to K0 pi and K0->Ks, so check for that, genGDau is the K0
	       for ( uint m = 0; m < genGDau->numberOfDaughters(); m++){
		 cout << " ===" << genGDau->daughter(m)->pdgId();
		 if (genGDau->daughter(m)->pdgId()==310) iks = m;
	       }
	       if (iks!=-1) {
		 const Candidate * ks = genGDau->daughter(iks);
		 int ipi1(-1), ipi2(-1);
		 for ( uint k = 0; k < ks->numberOfDaughters(); k++){
		   cout << " ====" << ks->daughter(k)->pdgId();
		   if (ks->daughter(k)->pdgId()==211) ipi1 = k;
		   if (ks->daughter(k)->pdgId()==-211) ipi2 = k; 
		   if ( abs(ks->daughter(k)->pdgId())!=211 && ks->daughter(k)->pdgId()!=22 )
		     wrong = true;
		 }
		 if (ipi1!=-1&&ipi2!=-1) {
		   ik0 = i;
		 }
	       }
	     }
	     if ( abs(genGDau->pdgId())==211 && abs(genDau->pdgId())==323 ) ipi = j;
	     if ( abs(genDau->pdgId())==311 && abs(genGDau->pdgId())!=211 && genGDau->pdgId()!=22 )
	       wrong = true;
	   }
	   if ( genDau->pdgId()!=443 && genDau->pdgId()!=323 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ik0!=-1 && ipi!=-1&&!wrong) ikstp = i;
	 }
	 if ( ipsi!=-1 && ikstp!=-1 &&!wrong) {
	   cout << " found genKstarpPsi";
	   genKstarpPsi =1;
	 }
	 cout << endl;
       } //check for genparticle with id = 521 for a B+
       
       // check for Lambda_b
       if (abs(BCand.pdgId())==5122) {
	 cout << "found Lambda_B";
	 int ipsi(-1), ilam(-1);
	 bool wrong = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check Lambda_b for psi and Lambda daughters
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1), ipi(-1), ip(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     cout << " ==" << genDau->daughter(j)->pdgId();
	     if (genDau->daughter(j)->pdgId()==13 && genDau->pdgId()==443) imu1 = j;
	     if (genDau->daughter(j)->pdgId()==-13 && genDau->pdgId()==443) imu2 = j; 
	     if ( abs(genDau->daughter(j)->pdgId())==211 && abs(genDau->pdgId())==3122) ipi = j;
	     if ( abs(genDau->daughter(j)->pdgId())==2212 && abs(genDau->pdgId())==3122) ip = j;
	     if ( genDau->pdgId()==443 && abs(genDau->daughter(j)->pdgId())!=13 && genDau->daughter(j)->pdgId()!=22 )
	       wrong = true;
	     if ( abs(genDau->pdgId())==3122 && abs(genDau->daughter(j)->pdgId())!=211 && abs(genDau->daughter(j)->pdgId())!=2212 && genDau->daughter(j)->pdgId()!=22 )
	       wrong = true;
	   }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=3122 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi!=-1&&ip!=-1&&!wrong) ilam = i;
	 }
	 if (ipsi!=-1&&ilam!=-1 &&!wrong) {
	   cout << " found genLambdaPsi";
	   genLambdaPsi = 1;
	 }
	 cout << endl;
       } // if (id==LambdaBID) 
       
       if ( abs(BCand.pdgId())==531 && abs(BCand.daughter(0)->pdgId())!=531 ) {
	 //only check after B_s oscilation
	 cout << "found B_s";
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     for ( uint k = 0; k < genGDau->numberOfDaughters(); k++){
	       cout << " ===" << genGDau->daughter(k)->pdgId();
	     }
	   }
	 }
	 cout << endl;
       }
       

       // check to see if JPsi is prompt
       bool isPrompt = true;
       const Candidate & PsiCand = (*genParticles)[ k ];
       if (abs(PsiCand.pdgId())==443) {
	 for ( uint i = 0; i < PsiCand.numberOfMothers(); i++){
	   const Candidate * psiMom = PsiCand.mother(i);
	   cout << "psi mother has ID = " << psiMom->pdgId() << endl;
	   if ( (abs(psiMom->pdgId())<600 && abs(psiMom->pdgId())>500) || (abs(psiMom->pdgId())<6000 && abs(psiMom->pdgId())>5000) ) {
	     isPrompt = false;
	     continue;
	   } else {
	     for ( uint i = 0; i < psiMom->numberOfMothers(); i++){
	       const Candidate * psiGMom = psiMom->mother(i);
	       cout << "psi grandmother has ID = " << psiGMom->pdgId() << endl;
	       if ( (abs(psiGMom->pdgId())<600 && abs(psiGMom->pdgId())>500) ||  (abs(psiGMom->pdgId())<6000 && abs(psiGMom->pdgId())>5000) ) {
		 isPrompt = false;
		 continue;
	       } else {
		 for ( uint i = 0; i < psiGMom->numberOfMothers(); i++){
		   const Candidate * psiGGMom = psiGMom->mother(i);
		   cout << "psi greatgrandmother has ID = " << psiGGMom->pdgId() << endl;
		   if ( (abs(psiGGMom->pdgId())<600 && abs(psiGGMom->pdgId())>500) ||  (abs(psiGGMom->pdgId())<6000 && abs(psiGGMom->pdgId())>5000) ) {
		     isPrompt = false;
		     continue;
		   }
		 }
	       }
	     }
	   }
	 }
	 if (!isPrompt) prompt = -1;
       }

       // check for JPsiKs feed up // currently from B+->JPsi K+
       if (abs(BCand.pdgId())==521) {
	 bool psidau = false;
	 bool kdau = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   if ( genDau->pdgId()==443 ) psidau = true;
	   if ( abs(genDau->pdgId())==321 ) kdau = true;
	 }
	 if ( psidau && kdau ) feedup = 1;
       }

       if (abs(BCand.pdgId())==5122) {  //Lambda_B
	 bool psidau = false;
	 bool lamdau = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   if ( genDau->pdgId()==443 ) psidau = true;
	   if ( abs(genDau->pdgId())==3122 ) lamdau = true;
	 }
	 if ( psidau && lamdau ) feedup = 2;
       }
       
       // check for JPsiKs feed down // currently from 
       if (abs(BCand.pdgId())==511) {
	 bool psidau = false;
	 bool ksdau = false;
	 bool psi2Sdau = false;
	 bool chic1dau = false;
	 bool kst0dau = false;
	 bool kst20dau = false;
	 bool pipGDau = false;
	 bool pizGDau = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   if ( genDau->pdgId()==443 ) psidau = true;
	   if ( genDau->pdgId()==310 || abs(genDau->pdgId())==311 ) ksdau = true;
	   if ( abs(genDau->pdgId())==100443 ) psi2Sdau = true;
	   if ( abs(genDau->pdgId())==20443 ) chic1dau = true;
	   if ( abs(genDau->pdgId())==313 ) {
	     kst0dau = true;
	     for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	       const Candidate * genGDau = genDau->daughter(j);
	       if ( abs(genGDau->pdgId())==211 ) pipGDau = true;
	       if ( abs(genGDau->pdgId())==111 ) pizGDau = true;
	     }
	   }
	   if ( abs(genDau->pdgId())==315 ) kst20dau = true;
	 }
	 if ( psi2Sdau && ksdau ) feeddown = 1;
	 if ( chic1dau && ksdau ) feeddown = 2;
	 if ( psidau && kst0dau && pizGDau ) feeddown = 3;
	 if ( psidau && kst0dau && pipGDau ) feeddown = 4;
	 if ( psidau && kst20dau ) feeddown = 5;
       }
       
       if (abs(BCand.pdgId())==521) {
	 bool psidau = false;
	 bool kstpdau = false;
	 bool kst1pdau = false;
	 bool pipGDau = false;
	 bool kpGDau = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   if ( genDau->pdgId()==443 ) psidau = true;
	   if ( abs(genDau->pdgId())==323 ) {
	     kstpdau = true;
	     for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	       const Candidate * genGDau = genDau->daughter(j);
	       if ( abs(genGDau->pdgId())==211 ) pipGDau = true;
	       if ( abs(genGDau->pdgId())==321 ) kpGDau = true;
	     }    
	   }
	   if ( abs(genDau->pdgId())==10323 ) kst1pdau = true;
	 }
	 if ( psidau && kstpdau && pipGDau ) feeddown = 6; 
	 if ( psidau && kstpdau && kpGDau ) feeddown = 7; 
	 if ( psidau && kst1pdau ) feeddown = 8; 
       }

       // check deeper truth match for KEVIN
       
       const Candidate *v0(0);
       const Candidate *psi(0);
       bool foundV0 = false;
       bool foundPsi = false;
       if ( ( abs(BCand.pdgId())==511 || abs(BCand.pdgId())==521 || 
	      abs(BCand.pdgId())==531 || abs(BCand.pdgId())==5122 ) &&
	    genKsPsi != 1 ) {
	 // loop through daughters to search for JPsi (443) or V0 (310, 3122)
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   const Candidate * genDau = BCand.daughter(i);
	   if ( genDau->pdgId()==443 ) {psi = genDau; foundPsi = true;}
	   if ( genDau->pdgId()==310||abs(genDau->pdgId())==5122 ) {v0 = genDau; foundV0 = true;}
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     if ( genGDau->pdgId()==443 ) {psi = genGDau; foundPsi = true;}
	     if ( genGDau->pdgId()==310||abs(genGDau->pdgId())==5122 ) {v0 = genGDau; foundV0 = true;}
	     for ( uint k = 0; k < genGDau->numberOfDaughters(); k++){
	       const Candidate * genGGDau = genGDau->daughter(k);
	       if ( genGGDau->pdgId()==443 ) {psi = genGGDau; foundPsi = true;}
	       if ( genGGDau->pdgId()==310||abs(genGGDau->pdgId())==5122 ) {v0 = genGGDau; foundV0 = true;}
	       for ( uint l = 0; l < genGGDau->numberOfDaughters(); l++){
		 const Candidate * genGGGDau = genGGDau->daughter(l);
		 if ( genGGGDau->pdgId()==443 ) {psi = genGGGDau; foundPsi = true;}
		 if ( genGGGDau->pdgId()==310||abs(genGGGDau->pdgId())==5122 ) {v0 = genGGGDau; foundV0 = true;}
	       }
	     }
	   }
	 }
       }
       if (foundV0&&foundPsi) {
	 fillV0(*v0);
	 fillPsi(*psi);
       }
     }
     
     if ( truthMatch.size()==0 ) { // if no truth match to signal has been found, fill with zeros
       for (uint i = 0; i<mupPx->size(); i++) {
	 truthPsi.push_back(0);
	 truthKs.push_back(0);
	 truthMatch.push_back(0);
       }
     }
     
     
   }

   //Check genParticles for every event to get number of generated signal and pT values
   if (doMC_) {
     genBPt = -1; genHLT_2mu0 = -1; genHLT_2mu3 = -1; genHLT_L12muOpen = -1;
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
         // only check for signal decay after possible B0 B0bar oscilation
         int ipsi(-1), iks(-1);
         for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   int imu1(-1), imu2(-1),  ipi1(-1), ipi2(-1);
	   bool wrong = false;
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     if ( genGDau->pdgId()==13 && abs(genDau->pdgId())==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && abs(genDau->pdgId())==443 ) imu2 = j;	       
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
	     for ( uint m = 0; m < genGDau->numberOfDaughters(); m++){
	       const Candidate * genGGDau = genGDau->daughter(m);
	       if ( genGGDau->pdgId()==211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi1 = m;
	       if ( genGGDau->pdgId()==-211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi2 = m;
	       if ( genGDau->pdgId()==310 && abs(genGGDau->pdgId())!=211 && genGGDau->pdgId()!=22 ) wrong = true; 
	     }
	   }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=311 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi1!=-1 && ipi2!=-1&&!wrong) iks = i;
         }
         if (ipsi!=-1&&iks!=-1) {
	   genBPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
           genHLT_2mu0 = hlt_2mu0;
	   genHLT_2mu3 = hlt_2mu3;
	   genHLT_L12muOpen = hlt_L12muOpen;
         }
       }
     }

     if (genBPt>0) treeMC_->Fill();
   }
   
   //fill the tree and clear the vectors
   if (nB > 0 ) {
     cout << "filling tree" << endl;
     tree_->Fill();
   }
   l1_mu3 = 0; l1_2mu3 = 0; l1_muOpen = 0; l1_mu0 = 0;
   hlt_mu3 = 0; hlt_mu5 = 0; hlt_mu7 = 0; hlt_mu9 = 0; hlt_2mu0 = 0; hlt_2mu3 = 0; hlt_2mu3JPsi = 0;
   hltBJPsiMuMu = 0; hlt_L1muOpen = 0; hlt_L12muOpen = 0;
   nB = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   priRfVtxZE->clear(); priRfVtxCL->clear(); priRfNTrkDif->clear();
   bMass->clear(); bVtxCL->clear(); bPx->clear(); bPy->clear(); bPz->clear(); 
   bPxE->clear(); bPyE->clear(); bPzE->clear(); 
   bctau->clear(); bctau2D->clear(); bctauBS->clear(); bctauMPV->clear(); 
   bctauRf->clear(); bctauMPVRf->clear(); bctauMPVBS->clear(); 
   bctauE->clear(); bctau2DE->clear(); bctauBSE->clear(); bctauMPVE->clear(); 
   bctauRfE->clear(); bctauMPVRfE->clear(); bctauMPVBSE->clear(); 
   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
   bResMass->clear(); bResVtxCL->clear(); bResPx->clear(); bResPy->clear(); bResPz->clear(); 
   bResDecayVtxX->clear(); bResDecayVtxY->clear(); bResDecayVtxZ->clear();
   bResDecayVtxXE->clear(); bResDecayVtxYE->clear(); bResDecayVtxZE->clear();
   VMass->clear(); VVtxCL->clear(); VPx->clear(); VPy->clear(); VPz->clear();
   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   JMass->clear(); JVtxCL->clear(); JPx->clear(); JPy->clear(); JPz->clear();
   JDecayVtxX->clear(); JDecayVtxY->clear(); JDecayVtxZ->clear();
   JDecayVtxXE->clear(); JDecayVtxYE->clear(); JDecayVtxZE->clear(); JmuOL->clear();
   mumPx->clear(); mumPy->clear(); mumPz->clear(); mumD0->clear(); mumD0E->clear(); mumC2->clear();
   mumCat->clear(); mumME1->clear(); mumAngT->clear(); mumNHits->clear(); mumNPHits->clear();
   mumTrigL1Open->clear(); mumTrig2mu0->clear(); mumTrig2mu3->clear(); mumTrigL1OpenD->clear(); mumTrig2mu0D->clear(); mumTrig2mu3D->clear();
   mupPx->clear(); mupPy->clear(); mupPz->clear(); mupD0->clear(); mupD0E->clear(); mupC2->clear();
   mupCat->clear(); mupME1->clear(); mupAngT->clear(); mupNHits->clear(); mupNPHits->clear();
   mupTrigL1Open->clear(); mupTrig2mu0->clear(); mupTrig2mu3->clear(); mupTrigL1OpenD->clear(); mupTrig2mu0D->clear(); mupTrig2mu3D->clear();
   VTrkpMass->clear(); VTrkpPx->clear(); VTrkpPy->clear(); VTrkpPz->clear(); 
   VTrkpD0->clear(); VTrkpD0E->clear();
   VTrkmMass->clear(); VTrkmPx->clear(); VTrkmPy->clear(); VTrkmPz->clear();
   VTrkmD0->clear(); VTrkmD0E->clear();
   bResTrkPx->clear(); bResTrkPy->clear(); bResTrkPz->clear();
   bResTrkD0->clear(); bResTrkD0E->clear(); bResTrkChg->clear();
   genKsPsi = 0; genKstarpPsi = 0; genLambdaPsi = 0; feedup = 0; feeddown = 0;
   truePriVtxX = 0; truePriVtxY = 0; truePriVtxZ = 0; trueBPx = 0; trueBPy = 0; trueBPz = 0;
   trueBDecayVtxX = 0; trueBDecayVtxY = 0; trueBDecayVtxZ = 0; trueBResPx = 0; trueBResPy = 0; trueBResPz = 0;
   trueBResDecayVtxX = 0; trueBResDecayVtxY = 0; trueBResDecayVtxZ = 0; 
   trueVPx = 0; trueVPy = 0; trueVPz = 0;
   trueVDecayVtxX = 0; trueVDecayVtxY = 0; trueVDecayVtxZ = 0; trueJPx = 0; trueJPy = 0; trueJPz = 0;
   trueJDecayVtxX = 0; trueJDecayVtxY = 0; trueJDecayVtxZ = 0;
   trueMumPx = 0; trueMumPy = 0; trueMumPz = 0; trueMumD0 = 0;
   trueMupPx = 0; trueMupPy = 0; trueMupPz = 0; trueMupD0 = 0;
   trueVTrkpPx = 0; trueVTrkpPy = 0; trueVTrkpPz = 0; trueVTrkpD0 = 0;
   trueVTrkmPx = 0; trueVTrkmPy = 0; trueVTrkmPz = 0; trueVTrkmD0 = 0;
   trueBResTrkPx = 0; trueBResTrkPy = 0; trueBResTrkPz = 0; trueBResTrkD0 = 0; trueBResTrkChg = 0;
   prompt = 0; truthMatch.clear(); truthKs.clear(); truthPsi.clear(); truthMatchPAT.clear(); truthKsPAT.clear(); truthPsiPAT.clear(); 
   
   bMass2->clear(); bMass3->clear(); bMass4->clear(); bMass5->clear(); 
   bPx2->clear(); bPx3->clear(); bPx4->clear(); bPx5->clear(); bDx2->clear(); bDx3->clear(); bDx4->clear(); bDx5->clear(); 
   bPy2->clear(); bPy3->clear(); bPy4->clear(); bPy5->clear(); bDy2->clear(); bDy3->clear(); bDy4->clear(); bDy5->clear(); 
   bPz2->clear(); bPz3->clear(); bPz4->clear(); bPz5->clear(); bDz2->clear(); bDz3->clear(); bDz4->clear(); bDz5->clear(); 
}

void JPsiPAT::fillPsi(const reco::Candidate& genpsi) {
  
  for (uint i=0; i<genpsi.numberOfDaughters(); i++) {
    if (genpsi.daughter(i)->pdgId()==13) { //13 is a mu-
      trueMumPx = genpsi.daughter(i)->px();
      trueMumPy = genpsi.daughter(i)->py();
      trueMumPz = genpsi.daughter(i)->pz();
    }
    if (genpsi.daughter(i)->pdgId()==-13) { //-13 is a mu+
      trueMupPx = genpsi.daughter(i)->px();
      trueMupPy = genpsi.daughter(i)->py();
      trueMupPz = genpsi.daughter(i)->pz();
    }
  }
}

void JPsiPAT::fillV0(const reco::Candidate& genv0) {
  
  for (uint i=0; i<genv0.numberOfDaughters(); i++) {
    if (genv0.daughter(i)->charge()>0 && genv0.numberOfDaughters()==2) {
      trueVTrkpPx = genv0.daughter(i)->px();
      trueVTrkpPy = genv0.daughter(i)->py();
      trueVTrkpPz = genv0.daughter(i)->pz();
      trueVDecayVtxX = genv0.daughter(i)->vx();
      trueVDecayVtxY = genv0.daughter(i)->vy();
      trueVDecayVtxZ = genv0.daughter(i)->vz();
    }
    if (genv0.daughter(i)->charge()<0 && genv0.numberOfDaughters()==2) {
      trueVTrkmPx = genv0.daughter(i)->px();
      trueVTrkmPy = genv0.daughter(i)->py();
      trueVTrkmPz = genv0.daughter(i)->pz();
    }
  }
}

bool const JPsiPAT::HasGoodME11(reco::Muon const& muon, double const dxdzCut) const{

  bool retVal = false;
  for(std::vector<reco::MuonChamberMatch>::const_iterator mcm = muon.matches().begin();
    mcm != muon.matches().end(); ++mcm) {
    DetId const& chamberId = mcm->id;
    if (chamberId.det() != DetId::Muon) continue;
    if (chamberId.subdetId() != MuonSubdetId::CSC) continue;
    CSCDetId id(chamberId.rawId());
    if (id.station() != 1) continue;
    if (fabs(mcm->dXdZ) > dxdzCut) continue;
    retVal = true;
  }
  return retVal;
}


// ------------ method called once each job just before starting event loop  ------------

void 
JPsiPAT::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuplePsi","JPsi to mumu ntuple");
  treeMC_ = fs->make<TTree>("mctruthPsi","jpsiks gen ntuple");

  treeMC_->Branch("genBPt",&genBPt,"genBPt/f");
  treeMC_->Branch("genHLT_2mu0",&genHLT_2mu0,"genHLT_2mu0/i");
  treeMC_->Branch("genHLT_2mu3",&genHLT_2mu3,"genHLT_2mu3/i");
  treeMC_->Branch("genHLT_L12muOpen",&genHLT_L12muOpen,"genHLT_L12muOpen/i");
  tree_->Branch("l1_mu3",&l1_mu3,"l1_mu3/i");
  tree_->Branch("l1_2mu3",&l1_2mu3,"l1_2mu3/i");
  tree_->Branch("l1_muOpen",&l1_muOpen,"l1_muOpen/i");
  tree_->Branch("l1_mu0",&l1_mu0,"l1_mu0/i");
  tree_->Branch("hlt_mu3",&hlt_mu3,"hlt_mu3/i");
  tree_->Branch("hlt_mu5",&hlt_mu5,"hlt_mu5/i");
  tree_->Branch("hlt_mu7",&hlt_mu7,"hlt_mu7/i");
  tree_->Branch("hlt_mu9",&hlt_mu9,"hlt_mu9/i");
  tree_->Branch("hlt_2mu0",&hlt_2mu0,"hlt_2mu0/i");
  tree_->Branch("hlt_2mu3",&hlt_2mu3,"hlt_2mu3/i");
  tree_->Branch("hlt_2mu3JPsi",&hlt_2mu3JPsi,"hlt_2mu3JPsi/i");
  tree_->Branch("hltBJPsiMuMu",&hltBJPsiMuMu,"hltBJPsiMuMu/i");
  tree_->Branch("hlt_L1muOpen",&hlt_L1muOpen,"hlt_L1muOpen/i");
  tree_->Branch("hlt_L12muOpen",&hlt_L12muOpen,"hlt_L12muOpen/i");
  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");
  tree_->Branch("priRfVtxX",&priRfVtxX);
  tree_->Branch("priRfVtxY",&priRfVtxY);
  tree_->Branch("priRfVtxZ",&priRfVtxZ);
  tree_->Branch("priRfVtxXE",&priRfVtxXE);
  tree_->Branch("priRfVtxYE",&priRfVtxYE);
  tree_->Branch("priRfVtxZE",&priRfVtxZE);
  tree_->Branch("priRfVtxCL",&priRfVtxCL);
  tree_->Branch("priRfNTrkDif",&priRfNTrkDif);
  tree_->Branch("bMass",&bMass);
  tree_->Branch("bVtxCL",&bVtxCL);
  tree_->Branch("bPx",&bPx);
  tree_->Branch("bPy",&bPy);
  tree_->Branch("bPz",&bPz);
  tree_->Branch("bPxE",&bPxE);
  tree_->Branch("bPyE",&bPyE);
  tree_->Branch("bPzE",&bPzE);
  tree_->Branch("bctau",&bctau);
  tree_->Branch("bctau2D",&bctau2D);
  tree_->Branch("bctauBS",&bctauBS);
  tree_->Branch("bctauMPV",&bctauMPV);
  tree_->Branch("bctauMPVBS",&bctauMPVBS);
  tree_->Branch("bctauRf",&bctauRf);
  tree_->Branch("bctauMPVRf",&bctauMPVRf);
  tree_->Branch("bctauE",&bctauE);
  tree_->Branch("bctau2DE",&bctau2DE);
  tree_->Branch("bctauBSE",&bctauBSE);
  tree_->Branch("bctauMPVE",&bctauMPVE);
  tree_->Branch("bctauMPVBSE",&bctauMPVBSE);
  tree_->Branch("bctauRfE",&bctauRfE);
  tree_->Branch("bctauMPVRfE",&bctauMPVRfE);
  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bResMass",&bResMass);
  tree_->Branch("bResVtxCL",&bResVtxCL);
  tree_->Branch("bResPx",&bResPx);
  tree_->Branch("bResPy",&bResPy);
  tree_->Branch("bResPz",&bResPz);
  tree_->Branch("bResDecayVtxX",&bResDecayVtxX);
  tree_->Branch("bResDecayVtxY",&bResDecayVtxY);
  tree_->Branch("bResDecayVtxZ",&bResDecayVtxZ);
  tree_->Branch("bResDecayVtxXE",&bResDecayVtxXE);
  tree_->Branch("bResDecayVtxYE",&bResDecayVtxYE);
  tree_->Branch("bResDecayVtxZE",&bResDecayVtxZE);
  tree_->Branch("VMass",&VMass);
  tree_->Branch("VVtxCL",&VVtxCL);
  tree_->Branch("VPx",&VPx);
  tree_->Branch("VPy",&VPy);
  tree_->Branch("VPz",&VPz);
  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("JMass",&JMass);
  tree_->Branch("JVtxCL",&JVtxCL);
  tree_->Branch("JPx",&JPx);
  tree_->Branch("JPy",&JPy);
  tree_->Branch("JPz",&JPz);
  tree_->Branch("JDecayVtxX",&JDecayVtxX);
  tree_->Branch("JDecayVtxY",&JDecayVtxY);
  tree_->Branch("JDecayVtxZ",&JDecayVtxZ);
  tree_->Branch("JDecayVtxXE",&JDecayVtxXE);
  tree_->Branch("JDecayVtxYE",&JDecayVtxYE);
  tree_->Branch("JDecayVtxZE",&JDecayVtxZE);
  tree_->Branch("JmuOL",&JmuOL);
  tree_->Branch("mumPx",&mumPx);
  tree_->Branch("mumPy",&mumPy);
  tree_->Branch("mumPz",&mumPz);
  tree_->Branch("mumD0",&mumD0);
  tree_->Branch("mumD0E",&mumD0E);
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumCat",&mumCat);
  tree_->Branch("mumME1",&mumME1);
  tree_->Branch("mumAngT",&mumAngT);
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mumTrigL1Open",&mumTrigL1Open);
  tree_->Branch("mumTrig2mu0",&mumTrig2mu0);  
  tree_->Branch("mumTrig2mu3",&mumTrig2mu3);
  tree_->Branch("mumTrigL1OpenD",&mumTrigL1OpenD);
  tree_->Branch("mumTrig2mu0D",&mumTrig2mu0D);  
  tree_->Branch("mumTrig2mu3D",&mumTrig2mu3D);
  tree_->Branch("mupTrigL1Open",&mupTrigL1Open);
  tree_->Branch("mupTrig2mu0",&mupTrig2mu0);  
  tree_->Branch("mupTrig2mu3",&mupTrig2mu3);
  tree_->Branch("mupTrigL1OpenD",&mupTrigL1OpenD);
  tree_->Branch("mupTrig2mu0D",&mupTrig2mu0D);  
  tree_->Branch("mupTrig2mu3D",&mupTrig2mu3D);
  tree_->Branch("mupPx",&mupPx);
  tree_->Branch("mupPy",&mupPy);
  tree_->Branch("mupPz",&mupPz);
  tree_->Branch("mupD0",&mupD0);
  tree_->Branch("mupD0E",&mupD0E);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupCat",&mupCat);
  tree_->Branch("mupME1",&mupME1);
  tree_->Branch("mupAngT",&mupAngT);
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("VTrkpTrk1Mass",&VTrkpMass);
  tree_->Branch("VTrkpPx",&VTrkpPx);
  tree_->Branch("VTrkpPy",&VTrkpPy);
  tree_->Branch("VTrkpPz",&VTrkpPz);
  tree_->Branch("VTrkpD0",&VTrkpD0);
  tree_->Branch("VTrkpD0E",&VTrkpD0E);
  tree_->Branch("VTrkmMass",&VTrkmMass);
  tree_->Branch("VTrkmPx",&VTrkmPx);
  tree_->Branch("VTrkmPy",&VTrkmPy);
  tree_->Branch("VTrkmPz",&VTrkmPz);
  tree_->Branch("VTrkmD0",&VTrkmD0);
  tree_->Branch("VTrkmD0E",&VTrkmD0E);
  tree_->Branch("bResTrkPx",&bResTrkPx);
  tree_->Branch("bResTrkPy",&bResTrkPy);
  tree_->Branch("bResTrkPz",&bResTrkPz);
  tree_->Branch("bResTrkD0",&bResTrkD0);
  tree_->Branch("bResTrkD0E",&bResTrkD0E);
  tree_->Branch("bResTrkChg",&bResTrkChg);
  tree_->Branch("genKsPsi", &genKsPsi, "genKsPsi/I");
  tree_->Branch("genKstarpPsi", &genKstarpPsi, "genKstarpPsi/I");
  tree_->Branch("genLambdaPsi", &genLambdaPsi, "genLambdaPsi/I");
  tree_->Branch("feedup", &feedup, "feedup/I");
  tree_->Branch("feeddown", &feeddown, "feeddown/I");

  tree_->Branch("bMass2",&bMass2);
  tree_->Branch("bMass3",&bMass3);
  tree_->Branch("bMass4",&bMass4);
  tree_->Branch("bMass5",&bMass5);
  tree_->Branch("bPx2",&bPx2);
  tree_->Branch("bPx3",&bPx3);
  tree_->Branch("bPx4",&bPx4);
  tree_->Branch("bPx5",&bPx5);
  tree_->Branch("bPy2",&bPy2);
  tree_->Branch("bPy3",&bPy3);
  tree_->Branch("bPy4",&bPy4);
  tree_->Branch("bPy5",&bPy5);
  tree_->Branch("bPz2",&bPz2);
  tree_->Branch("bPz3",&bPz3);
  tree_->Branch("bPz4",&bPz4);
  tree_->Branch("bPz5",&bPz5);
  tree_->Branch("bDx2",&bDx2);
  tree_->Branch("bDx3",&bDx3);
  tree_->Branch("bDx4",&bDx4);
  tree_->Branch("bDx5",&bDx5);
  tree_->Branch("bDy2",&bDy2);
  tree_->Branch("bDy3",&bDy3);
  tree_->Branch("bDy4",&bDy4);
  tree_->Branch("bDy5",&bDy5);
  tree_->Branch("bDz2",&bDz2);
  tree_->Branch("bDz3",&bDz3);
  tree_->Branch("bDz4",&bDz4);
  tree_->Branch("bDz5",&bDz5);


  // do branches for MC truth
  tree_->Branch("truePriVtxX", &truePriVtxX, "truePriVtxX/f");
  tree_->Branch("truePriVtxY", &truePriVtxY, "truePriVtxY/f");
  tree_->Branch("truePriVtxZ", &truePriVtxZ, "truePriVtxZ/f");
  tree_->Branch("trueBPx",&trueBPx, "trueBPx/f");
  tree_->Branch("trueBPy",&trueBPy, "trueBPy/f");
  tree_->Branch("trueBPz",&trueBPz, "trueBPz/f");
  tree_->Branch("trueBDecayVtxX",&trueBDecayVtxX, "trueBDecayVtxX/f");
  tree_->Branch("trueBDecayVtxY",&trueBDecayVtxY, "trueBDecayVtxY/f");
  tree_->Branch("trueBDecayVtxZ",&trueBDecayVtxZ, "trueBDecayVtxZ/f");
  tree_->Branch("trueBResPx",&trueBResPx, "trueBResPx/f");
  tree_->Branch("trueBResPy",&trueBResPy, "trueBResPy/f"); 
  tree_->Branch("trueBResPz",&trueBResPz, "trueBResPz/f");
  tree_->Branch("trueBResDecayVtxX",&trueBResDecayVtxX, "trueBResDecayVtxX/f");
  tree_->Branch("trueBResDecayVtxY",&trueBResDecayVtxY, "trueBResDecayVtxY/f");
  tree_->Branch("trueBResDecayVtxZ",&trueBResDecayVtxZ, "trueBResDecayVtxZ/f");
  tree_->Branch("trueVPx",&trueVPx, "trueVPx/f"); 
  tree_->Branch("trueVPy",&trueVPy, "trueVPy/f"); 
  tree_->Branch("trueVPz",&trueVPz, "trueVPz/f");
  tree_->Branch("trueVDecayVtxX",&trueVDecayVtxX, "trueVDecayVtxX/f"); 
  tree_->Branch("trueVDecayVtxY",&trueVDecayVtxY, "trueVDecayVtxY/f"); 
  tree_->Branch("trueVDecayVtxZ",&trueVDecayVtxZ, "trueVDecayVtxZ/f");
  tree_->Branch("trueJPx",&trueJPx, "trueJPx/f");
  tree_->Branch("trueJPy",&trueJPy, "trueJPy/f"); 
  tree_->Branch("trueJPz",&trueJPz, "trueJPz/f");
  tree_->Branch("trueJDecayVtxX",&trueJDecayVtxX, "trueJDecayVtxX/f"); 
  tree_->Branch("trueJDecayVtxY",&trueJDecayVtxY, "trueJDecayVtxY/f"); 
  tree_->Branch("trueJDecayVtxZ",&trueJDecayVtxZ, "trueJDecayVtxZ/f");
  tree_->Branch("trueMumPx",&trueMumPx, "trueMumPx/f"); 
  tree_->Branch("trueMumPy",&trueMumPy, "trueMumPy/f"); 
  tree_->Branch("trueMumPz",&trueMumPz, "trueMumPz/f"); 
  tree_->Branch("trueMumD0",&trueMumD0, "trueMumD0/f");
  tree_->Branch("trueMupPx",&trueMupPx, "trueMupPx/f"); 
  tree_->Branch("trueMupPy",&trueMupPy, "trueMupPy/f"); 
  tree_->Branch("trueMupPz",&trueMupPz, "trueMupPz/f"); 
  tree_->Branch("trueMupD0",&trueMupD0, "trueMupD0/f");
  tree_->Branch("trueVTrkpPx",&trueVTrkpPx, "trueVTrkpPx/f"); 
  tree_->Branch("trueVTrkpPy",&trueVTrkpPy, "trueVTrkpPy/f"); 
  tree_->Branch("trueVTrkpPz",&trueVTrkpPz, "trueVTrkpPz/f"); 
  tree_->Branch("trueVTrkpD0",&trueVTrkpD0, "trueVTrkpD0/f");
  tree_->Branch("trueVTrkmPx",&trueVTrkmPx, "trueVTrkmPx/f"); 
  tree_->Branch("trueVTrkmPy",&trueVTrkmPy, "trueVTrkmPy/f"); 
  tree_->Branch("trueVTrkmPz",&trueVTrkmPz, "trueVTrkmPz/f"); 
  tree_->Branch("trueVTrkmD0",&trueVTrkmD0, "trueVTrkmD0/f");
  tree_->Branch("trueBResTrkPx",&trueBResTrkPx, "trueBResTrkPx/f"); 
  tree_->Branch("trueBResTrkPy",&trueBResTrkPy, "trueBResTrkPy/f"); 
  tree_->Branch("trueBResTrkPz",&trueBResTrkPz, "trueBResTrkPz/f"); 
  tree_->Branch("trueBResTrkD0",&trueBResTrkD0, "trueBResTrkD0/f");
  tree_->Branch("trueBResTrkChg",&trueBResTrkChg, "trueBResTrkChg/I");
  tree_->Branch("prompt",&prompt, "prompt/I");
  tree_->Branch("truthMatch",&truthMatch);
  tree_->Branch("truthKs",&truthKs);
  tree_->Branch("truthPsi",&truthPsi);
  tree_->Branch("truthMatchPAT",&truthMatchPAT);
  tree_->Branch("truthKsPAT",&truthKsPAT);
  tree_->Branch("truthPsiPAT",&truthPsiPAT);
}


// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiPAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
  treeMC_->GetDirectory()->cd();
  treeMC_->Write();  
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiPAT);

