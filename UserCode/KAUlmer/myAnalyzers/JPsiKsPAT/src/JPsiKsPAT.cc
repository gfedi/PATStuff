// -*- C++ -*-
//
// Package:    JPsiKsPAT
// Class:      JPsiKsPAT
// 
/**\class JPsiKsPAT JPsiKsPAT.cc myAnalyzers/JPsiKsPAT/src/JPsiKsPAT.cc

 Description: <one line class summary>
Make rootTuple for b->s JPsi(mu+mu-) analyses

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer
//         Created:  Wed May  7 13:15:04 MDT 2008
// $Id: JPsiKsPAT.cc,v 1.33 2011/05/08 17:59:30 kaulmer Exp $
//
//


// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/interface/JPsiKsPAT.h"

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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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
JPsiKsPAT::JPsiKsPAT(const edm::ParameterSet& iConfig)
  :
  hlTriggerResults_ (iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) ),
  vtxSample ( iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices")) ),
  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  v0Collection_ ( iConfig.getUntrackedParameter<std::string>("V0Collection",std::string("generalV0Candidates")) ),
  muonType ( iConfig.getUntrackedParameter<std::string>("MuonType",std::string("cleanPatMuons")) ),
  muonTypeForPAT ( iConfig.getUntrackedParameter<std::string>("MuonTypeForPAT",std::string("muons")) ),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  onlyCount_ ( iConfig.getUntrackedParameter<bool>("onlyCount",false) ),

  tree_(0), treeMC_(0), treeMC2_(0), treeMC3_(0), treeMC4_(0), treeMC5_(0), l1_mu3(0), l1_2mu3(0), l1_muOpen(0), l1_mu0(0),
  hlt_mu3(0), hlt_mu5(0), hlt_mu7(0), hlt_mu9(0), hlt_2mu0(0), hlt_2mu0_quark(0), hlt_2mu3(0), hlt_2mu0L2(0), hlt_2mu3JPsi(0), hlt_BJPsiMuMu(0), 
  hlt_mu0trk0(0), hlt_mu3trk0(0), hlt_mu0trkmu0(0), hlt_mu3trkmu0(0), hlt_mu0trkmu0OST(0), hlt_mu3trkmu0OST(0), hlt_mu0trkmu0OST_tight(0),
  hlt_L1muOpen(0), hlt_L12muOpen(0), hlt_L12muOpenTight(0),
  nB(0), priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priRfVtxX(0), priRfVtxY(0), priRfVtxZ(0), priRfVtxXE(0), priRfVtxYE(0), priRfVtxZE(0), priRfVtxCL(0), 
  priRfNTrkDif(0),
  bMass(0), bVtxCL(0), bPx(0), bPy(0), bPz(0), bPxE(0), bPyE(0), bPzE(0), 
  bctau(0), bctau2D(0), bctauBS(0), bctauMPV(0), bctauMPV2D(0), bctauMPVBS(0), bctauRf(0), bctauMPVRf(0),
  bctauE(0), bctau2DE(0), bctauBSE(0), bctauMPVE(0), bctauMPV2DE(0), bctauMPVBSE(0), bctauRfE(0), bctauMPVRfE(0),
  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bResMass(0), bResVtxCL(0), bResPx(0), bResPy(0), bResPz(0),
  bResDecayVtxX(0), bResDecayVtxY(0), bResDecayVtxZ(0), bResDecayVtxXE(0), bResDecayVtxYE(0), bResDecayVtxZE(0),
  VMass(0), VVtxCL(0), VPx(0), VPy(0), VPz(0),
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0),
  VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  JMass(0), JVtxCL(0), JPx(0), JPy(0), JPz(0),
  JDecayVtxX(0), JDecayVtxY(0), JDecayVtxZ(0), JDecayVtxXE(0), JDecayVtxYE(0), JDecayVtxZE(0), JmuOL(0),
  mumPx(0), mumPy(0), mumPz(0), mumD0(0), mumD0E(0), mumC2(0), mumCat(0), mumME1(0), mumAngT(0), mumNHits(0), mumNPHits(0),
  mumTrigL1Open1mu(0), mumTrigL1Open2mu(0), mumTrigL1Open2muTight(0), mumTrig2mu0(0), mumTrig2mu3(0), mumTrig2mu0_quark(0),
  mumTrig2mu0L2(0), mumTrigmu0trk0(0), mumTrigmu3trk0(0), mumTrigmu0trkmu0(0), mumTrigmu3trkmu0(0), 
  mupPx(0), mupPy(0), mupPz(0), mupD0(0), mupD0E(0), mupC2(0), mupCat(0), mupME1(0), mupAngT(0), mupNHits(0), mupNPHits(0),
  mupTrigL1Open1mu(0), mupTrigL1Open2mu(0), mupTrigL1Open2muTight(0), mupTrig2mu0(0), mupTrig2mu3(0), mupTrig2mu0_quark(0),
  mupTrig2mu0L2(0), mupTrigmu0trk0(0), mupTrigmu3trk0(0), mupTrigmu0trkmu0(0), mupTrigmu3trkmu0(0), 
  VTrkpMass(0), VTrkpPx(0), VTrkpPy(0), VTrkpPz(0), 
  VTrkpD0(0), VTrkpD0E(0), 
  VTrkmMass(0), VTrkmPx(0), VTrkmPy(0), VTrkmPz(0), 
  VTrkmD0(0), VTrkmD0E(0), 
  bResTrkPx(0), bResTrkPy(0), bResTrkPz(0), 
  bResTrkD0(0), bResTrkD0E(0),bResTrkChg(0), 
  genKsPsi(0), genKsPsi2(0), genKstarpPsi(0), genLambdaPsi(0), feedup(0), feeddown(0),

  bMass2(0), bMass3(0), bMass4(0), bMass5(0), bPx2(0), bPx3(0), bPx4(0), bPx5(0), bPy2(0), bPy3(0), bPy4(0), bPy5(0), 
  bPz2(0), bPz3(0), bPz4(0), bPz5(0), 
  bDx2(0), bDx3(0), bDx4(0), bDx5(0), bDy2(0), bDy3(0), bDy4(0), bDy5(0), bDz2(0), bDz3(0), bDz4(0), bDz5(0),

  truePriVtxX(0), truePriVtxY(0), truePriVtxZ(0), trueBPx(0), trueBPy(0), trueBPz(0), trueBDecayVtxX(0), trueBDecayVtxY(0), trueBDecayVtxZ(0),
  trueBResPx(0), trueBResPy(0), trueBResPz(0), trueBResDecayVtxX(0), trueBResDecayVtxY(0), trueBResDecayVtxZ(0),
  trueVPx(0), trueVPy(0), trueVPz(0), trueVDecayVtxX(0), trueVDecayVtxY(0), trueVDecayVtxZ(0),
  trueJPx(0), trueJPy(0), trueJPz(0), trueJDecayVtxX(0), trueJDecayVtxY(0), trueJDecayVtxZ(0),
  trueMumPx(0), trueMumPy(0), trueMumPz(0), trueMumD0(0), trueMupPx(0), trueMupPy(0), trueMupPz(0), trueMupD0(0),
  genJPt(0), genJP(0),
  trueVTrkpPx(0), trueVTrkpPy(0), trueVTrkpPz(0), trueVTrkpD0(0),
  trueVTrkmPx(0), trueVTrkmPy(0), trueVTrkmPz(0), trueVTrkmD0(0),
  trueBResTrkPx(0), trueBResTrkPy(0), trueBResTrkPz(0), trueBResTrkD0(0), trueBResTrkChg(0),
  truthMatch(0), truthKs(0), truthPsi(0), truthMatchPAT(0), truthKsPAT(0), truthPsiPAT(0), prompt(0)

{
   //now do what ever initialization is needed
}


JPsiKsPAT::~JPsiKsPAT()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void


JPsiKsPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using std::vector;
   using namespace edm;
   using namespace reco;
   using namespace std;

   bool debug = false;

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
     std::string const &trig = std::string("TriggerResults::")+hlTriggerResults_;
     iEvent.getByLabel(edm::InputTag(trig),hltresults);
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
     const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
     
     for (int itrig=0; itrig< ntrigs; itrig++) {
       TString trigName = triggerNames_.triggerName(itrig);
       int hltflag = (*hltresults)[itrig].accept();
       //int wasrun  = (*hltresults)[itrig].wasrun();
       //cout << "Trigger " <<  trigName << " was passed = " <<  hltflag << endl;
       if (trigName=="HLT_DoubleMu3_BJPsi") hlt_BJPsiMuMu = hltflag;
       if (trigName=="HLT_Mu3") hlt_mu3 = hltflag;
       if (trigName=="HLT_Mu5") hlt_mu5 = hltflag;
       if (trigName=="HLT_Mu7") hlt_mu7 = hltflag;
       if (trigName=="HLT_Mu9") hlt_mu9 = hltflag;  
       if (trigName=="HLT_DoubleMu0") hlt_2mu0 = hltflag;
       if (trigName=="HLT_DoubleMu0_Quarkonium_v1") hlt_2mu0_quark = hltflag;
       if (trigName=="HLT_DoubleMu3") hlt_2mu3 = hltflag;
       if (trigName=="HLT_DoubleMu3_Jpsi") hlt_2mu3JPsi = hltflag;
       if (trigName=="HLT_Mu0_Track0_Jpsi") hlt_mu0trk0 = hltflag;
       if (trigName=="HLT_Mu3_Track0_Jpsi") hlt_mu3trk0 = hltflag;
       if (trigName=="HLT_Mu0_TkMu0_Jpsi") hlt_mu0trkmu0 = hltflag;
       if (trigName=="HLT_Mu3_TkMu0_Jpsi") hlt_mu3trkmu0 = hltflag;
       if (trigName=="HLT_Mu0_TkMu0_OST_Jpsi") hlt_mu0trkmu0OST = hltflag;
       if (trigName=="HLT_Mu3_TkMu0_OST_Jpsi") hlt_mu3trkmu0OST = hltflag;
       
       if (trigName=="HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1") hlt_mu0trkmu0OST_tight = hltflag;
       
       if (trigName=="HLT_L1MuOpen") hlt_L1muOpen = hltflag;
       if (trigName=="HLT_L1DoubleMuOpen") hlt_L12muOpen = hltflag;
       if (trigName=="HLT_L1DoubleMuOpen_Tight") hlt_L12muOpenTight = hltflag;
       if (trigName=="HLT_L2DoubleMu0 ") hlt_2mu0L2 = hltflag;
     }
   }

   //declare some things before the !onlyCount_ block
   	       
   float mb = 5.27953;	
   float pi = 3.14159265;
   Handle<vector<VertexCompositeCandidate> > theV0Handle;
   Handle< vector<pat::GenericParticle> > thePATTrackHandle;
   Handle< vector<pat::Muon> > thePATMuonHandle;
   reco::Vertex bestVtx;
   reco::Vertex bestVtxBS;
   if (!onlyCount_) {
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
     bestVtx = *(recVtxs->begin());

     priVtxX = bestVtx.x();
     priVtxY = bestVtx.y();
     priVtxZ = bestVtx.z();
     priVtxXE = bestVtx.xError();
     priVtxYE = bestVtx.yError();
     priVtxZE = bestVtx.zError();
     priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

     //get primary with beamspot constraint
     Handle<reco::VertexCollection> recVtxsBS;
     iEvent.getByLabel("offlinePrimaryVerticesWithBS", recVtxsBS);
   
     nVtxTrks = 0;
     bestVtxBS = *(recVtxsBS->begin());
   
     //get beamspot
     reco::BeamSpot beamSpot;
     edm::Handle<reco::BeamSpot> beamSpotHandle;
     iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
     if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
     else cout << "No beam spot available from EventSetup" << endl;
   
     //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     // try reconstruction without fitting modules
     //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
     iEvent.getByLabel(v0Collection_, "Kshort", theV0Handle);

     iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);

//     iEvent.getByLabel(muonType, thePATMuonHandle);
     iEvent.getByLabel("cleanPatMuonsTriggerMatch", thePATMuonHandle);


     Handle< vector<pat::Muon> > thePATMuonHandle2;
     iEvent.getByLabel("cleanPatMuonsTriggerMatch", thePATMuonHandle2);

     if (debug) cout << "event has " << theV0Handle->size() << "Ks and " << thePATMuonHandle->size() << "muons." << endl; 
   }
   if ( !onlyCount_ && theV0Handle->size()>0 && thePATMuonHandle->size()>=2 ) {
     for ( vector<VertexCompositeCandidate>::const_iterator iVee = theV0Handle->begin();
	   iVee != theV0Handle->end(); ++iVee ) {
       //get Ks tracks from V0 candidate
       vector<RecoChargedCandidate> v0daughters;
       vector<TrackRef> theDaughterTracks;
       v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
       			(iVee->daughter(0))) );
       v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
       			(iVee->daughter(1))) );

       for(unsigned int j = 0; j < v0daughters.size(); ++j) {
        theDaughterTracks.push_back(v0daughters[j].track());
       }

       pat::GenericParticle patTrack1;
       pat::GenericParticle patTrack2;
       //lets loop through the pat tracks to find the match for this reco::track
       for ( vector<pat::GenericParticle>::const_iterator iTrack = thePATTrackHandle->begin();
	     iTrack != thePATTrackHandle->end(); ++iTrack ) {
	 
	 // how to access the reco::Track object
	 TrackRef hadTrack = iTrack->track();
	 if ( hadTrack.isNull() ) {
	   cout << "continue due to no track ref" << endl;
	   continue;
	 }
	 if ( hadTrack==theDaughterTracks[0] )
	   patTrack1 = *iTrack;
	 if ( hadTrack==theDaughterTracks[1] )
	   patTrack2 = *iTrack;
       }

       for ( std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin();
	     iMuonP != thePATMuonHandle->end(); ++iMuonP ) {
	 //check for mu+ first
	 if (iMuonP->charge() == 1) {
	   const pat::Muon *patMuonP = &(*iMuonP);
	   //cout << "for muP: isGlobalMuon = " << iMuonP->isGlobalMuon() << "  isTrackerMuon = " << iMuonP->isTrackerMuon() << "  isStandAloneMuon = " << iMuonP->isStandAloneMuon() << "  isCaloMuon = " << iMuonP->isCaloMuon() << endl;

	   TrackRef glbTrackP = iMuonP->track();
	   //cout << "got track" << endl;
	   if ( glbTrackP.isNull() ) {
	     cout << "continue due to no track ref" << endl;
	     continue;
	   }
	   bool match = false;


	   /////////////////////////
	   //try PAT overlap check//
	   /////////////////////////

	     const reco::CandidatePtrVector & mu1P_overlaps = patTrack1.overlaps(muonTypeForPAT);
	     if ( mu1P_overlaps.size() > 0 ) cout << "patTrack1 overlaps with a muon." << endl;
	     for (size_t i = 0; i < mu1P_overlaps.size(); ++i) {
	       const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu1P_overlaps[i]);
	       if (mu) {
		 // check here that muon match isn't the same as a muon used in the reco...
		 if (mu==patMuonP) cout << "match between patTrack1 and patMuonP" << endl;
	       }
	     }

	     const reco::CandidatePtrVector & mu2P_overlaps = patTrack2.overlaps(muonTypeForPAT);
	     if ( mu2P_overlaps.size() > 0 ) cout << "patTrack2 overlaps with a muon." << endl;
	     for (size_t i = 0; i < mu2P_overlaps.size(); ++i) {
	       const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu2P_overlaps[i]);
	       if (mu) {
		 // check here that muon match isn't the same as a muon used in the reco...
		 if (mu==patMuonP) cout << "match between patTrack2 and patMuonP" << endl;
	       }
	     }
	     
	     ////////////////////////////


	   for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
	     if (glbTrackP->charge() == theDaughterTracks[j]->charge() && 
		 glbTrackP->momentum() == theDaughterTracks[j]->momentum() ) {
	       //std::cout << "Match found between muP and V0 track early with p = " << glbTrackP->momentum() << " and " << theDaughterTracks[j]->momentum() << endl;
	       match = true;
	     }
	     if (match) break;
	   }
	   if (match) { 
	     std::cout << "Match found between muP and V0 track" << endl;
	     match = false;
	     continue; 
	   } // Track is already used in making the V0

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

	       /////////////////////////
	       //try PAT overlap check//	   
	       ////////////////////////
		 
	       const reco::CandidatePtrVector & mu1M_overlaps = patTrack1.overlaps(muonTypeForPAT);
	       if ( mu1M_overlaps.size() > 0 ) cout << "patTrack1 overlaps with a muonM." << endl;
	       for (size_t i = 0; i < mu1M_overlaps.size(); ++i) {
		 const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu1M_overlaps[i]);
		 if (mu) {
		   // check here that muon match isn't the same as a muon used in the reco...
		   if (mu==patMuonM) cout << "match between patTrack1 and patMuonM" << endl;
		 }
	       }
		 
	       const reco::CandidatePtrVector & mu2M_overlaps = patTrack2.overlaps(muonTypeForPAT);
	       if ( mu2M_overlaps.size() > 0 ) cout << "patTrack2 overlaps with a muon." << endl;
	       for (size_t i = 0; i < mu2M_overlaps.size(); ++i) {
		 const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu2M_overlaps[i]);
		 if (mu) {
		   // check here that muon match isn't the same as a muon used in the reco...
		   if (mu==patMuonM) cout << "match between patTrack2 and patMuonM" << endl;
		 }
	       }
		 
	       ////////////////////////////
		 

	       //check that neither muon is the same track as one of the pion tracks
	       for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
		 if (glbTrackM->charge() == theDaughterTracks[j]->charge() && 
		     glbTrackM->momentum() == theDaughterTracks[j]->momentum() ) {
		   //std::cout << "Match found between muM and V0 track early with p = " << glbTrackM->momentum() << " and " << theDaughterTracks[j]->momentum() << endl;
		   match = true;
		 }
		 if (match) break;
	       }
	       if (match) { 
		 std::cout << "Match found between muM and V0 track" << endl;;
		 match = false;
		 continue; 
	       } // Track is already used in making the V0
	       //have two good oppositely charged muons and 2 pions. try to vertex them
	       if ( debug ) cout << "have 4 good tracks including good oppositely charged muons. " << endl;
	       
	       TransientTrack pion1TT(theDaughterTracks[0], &(*bFieldHandle) );
	       TransientTrack pion2TT(theDaughterTracks[1], &(*bFieldHandle) );
	       TransientTrack muon1TT(glbTrackP, &(*bFieldHandle) );
               TransientTrack muon2TT(glbTrackM, &(*bFieldHandle) );





               //Creating a KinematicParticleFactory
	       KinematicParticleFactoryFromTransientTrack pFactory;
	       
	       //The mass of a muon and the insignificant mass sigma 
	       //to avoid singularities in the covariance matrix.
	       ParticleMass muon_mass = 0.10565837; //pdg mass
	       ParticleMass pion_mass = 0.13957018;
	       ParticleMass ks_mass = 0.497614;
	       ParticleMass psi_mass = 3.096916;
	       float muon_sigma = muon_mass*1.e-6;
	       float pion_sigma = pion_mass*1.e-6;
	       float ks_sigma = 0.000024;
	       
	       //initial chi2 and ndf before kinematic fits.
	       float chi = 0.;
	       float ndf = 0.;
	       vector<RefCountedKinematicParticle> pionParticles;
	       vector<RefCountedKinematicParticle> muonParticles;
	       pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
	       pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
	       muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	       muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));



	       KinematicParticleVertexFitter fitter;   
	       RefCountedKinematicTree ksVertexFitTree;
	       ksVertexFitTree = fitter.fit(pionParticles); 
	       if (!ksVertexFitTree->isValid()) {
		 std::cout << "invalid vertex from the ks vertex fit" << std::endl;
		 continue; 
	       }
	       ksVertexFitTree->movePointerToTheTop();
	       
	       RefCountedKinematicParticle ks_vFit_noMC = ksVertexFitTree->currentParticle();
	       RefCountedKinematicVertex ks_vFit_vertex_noMC = ksVertexFitTree->currentDecayVertex();

	       if ( ks_vFit_vertex_noMC->chiSquared() < 0 ) cout << "negative chisq from ks fit" << endl;	       
	       ksVertexFitTree->movePointerToTheFirstChild();
	       RefCountedKinematicParticle ksPi1 = ksVertexFitTree->currentParticle();
	       ksVertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle ksPi2 = ksVertexFitTree->currentParticle();
	       ksVertexFitTree->movePointerToTheTop();

	       KinematicParameters ksPi1KP = ksPi1->currentState().kinematicParameters();
	       KinematicParameters ksPi2KP = ksPi2->currentState().kinematicParameters();
	       KinematicParameters ksPipKP;
	       KinematicParameters ksPimKP;

	       //this ksPi1KP momentum is defined at the ks fit vertex.

	       if ( ksPi1->currentState().particleCharge() > 0 ) ksPipKP = ksPi1KP;
	       if ( ksPi1->currentState().particleCharge() < 0 ) ksPimKP = ksPi1KP;
	       if ( ksPi2->currentState().particleCharge() > 0 ) ksPipKP = ksPi2KP;
	       if ( ksPi2->currentState().particleCharge() < 0 ) ksPimKP = ksPi2KP;

	       //(x,y,z,p_x,p_y,p_z,m)

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


	       // Do MC for Ks cand and do mass constrained vertex fit
	       // creating the constraint with a small sigma to put in the resulting covariance 
	       // matrix in order to avoid singularities
	       // JPsi mass constraint is applied in the final B fit

	       KinematicParticleFitter csFitterKs;
	       KinematicConstraint * ks_c = new MassKinematicConstraint(ks_mass,ks_sigma);
	       // add mass constraint to the ks fit to do a constrained fit:  
 
	       ksVertexFitTree = csFitterKs.fit(ks_c,ksVertexFitTree);
	       if (!ksVertexFitTree->isValid()){
		 std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		 continue; 
	       }
	       
	       ksVertexFitTree->movePointerToTheTop();
	       RefCountedKinematicParticle ks_vFit_withMC = ksVertexFitTree->currentParticle();

	       vector<RefCountedKinematicParticle> vFitMCParticles;
	       vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	       vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	       vFitMCParticles.push_back(ks_vFit_withMC);

	       MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
	       KinematicConstrainedVertexFitter kcvFitter;
	       RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
	       if (!vertexFitTree->isValid()) {
		 std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		 continue;
	       }

	       vertexFitTree->movePointerToTheTop();
	       RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
	       RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
	       if (!bDecayVertexMC->vertexIsValid()){
		 cout << "B MC fit vertex is not valid" << endl;
		 continue;
	       }
	       
	       if ( bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>1000 ) {
		 if ( bDecayVertexMC->chiSquared()<0 ) cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		 if (debug) cout << " continue from bad chi2 = " << bDecayVertexMC->chiSquared() << endl;
		 continue;
	       }
	       
	       
	       if ( bCandMC->currentState().mass() > 10 ) {
	         if (debug) cout << "continue from bmass > 10 = " << bCandMC->currentState().mass() << endl;
	         continue;
	       }


	       // get children from final B fit
	       vertexFitTree->movePointerToTheFirstChild();
	       RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
	       //cout << "mass of mu1 = " << mu1CandMC->currentState().mass() << " and charge = " << mu1CandMC->currentState().particleCharge() << endl;
	       vertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
	       //cout << "mass of mu2 = " << mu2CandMC->currentState().mass() << " and charge = " << mu2CandMC->currentState().particleCharge() << endl;
	       vertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle ksCandMC = vertexFitTree->currentParticle();
	       //cout << "mass of ks = " << ksCandMC->currentState().mass() << " and charge = " << ksCandMC->currentState().particleCharge() << endl;

	       // get mu+ and mu- momenta from final B fit
	       KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMupKP;
	       KinematicParameters psiMumKP;
	       
	       if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	       if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	       if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
	       if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 
	       
	       mupCategory = getMuCat( *iMuonP );
	       mumCategory = getMuCat( *iMuonM ); 
	       
               const reco::Muon *recoMuonM = patMuonM;
	       const reco::Muon *recoMuonP = patMuonP;

               mumME1Clean = HasGoodME11(*recoMuonM,2.);
               mupME1Clean = HasGoodME11(*recoMuonP,2.);

	       KinematicParameters VCandKP = ksCandMC->currentState().kinematicParameters();
	       
	       cout << "filling new candidate" << endl;
	       
	       // fill candidate variables now
	       bMass->push_back(bCandMC->currentState().mass());
	       bPx->push_back(bCandMC->currentState().globalMomentum().x());
	       bPy->push_back(bCandMC->currentState().globalMomentum().y());
	       bPz->push_back(bCandMC->currentState().globalMomentum().z());
	       
	       bPxE->push_back( sqrt( bCandMC->currentState().kinematicParametersError().matrix()(3,3) ) );
	       bPyE->push_back( sqrt( bCandMC->currentState().kinematicParametersError().matrix()(4,4) ) );
	       bPzE->push_back( sqrt( bCandMC->currentState().kinematicParametersError().matrix()(5,5) ) );

	       //ks_vFit_noMC->currentState().kinematicParametersError().matrix()(3,3);   //7x7 matrix: (x,y,z,p_x,p_y,p_z,m)
	       // the index starts at 0, so the position values are 0,1,2 the momentum values are 3,4,5 and the mass is 6
	       
	       bVtxCL->push_back( ChiSquaredProbability((double)(bDecayVertexMC->chiSquared()),(double)(bDecayVertexMC->degreesOfFreedom())) );
	       bDecayVtxX->push_back((*bDecayVertexMC).position().x());
	       bDecayVtxY->push_back((*bDecayVertexMC).position().y());
	       bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
	       bDecayVtxXE->push_back(sqrt( abs(bDecayVertexMC->error().cxx()) ));
	       bDecayVtxYE->push_back(sqrt( abs(bDecayVertexMC->error().cyy()) ));
	       bDecayVtxZE->push_back(sqrt( abs(bDecayVertexMC->error().czz()) ));

	       VMass->push_back( ks_vFit_noMC->currentState().mass() );
	       VVtxCL->push_back( ChiSquaredProbability((double)(ks_vFit_vertex_noMC->chiSquared()),(double)(ks_vFit_vertex_noMC->degreesOfFreedom())) );
               VDecayVtxX->push_back( ks_vFit_vertex_noMC->position().x() );
	       VDecayVtxY->push_back( ks_vFit_vertex_noMC->position().y() );
	       VDecayVtxZ->push_back( ks_vFit_vertex_noMC->position().z() );
               VDecayVtxXE->push_back( sqrt(ks_vFit_vertex_noMC->error().cxx()) );
	       VDecayVtxYE->push_back( sqrt(ks_vFit_vertex_noMC->error().cyy()) );
	       VDecayVtxZE->push_back( sqrt(ks_vFit_vertex_noMC->error().czz()) );
	       VPx->push_back( ks_vFit_noMC->currentState().globalMomentum().x() );
	       VPy->push_back( ks_vFit_noMC->currentState().globalMomentum().y() );
	       VPz->push_back( ks_vFit_noMC->currentState().globalMomentum().z() );   
	       
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

	       VDecayVtxX->push_back( ks_vFit_vertex_noMC->position().x() );
	       VDecayVtxY->push_back( ks_vFit_vertex_noMC->position().y() );
	       VDecayVtxZ->push_back( ks_vFit_vertex_noMC->position().z() );
	       
	       VDecayVtxXE->push_back( sqrt(ks_vFit_vertex_noMC->error().cxx()) );
	       VDecayVtxYE->push_back( sqrt(ks_vFit_vertex_noMC->error().cyy()) );
	       VDecayVtxZE->push_back( sqrt(ks_vFit_vertex_noMC->error().czz()) );
	       VVtxCL->push_back( ChiSquaredProbability((double)(ks_vFit_vertex_noMC->chiSquared()),(double)(ks_vFit_vertex_noMC->degreesOfFreedom())) );

	       VPx->push_back( VCandKP.momentum().x() );
	       VPy->push_back( VCandKP.momentum().y() );
	       VPz->push_back( VCandKP.momentum().z() );	

	       VTrkpPx->push_back(ksPipKP.momentum().x());
	       VTrkpPy->push_back(ksPipKP.momentum().y());
	       VTrkpPz->push_back(ksPipKP.momentum().z());
	       VTrkpMass->push_back(ksPipKP.mass());
	       if ( theDaughterTracks[0]->charge() > 0 ) {
		 VTrkpD0->push_back( theDaughterTracks[0]->d0() );
		 VTrkpD0E->push_back( theDaughterTracks[0]->d0Error() ); 
	       } else {
		 VTrkpD0->push_back( theDaughterTracks[1]->d0() );
		 VTrkpD0E->push_back( theDaughterTracks[1]->d0Error() ); 
	       }

	       VTrkmPx->push_back(ksPimKP.momentum().x());
	       VTrkmPy->push_back(ksPimKP.momentum().y());
	       VTrkmPz->push_back(ksPimKP.momentum().z());
	       VTrkmMass->push_back(ksPimKP.mass());
	       if ( theDaughterTracks[0]->charge() < 0 ) {
		 VTrkmD0->push_back( theDaughterTracks[0]->d0() );
		 VTrkmD0E->push_back( theDaughterTracks[0]->d0Error() ); 
	       } else {
		 VTrkmD0->push_back( theDaughterTracks[1]->d0() );
		 VTrkmD0E->push_back( theDaughterTracks[1]->d0Error() ); 
	       }

               /////////////////////////////////////////////////////////////////////
               // Try getting trigger match from PAT
	       
               mupTrigL1Open1mu->push_back(!patMuonP->triggerObjectMatchesByFilter("hltL1MuOpenL1Filtered0").empty());
               mumTrigL1Open1mu->push_back(!patMuonM->triggerObjectMatchesByFilter("hltL1MuOpenL1Filtered0").empty());
               mupTrigL1Open2mu->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty());
               mumTrigL1Open2mu->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty());
               mupTrigL1Open2muTight->push_back(!patMuonP->triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty());
               mumTrigL1Open2muTight->push_back(!patMuonM->triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty());
	       mupTrig2mu0->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty());
	       mumTrig2mu0->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty());
	       mupTrig2mu0_quark->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty());
	       mumTrig2mu0_quark->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty());
	       mupTrig2mu3->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty());
	       mumTrig2mu3->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered").empty());
	       mupTrig2mu0L2->push_back(!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty());
	       mumTrig2mu0L2->push_back(!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL2PreFiltered0").empty());
               
	       // next do check for mu+trk trigger. set 1 if matches trk, 2 if matches mu and 0 if no matches
               bool matchedMu = false, matchedTrack = false;
               pat::TriggerObjectStandAloneCollection mu0tkMatch = patMuonP->triggerObjectMatchesByFilter("hltMu0TrackJpsiTrackMassFiltered");
               for (unsigned k = 0; k < mu0tkMatch.size(); ++k) {
                 if (mu0tkMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_ ) matchedMu = true;
                 if (mu0tkMatch[k].collection() == string("hltMuTrackJpsiCtfTrackCands::")+hlTriggerResults_ ) matchedTrack = true;
               }
	       if (matchedMu&&matchedTrack) mupTrigmu0trk0->push_back( 3 ); 
               else if (matchedMu) mupTrigmu0trk0->push_back( 2 ); 
	       else if (matchedTrack) mupTrigmu0trk0->push_back( 1 ); 
	       else mupTrigmu0trk0->push_back( 0 );
	       
               matchedMu = false; matchedTrack = false;
	       mu0tkMatch = patMuonM->triggerObjectMatchesByFilter("hltMu0TrackJpsiTrackMassFiltered");
	       for (unsigned k = 0; k < mu0tkMatch.size(); ++k) {
                 if (mu0tkMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_ ) matchedMu = true;
                 if (mu0tkMatch[k].collection() == string("hltMuTrackJpsiCtfTrackCands::")+hlTriggerResults_ ) matchedTrack = true;
               }
	       if (matchedMu&&matchedTrack) mumTrigmu0trk0->push_back( 3 ); 
               else if (matchedMu) mumTrigmu0trk0->push_back( 2 ); 
	       else if (matchedTrack) mumTrigmu0trk0->push_back( 1 ); 
	       else mumTrigmu0trk0->push_back( 0 );
	       	       
               matchedMu = false; matchedTrack = false;
               pat::TriggerObjectStandAloneCollection mu3tkMatch = patMuonP->triggerObjectMatchesByFilter("hltMu3TrackJpsiTrackMassFiltered");
               for (unsigned k = 0; k < mu3tkMatch.size(); ++k) {
                 if (mu3tkMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_ ) matchedMu = true;
                 if (mu3tkMatch[k].collection() == string("hltMuTrackJpsiCtfTrackCands::")+hlTriggerResults_ ) matchedTrack = true;
               }
	       if (matchedMu&&matchedTrack) mupTrigmu3trk0->push_back( 3 );
               else if (matchedMu) mupTrigmu3trk0->push_back( 2 ); 
	       else if (matchedTrack) mupTrigmu3trk0->push_back( 1 ); 
	       else mupTrigmu3trk0->push_back( 0 );
	       
	       matchedMu = false; matchedTrack = false;
	       mu3tkMatch = patMuonM->triggerObjectMatchesByFilter("hltMu3TrackJpsiTrackMassFiltered");
	       for (unsigned k = 0; k < mu3tkMatch.size(); ++k) {
                 if (mu3tkMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_ ) matchedMu = true;
		 if (mu3tkMatch[k].collection() == string("hltMuTrackJpsiCtfTrackCands::")+hlTriggerResults_ ) matchedTrack = true;
               }
	       if (matchedMu&&matchedTrack) mumTrigmu3trk0->push_back( 3 ); 
               else if (matchedMu) mumTrigmu3trk0->push_back( 2 ); 
	       else if (matchedTrack) mumTrigmu3trk0->push_back( 1 ); 
	       else mumTrigmu3trk0->push_back( 0 );	 

	       matchedMu = false; matchedTrack = false;
               pat::TriggerObjectStandAloneCollection  mu0tkmuMatch = patMuonP->triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
               for (unsigned k = 0; k < mu0tkmuMatch.size(); ++k) {
                 if (mu0tkmuMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_) matchedMu = true;
                 if (mu0tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::")+hlTriggerResults_) matchedTrack = true;
               } 
	       if (matchedMu&&matchedTrack) mupTrigmu0trkmu0->push_back( 3 );
               else if (matchedMu) mupTrigmu0trkmu0->push_back( 2 ); 
	       else if (matchedTrack) mupTrigmu0trkmu0->push_back( 1 ); 
	       else mupTrigmu0trkmu0->push_back( 0 );	

	       matchedMu = false; matchedTrack = false;
               mu0tkmuMatch = patMuonM->triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
               for (unsigned k = 0; k < mu0tkmuMatch.size(); ++k) {
                 if (mu0tkmuMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_) matchedMu = true;
                 if (mu0tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::")+hlTriggerResults_) matchedTrack = true;
               } 
	       if (matchedMu&&matchedTrack) mumTrigmu0trkmu0->push_back( 3 ); 
               else if (matchedMu) mumTrigmu0trkmu0->push_back( 2 ); 
	       else if (matchedTrack) mumTrigmu0trkmu0->push_back( 1 ); 
	       else mumTrigmu0trkmu0->push_back( 0 );	
	       
	       matchedMu = false; matchedTrack = false;
               pat::TriggerObjectStandAloneCollection  mu3tkmuMatch = patMuonP->triggerObjectMatchesByFilter("hltMu3TkMuJpsiTkMuMassFiltered");
               for (unsigned k = 0; k < mu3tkmuMatch.size(); ++k) {
                 if (mu3tkmuMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_) matchedMu = true;
                 if (mu3tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::")+hlTriggerResults_) matchedTrack = true;
               } 
	       if (matchedMu&&matchedTrack) mupTrigmu3trkmu0->push_back( 3 ); 
               else if (matchedMu) mupTrigmu3trkmu0->push_back( 2 ); 
	       else if (matchedTrack) mupTrigmu3trkmu0->push_back( 1 ); 
	       else mupTrigmu3trkmu0->push_back( 0 );	

	       matchedMu = false; matchedTrack = false;
               mu3tkmuMatch = patMuonM->triggerObjectMatchesByFilter("hltMu3TkMuJpsiTkMuMassFiltered");
               for (unsigned k = 0; k < mu3tkmuMatch.size(); ++k) {
                 if (mu3tkmuMatch[k].collection() == string("hltL3MuonCandidates::")+hlTriggerResults_) matchedMu = true;
                 if (mu3tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::")+hlTriggerResults_) matchedTrack = true;
               } 
	       if (matchedMu&&matchedTrack) mumTrigmu3trkmu0->push_back( 3 );
               else if (matchedMu) mumTrigmu3trkmu0->push_back( 2 ); 
	       else if (matchedTrack) mumTrigmu3trkmu0->push_back( 1 ); 
	       else mumTrigmu3trkmu0->push_back( 0 );	



               //////////////////////////////////////////////////////////////////////////////////////////
               //calculate ctau and ctau2D with ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)

	       float betagamma = (bCandMC->currentState().globalMomentum().mag()/mb);
               float bPt = sqrt(bCandMC->currentState().globalMomentum().x()*bCandMC->currentState().globalMomentum().x()+
	                        bCandMC->currentState().globalMomentum().y()*bCandMC->currentState().globalMomentum().y() );
               float betagammaT = bPt/mb;

	       //bestVtx is the best primary vertex
	       //bestVtxBS is the best primary with beamspot constraint
	       
               GlobalPoint BVP = GlobalPoint( bDecayVertexMC->position() );
	       GlobalPoint PVP = GlobalPoint( bestVtx.position().x(), bestVtx.position().y(), bestVtx.position().z() );
	       GlobalVector sep3D = BVP-PVP;
			      
	       GlobalVector pBV = bCandMC->currentState().globalMomentum();	      
	       float bctau_temp = (mb*(sep3D.dot(pBV)))/(pBV.dot(pBV));
	       
	       bctau->push_back( bctau_temp );
               
	       // calculate ctau error. Momentum error is negligible compared to the vertex errors, so don't worry about it
	       GlobalError BVE = bDecayVertexMC->error();
	       GlobalError PVE = GlobalError( bestVtx.error() );
	       VertexDistance3D theVertexDistance3D; 
               Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVP, PVE) );
	       double myError = TheMeasurement.error();	       
	       
	       //  ctau is defined by the portion of the flight distance along the compoenent of the B momementum, so only
	       // consider the error of that component, too, which is accomplished by scaling by ((VB-VP)(dot)PB)/|VB-VP|*|PB|	       
               float scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );    	       
               float bctauE_temp =  (myError*scale)/betagamma;
	       bctauE->push_back( bctauE_temp );

	       // do ctau 2D calculation
	       GlobalVector pTBV = GlobalVector( bCandMC->currentState().globalMomentum().x(), bCandMC->currentState().globalMomentum().y(), 0);	       
	       float bctau2D_temp = (mb*(sep3D.dot(pTBV)))/(pTBV.dot(pTBV));
	       bctau2D->push_back(bctau2D_temp);

	       VertexDistanceXY theVertexDistanceXY; 
               Measurement1D TheMeasurementXY = theVertexDistanceXY.distance( VertexState(BVP, BVE), VertexState(PVP, PVE) );
	       double myErrorXY = TheMeasurementXY.error();	       
	       GlobalVector sep2D = GlobalVector( bDecayVertexMC->position().x()-bestVtx.position().x(), bDecayVertexMC->position().y()-bestVtx.position().y(), 0);
               float scaleXY = abs( (sep2D.dot(pTBV))/(sep2D.mag()*pTBV.mag()) );
	       bctau2DE->push_back( (myErrorXY*scaleXY)/betagammaT );
	       
	       // do ctau 3D calc with BS constraint
	       GlobalPoint PVBSP = GlobalPoint( bestVtxBS.position().x(), bestVtxBS.position().y(), bestVtxBS.position().z() );
               GlobalVector sep3DBS = BVP-PVBSP;
	       float bctauBS_temp = (mb*(sep3DBS.dot(pBV)))/(pBV.dot(pBV));
	       bctauBS->push_back(bctauBS_temp);

               GlobalError PVBSE = GlobalError( bestVtxBS.error() );
               TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVBSP, PVBSE) );
               myError = TheMeasurement.error();
               scale = abs( (sep3DBS.dot(pBV))/(sep3DBS.mag()*pBV.mag()) );    	 
	       bctauBSE->push_back( (myError*scale)/betagamma );
	       
	       //calculate most probable decay length
	       AlgebraicMatrix31 pB;
	       pB(0,0) = bCandMC->currentState().globalMomentum().x();
	       pB(1,0) = bCandMC->currentState().globalMomentum().y();
	       pB(2,0) = bCandMC->currentState().globalMomentum().z();

	       AlgebraicMatrix13 pBT;
	       pBT(0,0) = bCandMC->currentState().globalMomentum().x();
	       pBT(0,1) = bCandMC->currentState().globalMomentum().y();
	       pBT(0,2) = bCandMC->currentState().globalMomentum().z();

	       AlgebraicMatrix31 PV;
	       PV(0,0) = bestVtx.position().x();
	       PV(1,0) = bestVtx.position().y();
	       PV(2,0) = bestVtx.position().z();
	       AlgebraicMatrix31 BV;
	       BV(0,0) = bDecayVertexMC->position().x();
	       BV(1,0) = bDecayVertexMC->position().y();
	       BV(2,0) = bDecayVertexMC->position().z();
	       AlgebraicMatrix31 lxyz = BV-PV;
	       AlgebraicMatrix33 PVError(bestVtx.error());
	       AlgebraicMatrix33 BVError(bDecayVertexMC->error().matrix_new());
	       AlgebraicMatrix33 lxyzError = PVError + BVError;
	       lxyzError.Invert();

	       float bctauMPV_temp;
	       AlgebraicMatrix11 a = pBT * lxyzError * pB ;
	       AlgebraicMatrix11 b = pBT * lxyzError * lxyz;
	       double num(b(0,0));
	       double deno(a(0,0));
	       bctauMPV_temp = (num*mb)/(deno);
	       bctauMPV->push_back(bctauMPV_temp);
	       
	       // sigma=mb/sqrt(ptT*S*pt);
	       float bctauMPVE_temp = mb*(sqrt(1/(deno)));
	       
	       bctauMPVE->push_back( bctauMPVE_temp );

	       
	       // calculate MPV with beamspot constrained primary vertex
	       AlgebraicMatrix31 PVBS;
	       PVBS(0,0) = bestVtxBS.position().x();
	       PVBS(0,1) = bestVtxBS.position().y();
	       PVBS(0,2) = bestVtxBS.position().z();
	       AlgebraicMatrix31 lxyzBS = BV-PVBS;
	       AlgebraicMatrix33 PVBSError(bestVtxBS.error());
	       AlgebraicMatrix33 lxyzBSError = PVBSError + BVError;
	       lxyzBSError.Invert();

	       float bctauMPVBS_temp;
	       AlgebraicMatrix11 aBS = pBT * lxyzBSError * pB ;
	       AlgebraicMatrix11 bBS = pBT * lxyzBSError * lxyzBS;
	       double numBS(bBS(0,0));
	       double denoBS(aBS(0,0));
	       bctauMPVBS_temp = (numBS*mb)/(denoBS);
	       bctauMPVBS->push_back(bctauMPVBS_temp);
               bctauMPVBSE->push_back( mb*(sqrt(1/denoBS)) );

               //calculate 2D MPV
	       float bctauMPV2D_temp;
	       AlgebraicMatrix31 pB2D;
	       pB2D(0,0) = bCandMC->currentState().globalMomentum().x();
	       pB2D(1,0) = bCandMC->currentState().globalMomentum().y();
	       pB2D(2,0) = 0;
	       AlgebraicMatrix13 pBT2D;
	       pBT2D(0,0) = bCandMC->currentState().globalMomentum().x();
	       pBT2D(0,1) = bCandMC->currentState().globalMomentum().y();
	       pBT2D(0,2) = 0;

	       AlgebraicMatrix31 PV2D;
	       PV2D(0,0) = bestVtx.position().x();
	       PV2D(1,0) = bestVtx.position().y();
	       PV(2,0) = 0;
	       AlgebraicMatrix31 BV2D;
	       BV2D(0,0) = bDecayVertexMC->position().x();
	       BV2D(1,0) = bDecayVertexMC->position().y();
	       BV2D(2,0) = 0;
	       AlgebraicMatrix31 lxyz2D = BV2D-PV2D;
	       AlgebraicMatrix11 a2D = pBT2D * lxyzError * pB2D;
	       AlgebraicMatrix11 b2D = pBT2D * lxyzError * lxyz2D;
	       double num2D(b2D(0,0));
	       double deno2D(a2D(0,0));
	       bctauMPV2D_temp = (num2D*mb)/(deno2D);
	       bctauMPV2D->push_back(bctauMPV2D_temp);
	       
	       // sigma=mb/sqrt(ptT*S*pt);
	       float bctauMPV2DE_temp = mb*(sqrt(1/(deno2D)));
	       
	       bctauMPV2DE->push_back( bctauMPV2DE_temp );

	       // try refitting the primary without the tracks in the B reco candidate
	       
	       // first get tracks from the original primary
	       vector<reco::TransientTrack> vertexTracks;
	       
	       for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtx.tracks_begin();
		     iTrack != bestVtx.tracks_end(); ++iTrack) {
		 // compare primary tracks to check for matches with B cand
		 TrackRef trackRef = iTrack->castTo<TrackRef>();

		 // the 4 tracks in the B cand are theDaughterTracks[0] theDaughterTracks[1] glbTrackP glbTrackM
		 if (  !( (theDaughterTracks[0]==trackRef) || (theDaughterTracks[1]==trackRef) ||
		          (glbTrackP==trackRef) || (glbTrackM==trackRef) ) ) {
		   TransientTrack tt(trackRef, &(*bFieldHandle) );
		   vertexTracks.push_back(tt);
		 } //else { cout << "found track match with primary" << endl;}
	       }
	       
	       priRfNTrkDif->push_back( bestVtx.tracksSize() - vertexTracks.size() );
	       
	       // if no tracks in primary or no reco track included in primary then don't do anything
	       // if so, then update bctau_temp and bctauMPV_temp

               reco::Vertex bestVtxRf = bestVtx;

	       if (  vertexTracks.size()>0 && (bestVtx.tracksSize()!=vertexTracks.size()) ) {
		 
		 AdaptiveVertexFitter theFitter;
		 TransientVertex v = theFitter.vertex(vertexTracks);
		 if ( v.isValid() ) {
		   //calculate ctau with the new vertex to compare to the old one.
		   //1. standard 3D calculation
	           GlobalPoint PVRfP = GlobalPoint( v.position().x(), v.position().y(), v.position().z() );
	           GlobalVector sep3DRf = BVP-PVRfP;
                   bctau_temp = (mb*(sep3D.dot(pBV)))/(pBV.dot(pBV));

		   //cout << "bctauRf = " << bctau_temp << endl;

		   reco::Vertex recoV = (reco::Vertex)v;
	           GlobalError PVRfE = GlobalError( recoV.error() );
                   TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVRfP, PVRfE) );
	           myError = TheMeasurement.error();	       
                   scale = abs( (sep3DRf.dot(pBV))/(sep3DRf.mag()*pBV.mag()) );    	       
                   bctauE_temp =  (myError*scale)/betagamma;
		   
		   //2. 3D MPV calculation
		   AlgebraicMatrix31 PVnew;
		   PVnew(0,0) = v.position().x();
		   PVnew(0,1) = v.position().y();
		   PVnew(0,2) = v.position().z();
		   AlgebraicMatrix31 lxyznew = BV-PVnew;
		   AlgebraicMatrix33 PVnewError(recoV.error());
		   AlgebraicMatrix33 lxyznewError = PVnewError + BVError;
		   lxyznewError.Invert();
		   
		   AlgebraicMatrix11 anew = pBT * lxyznewError * pB ;
		   AlgebraicMatrix11 bnew = pBT * lxyznewError * lxyznew;
		   double numnew(bnew(0,0));
		   double denonew(anew(0,0));
		   bctauMPV_temp = (numnew*bCandMC->currentState().mass())/(denonew);   
		   //cout << "bctauMPVRf = " << bctauMPV_temp << endl;
		   
		   bctauMPVE_temp = mb*(sqrt(1/(denonew)));
		   
		   //set bestVtxRf as new best vertex to fill variables for ntuple
		   bestVtxRf = reco::Vertex(v);
		   
		 } 
	       } 
	       
	       bctauRf->push_back( bctau_temp );
               bctauRfE->push_back( bctauE_temp );
	       bctauMPVRf->push_back( bctauMPV_temp );
	       bctauMPVRfE->push_back( bctauMPVE_temp );
	       
	       priRfVtxX->push_back( bestVtxRf.x() );
	       priRfVtxY->push_back( bestVtxRf.y() );
	       priRfVtxZ->push_back( bestVtxRf.z() );
	       priRfVtxXE->push_back( bestVtxRf.xError() );
	       priRfVtxYE->push_back( bestVtxRf.yError() );
	       priRfVtxZE->push_back( bestVtxRf.zError() );
	       priRfVtxCL->push_back( ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );
     
	       nB++;
	       
	
	       /////////////////////////////////////////////////
	
               pionParticles.clear();
	       muonParticles.clear();
	       vFitMCParticles.clear();

/*
				 
	       //////////////////////////////
	       //Check PAT truth match here// PAT truth match doesn't work because the V0 tracks need to have the momentum defined at the V0 vertex
	       //////////////////////////////

	       GenParticleRef muPGP;
	       GenParticleRef muMGP;
	       GenParticleRef piPGP;
	       GenParticleRef piMGP;
	       bool psiMatch = false;
	       bool ksMatch = false;
	       bool bMatch = false;
	       
	       //cout << "reco mu+ has eta = " << iMuonP->eta() << " and phi = " << iMuonP->phi() << endl;
	       for ( size_t iGP = 0; iGP<iMuonP->genParticlesSize(); ++iGP) {
		 const GenParticleRef gpRef = iMuonP->genParticleRef(iGP);
		 if (!gpRef.isNull()) {
		   //cout << "genParticle " << iGP << " for mu+ has eta = " << gpRef->eta() << " and phi = " << gpRef->phi() << endl;
		   //cout << "deltaR for mup = " << sqrt(  (iMuonP->eta()-gpRef->eta())*(iMuonP->eta()-gpRef->eta()) + (iMuonP->phi()-gpRef->phi())*(iMuonP->phi()-gpRef->phi()) ) << endl;
		   muPGP = gpRef;
		 } else {
		   cout << "gpRef is null for mu+" << endl;
		 }
	       }
		 
	       //cout << "reco mu- has eta = " << iMuonM->eta() << " and phi = " << iMuonM->phi() << endl;
	       for ( size_t iGP = 0; iGP<iMuonM->genParticlesSize(); ++iGP) {
		 const GenParticleRef gpRef = iMuonM->genParticleRef(iGP);
		 if (!gpRef.isNull()) {
		   //cout << "genParticle " << iGP << " for mu- has eta = " << gpRef->eta() << " and phi = " << gpRef->phi() << endl;
		   //cout << "deltaR for mum = " << sqrt(  (iMuonM->eta()-gpRef->eta())*(iMuonM->eta()-gpRef->eta()) + (iMuonM->phi()-gpRef->phi())*(iMuonM->phi()-gpRef->phi()) ) << endl;
		   muMGP = gpRef;
		 } else {
		   cout << "gpRef is null for mu-" << endl;
		 }
	       }
	       
	       //look at V0 and track match
	       
	       //cout << "patTrack1 genParticleSize = " << patTrack1.genParticlesSize() << endl;
	       //cout << "reco pi+ has eta = " << patTrack1.eta() << " and phi = " << patTrack1.phi() << endl;
	       for ( size_t iGP = 0; iGP<patTrack1.genParticlesSize(); ++iGP) {
		 const GenParticleRef gpRef = patTrack1.genParticleRef(iGP);
		 //cout << "got gpRef" << endl;
		 if (!gpRef.isNull()) {
		   //cout << "genParticle " << iGP << " for pi+ has eta = " << gpRef->eta() << " and phi = " << gpRef->phi() << endl;
		   //cout << "deltaR for pip = " << sqrt(  (patTrack1.eta()-gpRef->eta())*(patTrack1.eta()-gpRef->eta()) + (patTrack1.phi()-gpRef->phi())*(patTrack1.phi()-gpRef->phi()) ) << endl;
		   piPGP = gpRef;
		 } else {
		   cout << "gpRef is null for pi+" << endl;
		 }
	       }

	       //cout << "reco pi- has eta = " << patTrack2.eta() << " and phi = " << patTrack2.phi() << endl;
	       for ( size_t iGP = 0; iGP<patTrack2.genParticlesSize(); ++iGP) {
		 const GenParticleRef gpRef = patTrack2.genParticleRef(iGP);
		 if (!gpRef.isNull()) {
		   //cout << "genParticle " << iGP << " for pi- has eta = " << gpRef->eta() << " and phi = " << gpRef->phi() << endl;
		   //cout << "number of mothers = " << gpRef->numberOfMothers() << endl;
		   //cout << "deltaR for pim = " << sqrt(  (patTrack2.eta()-gpRef->eta())*(patTrack2.eta()-gpRef->eta()) + (patTrack2.phi()-gpRef->phi())*(patTrack2.phi()-gpRef->phi()) ) << endl;
		   piMGP = gpRef;
		 } else {
		   cout << "gpRef is null for pi-" << endl;
		 }
	       }

	       cout << "moving to truth match check from PAT matching" << endl;

	       //have all 4 track genParticles. Check for truth match

	       if ( muMGP.isNonnull() && muPGP.isNonnull() ) {
		 if ( muMGP->numberOfMothers()==1 && muPGP->numberOfMothers()==1 ) {
		   if ( muMGP->mother() == muPGP->mother() ) {
		     //cout << "mu mothers match" << endl;
		     psiMatch = true;
		     //cout << "mu+ mother pdgid = " << muPGP->mother()->pdgId() << endl;
		   } 
		   //else cout << "mu mothers don't match" << endl;
		 } //else cout << "muGPs don't both have one mother: muP mothers = " << muPGP->numberOfMothers() << "  and muM mothers = " << muMGP->numberOfMothers() << endl;
	       } //else cout << "invalid ref in muons" << endl;
	       
	       if ( piMGP.isNonnull() && piPGP.isNonnull() ) {
		 if ( piMGP->numberOfMothers()==1 && piPGP->numberOfMothers()==1 ) {
		   if ( piMGP->mother() == piPGP->mother() ) {
		     cout << "pi mothers match with id = " << piMGP->mother()->pdgId() << "  and parent id = " << piMGP->mother()->mother()->pdgId() << endl;
		     ksMatch = true;
		   }  //else cout << "pi mothers don't match" << endl;
		 }  //else cout << "piGPs don't both have one mother: piP mothers = " << piPGP->numberOfMothers() << "  and piM mothers = " << piMGP->numberOfMothers() << endl;
	       } //else cout << "invalid ref in pions" << endl;
	       
	       if (psiMatch && ksMatch) {
		 //cout << "both ks and psi match, so check for common mother" << endl;
		 //cout << "psi mother is pdgid " << muMGP->mother()->mother()->pdgId() << " and ks mother is " << piMGP->mother()->mother()->pdgId() << " and ks gmother is " << piMGP->mother()->mother()->mother()->pdgId() << endl;
		 if ( piMGP->mother()->mother() == muMGP->mother()->mother() ) {
		   //cout << "ks and psi mothers match" << endl;
		   bMatch = true;
		 }
		 //cout << "psi mother pdgid = " << muMGP->mother()->mother()->pdgId() << endl;
		 //cout << "Ks mother pdgid = " << piMGP->mother()->mother()->pdgId() << endl;
	       }
	       
	       if(bMatch) truthMatchPAT.push_back(1);
	       else truthMatchPAT.push_back(-1);
	       
	       if(ksMatch) truthKsPAT.push_back(1);
	       else truthKsPAT.push_back(-1);
	       
	       if(psiMatch) truthPsiPAT.push_back(1);
	       else truthPsiPAT.push_back(-1);

*/
	       
	     }
	   }
	 }
       } 
     } 
   } // if V0Handle > 0 and muHandle > 1

   //////////////////////////////////////////////////////
   //////// get truth information from genParticles only for events with a B candidate
   if (nB > 0 && doMC_) {

     genKsPsi = -1; genKsPsi2 = -1; genKstarpPsi = -1; genLambdaPsi = -1; prompt = 0; feedup = -1; feeddown = -1;

     if (debug) cout << "Size of genParticle collection is " << genParticles->size() << endl;
     
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       // check if any of our signals were generated
       
       const Candidate & BCand = (*genParticles)[ k ];

       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
	 // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found B0";
	 int ipsi(-1), iks(-1), nGamma(0);
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
	   if (genDau->pdgId()==22) nGamma++;
	 }
	 if ( (BCand.numberOfDaughters()-nGamma) > 2) wrong = true;
	 if (ipsi!=-1&&iks!=-1&&!wrong)
	   cout << " found genKsPsi";
	 
	 cout << endl;
	 
	 if (ipsi!=-1&&iks!=-1&&!wrong) {
	   
	   genKsPsi = 1;
	   //write out info from daughters
	   const Candidate * genpsi =  BCand.daughter(ipsi);
	   const Candidate * genks =  BCand.daughter(iks)->daughter(0);  // must get daughter since Ks comes from K0->Ks
	   
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
              
       // check for JPsiKs where it doesn't go through K0
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
	 // only check for signal decay after possible B0 B0bar oscilation
	 int ipsi(-1), iks(-1), nGamma(0);
	 bool wrong = false;
	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   int imu1(-1), imu2(-1),  ipi1(-1), ipi2(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     if ( genGDau->pdgId()==13 && genDau->pdgId()==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && genDau->pdgId()==443 ) imu2 = j;
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
	     if ( genGDau->pdgId()==211 && genDau->pdgId()==310 ) ipi1 = j;
	     if ( genGDau->pdgId()==-211 && genDau->pdgId()==310 ) ipi2 = j;
	     if ( genDau->pdgId()==310 && abs(genGDau->pdgId())!=211 && genGDau->pdgId()!=22 ) wrong = true; 
           }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=310 && genDau->pdgId()!=22 ) wrong = true;
           if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi1!=-1 && ipi2!=-1&&!wrong) iks = i;
           if (genDau->pdgId()==22) nGamma++;
	 }
	 if ((BCand.numberOfDaughters()-nGamma)>2) wrong = true;
	 if (ipsi!=-1&&iks!=-1&&!wrong)
	   cout << " found genKsPsi2" << endl;
	 
	 if (ipsi!=-1&&iks!=-1&&!wrong)
	   genKsPsi2 = 1;
       }
       
       //check for B->JPsiK*+(kspi) decay   
       if ( abs(BCand.pdgId())==521 ) {
	 cout << "found B+";
	 int ipsi(-1), ikstp(-1), nGamma(0);
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
	   if (genDau->pdgId()==22) nGamma++;
	 }
	 if ((BCand.numberOfDaughters()-nGamma)>2) wrong = true;
	 if ( ipsi!=-1 && ikstp!=-1 &&!wrong) {
	   cout << " found genKstarpPsi";
	   genKstarpPsi =1;
	 }
	 cout << endl;
       } //check for genparticle with id = 521 for a B+
       
       // check for Lambda_b
       if (abs(BCand.pdgId())==5122) {
	 cout << "found Lambda_B";
	 int ipsi(-1), ilam(-1), nGamma(0);
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
	   if (genDau->pdgId()==22) nGamma++;
	 }
         if ( (BCand.numberOfDaughters() - nGamma) > 2) wrong = true;
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
	   //cout << "psi mother has ID = " << psiMom->pdgId() << endl;
	   if ( (abs(psiMom->pdgId())<600 && abs(psiMom->pdgId())>500) || (abs(psiMom->pdgId())<6000 && abs(psiMom->pdgId())>5000) ) {
	     isPrompt = false;
	     continue;
	   } else {
	     for ( uint i = 0; i < psiMom->numberOfMothers(); i++){
	       const Candidate * psiGMom = psiMom->mother(i);
	       //cout << "psi grandmother has ID = " << psiGMom->pdgId() << endl;
	       if ( (abs(psiGMom->pdgId())<600 && abs(psiGMom->pdgId())>500) ||  (abs(psiGMom->pdgId())<6000 && abs(psiGMom->pdgId())>5000) ) {
		 isPrompt = false;
		 continue;
	       } else {
		 for ( uint i = 0; i < psiGMom->numberOfMothers(); i++){
		   const Candidate * psiGGMom = psiGMom->mother(i);
		   //cout << "psi greatgrandmother has ID = " << psiGGMom->pdgId() << endl;
		   if ( (abs(psiGGMom->pdgId())<600 && abs(psiGGMom->pdgId())>500) ||  (abs(psiGGMom->pdgId())<6000 && abs(psiGGMom->pdgId())>5000) ) {
		     isPrompt = false;
		     continue;
		   }
		 }
	       }
	     }
	   }
	 }
	 if (!isPrompt) prompt = -1; else prompt = 1;
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
       if (foundPsi) {
         //fill the generated psi variables here
         genJP = psi->p();
	 genJPt = psi->pt();
       }
       if (foundV0) {
         //fill the generated ks variables here
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
   
   if (debug)  cout << "done with B reco. Now checking all generated particles for signal events" << endl;
   if (doMC_) {

     // first check for two reconstructed muons within acceptance
     bool accMup = false;
     bool accMum = false;
     bool trigMup = false;
     bool trigMum = false;
     if ( !onlyCount_ && thePATMuonHandle->size()>=2 ) {
     
       for ( std::vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
	     iMuon != thePATMuonHandle->end(); ++iMuon ) {
         if ( (iMuon->pt()>3.3&&abs(iMuon->eta())<1.3) ||
              (iMuon->p()>2.9&&abs(iMuon->eta())>1.3&&abs(iMuon->eta())<2.2) ||
              (iMuon->pt()>0.8&&abs(iMuon->eta())>2.2&&abs(iMuon->eta())<2.4 )) {
           if (iMuon->charge() == 1) {
	     accMup = true;
	     if (!iMuon->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty() ||
	         !iMuon->triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty() ) 
	       trigMup = true;
	   }
	   if (iMuon->charge() == -1) {
	     accMum = true;
	     if (!iMuon->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty() ||
	         !iMuon->triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty() ) 
	       trigMum = true;
	   }
         }
       }
     }       
            
     //Check genParticles for every event to get number of generated signal and pT values
     genBPt = -1; genBeta = -99.0; genBy = -99.0; genMupEta = -99.0; genMupPt = -99.0; genMupP = -99.0; genMupPhi = -99.0; genMumEta = -99.0; genMumPt = -99.0; genMumP = -99.0; genMumPhi = -99.0;
     muAcc = 0; muTrig = 0;
     genHLT_2mu0 = -1; genHLT_2mu0_quark = -1; genHLT_2mu3 = -1; genHLT_2mu0L2 = -1; genHLT_mu0trk0 = -1; genHLT_mu3trk0 = -1; 
     genHLT_mu0trkmu0 = -1; genHLT_mu3trkmu0 = -1;genHLT_mu0trkmu0OST = -1; genHLT_mu3trkmu0OST = -1; genHLT_mu0trkmu0OST_tight = -1;
     genHLT_L1muOpen = -1; genHLT_L12muOpen = -1; genHLT_L12muOpenTight = -1; 
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];       
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
         // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found gen1 B0";
         int ipsi(-1), iks(-1), nGamma(0);
	 bool wrong = false;
         for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1), ipi1(-1), ipi2(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     if ( genGDau->pdgId()==13 && abs(genDau->pdgId())==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && abs(genDau->pdgId())==443 ) imu2 = j;	       
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
	     for ( uint m = 0; m < genGDau->numberOfDaughters(); m++){
	       const Candidate * genGGDau = genGDau->daughter(m);
	       cout << " ===" << genGGDau->pdgId();
	       if ( genGGDau->pdgId()==211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi1 = m;
	       if ( genGGDau->pdgId()==-211 && abs(genGDau->pdgId())==310 && abs(genDau->pdgId())==311 ) ipi2 = m;
	       if ( genGDau->pdgId()==310 && abs(genGGDau->pdgId())!=211 && genGGDau->pdgId()!=22 ) wrong = true; 
	     }
	   }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=311 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi1!=-1 && ipi2!=-1&&!wrong) iks = i;
           if (genDau->pdgId() == 22) nGamma++;
         }
	 if ( (BCand.numberOfDaughters() - nGamma) > 2) wrong = true;
         if (ipsi!=-1&&iks!=-1&&!wrong) {
	   cout << " found genKsPsi";
	   genBPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
           genBeta = BCand.eta();
 	   float genBP = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) +  (BCand.pz()*BCand.pz()) ); 
	   float genBE = sqrt( mb*mb + genBP*genBP );
	   genBy = fabs( 0.5*log( (genBE+BCand.pz())/(genBE-BCand.pz()) ));
	   for (uint j=0; j<BCand.daughter(ipsi)->numberOfDaughters(); j++) {
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==13) { //13 is a mu-
	       genMumEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMumPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMumP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMumPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMumPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMumPhi += pi;
	     }
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==-13) { //-13 is a mu+
               genMupEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMupPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMupP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMupPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMupPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMupPhi += pi;
	     }
	   }
	   muAcc = (accMum&&accMup);
	   muTrig = (trigMum&&trigMup);
	   cout << "filling muTrig with " << muTrig << endl;
	   genHLT_2mu0 = hlt_2mu0;
	   genHLT_2mu0_quark = hlt_2mu0_quark;
	   genHLT_2mu3 = hlt_2mu3;
	   genHLT_2mu0L2 = hlt_2mu0L2;
           genHLT_L1muOpen = hlt_L1muOpen;
     	   genHLT_L12muOpen = hlt_L12muOpen;
     	   genHLT_L12muOpenTight = hlt_L12muOpenTight;
	   genHLT_mu0trk0 = hlt_mu0trk0;
	   genHLT_mu3trk0 = hlt_mu3trk0;
	   genHLT_mu0trkmu0 = hlt_mu0trkmu0;
	   genHLT_mu3trkmu0 = hlt_mu3trkmu0;
	   genHLT_mu0trkmu0OST = hlt_mu0trkmu0OST;
	   genHLT_mu3trkmu0OST = hlt_mu3trkmu0OST;
	   genHLT_mu0trkmu0OST_tight = hlt_mu0trkmu0OST_tight;
         }
	 cout << endl;
       }
     }

     if (genBPt>0) {
       treeMC_->Fill();
       cout << "filling tree for MC generated" << endl;
     }
   
     cout << " done with original generated truth check, moving to second generated truth check" << endl;
   
     genBPt = -1; genBeta = -99.0; genBy = -99.0;  genMupEta = -99.0; genMupPt = -99.0; genMupP = -99.0; genMupPhi = -99.0; genMumEta = -99.0; genMumPt = -99.0; genMumP = -99.0; genMumPhi = -99.0;
     muAcc = 0; muTrig = 0;
     genHLT_2mu0 = -1; genHLT_2mu0_quark = -1; genHLT_2mu3 = -1; genHLT_2mu0L2 = -1; genHLT_mu0trk0 = -1; genHLT_mu3trk0 = -1; 
     genHLT_mu0trkmu0 = -1; genHLT_mu3trkmu0 = -1;genHLT_mu0trkmu0OST = -1; genHLT_mu3trkmu0OST = -1; genHLT_mu0trkmu0OST_tight = -1;
     genHLT_L1muOpen = -1; genHLT_L12muOpen = -1; genHLT_L12muOpenTight = -1; 
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
         // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found gen2 B0";
         int ipsi(-1), iks(-1), nGamma(0);
	 bool wrong = false;
         for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1), ipi1(-1), ipi2(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     if ( genGDau->pdgId()==13 && abs(genDau->pdgId())==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && abs(genDau->pdgId())==443 ) imu2 = j;	       
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
	     if ( genGDau->pdgId()==211 && abs(genDau->pdgId())==310 ) ipi1 = j;
	     if ( genGDau->pdgId()==-211 && abs(genDau->pdgId())==310 ) ipi2 = j;	       
	     if ( genDau->pdgId()==310 && abs(genGDau->pdgId())!=211 && genGDau->pdgId()!=22 ) wrong = true; 
	   }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=310 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
	   if (ipi1!=-1 && ipi2!=-1&&!wrong) iks = i;
           if (genDau->pdgId() == 22) nGamma++;
         }
	 if ( (BCand.numberOfDaughters() - nGamma) > 2) wrong = true;
	 cout << " ipsi = " << ipsi << " iks = " << iks << " wrong = " << wrong << endl;
         if (ipsi!=-1&&iks!=-1&&!wrong) {
	   cout << " found genKsPsi";
	   genBPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
           genBeta = BCand.eta();
	   float genBP = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) +  (BCand.pz()*BCand.pz()) ); 
	   float genBE = sqrt( mb*mb + genBP*genBP );
	   genBy = fabs( 0.5*log( (genBE+BCand.pz())/(genBE-BCand.pz()) ));
	   for (uint j=0; j<BCand.daughter(ipsi)->numberOfDaughters(); j++) {
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==13) { //13 is a mu-
	       genMumEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMumPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMumP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMumPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMumPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMumPhi += pi;
	     }
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==-13) { //-13 is a mu+
               genMupEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMupPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMupP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMupPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMupPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMupPhi += pi;
	     }
	   }
	   muAcc = (accMum&&accMup);
	   muTrig = (trigMum&&trigMup);
	   genHLT_2mu0 = hlt_2mu0;
	   genHLT_2mu0_quark = hlt_2mu0_quark;
	   genHLT_2mu3 = hlt_2mu3;
	   genHLT_2mu0L2 = hlt_2mu0L2;
           genHLT_L1muOpen = hlt_L1muOpen;
     	   genHLT_L12muOpen = hlt_L12muOpen;
     	   genHLT_L12muOpenTight = hlt_L12muOpenTight;
	   genHLT_mu0trk0 = hlt_mu0trk0;
	   genHLT_mu3trk0 = hlt_mu3trk0;
	   genHLT_mu0trkmu0 = hlt_mu0trkmu0;
	   genHLT_mu3trkmu0 = hlt_mu3trkmu0;
	   genHLT_mu0trkmu0OST = hlt_mu0trkmu0OST;
	   genHLT_mu3trkmu0OST = hlt_mu3trkmu0OST;
	   genHLT_mu0trkmu0OST_tight = hlt_mu0trkmu0OST_tight;
         }
	 cout << endl;
       }
     }

     if (genBPt>0) {
       treeMC2_->Fill();
       cout << "filling tree for MC2 generated" << endl;
     }

     //finally check for generated JPsi->mumu and K0->Ks without the Ks->pipi (no simulation)
     genBPt = -1; genBeta = -99.0; genBy = -99.0; genMupEta = -99.0; genMupPt = -99.0; genMupP = -99.0; genMupPhi = -99.0; genMumEta = -99.0; genMumPt = -99.0; genMumP = -99.0; genMumPhi = -99.0;
     muAcc = 0; muTrig = 0; weight = 0;
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
         // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found gen3 B0";
         int ipsi(-1), iks(-1), nGamma(0);
	 bool wrong = false;
         for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
	   // check B0 for psi and ks daughters
	   const Candidate * genDau = BCand.daughter(i);
	   cout << " =" << genDau->pdgId();
	   int imu1(-1), imu2(-1);
	   for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
	     const Candidate * genGDau = genDau->daughter(j);
	     cout << " ==" << genGDau->pdgId();
	     if ( genGDau->pdgId()==13 && abs(genDau->pdgId())==443 ) imu1 = j;
	     if ( genGDau->pdgId()==-13 && abs(genDau->pdgId())==443 ) imu2 = j;	       
	     if ( genDau->pdgId()==443 && abs(genGDau->pdgId())!=13 && genGDau->pdgId()!=22 ) wrong = true; 
             if ( abs(genDau->pdgId())==311 && genGDau->pdgId()==310 ) iks = j; 
	   }
	   if ( genDau->pdgId()!=443 && abs(genDau->pdgId())!=311 && genDau->pdgId()!=22 ) 
	     wrong = true;
	   if (imu1!=-1&&imu2!=-1&&!wrong) ipsi = i;
           if (genDau->pdgId() == 22) nGamma++;
         }
	 if ( (BCand.numberOfDaughters() - nGamma) > 2) wrong = true;
         if (ipsi!=-1&&iks!=-1&&!wrong) {
	   cout << " found genKsPsi";
	   genBPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
           genBeta = BCand.eta();
	   float genBP = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) +  (BCand.pz()*BCand.pz()) ); 
	   float genBE = sqrt( mb*mb + genBP*genBP );
	   genBy = fabs( 0.5*log( (genBE+BCand.pz())/(genBE-BCand.pz()) ));
	   for (uint j=0; j<BCand.daughter(ipsi)->numberOfDaughters(); j++) {
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==13) { //13 is a mu-
	       genMumEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMumPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMumP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMumPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMumPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMumPhi += pi;
	     }
	     if (BCand.daughter(ipsi)->daughter(j)->pdgId()==-13) { //-13 is a mu+
               genMupEta = BCand.daughter(ipsi)->daughter(j)->eta();
	       genMupPt = BCand.daughter(ipsi)->daughter(j)->pt();
	       genMupP = BCand.daughter(ipsi)->daughter(j)->p();
	       genMupPhi = atan(BCand.daughter(ipsi)->daughter(j)->py()/BCand.daughter(ipsi)->daughter(j)->px());
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() < 0 ) genMupPhi -= pi;
	       if ( BCand.daughter(ipsi)->daughter(j)->px() < 0 && BCand.daughter(ipsi)->daughter(j)->py() > 0 ) genMupPhi += pi;
	     }
	   }
	   muAcc = (accMum&&accMup);
	   muTrig = (trigMum&&trigMup);
           int weight_sign = 0;
           try {
	     edm::Handle<GenEventInfoProduct> evt_info;
	     iEvent.getByType(evt_info);
             weight_sign = (evt_info->weight() > 0) ? 1 : -1;
           }
	   catch ( ... ) {}
	   weight = weight_sign;
         }
	 cout << endl;
       }
     }
     if (genBPt>0) {
       treeMC3_->Fill();
       cout << "filling tree for MC3 generated" << endl;
     }

     //next check for any generated JPsi->mumu
     genJPsiPt = -1; genJPsieta = -99.0;
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];
       if ( abs(BCand.pdgId())==443 ) {
         cout << "found gen JPsi";
   	 int imu1(-1), imu2(-1), nGamma(0);
         bool wrong = false;
   	 for ( uint i = 0; i < BCand.numberOfDaughters(); i++){
           // check JPsi for pi+ and pi- daughters
           const Candidate * genDau = BCand.daughter(i);
           cout << " =" << genDau->pdgId();
           if ( genDau->pdgId()==13 ) imu1 = i;
           if ( genDau->pdgId()==-13 ) imu2 = i;	     
           if ( abs(genDau->pdgId())!=13 && genDau->pdgId()!=22 ) wrong = true; 
           for ( uint j = 0; j < genDau->numberOfDaughters(); j++){
             const Candidate * genGDau = genDau->daughter(j);
             cout << " ==" << genGDau->pdgId();
           }
   	   if (genDau->pdgId() == 22) nGamma++;
   	 }
         if ( (BCand.numberOfDaughters() - nGamma) > 2) wrong = true;
         if (imu1!=-1&&imu2!=-1&&!wrong) {
           cout << " found genPsiMuMu";
           genJPsiPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
   	   genJPsieta = BCand.eta();
   	 }
         cout << endl;
       }
     }
     if (genJPsiPt>0) {
       treeMC4_->Fill();
       cout << "filling tree for MC4 generated" << endl;
     }     
     
     //lastly check all generated B0
     genBPt = -1; genBeta = -99.0; genBy = -99.0; weight = 0;
     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate & BCand = (*genParticles)[ k ];
       if ( abs(BCand.pdgId())==511 && abs(BCand.daughter(0)->pdgId())!=511 ) {
         // only check for signal decay after possible B0 B0bar oscilation
	 cout << "found gen5 B0" << endl;
	 genBPt = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) );
         genBeta = BCand.eta();
	 float genBP = sqrt( (BCand.px()*BCand.px()) + (BCand.py()*BCand.py()) +  (BCand.pz()*BCand.pz()) ); 
	 float genBE = sqrt( mb*mb + genBP*genBP );
	 genBy = fabs( 0.5*log( (genBE+BCand.pz())/(genBE-BCand.pz()) ));
         int weight_sign = 0;
         try {
	   edm::Handle<GenEventInfoProduct> evt_info;
	   iEvent.getByType(evt_info);
           weight_sign = (evt_info->weight() > 0) ? 1 : -1;
         }
	 catch ( ... ) {}
	 weight = weight_sign;
       }
     }
     if (genBPt>0) {
       treeMC5_->Fill();
       cout << "filling tree for MC5 generated" << endl;
     }
   }
   
   //fill the tree and clear the vectors
   if (nB > 0 ) {
     cout << "filling tree" << endl;
     tree_->Fill();
   }
   l1_mu3 = 0; l1_2mu3 = 0; l1_muOpen = 0; l1_mu0 = 0;
   hlt_mu3 = 0; hlt_mu5 = 0; hlt_mu7 = 0; hlt_mu9 = 0; hlt_2mu0 = 0; hlt_2mu0_quark = 0; hlt_2mu3 = 0; hlt_2mu0L2 = 0; hlt_2mu3JPsi = 0;
   hlt_BJPsiMuMu = 0; hlt_mu0trk0 = 0; hlt_mu3trk0 = 0; hlt_mu0trkmu0 = 0; hlt_mu3trkmu0 = 0; hlt_mu0trkmu0OST = 0; 
   hlt_mu3trkmu0OST = 0; hlt_mu0trkmu0OST_tight = 0; 
   hlt_L1muOpen = 0; hlt_L12muOpen = 0; hlt_L12muOpenTight = 0;
   nB = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   priRfVtxZE->clear(); priRfVtxCL->clear(); priRfNTrkDif->clear();
   bMass->clear(); bVtxCL->clear(); bPx->clear(); bPy->clear(); bPz->clear(); 
   bPxE->clear(); bPyE->clear(); bPzE->clear(); 
   bctau->clear(); bctau2D->clear(); bctauBS->clear(); bctauMPV->clear(); bctauMPV2D->clear(); 
   bctauRf->clear(); bctauMPVRf->clear(); bctauMPVBS->clear(); 
   bctauE->clear(); bctau2DE->clear(); bctauBSE->clear(); bctauMPVE->clear(); bctauMPV2DE->clear();
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
   mumTrigL1Open1mu->clear(); mumTrigL1Open2mu->clear(); mumTrigL1Open2muTight->clear(); mumTrig2mu0->clear(); mumTrig2mu3->clear();
   mumTrig2mu0_quark->clear();
   mumTrig2mu0L2->clear(); mumTrigmu0trk0->clear(); mumTrigmu3trk0->clear(); mumTrigmu0trkmu0->clear(); mumTrigmu3trkmu0->clear(); 
   mupPx->clear(); mupPy->clear(); mupPz->clear(); mupD0->clear(); mupD0E->clear(); mupC2->clear();
   mupCat->clear(); mupME1->clear(); mupAngT->clear(); mupNHits->clear(); mupNPHits->clear();
   mupTrigL1Open1mu->clear(); mupTrigL1Open2mu->clear(); mupTrigL1Open2muTight->clear(); mupTrig2mu0->clear(); mupTrig2mu3->clear();
   mupTrig2mu0_quark->clear();
   mupTrig2mu0L2->clear(); mupTrigmu0trk0->clear(); mupTrigmu3trk0->clear(); mupTrigmu0trkmu0->clear(); mupTrigmu3trkmu0->clear(); 
   VTrkpMass->clear(); VTrkpPx->clear(); VTrkpPy->clear(); VTrkpPz->clear(); 
   VTrkpD0->clear(); VTrkpD0E->clear();
   VTrkmMass->clear(); VTrkmPx->clear(); VTrkmPy->clear(); VTrkmPz->clear();
   VTrkmD0->clear(); VTrkmD0E->clear();
   bResTrkPx->clear(); bResTrkPy->clear(); bResTrkPz->clear();
   bResTrkD0->clear(); bResTrkD0E->clear(); bResTrkChg->clear();
   genKsPsi = 0; genKsPsi2 = 0; genKstarpPsi = 0; genLambdaPsi = 0; feedup = 0; feeddown = 0;
   truePriVtxX = 0; truePriVtxY = 0; truePriVtxZ = 0; trueBPx = 0; trueBPy = 0; trueBPz = 0;
   trueBDecayVtxX = 0; trueBDecayVtxY = 0; trueBDecayVtxZ = 0; trueBResPx = 0; trueBResPy = 0; trueBResPz = 0;
   trueBResDecayVtxX = 0; trueBResDecayVtxY = 0; trueBResDecayVtxZ = 0; 
   trueVPx = 0; trueVPy = 0; trueVPz = 0;
   trueVDecayVtxX = 0; trueVDecayVtxY = 0; trueVDecayVtxZ = 0; trueJPx = 0; trueJPy = 0; trueJPz = 0;
   trueJDecayVtxX = 0; trueJDecayVtxY = 0; trueJDecayVtxZ = 0;
   trueMumPx = 0; trueMumPy = 0; trueMumPz = 0; trueMumD0 = 0;
   trueMupPx = 0; trueMupPy = 0; trueMupPz = 0; trueMupD0 = 0;
   genJPt = 0; genJP = 0;
   trueVTrkpPx = 0; trueVTrkpPy = 0; trueVTrkpPz = 0; trueVTrkpD0 = 0;
   trueVTrkmPx = 0; trueVTrkmPy = 0; trueVTrkmPz = 0; trueVTrkmD0 = 0;
   trueBResTrkPx = 0; trueBResTrkPy = 0; trueBResTrkPz = 0; trueBResTrkD0 = 0; trueBResTrkChg = 0;
   prompt = 0; truthMatch.clear(); truthKs.clear(); truthPsi.clear(); truthMatchPAT.clear(); truthKsPAT.clear(); truthPsiPAT.clear(); 
   
   bMass2->clear(); bMass3->clear(); bMass4->clear(); bMass5->clear(); 
   bPx2->clear(); bPx3->clear(); bPx4->clear(); bPx5->clear(); bDx2->clear(); bDx3->clear(); bDx4->clear(); bDx5->clear(); 
   bPy2->clear(); bPy3->clear(); bPy4->clear(); bPy5->clear(); bDy2->clear(); bDy3->clear(); bDy4->clear(); bDy5->clear(); 
   bPz2->clear(); bPz3->clear(); bPz4->clear(); bPz5->clear(); bDz2->clear(); bDz3->clear(); bDz4->clear(); bDz5->clear(); 
}

void JPsiKsPAT::fillPsi(const reco::Candidate& genpsi) {
  
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

void JPsiKsPAT::fillV0(const reco::Candidate& genv0) {
  
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

int const JPsiKsPAT::getMuCat(reco::Muon const& muon) const{
  int muCat = 0;
  if (muon.isGlobalMuon()) {
    if (muon.isTrackerMuon()) muCat = 1;
    else muCat = 10;
  }
  else if (muon.isTrackerMuon()) {
    if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)) {
      if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)) muCat = 6;
      else if (muon::isGoodMuon(muon, muon::TMOneStationTight)) muCat = 5;
      else if (muon::isGoodMuon(muon, muon::TMOneStationLoose)) muCat = 4;
      else muCat = 3;
    } else muCat = 2;
  }
  else if (muon.isStandAloneMuon()) muCat = 7;
  else if (muon.isCaloMuon()) muCat = 8;
  else muCat = 9;
  
  if ( !(muon::isGoodMuon(muon, muon::TMOneStationLoose)) && muon::isGoodMuon(muon, muon::TMOneStationTight) )
    std::cout << "inconsistent muon cat 1" << std::endl;
  if ( !(muon::isGoodMuon(muon, muon::TMOneStationTight)) && muon::isGoodMuon(muon, muon::TMLastStationAngTight) )
    std::cout << "inconsistent muon cat 2" << std::endl;

  return muCat;
}



bool const JPsiKsPAT::HasGoodME11(reco::Muon const& muon, double const dxdzCut) const{
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
JPsiKsPAT::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","btosmumu ntuple");
  treeMC_ = fs->make<TTree>("mctruth","jpsi(mumu) ks(pipi) gen ntuple");
  treeMC2_ = fs->make<TTree>("mctruth2","jpsi(mumu) k0(ks(pipi)) gen ntuple");
  treeMC3_ = fs->make<TTree>("mctruth3","jpsiks(mumu) k0(ks) gen ntuple");
  treeMC4_ = fs->make<TTree>("mctruth4","jpsi(mumu) gen ntuple");
  treeMC5_ = fs->make<TTree>("mctruth5","B0 gen ntuple");

  treeMC_->Branch("genBPt",&genBPt,"genBPt/f");
  treeMC_->Branch("genBeta",&genBeta,"genBeta/f");
  treeMC_->Branch("genBy",&genBy,"genBy/f");
  treeMC_->Branch("genMupEta",&genMupEta,"genMupEta/f");
  treeMC_->Branch("genMupPt",&genMupPt,"genMupPt/f");
  treeMC_->Branch("genMupP",&genMupP,"genMupP/f");
  treeMC_->Branch("genMupPhi",&genMupPhi,"genMupPhi/f");
  treeMC_->Branch("genMumEta",&genMumEta,"genMumEta/f");
  treeMC_->Branch("genMumPt",&genMumPt,"genMumPt/f");
  treeMC_->Branch("genMumP",&genMumP,"genMumP/f");
  treeMC_->Branch("genMumPhi",&genMumPhi,"genMumPhi/f");  
  treeMC_->Branch("muAcc",&muAcc,"muAcc/I");
  treeMC_->Branch("muTrig",&muTrig,"muTrig/I");
  treeMC_->Branch("genHLT_2mu0",&genHLT_2mu0,"genHLT_2mu0/I");
  treeMC_->Branch("genHLT_2mu0_quark",&genHLT_2mu0_quark,"genHLT_2mu0_quark/I");
  treeMC_->Branch("genHLT_2mu3",&genHLT_2mu3,"genHLT_2mu3/I");
  treeMC_->Branch("genHLT_2mu0L2",&genHLT_2mu0L2,"genHLT_2mu0L2/I");  
  treeMC_->Branch("genHLT_mu0trk0",&genHLT_mu0trk0,"genHLT_mu0trk0/I");
  treeMC_->Branch("genHLT_mu3trk0",&genHLT_mu3trk0,"genHLT_mu3trk0/I");
  treeMC_->Branch("genHLT_mu0trkmu0",&genHLT_mu0trkmu0,"genHLT_mu0trkmu0/I");
  treeMC_->Branch("genHLT_mu3trkmu0",&genHLT_mu3trkmu0,"genHLT_mu3trkmu0/I");  
  treeMC_->Branch("genHLT_mu0trkmu0OST",&genHLT_mu0trkmu0OST,"genHLT_mu0trkmu0OST/I");
  treeMC_->Branch("genHLT_mu3trkmu0OST",&genHLT_mu3trkmu0OST,"genHLT_mu3trkmu0OST/I");    
  treeMC_->Branch("genHLT_mu0trkmu0OST_tight",&genHLT_mu0trkmu0OST_tight,"genHLT_mu0trkmu0OST_tight/I");
  treeMC_->Branch("genHLT_L12muOpen",&genHLT_L12muOpen,"genHLT_L12muOpen/I");
  treeMC_->Branch("genHLT_L12muOpenTight",&genHLT_L12muOpenTight,"genHLT_L12muOpenTight/I");
  treeMC_->Branch("genHLT_L1muOpen",&genHLT_L1muOpen,"genHLT_L1muOpen/I");
  
  treeMC2_->Branch("genBPt",&genBPt,"genBPt/f");
  treeMC2_->Branch("genBeta",&genBeta,"genBeta/f");
  treeMC2_->Branch("genBy",&genBy,"genBy/f");  
  treeMC2_->Branch("genMupEta",&genMupEta,"genMupEta/f");
  treeMC2_->Branch("genMupPt",&genMupPt,"genMupPt/f");
  treeMC2_->Branch("genMupP",&genMupP,"genMupP/f");
  treeMC2_->Branch("genMupPhi",&genMupPhi,"genMupPhi/f");
  treeMC2_->Branch("genMumEta",&genMumEta,"genMumEta/f");
  treeMC2_->Branch("genMumPt",&genMumPt,"genMumPt/f");
  treeMC2_->Branch("genMumP",&genMumP,"genMumP/f");
  treeMC2_->Branch("genMumPhi",&genMumPhi,"genMumPhi/f");
  treeMC2_->Branch("muAcc",&muAcc,"muAcc/I");
  treeMC2_->Branch("muTrig",&muTrig,"muTrig/I");
  treeMC2_->Branch("genHLT_2mu0",&genHLT_2mu0,"genHLT_2mu0/I");
  treeMC2_->Branch("genHLT_2mu0_quark",&genHLT_2mu0_quark,"genHLT_2mu0_quark/I");
  treeMC2_->Branch("genHLT_2mu3",&genHLT_2mu3,"genHLT_2mu3/I");
  treeMC2_->Branch("genHLT_2mu0L2",&genHLT_2mu0L2,"genHLT_2mu0L2/I");  
  treeMC2_->Branch("genHLT_mu0trk0",&genHLT_mu0trk0,"genHLT_mu0trk0/I");
  treeMC2_->Branch("genHLT_mu3trk0",&genHLT_mu3trk0,"genHLT_mu3trk0/I");
  treeMC2_->Branch("genHLT_mu0trkmu0",&genHLT_mu0trkmu0,"genHLT_mu0trkmu0/I");
  treeMC2_->Branch("genHLT_mu3trkmu0",&genHLT_mu3trkmu0,"genHLT_mu3trkmu0/I"); 
  treeMC2_->Branch("genHLT_mu0trkmu0OST",&genHLT_mu0trkmu0OST,"genHLT_mu0trkmu0OST/I");
  treeMC2_->Branch("genHLT_mu3trkmu0OST",&genHLT_mu3trkmu0OST,"genHLT_mu3trkmu0OST/I");    
  treeMC2_->Branch("genHLT_mu0trkmu0OST_tight",&genHLT_mu0trkmu0OST_tight,"genHLT_mu0trkmu0OST_tight/I");
  treeMC2_->Branch("genHLT_L12muOpen",&genHLT_L12muOpen,"genHLT_L12muOpen/I");
  treeMC2_->Branch("genHLT_L12muOpenTight",&genHLT_L12muOpenTight,"genHLT_L12muOpenTight/I");
  treeMC2_->Branch("genHLT_L1muOpen",&genHLT_L1muOpen,"genHLT_L1muOpen/I");
  
  treeMC3_->Branch("genBPt",&genBPt,"genBPt/f");
  treeMC3_->Branch("genBeta",&genBeta,"genBeta/f");
  treeMC3_->Branch("genBy",&genBy,"genBy/f");  
  treeMC3_->Branch("genMupEta",&genMupEta,"genMupEta/f");
  treeMC3_->Branch("genMupPt",&genMupPt,"genMupPt/f");
  treeMC3_->Branch("genMupP",&genMupP,"genMupP/f");
  treeMC3_->Branch("genMupPhi",&genMupPhi,"genMupPhi/f");
  treeMC3_->Branch("genMumEta",&genMumEta,"genMumEta/f");
  treeMC3_->Branch("genMumPt",&genMumPt,"genMumPt/f");
  treeMC3_->Branch("genMumP",&genMumP,"genMumP/f");
  treeMC3_->Branch("genMumPhi",&genMumPhi,"genMumPhi/f");
  treeMC3_->Branch("muAcc",&muAcc,"muAcc/I");
  treeMC3_->Branch("muTrig",&muTrig,"muTrig/I");
  treeMC3_->Branch("weight",&weight,"weight/I");

  treeMC4_->Branch("genJPsiPt",&genJPsiPt,"genJPsiPt/f");
  treeMC4_->Branch("genJPsieta",&genJPsieta,"genJPsieta/f");
  
  treeMC5_->Branch("genBPt",&genBPt,"genBPt/f");
  treeMC5_->Branch("genBeta",&genBeta,"genBeta/f");
  treeMC5_->Branch("genBy",&genBy,"genBy/f");
  treeMC5_->Branch("weight",&weight,"weight/I");
  
  tree_->Branch("l1_mu3",&l1_mu3,"l1_mu3/i");
  tree_->Branch("l1_2mu3",&l1_2mu3,"l1_2mu3/i");
  tree_->Branch("l1_muOpen",&l1_muOpen,"l1_muOpen/i");
  tree_->Branch("l1_mu0",&l1_mu0,"l1_mu0/i");
  tree_->Branch("hlt_mu3",&hlt_mu3,"hlt_mu3/i");
  tree_->Branch("hlt_mu5",&hlt_mu5,"hlt_mu5/i");
  tree_->Branch("hlt_mu7",&hlt_mu7,"hlt_mu7/i");
  tree_->Branch("hlt_mu9",&hlt_mu9,"hlt_mu9/i");
  tree_->Branch("hlt_2mu0",&hlt_2mu0,"hlt_2mu0/i");
  tree_->Branch("hlt_2mu0_quark",&hlt_2mu0_quark,"hlt_2mu0_quark/i");
  tree_->Branch("hlt_2mu3",&hlt_2mu3,"hlt_2mu3/i");
  tree_->Branch("hlt_2mu0L2",&hlt_2mu0L2,"hlt_2mu0L2/i");
  tree_->Branch("hlt_2mu3JPsi",&hlt_2mu3JPsi,"hlt_2mu3JPsi/i");
  tree_->Branch("hlt_BJPsiMuMu",&hlt_BJPsiMuMu,"hlt_BJPsiMuMu/i");
  tree_->Branch("hlt_mu0trk0",&hlt_mu0trk0,"hlt_mu0trk0/i");
  tree_->Branch("hlt_mu3trk0",&hlt_mu3trk0,"hlt_mu3trk0/i");  
  tree_->Branch("hlt_mu0trkmu0",&hlt_mu0trkmu0,"hlt_mu0trkmu0/i");
  tree_->Branch("hlt_mu3trkmu0",&hlt_mu3trkmu0,"hlt_mu3trkmu0/i");    
  tree_->Branch("hlt_mu0trkmu0OST",&hlt_mu0trkmu0OST,"hlt_mu0trkmu0OST/i");
  tree_->Branch("hlt_mu3trkmu0OST",&hlt_mu3trkmu0OST,"hlt_mu3trkmu0OST/i");
  tree_->Branch("hlt_mu0trkmu0OST_tight",&hlt_mu0trkmu0OST_tight,"hlt_mu0trkmu0OST_tight/i");
  tree_->Branch("hlt_L1muOpen",&hlt_L1muOpen,"hlt_L1muOpen/i");
  tree_->Branch("hlt_L12muOpen",&hlt_L12muOpen,"hlt_L12muOpen/i");
  tree_->Branch("hlt_L12muOpenTight",&hlt_L12muOpenTight,"hlt_L12muOpenTight/i");
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
  tree_->Branch("bctauMPV2D",&bctauMPV2D);
  tree_->Branch("bctauMPVBS",&bctauMPVBS);
  tree_->Branch("bctauRf",&bctauRf);
  tree_->Branch("bctauMPVRf",&bctauMPVRf);
  tree_->Branch("bctauE",&bctauE);
  tree_->Branch("bctau2DE",&bctau2DE);
  tree_->Branch("bctauBSE",&bctauBSE);
  tree_->Branch("bctauMPVE",&bctauMPVE);
  tree_->Branch("bctauMPV2DE",&bctauMPV2DE);
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
  tree_->Branch("mumTrigL1Open1mu",&mumTrigL1Open1mu);
  tree_->Branch("mumTrigL1Open2mu",&mumTrigL1Open2mu);
  tree_->Branch("mumTrigL1Open2muTight",&mumTrigL1Open2muTight);
  tree_->Branch("mumTrig2mu0",&mumTrig2mu0);  
  tree_->Branch("mumTrig2mu3",&mumTrig2mu3);
  tree_->Branch("mumTrig2mu0_quark",&mumTrig2mu0_quark);  
  tree_->Branch("mumTrig2mu0L2",&mumTrig2mu0L2);
  tree_->Branch("mumTrigmu0trk0",&mumTrigmu0trk0);  
  tree_->Branch("mumTrigmu3trk0",&mumTrigmu3trk0);
  tree_->Branch("mumTrigmu0trkmu0",&mumTrigmu0trkmu0);  
  tree_->Branch("mumTrigmu3trkmu0",&mumTrigmu3trkmu0);    
  tree_->Branch("mupTrigL1Open1mu",&mupTrigL1Open1mu);
  tree_->Branch("mupTrigL1Open2mu",&mupTrigL1Open2mu);
  tree_->Branch("mupTrigL1Open2muTight",&mupTrigL1Open2muTight);
  tree_->Branch("mupTrig2mu0",&mupTrig2mu0);  
  tree_->Branch("mupTrig2mu3",&mupTrig2mu3);
  tree_->Branch("mupTrig2mu0_quark",&mupTrig2mu0_quark);  
  tree_->Branch("mupTrig2mu0L2",&mupTrig2mu0L2);  
  tree_->Branch("mupTrigmu0trk0",&mupTrigmu0trk0);  
  tree_->Branch("mupTrigmu3trk0",&mupTrigmu3trk0);
  tree_->Branch("mupTrigmu0trkmu0",&mupTrigmu0trkmu0);  
  tree_->Branch("mupTrigmu3trkmu0",&mupTrigmu3trkmu0);  
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
  tree_->Branch("genKsPsi2", &genKsPsi2, "genKsPsi2/I"); 
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
  tree_->Branch("genJPt",&genJPt, "genJPt/f"); 
  tree_->Branch("genJP",&genJP, "genJP/f");   
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
JPsiKsPAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
  treeMC_->GetDirectory()->cd();
  treeMC_->Write();  
  treeMC2_->GetDirectory()->cd();
  treeMC2_->Write();  
  treeMC3_->GetDirectory()->cd();
  treeMC3_->Write();  
  treeMC4_->GetDirectory()->cd();
  treeMC4_->Write();  
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKsPAT);

