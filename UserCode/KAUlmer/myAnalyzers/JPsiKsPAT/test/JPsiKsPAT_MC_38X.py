import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(
	#'/store/mc/Spring10/B0ToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START3X_V26_S09-v1/0004/FE4E52C8-8C4D-DF11-8731-0026189438FA.root',
	'file:/bestman/storage/cms/store/mc/Summer10/B0ToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START36_V9_S09-v1/0001/04467E1D-2181-DF11-813C-0018F3D0967A.root',
	'file:/bestman/storage/cms/store/mc/Summer10/B0ToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START36_V9_S09-v1/0001/047FE937-2181-DF11-97AC-001A92971BDA.root',
	'file:/bestman/storage/cms/store/mc/Summer10/B0ToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/START36_V9_S09-v1/0001/04B709FE-2081-DF11-853A-003048678D52.root',
	))

#from myAnalyzers.JPsiKsPAT.RecoInput2_cfi import *

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START36_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")


#######################
# Muon momentum scale #
#######################
process.poolDBESSource = cms.ESSource("PoolDBESSource",
                                      BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
                                      DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(2),
#    authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
    ),
                                      timetype = cms.untracked.string('runnumber'),
                                      connect = cms.string('frontier://FrontierProd/CMS_COND_31X_PHYSICSTOOLS'),
#                                      connect = cms.string('oracle://cms_orcoff_prod/CMS_COND_31X_PHYSICSTOOLS'),
                                      toGet = cms.VPSet(cms.PSet(
    record = cms.string('MuScleFitDBobjectRcd'),
    tag = cms.string('MuScleFit_Scale_JPsi_1_3_invNb_innerTrack'))))

process.MuScleFitMuonProducer = cms.EDProducer(
    'MuScleFitMuonProducer',
    MuonLabel = cms.InputTag("muons"),
    PatMuons = cms.bool(False))

# Put either the original collection: "muons" or the one coming from the MuScleFitMuonProducer: "MuScleFitMuonProducer"
muonTypeForPAT = "MuScleFitMuonProducer"


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.src       = cms.InputTag(muonTypeForPAT)
process.muonMatch.matched   = cms.InputTag("genParticlesPlusSim")
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,		      # patAODTrackCands
	label='TrackCands',		      # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
	tracks=cms.InputTag('generalTracks'), # input track collection
	particleType='pi+',		      # particle type (for assigning a mass)
	preselection='pt > 0.1',	      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
	selection='pt > 0.1',		      # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
	isolation={},			      # Isolations to use ('source':deltaR; set to {} for None)
	isoDeposits=[],
	mcAs='muon'			      # Replicate MC match as the one used for Muons
	);				      # you can specify more than one collection for this

process.patTrackCandsMCMatch.mcPdgId = cms.vint32(211)
process.patTrackCandsMCMatch.mcStatus = cms.vint32(1, 8)
process.patTrackCandsMCMatch.matched = cms.InputTag("genParticlesPlusSim")
process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.02)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
	src	      = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
	setStatus     = cms.int32(8),
#	 particleTypes = cms.vstring(""),
	filter        = cms.vstring("pt > 0.0"),  # just for testing
	genParticles   = cms.InputTag("genParticles"))

# Redo V0Producer with looser cuts to try to increase efficiency
process.localV0Candidates = cms.EDProducer("V0Producer",
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    useSmoothing = cms.bool(True),
    storeSmoothedTracksInRecoVertex = cms.bool(False),
    doPostFitCuts = cms.bool(True),
    doTrackQualityCuts = cms.bool(True),
    # The next parameters are cut values
    # Track quality cuts
    trackQualities = cms.vstring('loose'),
    #   Normalized track Chi2:
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track:
    tkNhitsCut = cms.int32(6),

    # Vertex cuts
    vtxChi2Cut = cms.double(7.0),
    collinearityCut = cms.double(0.02),
    #  Setting this one to zero; significance cut is sufficient
    rVtxCut = cms.double(0.0),
#    vtxSignificanceCut = cms.double(15.0),
    vtxSignificance2DCut = cms.double(5.0),
    vtxSignificance3DCut = cms.double(0.0),
    kShortMassCut = cms.double(0.07),
    lambdaMassCut = cms.double(0.05),
    kShortNormalizedMassCut = cms.double(0.0),
    lambdaNormalizedMassCut = cms.double(0.0),
    impactParameterSigCut = cms.double(0.5),
    mPiPiCut = cms.double(1.),
    tkDCACut = cms.double(1.),

    # These parameters decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    vertexFitter = cms.InputTag('KalmanVertexFitter'),
    innerHitPosCut = cms.double(4.))

# B->JPsiKs reconstruction done in module with coded fit
process.mkcands = cms.EDAnalyzer("JPsiKsPAT",
  VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
#  HLTriggerResults = cms.untracked.InputTag('TriggerResults::HLT'),
#  HLTriggerResults = cms.untracked.InputTag('TriggerResults::REDIGI36X'),
  HLTriggerResults = cms.untracked.string('REDIGI36X'),
  GenParticles = cms.untracked.string('genParticlesPlusSim'),
  V0Collection = cms.untracked.string('localV0Candidates'),
  MuonType = cms.untracked.string("cleanPatMuons"),
  MuonTypeForPAT = cms.untracked.string(muonTypeForPAT),
  onlyCount = cms.untracked.bool(False),
  doMC = cms.untracked.bool(True))


###################################################################
###################################################################
# New (easier) Onia2mumu trigger matching
# Make PAT Muons

process.patMuons.muonSource = cms.InputTag(muonTypeForPAT)
process.patMuons.embedCaloMETMuonCorrs = False
process.patMuons.embedTcMETMuonCorrs   = False

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()

from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
# with some customization
addMCinfo(process)
# since we match inner tracks, keep the matching tight and make it one-to-one
process.muonMatch.maxDeltaR = 0.05
process.muonMatch.resolveByMatchQuality = True
#process.muonMatch.matched = "genMuons"
#changeRecoMuonInput(process, "mergedMuons")
useL1MatchingWindowForSinglets(process)
changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0


#################################################################
#################################################################

### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""))

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('breco_JPsiKsPAT_MC_ntpl.root'))

process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

# can I do a replace of patMuons with the sequence that includes the trigger matching?
process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)

process.pat = cms.Path( process.patDefaultSequence )

#print(process.pat)

#patAODTrackCandsUnfiltered*patAODTrackCands*electronMatch*(patTrackCandsMCMatch*patTrackCands+patElectrons)+
#muonMatch*patMuons+pfPileUp+pfNoPileUp*(pfAllNeutralHadrons+pfAllChargedHadrons+pfAllPhotons)*tauIsoDepositPFCandidates*
#tauIsoDepositPFChargedHadrons*tauIsoDepositPFNeutralHadrons*tauIsoDepositPFGammas*tauMatch*tauGenJets*
#tauGenJetsSelectorAllHadrons*tauGenJetMatch*patTaus+photonMatch*patPhotons+patCandidateSummary*selectedPatElectrons+
#selectedPatTrackCands+selectedPatMuons+selectedPatTaus+selectedPatPhotons+selectedPatCandidateSummary*
#cleanPatMuons*(cleanPatElectrons+cleanPatTrackCands)*cleanPatPhotons*cleanPatTaus*cleanPatCandidateSummary*
#countPatElectrons+countPatMuons+countPatTaus+countPatLeptons+countPatPhotons

process.momScale = cms.Path( process.MuScleFitMuonProducer )
process.ntup = cms.Path(process.localV0Candidates * process.mkcands )
process.filter = cms.Path(process.hltLevel1GTSeed * process.noScraping)
process.gp = cms.Path( process.genParticlesPlusSim )
process.schedule = cms.Schedule(process.filter, process.gp, process.momScale, process.pat, process.ntup )
