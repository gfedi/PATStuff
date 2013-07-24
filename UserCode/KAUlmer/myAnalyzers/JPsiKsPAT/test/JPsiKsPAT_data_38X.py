import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(
	#'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/08FA1386-5762-DF11-9398-003048D4770E.root',  
	#'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/1A20B9F5-5A62-DF11-959A-003048D460FC.root',  
	#'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/6E6495CD-5962-DF11-B0E4-003048D476B8.root',  
        #'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/0CE8183F-5C62-DF11-85FF-0025B3E05C7E.root',  
        #'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/1AE9F905-5562-DF11-8C18-003048635E32.root',  
        #'/store/data/Commissioning10/MinimumBias/USER/v9/000/135/528/BACBC908-5562-DF11-B94C-003048635FA8.root'
        '/store/data/Run2010A/MuOnia/RAW-RECO/v6/000/142/191/8E4CF75C-05A4-DF11-85D3-003048D293B4.root'
	)
)

#from myAnalyzers.JPsiKsPAT.RecoInput2_cfi import *

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_36X_V12A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,					    #	      patAODTrackCands
	label='TrackCands',		      # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
	tracks=cms.InputTag('generalTracks'), # input track collection
	particleType='pi+',		      # particle type (for assigning a mass)
	preselection='pt > 0.1',	      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
	selection='pt > 0.1',		      # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
	isolation={},			      # Isolations to use ('source':deltaR; set to {} for None)
	isoDeposits=[],
	mcAs = None,			      # Replicate MC match as the one used for Muons
	);				      #  you can specify more than one collection for this

removeMCMatching(process, ['All'])

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
    innerHitPosCut = cms.double(4.)


)

# B->JPsiKs reconstruction done in module with coded fit
process.mkcands = cms.EDAnalyzer("JPsiKsPAT",
  VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
  HLTriggerResults = cms.untracked.string('HLT'),
#  HLTriggerResults = cms.untracked.string('REDIGI36X'),
  GenParticles = cms.untracked.string('genParticlesPlusSim'),
  V0Collection = cms.untracked.string('localV0Candidates'),
  onlyCount = cms.untracked.bool(False),
  doMC = cms.untracked.bool(False)
)







###################################################################
###################################################################
# New (easier) Onia2mumu trigger matching
#
#    # Make PAT Muons

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()


from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

#
#################################################################
#################################################################


### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""), 
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)



process.TFileService = cms.Service("TFileService",
    fileName = cms.string('breco_JPsiKsPAT_data_ntpl.root')
)

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


process.ntup = cms.Path(process.localV0Candidates * process.mkcands )
process.filter = cms.Path(process.hltLevel1GTSeed * process.noScraping)
process.schedule = cms.Schedule(process.filter, process.pat, process.ntup )
