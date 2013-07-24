runMC = True

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

if (runMC):
   process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(
   # one file from: /RelValQCD_Pt_80_120/CMSSW_4_1_5-START311_V2-v1/GEN-SIM-RECO
   #	'file:EC43622A-746F-E011-8800-0018F3D09664.root'
   # one file from /RelValJpsiMM/CMSSW_4_2_2-START42_V11-v1/GEN-SIM-RECO
        'file:/nfs/data36/cms/ulmerk/CMSSW_4_1_5/src/109D8F99-B16D-E011-A69A-0018F3D0968C.root'
   #one file from Brian's B0->KsMuMu production
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_990_1_bD4.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_991_1_jJO.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_992_1_IhK.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_993_1_iU4.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_994_1_Oe4.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_995_1_3YU.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_996_1_cci.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_997_1_1WM.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_998_1_YGm.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_999_1_XTB.root',
   #'file:/nfs/data36/cms/drell/reco_BKsMuMu_new/reco_PYTHIA6_B0toKsmumu_7TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1000_1_CUD.root'
	)
   )
   process.GlobalTag.globaltag = cms.string('START311_V2::All')
   #process.GlobalTag.globaltag = cms.string('START42_V11::All')

else:
   process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(
	'/store/data/Run2011A/SingleMu/RECO/PromptReco-v1/000/161/008/0E573388-C555-E011-B956-001D09F24FBA.root'
	)
   )
   process.GlobalTag.globaltag = cms.string('GR_P_V17::All')

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.matched = cms.InputTag("genParticlesPlusSim")
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

## add track candidates
from PhysicsTools.PatAlgos.tools.trackTools import *

if (runMC):
   makeTrackCandidates(process,
       label        = 'TrackCands',                  
       tracks       = cms.InputTag('generalTracks'), 
       particleType = 'pi+',                         
       preselection = 'pt > 0.1',                     
       selection    = 'pt > 0.1',                     
       isolation    = {},                            
       isoDeposits  = [],                            
       mcAs         = 'muon'           
   )                      
   process.patTrackCandsMCMatch.mcPdgId = cms.vint32(211)
   process.patTrackCandsMCMatch.mcStatus = cms.vint32(1, 8)
   process.patTrackCandsMCMatch.matched = cms.InputTag("genParticlesPlusSim")
   process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
   process.patTrackCandsMCMatch.maxDeltaR = cms.double(0.02)   


else:
    makeTrackCandidates(process,
       label        = 'TrackCands',                  
       tracks       = cms.InputTag('generalTracks'), 
       particleType = 'pi+',                         
       preselection = 'pt > 0.1',                     
       selection    = 'pt > 0.1',                     
       isolation    = {},                            
       isoDeposits  = [],                            
       mcAs         = None          
   )    
             


if not(runMC):
   removeMCMatching(process, ['All'])

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
	src	      = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
	setStatus     = cms.int32(8),
	filter        = cms.vstring("pt > 0.0"),  
	genParticles   = cms.InputTag("genParticles")
)

# Redo V0Producer with looser cuts to increase efficiency
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
if (runMC):
   process.mkcands = cms.EDAnalyzer("JPsiKsPAT",
     VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
     HLTriggerResults = cms.untracked.string('HLT'),
     #  HLTriggerResults = cms.untracked.string('REDIGI36X'),
     GenParticles = cms.untracked.string('genParticlesPlusSim'),
     V0Collection = cms.untracked.string('localV0Candidates'),
     onlyCount = cms.untracked.bool(False),
     doMC = cms.untracked.bool(True)
   )

else:
   process.mkcands = cms.EDAnalyzer("JPsiKsPAT",
     VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
     HLTriggerResults = cms.untracked.string('HLT'),
     #  HLTriggerResults = cms.untracked.string('REDIGI36X'),
     GenParticles = cms.untracked.string('genParticlesPlusSim'),
     V0Collection = cms.untracked.string('localV0Candidates'),
     onlyCount = cms.untracked.bool(False),
     doMC = cms.untracked.bool(False)
   )


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('breco_KsMuMuPAT_MC_ntpl_test.root')
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


# do trigger matching between PAT muons and HLT muons
process.cleanMuonTriggerMatchHLT = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"		  # match by DeltaR only, best match by DeltaR
, src	  = cms.InputTag( "cleanPatMuons" )
, matched = cms.InputTag( "patTrigger" )	  # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
, matchedCuts = cms.string( 'type( "TriggerMuon" )' )
, maxDeltaR = cms.double( 0.10 )
, resolveAmbiguities	= cms.bool( False )	  # only one match per trigger object
, resolveByMatchQuality = cms.bool( False )	  # take best match found per reco object: by DeltaR here (s. above)
)

# Embedding in muons
process.cleanPatMuonsTriggerMatch = cms.EDProducer(
  "PATTriggerMatchMuonEmbedder"
, src     = cms.InputTag( "cleanPatMuons" )
, matches = cms.VInputTag(
   'cleanMuonTriggerMatchHLT'
  )
)
process.trigMatch = cms.Sequence(process.cleanMuonTriggerMatchHLT *process.cleanPatMuonsTriggerMatch )


## let it run
process.pat = cms.Path(
    process.patDefaultSequence
)

# switch on PAT trigger info
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )

#print(process.pat)


#process.load( "HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi" )
#process.hltEventAnalyzerAOD.triggerName = cms.string( '@' ) # <-- wild-card for all HLT paths
#process.p = cms.Path( process.hltEventAnalyzerAOD )


# Schedule
process.ntup = cms.Path(process.localV0Candidates * process.mkcands )
process.gp = cms.Path( process.genParticlesPlusSim )
process.trigMatchPath = cms.Path( process.trigMatch )
if (runMC):
   process.schedule = cms.Schedule(process.gp, process.pat, process.trigMatchPath, process.ntup )
else:
   process.schedule = cms.Schedule(process.pat, process.trigMatchPath, process.ntup )
