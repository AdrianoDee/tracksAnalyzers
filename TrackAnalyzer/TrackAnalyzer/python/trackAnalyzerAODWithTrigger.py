import FWCore.ParameterSet.Config as cms
process = cms.Process('2mu2k')

from FWCore.ParameterSet.VarParsing import VarParsing

import sys
sys.path.append("mclists/")

par = VarParsing ('analysis')

samplefile = "/store/data/Run2018D/Charmonium/AOD/12Nov2019_UL2018-v1/280003/F6E94749-EEF3-B349-8D56-2EC1E585B9EC.root"

par.register ('f',
                                  samplefile,
                                  VarParsing.multiplicity.list,
                                  VarParsing.varType.string,
                                  "file")

par.register ('T',
                                  4,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.int,
                                  "Threads & Streams")

par.register ('e',
                                  -1,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.int,
                                  "Max Events")

par.parseArguments()

files = par.f

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
# process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v26', '')


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files),
    #eventsToProcess = cms.untracked.VEventRange(eventsToProcess),
)

from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
    label        = 'TrackCands',                  # output collection
    tracks       = cms.InputTag('generalTracks'), # input track collection
    particleType = 'pi+',                         # particle type (for assigning a mass)
    preselection = 'pt > 0.7',                    # preselection cut on candidates
    selection    = 'pt > 0.7',                    # selection on PAT Layer 1 objects
    isolation    = {},                            # isolations to use (set to {} for None)
    isoDeposits  = [],
    mcAs         = None                           # replicate MC match as the one used for Muons
)
process.patTrackCands.embedTrack = True
process.patTrackCands.src = cms.InputTag("patAODTrackCandsUnfiltered")

process.trackHLTMatcher = cms.EDProducer(
  # matching in DeltaR, sorting by best DeltaR
  "PATTriggerMatcherDRLessByR"
  # matcher input collections
, src     = cms.InputTag( 'patTrackCands' )
, matched = cms.InputTag( 'patTrigger' )
  # selections of trigger objects
, matchedCuts = cms.string( 'path( "HLT_DoubleMu4_JpsiTrk_Displaced_v*", 1, 0 )' ) # input does not yet have the 'saveTags' parameter in HLT
  # selection of matches
, maxDPtRel   = cms.double( 0.5 ) # no effect here
, maxDeltaR   = cms.double( 0.5 )
, maxDeltaEta = cms.double( 0.2 ) # no effect here
  # definition of matcher output
, resolveAmbiguities    = cms.bool( True )
, resolveByMatchQuality = cms.bool( True )
)



from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, None, None, None, None, '' )
switchOnTriggerMatching( process, triggerMatchers = [ 'trackHLTMatcher' ], outputModule = '' )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(par.e))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('aod_track_analyzer.root'),
)


process.TrackAnalyzer = cms.EDAnalyzer('TrackAnayzerAODWithTrigger',
    Tracks     = cms.InputTag('generalTracks'),
    TriggerFilters = cms.vstring(["hltDoubleMu2JpsiDoubleTrkL3Filtered"]),
    PrimaryVertex    = cms.InputTag("offlinePrimaryVertices"),
    TriggerEvent = cms.InputTag("hltTriggerSummaryAOD"),
    TreeName         = cms.string('Tracks')
)

# End of customisation functions
process.options.numberOfThreads=cms.untracked.uint32(par.T)
process.options.numberOfStreams=cms.untracked.uint32(par.T)
process.options.numberOfConcurrentLuminosityBlocks = 1

process.p = cms.Path(process.patAODTrackCandsUnfiltered+process.patTrackCands+process.TrackAnalyzer)
process.schedule = cms.Schedule([process.p])
