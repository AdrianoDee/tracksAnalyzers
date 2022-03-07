import FWCore.ParameterSet.Config as cms
process = cms.Process('trackAnalyzerAOD')

from FWCore.ParameterSet.VarParsing import VarParsing

par = VarParsing ('analysis')

par.register ('f',
                                  "file",
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
process.GlobalTag = GlobalTag(process.GlobalTag, '112X_mcRun3_2021_realistic_v7', '')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files),
    #eventsToProcess = cms.untracked.VEventRange(eventsToProcess),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(par.e))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('aod_track_analyzer.root'),
)


process.TrackAnalyzer = cms.EDAnalyzer('TrackAnayzerAOD',
    Tracks     = cms.InputTag('generalTracks'),
    PrimaryVertex    = cms.InputTag("offlinePrimaryVertices"),
    TreeName         = cms.string('Tracks')
)

# End of customisation functions
process.options.numberOfThreads=cms.untracked.uint32(par.T)
process.options.numberOfStreams=cms.untracked.uint32(par.T)
process.options.numberOfConcurrentLuminosityBlocks = 1

process.p = cms.Path(process.TrackAnalyzer)
process.schedule = cms.Schedule([process.p])
