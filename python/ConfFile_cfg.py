import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/tmp/soffi/003757A2-78B3-E711-93D5-0242AC130002.root'
        'file:'
    )
)
process.TFileService = cms.Service("TFileService", 
                                                      fileName = cms.string("output.root")
)
process.demo = cms.EDAnalyzer('DisplacedJetsAnalyzer',
#                              genparts = cms.untracked.InputTag("genParticles","","HLT"),
#                              genjets  = cms.untracked.InputTag("ak4GenJets","","HLT"),
                              jets = cms.untracked.InputTag("ak4PFJets", "", "RECO"),
)


process.p = cms.Path(process.demo)
