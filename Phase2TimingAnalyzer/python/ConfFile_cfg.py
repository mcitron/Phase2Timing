import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                'file:reco_8.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_phase2timing.root")
)

process.demo = cms.EDAnalyzer('Phase2TimingAnalyzer',
                              genParticles    = cms.InputTag("genParticles", "", "HLT"),
                              recoak4PFJets    = cms.InputTag("ak4PFJets", "", "RECO"),
                              ebRecHitsColl = cms.InputTag( 'ecalRecHit','EcalRecHitsEB',"RECO" ),
)


process.p = cms.Path(process.demo)
