import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonID")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring()
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("MuonID_upgrade2023.root"),
      #closeFileFast = cms.untracked.bool(True)
)

process.muonid = cms.EDAnalyzer('Identification_Reco',
                              bsTag = cms.InputTag("offlineBeamSpot",""),
                              pvTag = cms.InputTag("offlinePrimaryVertices",""),
                              muTag = cms.InputTag("muons",""),
                              eTag = cms.InputTag("gsfElectrons",""),
                              jtTag = cms.InputTag("ak5PFJetsCHS",""),
                              RhoTag = cms.InputTag("kt6PFJetsCentralNeutral","rho"),
                              PFJetID = cms.PSet(
                                                 version = cms.string('FIRSTDATA'),
                                                 quality = cms.string('LOOSE') 
                                                ) 
)


process.p = cms.Path(process.muonid)
