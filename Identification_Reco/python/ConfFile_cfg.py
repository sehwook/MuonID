import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonID")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/18BBFD9C-12AE-E311-983E-002618943833.root',
        #'file:/afs/cern.ch/work/s/seh/private/MC_samples/upgrade2023/18BBFD9C-12AE-E311-983E-002618943833.root',
        'file:/afs/cern.ch/work/s/seh/private/MC_samples/upgrade2023/1E87BA41-FDAD-E311-B614-0025905A6070.root'

        #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/18BBFD9C-12AE-E311-983E-002618943833.root',
        #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/1E87BA41-FDAD-E311-B614-0025905A6070.root'
       #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValFourMuPt1_200/GEN-SIM-RECO/DES19_62_V8_UPG2019-v1/00000/E60624C9-19AE-E311-A962-0025905A60E4.root' 
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("MuonID.root"),
      #closeFileFast = cms.untracked.bool(True)
)

process.muonid = cms.EDAnalyzer('Identification_Reco',
                              pvTag = cms.InputTag("offlinePrimaryVertices",""),
                              muTag = cms.InputTag("muons",""),
                              jtTag = cms.InputTag("ak5PFJetsCHS",""),
                              RhoTag = cms.InputTag("kt6PFJetsCentralNeutral","rho")
)


process.p = cms.Path(process.muonid)
