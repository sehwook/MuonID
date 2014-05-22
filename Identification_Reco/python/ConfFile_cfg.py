import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonID")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/18BBFD9C-12AE-E311-983E-002618943833.root',

        'file:/afs/cern.ch/work/s/seh/private/MC_samples/upgrade2023/18BBFD9C-12AE-E311-983E-002618943833.root',
        #'file:/afs/cern.ch/work/s/seh/private/MC_samples/upgrade2023/1E87BA41-FDAD-E311-B614-0025905A6070.root'

        #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/18BBFD9C-12AE-E311-983E-002618943833.root',
        #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/1E87BA41-FDAD-E311-B614-0025905A6070.root'
        #'root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC9/RelValFourMuPt1_200/GEN-SIM-RECO/DES19_62_V8_UPG2019-v1/00000/E60624C9-19AE-E311-A962-0025905A60E4.root' 
        #'root://xrootd.unl.edu//store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/2668E6A9-E92C-E211-BA7D-003048D37666.root' 
        #'root://xrootd.unl.edu//store/user/calabria/calabria_SingleMuPt5_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC9_DIGIv7_2023_TeVMuon_NewV2_Case4/calabria_SingleMuPt5_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC9_DIGIv7_2023_TeVMuon_NewV2_Case4/baf4501773a396a3ecf58e228df69850/out_reco_100_1_O6B.root'
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("MuonID.root"),
      #closeFileFast = cms.untracked.bool(True)
)

process.muonid = cms.EDAnalyzer('Identification_Reco',
                              pileupTag = cms.InputTag("addPileupInfo",""),
                              genParTag = cms.InputTag("genParticles",""),
                              bsTag = cms.InputTag("offlineBeamSpot",""),
                              pvTag = cms.InputTag("offlinePrimaryVertices",""),
                              muTag = cms.InputTag("muons",""),
                              eTag = cms.InputTag("gsfElectrons",""),
                              jtTag = cms.InputTag("ak5PFJetsCHS",""),
                              #jtTag = cms.InputTag("ak5PFJets",""),
                              RhoTag = cms.InputTag("kt6PFJetsCentralNeutral","rho"),
                              PFJetID = cms.PSet(
                                                 version = cms.string('FIRSTDATA'),
                                                 quality = cms.string('LOOSE') 
                                                ) 
)


process.p = cms.Path(process.muonid)
