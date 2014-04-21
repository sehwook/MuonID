import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonID")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/18BBFD9C-12AE-E311-983E-002618943833.root',
        #'/store/relval/CMSSW_6_2_0_SLHC9/RelValZMM_14TeV/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/1E87BA41-FDAD-E311-B614-0025905A6070.root'
        #'file:/afs/cern.ch/work/s/seh/private/PAT_Production/CMSSW_5_3_11/src/PhysicsTools/PatAlgos/patTuple/patTuple_addVertexInfo.root'
        #'file:/afs/cern.ch/work/s/seh/private/PAT_Production/CMSSW_5_3_11/src/PhysicsTools/PatAlgos/patTuple/patTuple_standard.root'
         #'file:/afs/cern.ch/work/s/seh/private/PAT_Production/CMSSW_5_3_11/src/PhysicsTools/PatAlgos/patTuple/patTuple_addJets.root'
         'file:/afs/cern.ch/work/s/seh/private/PAT_Production/CMSSW_5_3_11/src/PhysicsTools/PatAlgos/patTuple/patTuple_PF2PAT.root'
        # '/store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/2668E6A9-E92C-E211-BA7D-003048D37666.root',
        # '/store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/2E11C1DD-E82C-E211-8FC5-003048D37694.root',
        # '/store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/8E233104-EA2C-E211-A65B-0030486730C6.root'
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("MuonID.root"),
      #closeFileFast = cms.untracked.bool(True)
)

process.muonid = cms.EDAnalyzer('Identification_PAT',
                              pvTag = cms.InputTag("offlinePrimaryVertices",""),#bestVertex","")
                              muTag = cms.InputTag("selectedPatMuonsPFlow",""),#cleanPatMuons",""),
                              jtTag = cms.InputTag("selectedPatJetsPFlow","")#cleanPatJets","")
)


process.p = cms.Path(process.muonid)
