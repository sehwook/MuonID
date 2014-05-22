// -*- C++ -*-
//
// Package:    Identification_Reco
// Class:      Identification_Reco
// 
/**\class Identification_Reco Identification_Reco.cc MuonID/Identification_Reco/plugins/Identification_Reco.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seh Wook Lee
//         Created:  Fri, 18 Apr 2014 01:47:41 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>     // std::cout
#include <algorithm>    // std::max

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
         
// needed for TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

/// needed for pat::Muon, pat::Jet
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

/// ROOT
#include "TTree.h"

//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"    
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// muon ID
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

/// Jet ID
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

//
// class declaration
//

using namespace std;
using namespace reco;
using namespace muon;

class Identification_Reco : public edm::EDAnalyzer {

   typedef std::vector<reco::GenParticle> GenParticleCollection;
   typedef std::vector<reco::Vertex> VertexCollection;
   typedef std::vector<reco::Muon> MuonCollection;
   typedef std::vector<reco::GsfElectron> ElectronCollection;
   //typedef std::vector<reco::Jet> JetCollection;
   typedef std::vector<reco::PFJet> JetCollection;
   typedef std::vector<PileupSummaryInfo> PUInfo;
   typedef reco::BeamSpot beamSpot;

   public:
      explicit Identification_Reco(const edm::ParameterSet&);
      ~Identification_Reco();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      void Book();
      void Reset();
     
      //rho correction
      void setEffectiveArea();
      double getEffectiveArea(double eta, double dr);
      
      // Jet
      bool recoPFJetID(const reco::PFJet & jet, string quality);

      // TFileService
      edm::Service<TFileService> muisofs;
      //Input Tag
      edm::InputTag genParInfoTag;
      edm::InputTag beamspotTag;
      edm::InputTag vertexTag;
      edm::InputTag muonTag;
      edm::InputTag electronTag;
      edm::InputTag jetTag;
      edm::InputTag puTag;
      edm::InputTag rhoTag;
      edm::ParameterSet pfJetIDparam;

    private:

      TTree* MuonID;

      // Generator information
      std::vector<int> genPar_pdgId;
      std::vector<double> genPar_pt;
      std::vector<double> genPar_eta;
      std::vector<double> genPar_phi;
      std::vector<double> genPar_px;
      std::vector<double> genPar_py;
      std::vector<double> genPar_pz;
      std::vector<double> genPar_energy;
      std::vector<int> genPar_charge;
      std::vector<int> genPar_status;
      std::vector<int> genPar_mother_pdgId;
      std::vector<bool> genPar_NoMother;
      
      // muon variable containers
      std::vector<double> muons_Px;
      std::vector<double> muons_Py;
      std::vector<double> muons_Pz;
      std::vector<double> muons_Pt;
      std::vector<double> muons_eta;
      std::vector<double> muons_phi;
      std::vector<double> muons_pdgid;
      std::vector<double> muons_et;
      std::vector<double> muons_energy;
      std::vector<double> muons_energy_dep_em;
      std::vector<double> muons_energy_dep_had;
      std::vector<double> muons_energy_dep_ho;
      std::vector<double> muons_charge;
      std::vector<double> muons_trkIso03;
      std::vector<double> muons_caloIso03;
      std::vector<double> muons_relIso03;
      std::vector<double> muons_trkIso05;
      std::vector<double> muons_caloIso05;
      std::vector<double> muons_relIso05;
      std::vector<double> muons_pfiso03_sumChargedHadronPt;
      std::vector<double> muons_pfiso03_sumChargedParticlePt;
      std::vector<double> muons_pfiso03_sumNeutralHadronEt;
      std::vector<double> muons_pfiso03_sumPhotonEt;
      std::vector<double> muons_pfiso03_sumPUPt;
      std::vector<double> muons_PFIsodbeta03;
      std::vector<double> muons_PFIsorho03;
      std::vector<double> muons_pfiso04_sumChargedHadronPt;
      std::vector<double> muons_pfiso04_sumChargedParticlePt;
      std::vector<double> muons_pfiso04_sumNeutralHadronEt;
      std::vector<double> muons_pfiso04_sumPhotonEt;
      std::vector<double> muons_pfiso04_sumPUPt;
      std::vector<double> muons_PFIsodbeta04;
      std::vector<double> muons_PFIsorho04;
      std::vector<bool> muons_isLoose;
      std::vector<bool> muons_isSoft;
      std::vector<bool> muons_isTight;
      std::vector<bool> muons_isHighPt;
      std::vector<bool> muons_isFalseSoft;
      std::vector<bool> muons_isFalseTight;
      // muon id
      std::vector<bool> muons_isPF;
      std::vector<bool> muons_isGlobal;
      std::vector<bool> muons_isTracker;
      std::vector<bool> muons_isStandAlone;
      std::vector<bool> muons_isThereGlobalTrk;
      std::vector<bool> muons_isThereMuonBestTrk;
      std::vector<bool> muons_isThereInnerTrk;
      std::vector<bool> muons_isThereOuterTrk;
      std::vector<double> muons_globalTrack_normChi2;
      std::vector<int> muons_globalTrack_NValidMuonHits;
      std::vector<int> muons_NMatchedStations;
      std::vector<int> muons_innerTrack_trackerLayersWithMeasurement;
      std::vector<int> muons_innerTrack_NValidPixelHits;
      std::vector<double> muons_BestTrack_dxy;
      std::vector<double> muons_BestTrack_dz;
      std::vector<double> muons_BestTrack_Pt;
      std::vector<double> muons_BestTrack_PtErr;
      std::vector<double> muons_InnerTrk_pt;
      std::vector<double> muons_InnerTrk_eta;
      std::vector<double> muons_InnerTrk_phi;
      std::vector<double> muons_OuterTrk_pt;
      std::vector<double> muons_OuterTrk_eta;
      std::vector<double> muons_OuterTrk_phi;
      //To sort objects in Pt order.
      std::map<double, double> map_muons_Px;
      std::map<double, double> map_muons_Py;
      std::map<double, double> map_muons_Pz;
      std::map<double, double> map_muons_eta;
      std::map<double, double> map_muons_phi;
      std::map<double, double> map_muons_pdgid;
      std::map<double, double> map_muons_et;
      std::map<double, double> map_muons_energy;
      std::map<double, double> map_muons_energy_dep_em;
      std::map<double, double> map_muons_energy_dep_had;
      std::map<double, double> map_muons_energy_dep_ho;
      std::map<double, double> map_muons_charge;
      std::map<double, double> map_muons_trkIso03;
      std::map<double, double> map_muons_caloIso03;
      std::map<double, double> map_muons_relIso03;
      std::map<double, double> map_muons_trkIso05;
      std::map<double, double> map_muons_caloIso05;
      std::map<double, double> map_muons_relIso05;
      std::map<double, double> map_muons_pfiso03_sumChargedHadronPt;
      std::map<double, double> map_muons_pfiso03_sumChargedParticlePt;
      std::map<double, double> map_muons_pfiso03_sumNeutralHadronEt;
      std::map<double, double> map_muons_pfiso03_sumPhotonEt;
      std::map<double, double> map_muons_pfiso03_sumPUPt;
      std::map<double, double> map_muons_PFIsodbeta03;
      std::map<double, double> map_muons_PFIsorho03;
      std::map<double, double> map_muons_pfiso04_sumChargedHadronPt;
      std::map<double, double> map_muons_pfiso04_sumChargedParticlePt;
      std::map<double, double> map_muons_pfiso04_sumNeutralHadronEt;
      std::map<double, double> map_muons_pfiso04_sumPhotonEt;
      std::map<double, double> map_muons_pfiso04_sumPUPt;
      std::map<double, double> map_muons_PFIsodbeta04;
      std::map<double, double> map_muons_PFIsorho04;
      std::map<double, bool> map_muons_isLoose;
      std::map<double, bool> map_muons_isSoft;
      std::map<double, bool> map_muons_isTight;
      std::map<double, bool> map_muons_isHighPt;
      std::map<double, bool> map_muons_isFalseSoft;
      std::map<double, bool> map_muons_isFalseTight;
      // muon id
      std::map<double, bool> map_muons_isPF;
      std::map<double, bool> map_muons_isGlobal;
      std::map<double, bool> map_muons_isTracker;
      std::map<double, bool> map_muons_isStandAlone;
      std::map<double, bool> map_muons_isThereGlobalTrk;
      std::map<double, bool> map_muons_isThereMuonBestTrk;
      std::map<double, bool> map_muons_isThereInnerTrk;
      std::map<double, bool> map_muons_isThereOuterTrk;
      std::map<double, double> map_muons_globalTrack_normChi2;
      std::map<double, int> map_muons_globalTrack_NValidMuonHits;
      std::map<double, int> map_muons_NMatchedStations;
      std::map<double, int> map_muons_innerTrack_trackerLayersWithMeasurement;
      std::map<double, int> map_muons_innerTrack_NValidPixelHits;
      std::map<double, double> map_muons_BestTrack_dxy;
      std::map<double, double> map_muons_BestTrack_dz;
      std::map<double, double> map_muons_BestTrack_Pt;
      std::map<double, double> map_muons_BestTrack_PtErr;
      std::map<double, double> map_muons_InnerTrk_pt;
      std::map<double, double> map_muons_InnerTrk_eta;
      std::map<double, double> map_muons_InnerTrk_phi;
      std::map<double, double> map_muons_OuterTrk_pt;
      std::map<double, double> map_muons_OuterTrk_eta;
      std::map<double, double> map_muons_OuterTrk_phi;

      // Electron variable containers
      std::vector<double> electrons_Pt;
      std::vector<double> electrons_eta;
      std::vector<double> electrons_phi;
      std::vector<double> electrons_energy;
      std::vector<double> electrons_et;
      // To sort object in Pt descending order
      std::map<double, double> map_electrons_eta;
      std::map<double, double> map_electrons_phi;
      std::map<double, double> map_electrons_energy;
      std::map<double, double> map_electrons_et;
   
      // Jet variable containers
      std::vector<double> jets_Px;
      std::vector<double> jets_Py;
      std::vector<double> jets_Pz;
      std::vector<double> jets_Pt;
      std::vector<double> jets_eta;
      std::vector<double> jets_phi;
      std::vector<double> jets_et;
      std::vector<double> jets_energy;
      std::vector<double> jets_charge;
      std::vector<double> jets_chargedMultiplicity;
      std::vector<double> jets_isPFJet;
      std::vector<double> jets_pdgid;
      std::vector<double> jets_isJet;
      std::vector<double> diff_jetE_muE;
      std::vector<double> frac_muE_jetE;
      std::vector<double> muPt_inJet;
      std::vector<double> frac_eleE_jetE;
      //To sort object in Pt order
      std::map<double, double> map_jets_Px;
      std::map<double, double> map_jets_Py;
      std::map<double, double> map_jets_Pz;
      std::map<double, double> map_jets_eta;
      std::map<double, double> map_jets_phi;
      std::map<double, double> map_jets_et;
      std::map<double, double> map_jets_energy;
      std::map<double, double> map_jets_charge;
      std::map<double, double> map_jets_chargedMultiplicity;
      std::map<double, double> map_jets_isPFJet;
      std::map<double, double> map_jets_pdgid;
      std::map<double, double> map_jets_isJet;

      ///Vertex
      std::vector<double> vtx_sumptsquare;
      std::vector<double> vtx_rho;
      std::vector<double> vtx_z;
      std::map<double, reco::Vertex> map_vtx_sumptsquare;
      std::map<double, double> map_vtx_rho;
      std::map<double, double> map_vtx_z;

      /// Pile Up Info
      std::vector<int> bunchX;
      std::vector<int> N_PU;
      std::vector<double> N_trueInt;
                         
      // Effective Area
      std::map<unsigned int, double> EffectiveArea03;
      std::map<unsigned int, double> EffectiveArea04;

      //
      ////////////////////////////////////////////////////////
      //
      // Event info.
      int Event;
      int Run;
      int Lumi;
      bool isData;

      /// Generator information
      int N_genPar;

      /// Muon
      unsigned int NPV;
      unsigned int N_muons;
      bool isLoose;
      bool isSoft;
      bool isTight;
      bool isHighPt;
      bool isFalseSoft;
      bool isFalseTight;
      double trkIso03;
      double caloIso03;
      double relIso03;
      double trkIso05;
      double caloIso05;
      double relIso05;
      double dimuon_inv_mass;
      double PFIsodbeta03;
      double PFIsodbeta04;
      double PFIsorho03;
      double PFIsorho04;
      double effA03;
      double effA04;
      double mu_vtx_dxy;
      double mu_vtx_dz;
      //muon id
      bool _isGlobal;
      bool _isPF;
      bool _isTracker;
      bool _isStandAlone;
      bool _isThereGlobalTrk;
      bool _isThereMuonBestTrk;
      bool _isThereInnerTrk;
      bool _isThereOuterTrk;
      double _globalTrack_normChi2;
      int _globalTrack_NValidMuonHits;
      int _NMatchedStations;
      int _innerTrack_trackerLayersWithMeasurement;
      int _innerTrack_NValidPixelHits;
      double _BestTrack_dxy;
      double _BestTrack_dz;
      double _BestTrack_Pt;
      double _BestTrack_PtErr;
      double _inner_trk_pt;
      double _inner_trk_eta;
      double _inner_trk_phi;
      double _outer_trk_pt;
      double _outer_trk_eta;
      double _outer_trk_phi;

      // Electron
      unsigned int N_electrons;

      /// Jet
      bool jetid_pass;
      bool jet_isMuon;
      bool jet_iselectron;
      unsigned int NJets;
      double rho;
      double dR;

      // Beam Spot
      double beamSpot_x;
      double beamSpot_y;
      double beamSpot_z;

      //Vertex
      double sumtrackptsquare;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Identification_Reco::Identification_Reco(const edm::ParameterSet& iConfig)
:genParInfoTag( iConfig.getParameter<edm::InputTag>("genParTag") ),
 beamspotTag( iConfig.getParameter<edm::InputTag>("bsTag") ),
 vertexTag( iConfig.getParameter<edm::InputTag>("pvTag") ),
 muonTag( iConfig.getParameter<edm::InputTag>("muTag") ),
 electronTag( iConfig.getParameter<edm::InputTag>("eTag") ),
 jetTag( iConfig.getParameter<edm::InputTag>("jtTag") ),
 puTag( iConfig.getParameter<edm::InputTag>("pileupTag") ),
 rhoTag( iConfig.getParameter<edm::InputTag>("RhoTag") ),
 pfJetIDparam( iConfig.getParameter<edm::ParameterSet> ("PFJetID") )
{
   //now do what ever initialization is needed
   setEffectiveArea();
}


Identification_Reco::~Identification_Reco()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Identification_Reco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Reset variables
   Reset();

   Event = iEvent.id().event();
   Run = iEvent.id().run();
   Lumi = iEvent.id().luminosityBlock();
   isData = iEvent.isRealData();
   //cout << isData << "   " << EVENT << "   " << RUN << "   " << LUMI << endl;

   // Generator Particle Information
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel(genParInfoTag, genParticles);
   for ( GenParticleCollection::const_iterator itGenPar = genParticles->begin(); itGenPar != genParticles->end(); itGenPar++ )
   {
      /*
      cout << (*itGenPar).pdgId() << "   " << (*itGenPar).pt() << "   " << (*itGenPar).eta()  
           << (*itGenPar).phi() << "   " << (*itGenPar).px() << "   " << (*itGenPar).py() << "   "  
           << (*itGenPar).pz() << "   " << (*itGenPar).charge() << "    " << (*itGenPar).status() << "   "
           << (*itGenPar).energy() << "   " << (*itGenPar).numberOfMothers() << endl;
      */
      /*
      int n_mothers = (*itGenPar).numberOfMothers();
      for(int j = 0; j < n_mothers; j++) 
      {
         const Candidate * mom = (*itGenPar).mother( j );
         int momId = mom->pdgId();
         cout << momId << endl;
      }
      */

      genPar_pdgId.push_back( (*itGenPar).pdgId() );
      genPar_pt.push_back( (*itGenPar).pt() );
      genPar_eta.push_back( (*itGenPar).eta() );
      genPar_phi.push_back( (*itGenPar).phi() );
      genPar_px.push_back( (*itGenPar).px() );
      genPar_py.push_back( (*itGenPar).py() );
      genPar_pz.push_back( (*itGenPar).pz() );
      genPar_energy.push_back( (*itGenPar).energy() );
      genPar_charge.push_back( (*itGenPar).charge() );
      genPar_status.push_back( (*itGenPar).status() );

      if ( (*itGenPar).numberOfMothers() == 0)
      {
         genPar_mother_pdgId.push_back( (*itGenPar).pdgId() );
         genPar_NoMother.push_back(true);

      } else {
         const Candidate *mom = (*itGenPar).mother(0);
         genPar_mother_pdgId.push_back( mom->pdgId() );
         genPar_NoMother.push_back(false);

      }

      N_genPar++;
   }

   // Beam Spot
   Handle<beamSpot> BS;
   iEvent.getByLabel(beamspotTag, BS);
   //cout << (*BS).position() << "   " << (*BS).position().x() << "   " << (*BS).x0() << "   " << (*BS).y0() << "   " << (*BS).z0() << endl;
   beamSpot_x = (*BS).position().x();
   beamSpot_y = (*BS).position().y();
   beamSpot_z = (*BS).position().z();

   /// PileUp information
   Handle<PUInfo> pPUInfo;
   iEvent.getByLabel(puTag, pPUInfo);
   for (PUInfo::const_iterator itPU = pPUInfo->begin(); itPU != pPUInfo->end(); itPU++)
   {
      //cout << itPU->getBunchCrossing() << "   " << itPU->getPU_NumInteractions() << "   " << itPU->getTrueNumInteractions() << "   " 
      //     << itPU->getPU_zpositions()[0] << endl;

      // -1: previous BX, 0: current BX,  1: next BX
      bunchX.push_back(itPU->getBunchCrossing());
      // # of PU vertices
      N_PU.push_back(itPU->getPU_NumInteractions());
      // True # of interactions
      N_trueInt.push_back(itPU->getPU_NumInteractions());
   }

   //Find Primary Vertices
   Handle<VertexCollection> vertices;
   iEvent.getByLabel(vertexTag, vertices);
   for (VertexCollection::const_iterator itPV = vertices->begin(); itPV != vertices->end(); itPV++)
   {
      //if ( !itPV->isFake() && itPV->ndof() > 4.0 && itPV->position().Rho() < 2. && abs(itPV->z()) < 24. )
      if ( !itPV->isFake() && itPV->ndof() > 4.0 )
      {
         sumtrackptsquare = 0;
         for(reco::Vertex::trackRef_iterator trkItr = itPV->tracks_begin(); trkItr != itPV->tracks_end(); ++trkItr)
         {
            sumtrackptsquare += (*trkItr)->pt()*(*trkItr)->pt();
            //cout << (*trkItr)->pt() << endl; 
         }

         vtx_sumptsquare.push_back(sumtrackptsquare);
         map_vtx_sumptsquare[sumtrackptsquare] = *itPV;
         map_vtx_rho[sumtrackptsquare] = (*itPV).position().Rho();
         map_vtx_z[sumtrackptsquare] = (*itPV).z();

         NPV++;
         //cout << itPV->position().Rho() << endl;
      }
   }

   //cout << NPV << endl;
   std::sort(vtx_sumptsquare.begin(), vtx_sumptsquare.end(), std::greater<double>());
   
   for (unsigned int i=0; i < vtx_sumptsquare.size(); i++)
   {
      //cout << map_vtx_sumptsquare.find(vtx_sumptsquare[i])->second.ndof() << "   " << vtx_sumptsquare[i] << endl;
      vtx_rho.push_back( map_vtx_rho.find(vtx_sumptsquare[i])->second );
      vtx_z.push_back( map_vtx_z.find(vtx_sumptsquare[i])->second );
   }
   

   //////////////////////
   ///       Rho      ///
   //////////////////////
   Handle<double> rhoHandle;
   iEvent.getByLabel(rhoTag, rhoHandle);

   if(rhoHandle.isValid()) 
   {
      rho = *(rhoHandle.product());
      //cout << rho << "   " << *rhoHandle << endl;
   }

   ////////////////////////////
   /// Let's find muons     ///
   ////////////////////////////
   Handle<MuonCollection> muons;
   iEvent.getByLabel(muonTag, muons);

   // fill muon branches
   for (MuonCollection::const_iterator itMu = muons->begin() ; itMu != muons->end(); itMu++)
   {
      isLoose = false;
      isSoft = false;
      isTight = false;
      isHighPt = false;
      isFalseSoft = false;
      isFalseTight = false;

      _isGlobal = false;
      _isPF = false;   
      _isTracker = false;
      _isStandAlone = false;
      _isThereGlobalTrk = false;
      _isThereMuonBestTrk = false;
      _isThereInnerTrk = false;
      _isThereOuterTrk = false;
      _globalTrack_normChi2 = -333.;
      _globalTrack_NValidMuonHits = -333;
      _NMatchedStations = -333;
      _innerTrack_trackerLayersWithMeasurement = -333;
      _innerTrack_NValidPixelHits = -333;
      _BestTrack_dxy = -333333.;
      _BestTrack_dz = -333333.;
      _BestTrack_Pt = -333333.;
      _BestTrack_PtErr = -333333.;
      _inner_trk_pt = -333333.;
      _inner_trk_eta = -333333.;
      _inner_trk_phi = -333333.;
      _outer_trk_pt = -333333.;
      _outer_trk_eta = -333333.;
      _outer_trk_phi = -333333.;

      mu_vtx_dxy = -333;
      mu_vtx_dz = -333;

      // isLooseMuon ?
      /*
      if ( (*itMu).isPFMuon() && ( (*itMu).isGlobalMuon() || (*itMu).isTrackerMuon() ) )
      {
         isLoose = 1;
      } else {
         isLoose = 0;
      }
      */
      isLoose = muon::isLooseMuon( (*itMu) );
      //cout << "isLoose " << isLoose << "   " << muon::isLooseMuon((*itMu)) << endl;
   
      if (NPV == 0)
      {   
         // isSoftMuon ?
         isSoft = muon::isSoftMuonWithBS( (*itMu), (*BS) );

         // isTightMuon ?
         isTight = muon::isTightMuonWithBS( (*itMu), (*BS) );
     
         // isHighPtMuon ?
         isHighPt = muon::isHighPtMuonWithBS( (*itMu), (*BS) );
          
      } else {
         // isSoftMuon ?
         isSoft = muon::isSoftMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
         //cout << "isSoft " << isSoft << endl;

         // isTightMuon ?
         isTight = muon::isTightMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
         //cout << "isTight " << isTight << endl;
     
         // isHighPtMuon ?
         isHighPt = muon::isHighPtMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
         //cout << "isHighPt " << isHighPt << endl;
      }

      // For sanity check
      if (isLoose == 0 && isSoft == 1)
      {
         isFalseSoft = true;
         cout << "Soft muon was found despite it didn't pass Loose Muon requirements!" << endl;
         //break;
      } 

      if (isLoose == 0 && isTight == 1) 
      {
         isFalseTight = true;
         cout << "Tight muon was found despite it didn't pass Loose Muon requirements!" << endl;
         //break;
      }

      //cout << muon::isGoodMuon((*itMu), TMOneStationTight) << "    " << TMOneStationTight << endl;

      // Muon ID
      //cout << (*itMu).isGlobalMuon() << "   " << (*itMu).isPFMuon() << endl;
      /*
      if (!(*itMu).globalTrack().isNull())
      {
         cout << (*itMu).globalTrack()->normalizedChi2() << "   " << (*itMu).globalTrack()->hitPattern().numberOfValidMuonHits() << endl;
      }
      */
      //cout << (*itMu).numberOfMatchedStations() << "   " << (*itMu).innerTrack()->hitPattern().trackerLayersWithMeasurement() << "   "
      //     << (*itMu).innerTrack()->hitPattern().numberOfValidPixelHits() << endl;
      //cout << (*itMu).muonBestTrack()->dxy(map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second.position()) << "   " 
      //     << (*itMu).muonBestTrack()->dz(map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second.position()) << endl;
      //cout << (*itMu).muonBestTrack()->dxy((*BS).position()) << "   " 
      //     << (*itMu).muonBestTrack()->dz((*BS).position()) << endl;
      _isGlobal = (*itMu).isGlobalMuon();
      _isPF = (*itMu).isPFMuon();
      _isTracker = (*itMu).isTrackerMuon();
      _isStandAlone = (*itMu).isStandAloneMuon() ;
      _NMatchedStations = (*itMu).numberOfMatchedStations();
      if (!(*itMu).globalTrack().isNull())
      {
         _isThereGlobalTrk = true;
         _globalTrack_normChi2 = (*itMu).globalTrack()->normalizedChi2();
         _globalTrack_NValidMuonHits = (*itMu).globalTrack()->hitPattern().numberOfValidMuonHits();
      }
      if (!(*itMu).muonBestTrack().isNull())
      {
         _isThereMuonBestTrk = true;
         _BestTrack_PtErr = (*itMu).muonBestTrack()->ptError();
         _BestTrack_Pt = (*itMu).muonBestTrack()->pt();
         if (NPV == 0)
         {
            _BestTrack_dxy = (*itMu).muonBestTrack()->dxy((*BS).position());
            _BestTrack_dz = (*itMu).muonBestTrack()->dz((*BS).position());

         } else {
            _BestTrack_dxy = (*itMu).muonBestTrack()->dxy(map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second.position());
            _BestTrack_dz = (*itMu).muonBestTrack()->dz(map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second.position());

         }
      }
      if (!(*itMu).innerTrack().isNull())
      {
         _isThereInnerTrk = true;
         _innerTrack_trackerLayersWithMeasurement = (*itMu).innerTrack()->hitPattern().trackerLayersWithMeasurement();
         _innerTrack_NValidPixelHits = (*itMu).innerTrack()->hitPattern().numberOfValidPixelHits();
         _inner_trk_pt = (*itMu).innerTrack()->pt();
         _inner_trk_eta = (*itMu).innerTrack()->eta();
         _inner_trk_phi = (*itMu).innerTrack()->phi();
      }
      if (!(*itMu).outerTrack().isNull())
      {
         _isThereOuterTrk = true;
         _outer_trk_pt = (*itMu).outerTrack()->pt();
         _outer_trk_eta = (*itMu).outerTrack()->eta();
         _outer_trk_phi = (*itMu).outerTrack()->phi();
      }
      map_muons_isGlobal[(*itMu).pt()] = _isGlobal;
      map_muons_isPF[(*itMu).pt()] = _isPF;
      map_muons_isTracker[(*itMu).pt()] = _isTracker;
      map_muons_isStandAlone[(*itMu).pt()] = _isStandAlone; 
      map_muons_isThereGlobalTrk[(*itMu).pt()] = _isThereGlobalTrk;
      map_muons_isThereMuonBestTrk[(*itMu).pt()] = _isThereMuonBestTrk;
      map_muons_isThereInnerTrk[(*itMu).pt()] = _isThereInnerTrk;
      map_muons_isThereOuterTrk[(*itMu).pt()] = _isThereOuterTrk;
      map_muons_globalTrack_normChi2[(*itMu).pt()] = _globalTrack_normChi2;
      map_muons_globalTrack_NValidMuonHits[(*itMu).pt()] = _globalTrack_NValidMuonHits;
      map_muons_NMatchedStations[(*itMu).pt()] = _NMatchedStations;
      map_muons_innerTrack_trackerLayersWithMeasurement[(*itMu).pt()] = _innerTrack_trackerLayersWithMeasurement;
      map_muons_innerTrack_NValidPixelHits[(*itMu).pt()] = _innerTrack_NValidPixelHits;
      map_muons_BestTrack_dxy[(*itMu).pt()] = _BestTrack_dxy;
      map_muons_BestTrack_dz[(*itMu).pt()] = _BestTrack_dz;
      map_muons_BestTrack_Pt[(*itMu).pt()] = _BestTrack_Pt;
      map_muons_BestTrack_PtErr[(*itMu).pt()] = _BestTrack_PtErr;
      map_muons_InnerTrk_pt[(*itMu).pt()] = _inner_trk_pt;
      map_muons_InnerTrk_eta[(*itMu).pt()] = _inner_trk_eta;
      map_muons_InnerTrk_phi[(*itMu).pt()] = _inner_trk_phi;
      map_muons_OuterTrk_pt[(*itMu).pt()] = _outer_trk_pt;
      map_muons_OuterTrk_eta[(*itMu).pt()] = _outer_trk_eta;
      map_muons_OuterTrk_phi[(*itMu).pt()] = _outer_trk_phi;

      map_muons_isLoose[(*itMu).pt()] = isLoose;
      map_muons_isSoft[(*itMu).pt()] = isSoft;
      map_muons_isTight[(*itMu).pt()] = isTight;
      map_muons_isHighPt[(*itMu).pt()] = isHighPt;
      map_muons_isFalseSoft[(*itMu).pt()] = isFalseSoft;
      map_muons_isFalseTight[(*itMu).pt()] = isFalseTight;

      trkIso03 = (*itMu).isolationR03().sumPt;
      caloIso03 = (*itMu).isolationR03().emEt + (*itMu).isolationR03().hadEt;
      relIso03 = ( trkIso03 + caloIso03 ) / (*itMu).pt();

      trkIso05 = (*itMu).isolationR05().sumPt;
      caloIso05 = (*itMu).isolationR05().emEt + (*itMu).isolationR05().hadEt;
      relIso05 = ( trkIso05 + caloIso05 ) / (*itMu).pt();

      muons_Pt.push_back((*itMu).pt());

      map_muons_Px[(*itMu).pt()] = (*itMu).px();
      map_muons_Py[(*itMu).pt()] = (*itMu).py();
      map_muons_Pz[(*itMu).pt()] = (*itMu).pz();
      map_muons_eta[(*itMu).pt()] = (*itMu).eta();
      map_muons_phi[(*itMu).pt()] = (*itMu).phi();
      map_muons_energy[(*itMu).pt()] = (*itMu).energy();
      map_muons_et[(*itMu).pt()] = (*itMu).et();

      map_muons_pdgid[(*itMu).pt()] = (*itMu).pdgId();
      map_muons_charge[(*itMu).pt()] = (*itMu).charge();

      map_muons_trkIso03[(*itMu).pt()] = trkIso03;
      map_muons_caloIso03[(*itMu).pt()] = caloIso03;
      map_muons_relIso03[(*itMu).pt()] = relIso03;

      map_muons_trkIso05[(*itMu).pt()] = trkIso05;
      map_muons_caloIso05[(*itMu).pt()] = caloIso05;
      map_muons_relIso05[(*itMu).pt()] = relIso05;

      map_muons_pfiso03_sumChargedHadronPt[(*itMu).pt()] = (*itMu).pfIsolationR03().sumChargedHadronPt;
      map_muons_pfiso03_sumChargedParticlePt[(*itMu).pt()] = (*itMu).pfIsolationR03().sumChargedParticlePt;
      map_muons_pfiso03_sumNeutralHadronEt[(*itMu).pt()] = (*itMu).pfIsolationR03().sumNeutralHadronEt;
      map_muons_pfiso03_sumPhotonEt[(*itMu).pt()] = (*itMu).pfIsolationR03().sumPhotonEt;
      map_muons_pfiso03_sumPUPt[(*itMu).pt()] = (*itMu).pfIsolationR03().sumPUPt;

      PFIsodbeta03 = ( (*itMu).pfIsolationR03().sumChargedHadronPt +
              std::max( 0., (*itMu).pfIsolationR03().sumNeutralHadronEt + (*itMu).pfIsolationR03().sumPhotonEt - 0.5*(*itMu).pfIsolationR03().sumPUPt ) )
              / (*itMu).pt();
      map_muons_PFIsodbeta03[(*itMu).pt()] = PFIsodbeta03;

      effA03 = getEffectiveArea( fabs( (*itMu).eta() ), 0.3 ); //cout << effA03 << endl;
      PFIsorho03 = ( (*itMu).pfIsolationR03().sumChargedHadronPt +
              std::max( 0., (*itMu).pfIsolationR03().sumNeutralHadronEt + (*itMu).pfIsolationR03().sumPhotonEt - max(0.0, rho)*effA03 ) )
              / (*itMu).pt();
      map_muons_PFIsorho03[(*itMu).pt()] = PFIsorho03;

      map_muons_pfiso04_sumChargedHadronPt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumChargedHadronPt;
      map_muons_pfiso04_sumChargedParticlePt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumChargedParticlePt;
      map_muons_pfiso04_sumNeutralHadronEt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumNeutralHadronEt;
      map_muons_pfiso04_sumPhotonEt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumPhotonEt;
      map_muons_pfiso04_sumPUPt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumPUPt;

      PFIsodbeta04 = ( (*itMu).pfIsolationR04().sumChargedHadronPt +
              std::max( 0., (*itMu).pfIsolationR04().sumNeutralHadronEt + (*itMu).pfIsolationR04().sumPhotonEt - 0.5*(*itMu).pfIsolationR04().sumPUPt ) )
              / (*itMu).pt();
      map_muons_PFIsodbeta04[(*itMu).pt()] = PFIsodbeta04;

      effA04 = getEffectiveArea( fabs( (*itMu).eta() ), 0.4 ); //cout << effA04 << endl;
      PFIsorho04 = ( (*itMu).pfIsolationR04().sumChargedHadronPt +
              std::max( 0., (*itMu).pfIsolationR04().sumNeutralHadronEt + (*itMu).pfIsolationR04().sumPhotonEt - max(0.0, rho)*effA04 ) )
              / (*itMu).pt();
      map_muons_PFIsorho04[(*itMu).pt()] = PFIsorho04;

      //cout << (*itMu).calEnergy().em << "   " << (*itMu).calEnergy().had << "   " << (*itMu).calEnergy().ho << endl;
      map_muons_energy_dep_em[(*itMu).pt()] = (*itMu).calEnergy().em;
      map_muons_energy_dep_had[(*itMu).pt()] = (*itMu).calEnergy().had;
      map_muons_energy_dep_ho[(*itMu).pt()] = (*itMu).calEnergy().ho;

      N_muons++;
   }

   // Sort Muon Pt in descending order
   std::sort(muons_Pt.begin(), muons_Pt.end(), std::greater<double>());
   for (unsigned int i=0; i < muons_Pt.size(); i++)
   {
      muons_Px.push_back( map_muons_Px.find(muons_Pt[i])->second );
      muons_Py.push_back( map_muons_Py.find(muons_Pt[i])->second );
      muons_Pz.push_back( map_muons_Pz.find(muons_Pt[i])->second );
      muons_energy.push_back( map_muons_energy.find(muons_Pt[i])->second );
      muons_energy_dep_em.push_back( map_muons_energy_dep_em.find(muons_Pt[i])->second );
      muons_energy_dep_had.push_back( map_muons_energy_dep_had.find(muons_Pt[i])->second );
      muons_energy_dep_ho.push_back( map_muons_energy_dep_ho.find(muons_Pt[i])->second );
      muons_et.push_back( map_muons_et.find(muons_Pt[i])->second );
      muons_eta.push_back( map_muons_eta.find(muons_Pt[i])->second );
      muons_phi.push_back( map_muons_phi.find(muons_Pt[i])->second );

      muons_pdgid.push_back( map_muons_pdgid.find(muons_Pt[i])->second );
      muons_charge.push_back( map_muons_charge.find(muons_Pt[i])->second );

      muons_trkIso03.push_back( map_muons_trkIso03.find(muons_Pt[i])->second );
      muons_caloIso03.push_back( map_muons_caloIso03.find(muons_Pt[i])->second );
      muons_relIso03.push_back( map_muons_relIso03.find(muons_Pt[i])->second );

      muons_trkIso05.push_back( map_muons_trkIso05.find(muons_Pt[i])->second );
      muons_caloIso05.push_back( map_muons_caloIso05.find(muons_Pt[i])->second );
      muons_relIso05.push_back( map_muons_relIso05.find(muons_Pt[i])->second );

      muons_pfiso03_sumChargedHadronPt.push_back( map_muons_pfiso03_sumChargedHadronPt.find(muons_Pt[i])->second );
      muons_pfiso03_sumChargedParticlePt.push_back( map_muons_pfiso03_sumChargedParticlePt.find(muons_Pt[i])->second );
      muons_pfiso03_sumNeutralHadronEt.push_back( map_muons_pfiso03_sumNeutralHadronEt.find(muons_Pt[i])->second );
      muons_pfiso03_sumPhotonEt.push_back( map_muons_pfiso03_sumPhotonEt.find(muons_Pt[i])->second );
      muons_pfiso03_sumPUPt.push_back( map_muons_pfiso03_sumPUPt.find(muons_Pt[i])->second );
      muons_PFIsodbeta03.push_back( map_muons_PFIsodbeta03.find(muons_Pt[i])->second );
      muons_PFIsorho03.push_back( map_muons_PFIsorho03.find(muons_Pt[i])->second );

      muons_pfiso04_sumChargedHadronPt.push_back( map_muons_pfiso04_sumChargedHadronPt.find(muons_Pt[i])->second );
      muons_pfiso04_sumChargedParticlePt.push_back( map_muons_pfiso04_sumChargedParticlePt.find(muons_Pt[i])->second );
      muons_pfiso04_sumNeutralHadronEt.push_back( map_muons_pfiso04_sumNeutralHadronEt.find(muons_Pt[i])->second );
      muons_pfiso04_sumPhotonEt.push_back( map_muons_pfiso04_sumPhotonEt.find(muons_Pt[i])->second );
      muons_pfiso04_sumPUPt.push_back( map_muons_pfiso04_sumPUPt.find(muons_Pt[i])->second );
      muons_PFIsodbeta04.push_back( map_muons_PFIsodbeta04.find(muons_Pt[i])->second );
      muons_PFIsorho04.push_back( map_muons_PFIsorho04.find(muons_Pt[i])->second );

      muons_isLoose.push_back( map_muons_isLoose.find(muons_Pt[i])->second );
      muons_isSoft.push_back( map_muons_isSoft.find(muons_Pt[i])->second );
      muons_isTight.push_back( map_muons_isTight.find(muons_Pt[i])->second );
      muons_isHighPt.push_back( map_muons_isHighPt.find(muons_Pt[i])->second );
      muons_isFalseSoft.push_back( map_muons_isFalseSoft.find(muons_Pt[i])->second );
      muons_isFalseTight.push_back( map_muons_isFalseTight.find(muons_Pt[i])->second );

      // Muon ID
      muons_isGlobal.push_back( map_muons_isGlobal.find(muons_Pt[i])->second );
      muons_isTracker.push_back( map_muons_isTracker.find(muons_Pt[i])->second );
      muons_isStandAlone.push_back( map_muons_isStandAlone.find(muons_Pt[i])->second );
      muons_isPF.push_back( map_muons_isPF.find(muons_Pt[i])->second );
      muons_isThereGlobalTrk.push_back( map_muons_isThereGlobalTrk.find(muons_Pt[i])->second );
      muons_isThereMuonBestTrk.push_back( map_muons_isThereMuonBestTrk.find(muons_Pt[i])->second );
      muons_isThereInnerTrk.push_back( map_muons_isThereInnerTrk.find(muons_Pt[i])->second );
      muons_isThereOuterTrk.push_back( map_muons_isThereOuterTrk.find(muons_Pt[i])->second );
      muons_globalTrack_normChi2.push_back( map_muons_globalTrack_normChi2.find(muons_Pt[i])->second );
      muons_globalTrack_NValidMuonHits.push_back( map_muons_globalTrack_NValidMuonHits.find(muons_Pt[i])->second );
      muons_NMatchedStations.push_back( map_muons_NMatchedStations.find(muons_Pt[i])->second );
      muons_innerTrack_trackerLayersWithMeasurement.push_back( map_muons_innerTrack_trackerLayersWithMeasurement.find(muons_Pt[i])->second );
      muons_innerTrack_NValidPixelHits.push_back( map_muons_innerTrack_NValidPixelHits.find(muons_Pt[i])->second );
      muons_BestTrack_dxy.push_back( map_muons_BestTrack_dxy.find(muons_Pt[i])->second );
      muons_BestTrack_dz.push_back( map_muons_BestTrack_dz.find(muons_Pt[i])->second );
      muons_BestTrack_Pt.push_back( map_muons_BestTrack_Pt.find(muons_Pt[i])->second );
      muons_BestTrack_PtErr.push_back( map_muons_BestTrack_PtErr.find(muons_Pt[i])->second );
      muons_InnerTrk_pt.push_back( map_muons_InnerTrk_pt.find(muons_Pt[i])->second );
      muons_InnerTrk_eta.push_back( map_muons_InnerTrk_eta.find(muons_Pt[i])->second );
      muons_InnerTrk_phi.push_back( map_muons_InnerTrk_phi.find(muons_Pt[i])->second );
      muons_OuterTrk_pt.push_back( map_muons_OuterTrk_pt.find(muons_Pt[i])->second );
      muons_OuterTrk_eta.push_back( map_muons_OuterTrk_eta.find(muons_Pt[i])->second );
      muons_OuterTrk_phi.push_back( map_muons_OuterTrk_phi.find(muons_Pt[i])->second );

      //cout << i << "   " << muons_Pt[i] << "   " << muons_Px[i] 
      //     << map_muons_isTight.find(muons_Pt[i])->second <<"   "
      //     << endl;
   }

   if ( N_muons >= 2 )
   {
      dimuon_inv_mass = sqrt( (muons_energy[0]+muons_energy[1])*(muons_energy[0]+muons_energy[1])
                             -(muons_Px[0]+muons_Px[1])*(muons_Px[0]+muons_Px[1])
                             -(muons_Py[0]+muons_Py[1])*(muons_Py[0]+muons_Py[1])
                             -(muons_Pz[0]+muons_Pz[1])*(muons_Pz[0]+muons_Pz[1])
                            );
   }


   ////////////////////////////
   /// Let's find electrons ///
   ////////////////////////////
   Handle<ElectronCollection> electrons;
   iEvent.getByLabel(electronTag, electrons);
   
   // fill muon branches
   for (ElectronCollection::const_iterator itele = electrons->begin() ; itele != electrons->end(); itele++)
   {
      electrons_Pt.push_back((*itele).pt());

      map_electrons_eta[(*itele).pt()] = (*itele).eta();
      map_electrons_phi[(*itele).pt()] = (*itele).phi();
      map_electrons_energy[(*itele).pt()] = (*itele).energy();   
      map_electrons_et[(*itele).pt()] = (*itele).et();   

      N_electrons++;
   }
                        
   // Sort Muon Pt in descending order
   std::sort(electrons_Pt.begin(), electrons_Pt.end(), std::greater<double>());
   for (unsigned int i=0; i < electrons_Pt.size(); i++)
   {
      electrons_energy.push_back( map_electrons_energy.find(electrons_Pt[i])->second );
      electrons_et.push_back( map_electrons_et.find(electrons_Pt[i])->second );
      electrons_eta.push_back( map_electrons_eta.find(electrons_Pt[i])->second );
      electrons_phi.push_back( map_electrons_phi.find(electrons_Pt[i])->second );
   }



   ////////////////////////////
   /// Let's find jets      ///
   ////////////////////////////


   // Utility for Jet ID
   //PFJetIDSelectionFunctor JetID(pfJetIDparam);
   //pat::strbitset looseJetIdSel = JetID.getBitTemplate();

   Handle<JetCollection> jets;
   iEvent.getByLabel(jetTag, jets);
   //
   // fill jet branches
   for (JetCollection::const_iterator itJet = jets->begin() ; itJet != jets->end(); itJet++)
   {

      ////////////////////////////////////
      /// selection of resaonable jets ///
      ////////////////////////////////////
      //
      
      if ( abs( (*itJet).eta() ) > 2.5 ) {continue;}

      if ( (*itJet).pt() < 20 ) {continue;}

      //
      jetid_pass = false;
      jet_isMuon = false;
      jet_iselectron = false;

      ///For PFJetIDSelectionFunctor
      //looseJetIdSel.set(false);
      //jetid_pass = JetID(*itJet, looseJetIdSel);
 
      // Jet ID using my function
      jetid_pass = recoPFJetID(*itJet, "Loose");

      //cout << jetid_pass << "   " << recoPFJetID(*itJet, "Loose") << endl;

      if (!jetid_pass) {continue;}
     
      // let's find funcking muons and electrons which were stored as jets
      // Muon
      for (unsigned int j=0;j<N_muons;j++)
      {
         dR = deltaR(muons_eta[j], muons_phi[j], (*itJet).eta(), (*itJet).phi());
         //sqrt( (muons_eta[j]-jets_eta[i])*(muons_eta[j]-jets_eta[i]) + (muons_phi[j]-jets_phi[i])*(muons_phi[j]-jets_phi[i]) );
         if (dR < 0.5)
         {
            //cout << "muon " << muons_eta[j] << "   " << muons_phi[j] << "   " << muons_Pt[j] << endl;
            //cout << "jet " << (*itJet).eta() << "   " << (*itJet).phi() << "   " << (*itJet).pt() << endl;
            //if ( (*itJet).energy() - muons_energy[j] < 0 ){
            //cout << "muon " << muons_phi[j] << "   " << muons_eta[j] << "   " << muons_Pt[j] << "   " << muons_energy[j] << endl; 
            //cout << "jet " << (*itJet).phi() << "   " << (*itJet).eta() << "   " << (*itJet).pt() << "   " << (*itJet).energy() << endl; 
            //cout << "jet - muon " << (*itJet).energy() - muons_energy[j] << endl; 
            //}
            diff_jetE_muE.push_back( (*itJet).energy() - muons_energy[j] );
            frac_muE_jetE.push_back( muons_energy[j]/(*itJet).energy() );
            muPt_inJet.push_back( muons_Pt[j] );
            if ( muons_energy[j]/(*itJet).energy() > 0.8 )
            {
               jet_isMuon = true; 
               break;
            }
         }
      }
   
      if (jet_isMuon) {continue;}

      // electron
      for (unsigned int k=0;k<N_electrons;k++)
      {
         dR = deltaR(electrons_eta[k], electrons_phi[k], (*itJet).eta(), (*itJet).phi());
        
         if (dR < 0.5)
         {
            //cout << "electron " << electrons_phi[k] << "   " << electrons_eta[k] << "   " << electrons_Pt[k] << "   " << electrons_energy[k] << endl;
            //cout << "jet " << (*itJet).phi() << "   " << (*itJet).eta() << "   " << (*itJet).pt() << "   " << (*itJet).energy() << endl;
            //cout << "electron energy fraction " << (*itJet).electronEnergyFraction() << "   " << (*itJet).chargedEmEnergy() << "   " << electrons_Pt[k]/(*itJet).et() << endl; 
           // frac_eleE_jetE.push_back( electrons_energy[k]/(*itJet).energy() );
            if ( electrons_energy[k]/(*itJet).energy() > 0.95 )
            {
               jet_iselectron = true; 
               break;
            }
         }
      }

      if (jet_iselectron) {continue;}
      ////////////////////////////
      /// End of jet selection ///
      ////////////////////////////

      /// put jet information into ntuple
      jets_Pt.push_back( (*itJet).pt() );

      map_jets_Px[(*itJet).pt()] = (*itJet).px();
      map_jets_Py[(*itJet).pt()] = (*itJet).py();
      map_jets_Pz[(*itJet).pt()] = (*itJet).pz();
      map_jets_eta[(*itJet).pt()] = (*itJet).eta();
      map_jets_phi[(*itJet).pt()] = (*itJet).phi();
      map_jets_et[(*itJet).pt()] = (*itJet).et();
      map_jets_energy[(*itJet).pt()] = (*itJet).energy();
      map_jets_charge[(*itJet).pt()] = (*itJet).charge();
      /// works only for JPT or PF jet
      //map_jets_chargedMultiplicity[(*itJet).pt()] = (*itJet).chargedMultiplicity();
      //map_jets_isPFJet[(*itJet).pt()] = (*itJet).isPFJet();
      map_jets_pdgid[(*itJet).pt()] = (*itJet).pdgId();
      map_jets_isJet[(*itJet).pt()] = (*itJet).isJet();

      NJets++;
   }

   // Sort Muon Pt in descending order
   std::sort(jets_Pt.begin(), jets_Pt.end(), std::greater<double>());
   for (unsigned int i=0;i<jets_Pt.size();i++)
   {
      jets_Px.push_back( map_jets_Px.find(jets_Pt[i])->second );
      jets_Py.push_back( map_jets_Py.find(jets_Pt[i])->second );
      jets_Pz.push_back( map_jets_Pz.find(jets_Pt[i])->second );
      jets_eta.push_back( map_jets_eta.find(jets_Pt[i])->second );
      jets_phi.push_back( map_jets_phi.find(jets_Pt[i])->second );
      jets_et.push_back( map_jets_et.find(jets_Pt[i])->second );
      jets_energy.push_back( map_jets_energy.find(jets_Pt[i])->second );
      jets_charge.push_back( map_jets_charge.find(jets_Pt[i])->second );
      jets_isPFJet.push_back( map_jets_isPFJet.find(jets_Pt[i])->second );
      jets_pdgid.push_back( map_jets_pdgid.find(jets_Pt[i])->second );
      jets_isJet.push_back( map_jets_isJet.find(jets_Pt[i])->second );

      //cout << i << "   " << jets_Pt[i] << "   " << jets_Px[i] 
      //     << endl;
   }

   /*
   double dR;
   /// Jet and muon
   for (unsigned int i=0;i<NJets;i++)
   {
      for (unsigned int j=0;j<N_muons;j++)
      {
         dR = deltaR(muons_eta[j], muons_phi[j], jets_eta[i], jets_phi[i]);
         //sqrt( (muons_eta[j]-jets_eta[i])*(muons_eta[j]-jets_eta[i]) + (muons_phi[j]-jets_phi[i])*(muons_phi[j]-jets_phi[i]) );
         if (dR < 0.05 && muons_Pt[j] < 15)
         {
            cout << "muon " << muons_eta[j] << "   " << muons_phi[j] << "   " << muons_Pt[j] << endl;
            cout << "jet " << jets_eta[j] << "   " << jets_phi[j] << "   " << jets_Pt[j] << endl;
         }
      }
   }
   */

   /////////////////
   ///Fill ntuple///
   /////////////////
   MuonID->Fill();


      
}


// ------------ method called once each job just before starting event loop  ------------
void 
Identification_Reco::beginJob()
{
   MuonID = muisofs->make<TTree>("MuID", "Muon Identification");
   Book();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Identification_Reco::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Identification_Reco::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Identification_Reco::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Identification_Reco::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Identification_Reco::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Identification_Reco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Identification_Reco::Book()
{

   // Event info.
   MuonID->Branch("Event", &Event, "Event/I");
   MuonID->Branch("Run", &Run, "Run/I");
   MuonID->Branch("Lumi", &Lumi, "Lumi/I");
   MuonID->Branch("isData", &isData, "isData/O");

   // Generator information
   MuonID->Branch("N_genPar", &N_genPar, "N_genPar/I");
   MuonID->Branch("genPar_pdgId", &genPar_pdgId);
   MuonID->Branch("genPar_pt", &genPar_pt);
   MuonID->Branch("genPar_eta", &genPar_eta);
   MuonID->Branch("genPar_phi", &genPar_phi);
   MuonID->Branch("genPar_px", &genPar_px);
   MuonID->Branch("genPar_py", &genPar_py);
   MuonID->Branch("genPar_pz", &genPar_pz);
   MuonID->Branch("genPar_energy", &genPar_energy);
   MuonID->Branch("genPar_charge", &genPar_charge);
   MuonID->Branch("genPar_mother_pdgId", &genPar_mother_pdgId);
   MuonID->Branch("genPar_NoMother", &genPar_NoMother);
   MuonID->Branch("genPar_status", &genPar_status);

   // Beam Position
   MuonID->Branch("beamSpot_x", &beamSpot_x,"beamSpot_x/D");
   MuonID->Branch("beamSpot_y", &beamSpot_y,"beamSpot_y/D");
   MuonID->Branch("beamSpot_z", &beamSpot_z,"beamSpot_z/D");

   // Pile UP
   MuonID->Branch("bunchX", &bunchX);
   MuonID->Branch("N_PU", &N_PU);
   MuonID->Branch("N_trueInt", &N_trueInt);

   // Vertex
   MuonID->Branch("NPV", &NPV, "NPV/I");
   MuonID->Branch("vtx_sumptsquare", &vtx_sumptsquare);
   MuonID->Branch("vtx_rho", &vtx_rho);
   MuonID->Branch("vtx_z", &vtx_z);

   // Muon
   MuonID->Branch("N_muons", &N_muons, "N_muons/I");
   MuonID->Branch("dimuon_inv_mass", &dimuon_inv_mass, "dimuon_inv_mass/D");
   MuonID->Branch("muons_Px", &muons_Px);
   MuonID->Branch("muons_Py", &muons_Py);
   MuonID->Branch("muons_Pz", &muons_Pz);
   MuonID->Branch("muons_Pt", &muons_Pt);
   MuonID->Branch("muons_eta", &muons_eta);
   MuonID->Branch("muons_et", &muons_et);
   MuonID->Branch("muons_phi", &muons_phi);
   MuonID->Branch("muons_pdgid", &muons_pdgid);
   MuonID->Branch("muons_energy", &muons_energy);
   MuonID->Branch("muons_energy_dep_em", &muons_energy_dep_em);
   MuonID->Branch("muons_energy_dep_had", &muons_energy_dep_had);
   MuonID->Branch("muons_energy_dep_ho", &muons_energy_dep_ho);
   MuonID->Branch("muons_charge", &muons_charge);
   MuonID->Branch("muons_trkIso03", &muons_trkIso03);
   MuonID->Branch("muons_caloIso03", &muons_caloIso03);
   MuonID->Branch("muons_relIso03", &muons_relIso03);
   MuonID->Branch("muons_trkIso05", &muons_trkIso05);
   MuonID->Branch("muons_caloIso05", &muons_caloIso05);
   MuonID->Branch("muons_relIso05", &muons_relIso05);
   MuonID->Branch("muons_pfiso03_sumChargedHadronPt", &muons_pfiso03_sumChargedHadronPt);
   MuonID->Branch("muons_pfiso03_sumChargedParticlePt", &muons_pfiso03_sumChargedParticlePt);
   MuonID->Branch("muons_pfiso03_sumNeutralHadronEt", &muons_pfiso03_sumNeutralHadronEt);
   MuonID->Branch("muons_pfiso03_sumPhotonEt", &muons_pfiso03_sumPhotonEt);
   MuonID->Branch("muons_pfiso03_sumPUPt", &muons_pfiso03_sumPUPt);
   MuonID->Branch("muons_PFIsodbeta03", &muons_PFIsodbeta03);
   MuonID->Branch("muons_PFIsorho03", &muons_PFIsorho03);
   MuonID->Branch("muons_pfiso04_sumChargedHadronPt", &muons_pfiso04_sumChargedHadronPt);
   MuonID->Branch("muons_pfiso04_sumChargedParticlePt", &muons_pfiso04_sumChargedParticlePt);
   MuonID->Branch("muons_pfiso04_sumNeutralHadronEt", &muons_pfiso04_sumNeutralHadronEt);
   MuonID->Branch("muons_pfiso04_sumPhotonEt", &muons_pfiso04_sumPhotonEt);
   MuonID->Branch("muons_pfiso04_sumPUPt", &muons_pfiso04_sumPUPt);
   MuonID->Branch("muons_PFIsodbeta04", &muons_PFIsodbeta04);
   MuonID->Branch("muons_PFIsorho04", &muons_PFIsorho04);
   MuonID->Branch("muons_isLoose", &muons_isLoose);
   MuonID->Branch("muons_isSoft", &muons_isSoft);
   MuonID->Branch("muons_isTight", &muons_isTight);
   MuonID->Branch("muons_isHighPt", &muons_isHighPt);
   MuonID->Branch("muons_isFalseSoft", &muons_isFalseSoft);
   MuonID->Branch("muons_isFalseTight", &muons_isFalseTight);
   MuonID->Branch("mu_vtx_dxy", &mu_vtx_dxy, "mu_vtx_dxy/D");
   MuonID->Branch("mu_vtx_dz", &mu_vtx_dz, "mu_vtx_dz/D");
   // Muon ID
   MuonID->Branch("muons_isGlobal", &muons_isGlobal);
   MuonID->Branch("muons_isTracker", &muons_isTracker);
   MuonID->Branch("muons_isStandAlone", &muons_isStandAlone);
   MuonID->Branch("muons_isPF", &muons_isPF);
   MuonID->Branch("muons_isThereGlobalTrk", &muons_isThereGlobalTrk);
   MuonID->Branch("muons_isThereMuonBestTrk", &muons_isThereMuonBestTrk);
   MuonID->Branch("muons_isThereInnerTrk", &muons_isThereInnerTrk);
   MuonID->Branch("muons_isThereOuterTrk", &muons_isThereOuterTrk);
   MuonID->Branch("muons_globalTrack_normChi2", &muons_globalTrack_normChi2);
   MuonID->Branch("muons_globalTrack_NValidMuonHits", &muons_globalTrack_NValidMuonHits);
   MuonID->Branch("muons_NMatchedStations", &muons_NMatchedStations);
   MuonID->Branch("muons_innerTrack_trackerLayersWithMeasurement", &muons_innerTrack_trackerLayersWithMeasurement);
   MuonID->Branch("muons_innerTrack_NValidPixelHits", &muons_innerTrack_NValidPixelHits);
   MuonID->Branch("muons_BestTrack_dxy", &muons_BestTrack_dxy);
   MuonID->Branch("muons_BestTrack_dz", &muons_BestTrack_dz);
   MuonID->Branch("muons_BestTrack_Pt", &muons_BestTrack_Pt);
   MuonID->Branch("muons_BestTrack_PtErr", &muons_BestTrack_PtErr);
   MuonID->Branch("muons_InnerTrk_pt", &muons_InnerTrk_pt);
   MuonID->Branch("muons_InnerTrk_eta", &muons_InnerTrk_eta);
   MuonID->Branch("muons_InnerTrk_phi", &muons_InnerTrk_phi);
   MuonID->Branch("muons_OuterTrk_pt", &muons_OuterTrk_pt);
   MuonID->Branch("muons_OuterTrk_eta", &muons_OuterTrk_eta);
   MuonID->Branch("muons_OuterTrk_phi", &muons_OuterTrk_phi);

   // Electron
   MuonID->Branch("N_electrons", &N_electrons, "N_electrons/I");
   MuonID->Branch("electrons_Pt", &electrons_Pt);
   MuonID->Branch("electrons_eta", &electrons_eta);
   MuonID->Branch("electrons_phi", &electrons_phi);
   MuonID->Branch("electrons_energy", &electrons_energy);
   MuonID->Branch("electrons_et", &electrons_et);

   // Jet
   MuonID->Branch("NJets", &NJets, "NJets/I");
   MuonID->Branch("jets_Px", &jets_Px);
   MuonID->Branch("jets_Py", &jets_Py);
   MuonID->Branch("jets_Pz", &jets_Pz);
   MuonID->Branch("jets_Pt", &jets_Pt);
   MuonID->Branch("jets_eta", &jets_eta);
   MuonID->Branch("jets_phi", &jets_phi);
   MuonID->Branch("jets_et", &jets_et);
   MuonID->Branch("jets_energy", &jets_energy);
   MuonID->Branch("jets_charge", &jets_charge);
   MuonID->Branch("jets_chargedMultiplicity", &jets_chargedMultiplicity);
   MuonID->Branch("jets_isPFJet", &jets_isPFJet);
   MuonID->Branch("jets_pdgid", &jets_pdgid);
   MuonID->Branch("jets_isJet", &jets_isJet);
   MuonID->Branch("diff_jetE_muE", &diff_jetE_muE);
   MuonID->Branch("frac_muE_jetE", &frac_muE_jetE);
   MuonID->Branch("muPt_inJet", &muPt_inJet);
   MuonID->Branch("frac_eleE_jetE", &frac_eleE_jetE);

   // Rho
   MuonID->Branch("rho", &rho, "rho/D");

   cout << "Tree was booked!" << endl;
}

void Identification_Reco::Reset()
{
   /// Basic event info.
   Event = -333;
   Run = -333;
   Lumi = -333;
   isData = 0;

   /// Generator information
   N_genPar = 0;
   // vector for generator particles
   genPar_pdgId.clear();
   genPar_pt.clear();
   genPar_eta.clear();
   genPar_phi.clear();
   genPar_px.clear();
   genPar_py.clear();
   genPar_pz.clear();
   genPar_energy.clear();
   genPar_charge.clear();
   genPar_mother_pdgId.clear();
   genPar_status.clear();
   genPar_NoMother.clear();

   /// Beam Spot
   beamSpot_x = -333333.;
   beamSpot_y = -333333.;
   beamSpot_z = -333333.;

   /// Primary Vertex
   NPV = 0;

   ///Muon
   N_muons = 0;
   dimuon_inv_mass = -333;
   // vector for muon
   muons_Pt.clear();
   muons_Px.clear();
   muons_Py.clear();
   muons_Pz.clear();
   muons_eta.clear();
   muons_et.clear();
   muons_phi.clear();
   muons_energy.clear();
   muons_energy_dep_em.clear();
   muons_energy_dep_had.clear();
   muons_energy_dep_ho.clear();

   muons_pdgid.clear();
   muons_charge.clear();

   muons_trkIso03.clear();
   muons_caloIso03.clear();
   muons_relIso03.clear();

   muons_trkIso05.clear();
   muons_caloIso05.clear();
   muons_relIso05.clear();

   muons_pfiso03_sumChargedHadronPt.clear();
   muons_pfiso03_sumChargedParticlePt.clear();
   muons_pfiso03_sumNeutralHadronEt.clear();
   muons_pfiso03_sumPhotonEt.clear();
   muons_pfiso03_sumPUPt.clear();
   muons_PFIsodbeta03.clear();
   muons_PFIsorho03.clear();

   muons_pfiso04_sumChargedHadronPt.clear();
   muons_pfiso04_sumChargedParticlePt.clear();
   muons_pfiso04_sumNeutralHadronEt.clear();
   muons_pfiso04_sumPhotonEt.clear();
   muons_pfiso04_sumPUPt.clear();
   muons_PFIsodbeta04.clear();
   muons_PFIsorho04.clear();

   muons_isLoose.clear();
   muons_isSoft.clear();
   muons_isTight.clear();
   muons_isHighPt.clear();
   muons_isFalseSoft.clear();
   muons_isFalseTight.clear();

   muons_isGlobal.clear();
   muons_isTracker.clear();
   muons_isStandAlone.clear();
   muons_isPF.clear();
   muons_isThereGlobalTrk.clear();
   muons_isThereMuonBestTrk.clear();
   muons_isThereInnerTrk.clear();
   muons_isThereOuterTrk.clear();
   muons_globalTrack_normChi2.clear();
   muons_globalTrack_NValidMuonHits.clear();
   muons_NMatchedStations.clear();
   muons_innerTrack_trackerLayersWithMeasurement.clear();
   muons_innerTrack_NValidPixelHits.clear();
   muons_BestTrack_dxy.clear();
   muons_BestTrack_dz.clear();
   muons_BestTrack_Pt.clear();
   muons_BestTrack_PtErr.clear();
   muons_InnerTrk_pt.clear();
   muons_InnerTrk_eta.clear();
   muons_InnerTrk_phi.clear();
   muons_OuterTrk_pt.clear();
   muons_OuterTrk_eta.clear();
   muons_OuterTrk_phi.clear();

   // map for muon
   map_muons_Px.clear();
   map_muons_Py.clear();
   map_muons_Pz.clear();
   map_muons_eta.clear();
   map_muons_et.clear();
   map_muons_phi.clear();
   map_muons_pdgid.clear();
   map_muons_energy.clear();
   map_muons_energy_dep_em.clear();
   map_muons_energy_dep_had.clear();
   map_muons_energy_dep_ho.clear();
   map_muons_charge.clear();

   map_muons_trkIso03.clear();
   map_muons_caloIso03.clear();
   map_muons_relIso03.clear();

   map_muons_trkIso05.clear();
   map_muons_caloIso05.clear();
   map_muons_relIso05.clear();

   map_muons_pfiso03_sumChargedHadronPt.clear();
   map_muons_pfiso03_sumChargedParticlePt.clear();
   map_muons_pfiso03_sumNeutralHadronEt.clear();
   map_muons_pfiso03_sumPhotonEt.clear();
   map_muons_pfiso03_sumPUPt.clear();
   map_muons_PFIsodbeta03.clear();
   map_muons_PFIsorho03.clear();

   map_muons_pfiso04_sumChargedHadronPt.clear();
   map_muons_pfiso04_sumChargedParticlePt.clear();
   map_muons_pfiso04_sumNeutralHadronEt.clear();
   map_muons_pfiso04_sumPhotonEt.clear();
   map_muons_pfiso04_sumPUPt.clear();
   map_muons_PFIsodbeta04.clear();
   map_muons_PFIsorho04.clear();

   map_muons_isLoose.clear();
   map_muons_isSoft.clear();
   map_muons_isTight.clear();
   map_muons_isHighPt.clear();
   map_muons_isFalseSoft.clear();
   map_muons_isFalseTight.clear();
 
   map_muons_isGlobal.clear();
   map_muons_isTracker.clear();
   map_muons_isStandAlone.clear();
   map_muons_isPF.clear();
   map_muons_isThereGlobalTrk.clear();
   map_muons_isThereMuonBestTrk.clear();
   map_muons_isThereInnerTrk.clear();
   map_muons_isThereOuterTrk.clear();
   map_muons_globalTrack_normChi2.clear();
   map_muons_globalTrack_NValidMuonHits.clear();
   map_muons_NMatchedStations.clear();
   map_muons_innerTrack_trackerLayersWithMeasurement.clear();
   map_muons_innerTrack_NValidPixelHits.clear();   
   map_muons_BestTrack_dxy.clear();
   map_muons_BestTrack_dz.clear();
   map_muons_BestTrack_Pt.clear();
   map_muons_BestTrack_PtErr.clear();
   map_muons_InnerTrk_pt.clear();
   map_muons_InnerTrk_eta.clear();
   map_muons_InnerTrk_phi.clear();
   map_muons_OuterTrk_pt.clear();
   map_muons_OuterTrk_eta.clear();
   map_muons_OuterTrk_phi.clear();

   // electron
   N_electrons = 0;
   electrons_Pt.clear();
   electrons_eta.clear();
   electrons_phi.clear();
   electrons_energy.clear();
   electrons_et.clear();
   map_electrons_eta.clear();
   map_electrons_phi.clear();
   map_electrons_energy.clear();
   map_electrons_et.clear();
   
   /// Jet                      
   NJets = 0;
   rho = 0.;
   //vector for jet
   jets_Px.clear();
   jets_Py.clear();
   jets_Pz.clear();
   jets_Pt.clear();
   jets_eta.clear();
   jets_phi.clear();
   jets_et.clear();
   jets_energy.clear();
   jets_charge.clear();
   jets_chargedMultiplicity.clear();
   jets_isPFJet.clear();
   jets_pdgid.clear();
   jets_isJet.clear();
   diff_jetE_muE.clear();
   frac_muE_jetE.clear();
   frac_eleE_jetE.clear();
   muPt_inJet.clear();
   //map for jet
   map_jets_Px.clear();
   map_jets_Py.clear();
   map_jets_Pz.clear();
   map_jets_eta.clear();
   map_jets_phi.clear();
   map_jets_et.clear();
   map_jets_energy.clear();
   map_jets_charge.clear();
   map_jets_chargedMultiplicity.clear();
   map_jets_isPFJet.clear();
   map_jets_pdgid.clear();
   map_jets_isJet.clear();

   // Pile Up
   bunchX.clear();
   N_PU.clear();
   N_trueInt.clear();

   //Vertex
   vtx_sumptsquare.clear();
   vtx_rho.clear();
   vtx_z.clear();
   map_vtx_sumptsquare.clear();
   map_vtx_rho.clear();
   map_vtx_z.clear();
}


void Identification_Reco::setEffectiveArea()
{
   EffectiveArea03[10] = .01;  
   EffectiveArea03[1015] = .01;  
   EffectiveArea03[1520] = .01;  
   EffectiveArea03[2022] = .01;  
   EffectiveArea03[2223] = .01;  
   EffectiveArea03[2324] = .01;  

   EffectiveArea04[10] = .01;
   EffectiveArea04[1015] = .01;
   EffectiveArea04[1520] = .01;
   EffectiveArea04[2022] = .01;
   EffectiveArea04[2223] = .01;
   EffectiveArea04[2324] = .01;
}



double Identification_Reco::getEffectiveArea(double eta, double dr)
{
   double EA = 0;
  
   if (dr == 0.3)
   { 
      if (eta >= 0.0 && eta < 1.0)
      {
         EA = EffectiveArea03.find(10)->second;

      } else if (eta >= 1.0 && eta < 1.5) {
         EA = EffectiveArea03.find(1015)->second;

      } else if (eta >= 1.5 && eta < 2.0) {
         EA = EffectiveArea03.find(1520)->second;

      } else if (eta >= 2.0 && eta < 2.2) {
         EA = EffectiveArea03.find(2022)->second;

      } else if (eta >= 2.2 && eta < 2.3) {
         EA = EffectiveArea03.find(2223)->second;

      } else if (eta >= 2.3 && eta < 2.4) {
         EA = EffectiveArea03.find(2324)->second;

      } else { EA = 0; }

   } else if (dr == 0.4) {
      if (eta >= 0.0 && eta < 1.0)
      {
         EA = EffectiveArea04.find(10)->second;

      } else if (eta >= 1.0 && eta < 1.5) {
         EA = EffectiveArea04.find(1015)->second;

      } else if (eta >= 1.5 && eta < 2.0) {
         EA = EffectiveArea04.find(1520)->second;

      } else if (eta >= 2.0 && eta < 2.2) {
         EA = EffectiveArea04.find(2022)->second;

      } else if (eta >= 2.2 && eta < 2.3) {
         EA = EffectiveArea04.find(2223)->second;

      } else if (eta >= 2.3 && eta < 2.4) {
         EA = EffectiveArea04.find(2324)->second;

      } else { EA = 0; }
   }


   return EA;
}



bool Identification_Reco::recoPFJetID(const reco::PFJet & Jet, string Quality)
{
   bool pass = false;

   if (Quality == "Loose")
   {
     if ( Jet.numberOfDaughters() > 1 
         && ( ( Jet.neutralHadronEnergy() + Jet.HFHadronEnergy() ) / Jet.energy() ) < 0.99
         && Jet.neutralEmEnergyFraction() < 0.99
         && ( Jet.chargedEmEnergyFraction() < 0.99 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedHadronEnergyFraction() > 0 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedMultiplicity() > 0 || abs(Jet.eta()) > 2.4 )
        ) {pass = true;}

   } else if (Quality == "Medium") {
      if ( Jet.numberOfDaughters() > 1
         && ( ( Jet.neutralHadronEnergy() + Jet.HFHadronEnergy() ) / Jet.energy() ) < 0.95
         && Jet.neutralEmEnergyFraction() < 0.95
         && ( Jet.chargedEmEnergyFraction() < 0.99 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedHadronEnergyFraction() > 0 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedMultiplicity() > 0 || abs(Jet.eta()) > 2.4 )
        ) {pass = true;}

   } else if (Quality == "Tight") {
      if ( Jet.numberOfDaughters() > 1
         && ( ( Jet.neutralHadronEnergy() + Jet.HFHadronEnergy() ) / Jet.energy() ) < 0.90
         && Jet.neutralEmEnergyFraction() < 0.90
         && ( Jet.chargedEmEnergyFraction() < 0.99 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedHadronEnergyFraction() > 0 || abs(Jet.eta()) > 2.4 )
         && ( Jet.chargedMultiplicity() > 0 || abs(Jet.eta()) > 2.4 )
        ) {pass = true;}

   }

   return pass;
}



//define this as a plug-in
DEFINE_FWK_MODULE(Identification_Reco);
