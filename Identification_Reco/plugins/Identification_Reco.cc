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
    
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"


//
// class declaration
//

using namespace std;
using namespace reco;
using namespace muon;

class Identification_Reco : public edm::EDAnalyzer {

   typedef std::vector<reco::Vertex> VertexCollection;
   typedef std::vector<reco::Muon> MuonCollection;
   //typedef std::vector<reco::Jet> JetCollection;
   typedef std::vector<reco::PFJet> JetCollection;

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
      void setEffectiveArea();
      double getEffectiveArea(double eta, double dr);
      
      // TFileService
      edm::Service<TFileService> muisofs;
      //Input Tag
      edm::InputTag vertexTag;
      edm::InputTag muonTag;
      edm::InputTag jetTag;
      edm::InputTag rhoTag;


    private:

      TTree* MuonID;

      std::vector<double> muons_Px;
      std::vector<double> muons_Py;
      std::vector<double> muons_Pz;
      std::vector<double> muons_Pt;
      std::vector<double> muons_eta;
      std::vector<double> muons_phi;
      std::vector<double> muons_pdgid;
      std::vector<double> muons_energy;
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
      //To sort objects in Pt order.
      std::map<double, double> map_muons_Px;
      std::map<double, double> map_muons_Py;
      std::map<double, double> map_muons_Pz;
      std::map<double, double> map_muons_eta;
      std::map<double, double> map_muons_phi;
      std::map<double, double> map_muons_pdgid;
      std::map<double, double> map_muons_energy;
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
      std::map<double, reco::Vertex> map_vtx_sumptsquare;

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

      /// Jet
      unsigned int NJets;
      double rho;

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
:vertexTag( iConfig.getParameter<edm::InputTag>("pvTag") ),
 muonTag( iConfig.getParameter<edm::InputTag>("muTag") ),
 jetTag( iConfig.getParameter<edm::InputTag>("jtTag") ),
 rhoTag( iConfig.getParameter<edm::InputTag>("RhoTag") )
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

   //To read Reco format
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

         NPV++;
         //cout << itPV->position().Rho() << endl;
      }
   }

   //cout << NPV << endl;
   std::sort(vtx_sumptsquare.begin(), vtx_sumptsquare.end(), std::greater<double>());
   /*
   for (unsigned int i=0; i < vtx_sumptsquare.size(); i++)
   {
      cout << map_vtx_sumptsquare.find(vtx_sumptsquare[i])->second.ndof() << "   " << vtx_sumptsquare[i] << endl;
   }
   */

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

   // fill electron branches
   for (MuonCollection::const_iterator itMu = muons->begin() ; itMu != muons->end(); itMu++)
   {
      isLoose = false;
      isSoft = false;
      isTight = false;
      isHighPt = false;
      isFalseSoft = false;
      isFalseTight = false;

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
      
      // isSoftMuon ?
      isSoft = muon::isSoftMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
      //cout << "isSoft " << isSoft << endl;

      // isTightMuon ?
      isTight = muon::isTightMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
      //cout << "isTight " << isTight << endl;
     
      // isHighPtMuon ?
      isHighPt = muon::isHighPtMuon( (*itMu), map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
      //cout << "isHighPt " << isHighPt << endl;

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
   /// Let's find jets      ///
   ////////////////////////////
   Handle<JetCollection> jets;
   iEvent.getByLabel(jetTag, jets);
   //
   // fill jet branches
   for (JetCollection::const_iterator itJet = jets->begin() ; itJet != jets->end(); itJet++)
   {
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

   // Vertex
   MuonID->Branch("NPV", &NPV, "NPV/I");
   MuonID->Branch("vtx_sumptsquare", &vtx_sumptsquare);

   // Muon
   MuonID->Branch("N_muons", &N_muons, "N_muons/I");
   MuonID->Branch("dimuon_inv_mass", &dimuon_inv_mass, "dimuon_inv_mass/D");
   MuonID->Branch("muons_Px", &muons_Px);
   MuonID->Branch("muons_Py", &muons_Py);
   MuonID->Branch("muons_Pz", &muons_Pz);
   MuonID->Branch("muons_Pt", &muons_Pt);
   MuonID->Branch("muons_eta", &muons_eta);
   MuonID->Branch("muons_phi", &muons_phi);
   MuonID->Branch("muons_pdgid", &muons_pdgid);
   MuonID->Branch("muons_energy", &muons_energy);
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

   // Rho
   MuonID->Branch("rho", &rho, "rho/D");

   cout << "Tree was booked!" << endl;
}

void Identification_Reco::Reset()
{
   Event = -333;
   Run = -333;
   Lumi = -333;
   isData = 0;

   NPV = 0;
   N_muons = 0;
   dimuon_inv_mass = -333;
   // vector for muon
   muons_Pt.clear();
   muons_Px.clear();
   muons_Py.clear();
   muons_Pz.clear();
   muons_eta.clear();
   muons_phi.clear();
   muons_energy.clear();

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

   // map for muon
   map_muons_Px.clear();
   map_muons_Py.clear();
   map_muons_Pz.clear();
   map_muons_eta.clear();
   map_muons_phi.clear();
   map_muons_pdgid.clear();
   map_muons_energy.clear();
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

   //Vertex
   vtx_sumptsquare.clear();
   map_vtx_sumptsquare.clear();
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

//define this as a plug-in
DEFINE_FWK_MODULE(Identification_Reco);