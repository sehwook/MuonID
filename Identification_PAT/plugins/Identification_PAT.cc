// -*- C++ -*-
//
// Package:    Identification_PAT
// Class:      Identification_PAT
// 
/**\class Identification_PAT Identification_PAT.cc MuonID/Identification_PAT/plugins/Identification_PAT.cc

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

class Identification_PAT : public edm::EDAnalyzer {

   typedef std::vector<reco::Vertex> VertexCollection;
   typedef std::vector<pat::Muon> MuonCollection;
   typedef std::vector<pat::Jet> JetCollection;
   //typedef std::vector<reco::PFCandidate> JetCollection;

   public:
      explicit Identification_PAT(const edm::ParameterSet&);
      ~Identification_PAT();

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

      // TFileService
      edm::Service<TFileService> muisofs;
      //Input Tag
      edm::InputTag vertexTag;
      edm::InputTag muonTag;
      edm::InputTag jetTag;


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
      std::vector<double> muons_pfiso04_sumChargedHadronPt;
      std::vector<double> muons_pfiso04_sumChargedParticlePt;
      std::vector<double> muons_pfiso04_sumNeutralHadronEt;
      std::vector<double> muons_pfiso04_sumPhotonEt;
      std::vector<double> muons_pfiso04_sumPUPt;
      std::vector<double> muons_PFIsodbeta04;
      std::vector<double> muons_puChargedHadronIso;
      std::vector<bool> muons_isLoose;
      std::vector<bool> muons_isSoft;
      std::vector<bool> muons_isTight;
      std::vector<bool> muons_isHighPt;
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
      std::map<double, double> map_muons_pfiso04_sumChargedHadronPt;
      std::map<double, double> map_muons_pfiso04_sumChargedParticlePt;
      std::map<double, double> map_muons_pfiso04_sumNeutralHadronEt;
      std::map<double, double> map_muons_pfiso04_sumPhotonEt;
      std::map<double, double> map_muons_pfiso04_sumPUPt;
      std::map<double, double> map_muons_PFIsodbeta04;
      std::map<double, double> map_muons_puChargedHadronIso;
      std::map<double, bool> map_muons_isLoose;
      std::map<double, bool> map_muons_isSoft;
      std::map<double, bool> map_muons_isTight;
      std::map<double, bool> map_muons_isHighPt;

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
      
      /// Muon
      unsigned int NPV;
      unsigned int N_muons;
      bool isLoose;
      bool isSoft;
      bool isTight;
      bool isHighPt;
      double trkIso03;
      double caloIso03;
      double relIso03;
      double trkIso05;
      double caloIso05;
      double relIso05;
      double dimuon_inv_mass;
      double PFIsodbeta03;
      double PFIsodbeta04;

      /// Jet
      unsigned int NJets;

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
Identification_PAT::Identification_PAT(const edm::ParameterSet& iConfig)
:vertexTag( iConfig.getParameter<edm::InputTag>("pvTag") ),
 muonTag( iConfig.getParameter<edm::InputTag>("muTag") ),
 jetTag( iConfig.getParameter<edm::InputTag>("jtTag") )
{
   //now do what ever initialization is needed

}


Identification_PAT::~Identification_PAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Identification_PAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Reset variables
   Reset();

   //To read reco format
   //NPV = MuonIso::getNPV(iEvent,iSetup);
   //cout << NPV << endl;

   //To read PAT format
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


   ////////////////////////////
   /// Let's find muons     ///
   ////////////////////////////
   Handle<MuonCollection> muons;
   iEvent.getByLabel(muonTag, muons);

   // fill electron branches
   for (MuonCollection::const_iterator itMu = muons->begin() ; itMu != muons->end(); itMu++)
   {
      /// For pat::Muon
      //
      //
      /*
      cout << (*itMu).pfIsolationR03().sumChargedHadronPt << "   " 
           << (*itMu).pfIsolationR03().sumChargedParticlePt << "  " << (*itMu).pfIsolationR03().sumNeutralHadronEt << "   "
           << (*itMu).pfIsolationR03().sumPhotonEt << "   " << (*itMu).pfIsolationR03().sumPUPt << "   " << (*itMu).puChargedHadronIso() << "   "
           << (*itMu).caloIso() << "   " << (*itMu).isolationR03().emEt + (*itMu).isolationR03().hadEt + (*itMu).isolationR03().hoEt << "   "
           << (*itMu).trackIso() << "   " << (*itMu).isolationR03().sumPt << "   " 
           << (*itMu).pt() << "   " << (*itMu).eta() << "   " << (*itMu).phi() << "   " << (*itMu).pdgId() << "   " << (*itMu).energy() << "   "
           << (*itMu).charge() 
           << endl;
      */

      //cout << "loose " << (*itMu).isLooseMuon() << endl;
      //cout << "soft " << (*itMu).isSoftMuon(*itPV) << endl;
      //cout << "tight " << (*itMu).isTightMuon(*itPV) << endl;
      //cout << "high pt " << (*itMu).isHighPtMuon(*itPV) << endl;

      isLoose = (*itMu).isLooseMuon();
      //if (! isLoose) {continue;}

      isSoft = (*itMu).isSoftMuon( map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
      isTight = (*itMu).isTightMuon( map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );
      isHighPt = (*itMu).isHighPtMuon( map_vtx_sumptsquare.find(vtx_sumptsquare[0])->second );

      //cout << "loose " << isLoose << endl;
      //cout << "soft " << isSoft << endl;
      //cout << "tight " << isTight << endl;
      //cout << "high pt " << isHighPt << endl;

      map_muons_isLoose[(*itMu).pt()] = isLoose;
      map_muons_isSoft[(*itMu).pt()] = isSoft;
      map_muons_isTight[(*itMu).pt()] = isTight;
      map_muons_isHighPt[(*itMu).pt()] = isHighPt;

      if ((isLoose == 0 && isTight == 1) || (isLoose == 0 && isSoft == 1))
      {
         cout << "Tight or Soft muon was found despite it didn't pass Loose Muon requirements!" << endl;
         break;
      }


      trkIso03 = (*itMu).trackIso();
      caloIso03 = (*itMu).caloIso();
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

      map_muons_pfiso04_sumChargedHadronPt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumChargedHadronPt;
      map_muons_pfiso04_sumChargedParticlePt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumChargedParticlePt;
      map_muons_pfiso04_sumNeutralHadronEt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumNeutralHadronEt;
      map_muons_pfiso04_sumPhotonEt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumPhotonEt;
      map_muons_pfiso04_sumPUPt[(*itMu).pt()] = (*itMu).pfIsolationR04().sumPUPt;

      PFIsodbeta04 = ( (*itMu).pfIsolationR04().sumChargedHadronPt +
              std::max( 0., (*itMu).pfIsolationR04().sumNeutralHadronEt + (*itMu).pfIsolationR04().sumPhotonEt - 0.5*(*itMu).pfIsolationR04().sumPUPt ) )
              / (*itMu).pt();
      map_muons_PFIsodbeta04[(*itMu).pt()] = PFIsodbeta04;

      map_muons_puChargedHadronIso[(*itMu).pt()] = (*itMu).puChargedHadronIso();

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

      muons_pfiso04_sumChargedHadronPt.push_back( map_muons_pfiso04_sumChargedHadronPt.find(muons_Pt[i])->second );
      muons_pfiso04_sumChargedParticlePt.push_back( map_muons_pfiso04_sumChargedParticlePt.find(muons_Pt[i])->second );
      muons_pfiso04_sumNeutralHadronEt.push_back( map_muons_pfiso04_sumNeutralHadronEt.find(muons_Pt[i])->second );
      muons_pfiso04_sumPhotonEt.push_back( map_muons_pfiso04_sumPhotonEt.find(muons_Pt[i])->second );
      muons_pfiso04_sumPUPt.push_back( map_muons_pfiso04_sumPUPt.find(muons_Pt[i])->second );
      muons_PFIsodbeta04.push_back( map_muons_PFIsodbeta04.find(muons_Pt[i])->second );

      muons_puChargedHadronIso.push_back( map_muons_puChargedHadronIso.find(muons_Pt[i])->second );

      muons_isLoose.push_back( map_muons_isLoose.find(muons_Pt[i])->second );
      muons_isSoft.push_back( map_muons_isSoft.find(muons_Pt[i])->second );
      muons_isTight.push_back( map_muons_isTight.find(muons_Pt[i])->second );
      muons_isHighPt.push_back( map_muons_isHighPt.find(muons_Pt[i])->second );

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
      map_jets_isPFJet[(*itJet).pt()] = (*itJet).isPFJet();
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

      

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Identification_PAT::beginJob()
{
   MuonID = muisofs->make<TTree>("MuID", "Muon Identification");
   Book();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Identification_PAT::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Identification_PAT::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Identification_PAT::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Identification_PAT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Identification_PAT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Identification_PAT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Identification_PAT::Book()
{

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
   MuonID->Branch("muons_pfiso04_sumChargedHadronPt", &muons_pfiso04_sumChargedHadronPt);
   MuonID->Branch("muons_pfiso04_sumChargedParticlePt", &muons_pfiso04_sumChargedParticlePt);
   MuonID->Branch("muons_pfiso04_sumNeutralHadronEt", &muons_pfiso04_sumNeutralHadronEt);
   MuonID->Branch("muons_pfiso04_sumPhotonEt", &muons_pfiso04_sumPhotonEt);
   MuonID->Branch("muons_pfiso04_sumPUPt", &muons_pfiso04_sumPUPt);
   MuonID->Branch("muons_PFIsodbeta04", &muons_PFIsodbeta04);
   MuonID->Branch("muons_puChargedHadronIso", &muons_puChargedHadronIso);
   MuonID->Branch("muons_isLoose", &muons_isLoose);
   MuonID->Branch("muons_isSoft", &muons_isSoft);
   MuonID->Branch("muons_isTight", &muons_isTight);
   MuonID->Branch("muons_isHighPt", &muons_isHighPt);

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

   cout << "Tree was booked!" << endl;
}

void Identification_PAT::Reset()
{
   NPV = 0;
   N_muons = 0;
   isLoose = 0;
   isSoft = 0;
   isTight = 0;
   isHighPt = 0;
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

   muons_pfiso04_sumChargedHadronPt.clear();
   muons_pfiso04_sumChargedParticlePt.clear();
   muons_pfiso04_sumNeutralHadronEt.clear();
   muons_pfiso04_sumPhotonEt.clear();
   muons_pfiso04_sumPUPt.clear();
   muons_PFIsodbeta04.clear();

   muons_puChargedHadronIso.clear();

   muons_isLoose.clear();
   muons_isSoft.clear();
   muons_isTight.clear();
   muons_isHighPt.clear();

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

   map_muons_pfiso04_sumChargedHadronPt.clear();
   map_muons_pfiso04_sumChargedParticlePt.clear();
   map_muons_pfiso04_sumNeutralHadronEt.clear();
   map_muons_pfiso04_sumPhotonEt.clear();
   map_muons_pfiso04_sumPUPt.clear();
   map_muons_PFIsodbeta04.clear();

   map_muons_puChargedHadronIso.clear();

   map_muons_isLoose.clear();
   map_muons_isSoft.clear();
   map_muons_isTight.clear();
                                                
   map_muons_isHighPt.clear();

   NJets = 0;
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

//define this as a plug-in
DEFINE_FWK_MODULE(Identification_PAT);
