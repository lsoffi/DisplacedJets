
// -*- C++ -*-
//
// Package:    Analyzer/DisplacedJetsAnalyzer
// Class:      DisplacedJetsAnalyzer
// 
/**\class DisplacedJetsAnalyzer DisplacedJetsAnalyzer.cc Analyzer/DisplacedJetsAnalyzer/plugins/DisplacedJetsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Wed, 11 Apr 2018 19:18:43 GMT
//
//


// system include files
#include <memory>

#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"                                                                                                                                   
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
// user include files                                                                                                                                                                                 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/src/one/implementorsMethods.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
//#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"//
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/JetReco/interface/PFJet.h"



// ROOT                                                                                                                                                                                                
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
//
// class declaration
//
struct Tree_struc_{
  int event;
  int lumi;
  int run;
  int ngenjets;
  std::vector<float> genjet_e;
  std::vector<float> genjet_pt;
  std::vector<float> genjet_eta;
  std::vector<float> genjet_phi;
  std::vector<float> genjet_numdoug;
  std::vector<float> genjet_numconst;
  std::vector<int>   genjet_electrons;
  std::vector<int>   genjet_muons;
  std::vector<int>   genjet_taus;
  std::vector<int>   genjet_quarks;
  std::vector<int>   genjet_gluons;
  std::vector<int>   genjet_photons;
  std::vector<int>   genjet_pizeros;
  std::vector<int>   genjet_pipm;
  std::vector<int>   genjet_kzeros;
  std::vector<int>   genjet_kpm;
  std::vector<int>   genjet_protons;
  std::vector<int>   genjet_neutrons;
  std::vector<int>   genjet_etas;
  std::vector<int>   genjet_lambdas;
  std::vector<int>   genjet_chHad;
  int njets;
  std::vector<float> jet_e;
  std::vector<float> jet_pt;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;

};

class DisplacedJetsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DisplacedJetsAnalyzer(const edm::ParameterSet&);
      ~DisplacedJetsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void initTreeStructure();  
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;
  //  edm::EDGetTokenT<View<reco::GenParticle> > genPartToken_;                                                                                                                                         
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genPartToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genJetToken_;
  edm::EDGetTokenT<std::vector<reco::PFJet> > jetToken_;    
  // event info                                                                                                                                                                                         
      unsigned long int event;
      int run, lumi;
      int ngenjet;
      std::vector<float> genjet_e;
      std::vector<float> genjet_pt;
      std::vector<float> genjet_eta;
      std::vector<float> genjet_phi;
      std::vector<float> genjet_numdoug;
      std::vector<float> genjet_numconst;

  std::vector<int>   genjet_electrons;
  std::vector<int>   genjet_muons;
  std::vector<int>   genjet_taus;
  std::vector<int>   genjet_quarks;
  std::vector<int>   genjet_gluons;
  std::vector<int>   genjet_photons;
  std::vector<int>   genjet_pizeros;
  std::vector<int>   genjet_pipm;
  std::vector<int>   genjet_kzeros;
  std::vector<int>   genjet_kpm;
  std::vector<int>   genjet_protons;
  std::vector<int>   genjet_neutrons;
  std::vector<int>   genjet_etas;
  std::vector<int>   genjet_lambdas;
  std::vector<int>   genjet_chHad;

      int njet;
      std::vector<float> jet_e;
      std::vector<float> jet_pt;
      std::vector<float> jet_eta;
      std::vector<float> jet_phi;




      TTree* tree; 
      Tree_struc_ tree_;
  std::string SimHitLabel;
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
DisplacedJetsAnalyzer::DisplacedJetsAnalyzer(const edm::ParameterSet& iConfig)
{

  //  genPartToken_   = consumes<std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparts"));    
  //genJetToken_   = consumes<std::vector<reco::GenJet> > (iConfig.getParameter<edm::InputTag>("genjets"));    
  jetToken_ = consumes<std::vector<reco::PFJet> >(iConfig.getUntrackedParameter<edm::InputTag>("jets"));
}


DisplacedJetsAnalyzer::~DisplacedJetsAnalyzer()
{
  
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DisplacedJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   /*   edm::Handle<std::vector<reco::GenParticle> > genParticlesH;
   iEvent.getByToken( genPartToken_, genParticlesH );

   edm::Handle<std::vector<reco::GenJet> > genJetsH;
   iEvent.getByToken(genJetToken_, genJetsH);
   */
   edm::Handle<std::vector<reco::PFJet> > jetsH;
   iEvent.getByToken(jetToken_, jetsH);


   int   run_   = iEvent.id().run();
   int lumi  = iEvent.luminosityBlock();
   int  event = iEvent.id().event();

   tree_.run = run_;
   tree_.lumi=lumi;
   tree_.event=event;

   //   int ngenjets=0;
   std::vector<float> genjet_e;
   std::vector<float> genjet_pt;
   std::vector<float> genjet_eta;
   std::vector<float> genjet_phi;
   std::vector<float> genjet_numdoug;
   std::vector<float> genjet_numconst;

   std::vector<int>   genjet_electrons;
   std::vector<int>   genjet_muons;
   std::vector<int>   genjet_taus;
   std::vector<int>   genjet_quarks;
   std::vector<int>   genjet_gluons;
   std::vector<int>   genjet_photons;
   std::vector<int>   genjet_pizeros;
   std::vector<int>   genjet_pipm;
   std::vector<int>   genjet_kzeros;
   std::vector<int>   genjet_kpm;
   std::vector<int>   genjet_protons;
   std::vector<int>   genjet_neutrons;
   std::vector<int>   genjet_etas;
   std::vector<int>   genjet_lambdas;
   std::vector<int>   genjet_chHad;


   int njets=0;
   std::vector<float> jet_e;
   std::vector<float> jet_pt;
   std::vector<float> jet_eta;
   std::vector<float> jet_phi;
   std::vector<float> jet_numdoug;
   std::vector<float> jet_numconst;

   /*
   //accessing genjets
   if (genJetsH.isValid()) // make sure gen particles exist
     {
       ngenjets = genJetsH->size();
       
       for (const auto & gjetiter : *genJetsH) // loop over gen jets
	 {

       int electrons=0;
       int muons=0;
       int taus=0;
       int gluons=0;
       int photons=0;
       int pizeros=0;
       int pipm=0;
       int kzeros=0;
       int kpm=0;
       int protons=0;
       int neutrons=0;
       int quarks=0;
       int etas=0;
       int lambdas=0;
       
       genjet_e.push_back(gjetiter.emEnergy());
       genjet_pt.push_back(gjetiter.pt());
       genjet_phi.push_back(gjetiter.phi());
       genjet_eta.push_back(gjetiter.eta());
       genjet_numdoug.push_back(gjetiter.numberOfDaughters());
       
       // Get the constituent particles
       std::vector<const reco::GenParticle*> particles = gjetiter.getGenConstituents();
       genjet_numconst.push_back(particles.size());
       
       //       std::cout<<" ----------- Look at its gen constituents particles: -----------------"<<std::endl;
       //investigsate nature of gen jet constituents
       for(unsigned int ip = 0; ip< particles.size(); ip++){
	 
	 int  genPdgId= particles[ip]->pdgId();
	 //	 float  genVtxX = particles[ip]->vertex().x();                                                                                                                       
	 //float genVtxY = particles[ip]->vertex().y();  
	 //float  genVtxZ = particles[ip]->vertex().z(); 
	 //std::cout<<genPdgId<< " "<<genVtxX<< " "<<genVtxY<<" "<<genVtxZ<<std::endl;    
	 if(abs(genPdgId) < 9) quarks++;
	 else if(abs(genPdgId) ==11) electrons++;
	 else if(abs(genPdgId) ==13) muons++;
	 else if(abs(genPdgId) ==15) taus++;
	 else if(abs(genPdgId) ==21) gluons++;
	 else if(abs(genPdgId) ==22) photons++;
	 else if(abs(genPdgId) ==111) pizeros++;
	 else if(abs(genPdgId) ==211) pipm++;
	 else if(abs(genPdgId) ==130 || abs(genPdgId) ==310) kzeros++;
	 else if(abs(genPdgId) ==321) kpm++;
	 else if(abs(genPdgId) ==2212) protons++;
	 else if(abs(genPdgId) ==2112) neutrons++;
	 else if(abs(genPdgId) ==221) etas++;
	 else if(abs(genPdgId) ==3122) lambdas++;
	 
   }


	   genjet_electrons.push_back(electrons);
	   genjet_muons.push_back(muons);
	   genjet_taus.push_back(taus);
	   genjet_quarks.push_back(quarks);
	   genjet_gluons.push_back(gluons);
	   genjet_photons.push_back(photons);
	   genjet_pizeros.push_back(pizeros);
	   genjet_pipm.push_back(pipm);
	   genjet_kzeros.push_back(kzeros);
	   genjet_kpm.push_back(kpm);
	   genjet_protons.push_back(protons);
	   genjet_neutrons.push_back(neutrons);
	   genjet_etas.push_back(etas);
	   genjet_lambdas.push_back(lambdas);

	   genjet_chHad.push_back(pipm+kpm+protons);


	 }

	   //check gensim hit inside the cone of the gen jet


	   */
       ///checking RECO jets
   if (jetsH.isValid()) // make sure gen particles exist                                                                                                                                               
     {
       njets = jetsH->size();
       
       for (const auto & jetiter : *jetsH) // loop over gen jets                                                                                                                                      
	 {
	   
	   if(fabs(jetiter.eta())>1.5 || jetiter.pt()<30 || jetiter.pt()>=1000)
		 {
		   jet_pt.push_back(jetiter.pt());
		   jet_eta.push_back(jetiter.eta());
		   jet_phi.push_back(jetiter.phi());
		   jet_e.push_back(jetiter.energy());
		   
		   
		 }
	 }
     }
   

//   tree_.ngenjets=ngenjets;   

/*
   for(int ij = 0; ij<ngenjets; ij++){
     tree_.genjet_e.push_back(genjet_e[ij]);
     tree_.genjet_pt.push_back(genjet_pt[ij]);
     tree_.genjet_eta.push_back(genjet_eta[ij]);
     tree_.genjet_phi.push_back(genjet_phi[ij]);
     tree_.genjet_numdoug.push_back(genjet_numdoug[ij]);
     tree_.genjet_numconst.push_back(genjet_numconst[ij]);


     tree_.genjet_electrons.push_back(genjet_electrons[ij]);
     tree_.genjet_muons.push_back(genjet_muons[ij]);
     tree_.genjet_taus.push_back(genjet_taus[ij]);
     tree_.genjet_quarks.push_back(genjet_quarks[ij]);
     tree_.genjet_gluons.push_back(genjet_gluons[ij]);
     tree_.genjet_photons.push_back(genjet_photons[ij]);
     tree_.genjet_pizeros.push_back(genjet_pizeros[ij]);
     tree_.genjet_pipm.push_back(genjet_pipm[ij]);
     tree_.genjet_kzeros.push_back(genjet_kzeros[ij]);
     tree_.genjet_kpm.push_back(genjet_kpm[ij]);
     tree_.genjet_protons.push_back(genjet_protons[ij]);
     tree_.genjet_neutrons.push_back(genjet_neutrons[ij]);
     tree_.genjet_etas.push_back(genjet_etas[ij]);
     tree_.genjet_lambdas.push_back(genjet_lambdas[ij]);
     tree_.genjet_chHad.push_back(genjet_chHad[ij]);

   }

*/
   tree_.njets=njets;   
   
   for(int ij = 0; ij<njets; ij++){
     tree_.jet_e.push_back(jet_e[ij]);
     tree_.jet_pt.push_back(jet_pt[ij]);
     tree_.jet_eta.push_back(jet_eta[ij]);
     tree_.jet_phi.push_back(jet_phi[ij]);
   }



      tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedJetsAnalyzer::beginJob()
{
  std::cout<<"begin job" <<std::endl;
  tree = fs->make<TTree>("tree","tree");
  
  tree->Branch("run", &tree_.run, "run/i");
  tree->Branch("lumi", &tree_.lumi, "lumi/i");
  tree->Branch("event", &tree_.event, "event/l");
  tree->Branch("ngenjets", &tree_.ngenjets, "ngenjets/i");
  tree->Branch("genjet_e",&tree_.genjet_e); 
  tree->Branch("genjet_pt",&tree_.genjet_pt); 
  tree->Branch("genjet_eta",&tree_.genjet_eta); 
  tree->Branch("genjet_phi",&tree_.genjet_phi); 
  tree->Branch("genjet_numdoug",&tree_.genjet_numdoug); 
  tree->Branch("genjet_numconst",&tree_.genjet_numconst); 
  
  
  tree->Branch("genjet_electrons", &tree_.genjet_electrons);
  tree->Branch("genjet_muons", &tree_.genjet_muons);	  
  tree->Branch("genjet_taus", &tree_.genjet_taus);	  
  tree->Branch("genjet_quarks", &tree_.genjet_quarks);
  tree->Branch("genjet_gluons", &tree_.genjet_gluons);	  
  tree->Branch("genjet_photons", &tree_.genjet_photons);  
  tree->Branch("genjet_pizeros", &tree_.genjet_pizeros);
  tree->Branch("genjet_pipm", &tree_.genjet_pipm);	  
  tree->Branch("genjet_kzeros", &tree_.genjet_kzeros);	  
  tree->Branch("genjet_kpm", &tree_.genjet_kpm);
  tree->Branch("genjet_protons", &tree_.genjet_protons);  
  tree->Branch("genjet_neutrons", &tree_.genjet_neutrons); 
  tree->Branch("genjet_etas", &tree_.genjet_etas); 
  tree->Branch("genjet_lambdas", &tree_.genjet_lambdas); 
  tree->Branch("genjet_chHad", &tree_.genjet_chHad);     
  

  tree->Branch("njets", &tree_.njets, "njets/i");
  tree->Branch("jet_e",&tree_.jet_e); 
  tree->Branch("jet_pt",&tree_.jet_pt); 
  tree->Branch("jet_eta",&tree_.jet_eta); 
  tree->Branch("jet_phi",&tree_.jet_phi); 
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedJetsAnalyzer::endJob() 
{
}

void DisplacedJetsAnalyzer::initTreeStructure() {
  tree_.run = 0.;

}


// ------------ method called when starting to processes a run  ------------
/*
void 
DisplacedJetsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DisplacedJetsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DisplacedJetsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DisplacedJetsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedJetsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedJetsAnalyzer);
