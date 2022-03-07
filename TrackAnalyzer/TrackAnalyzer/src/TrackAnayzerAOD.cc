/*
   Package:    TrackAnayzerAOD
   Class:      TrackAnayzerAOD

   Original Author: Adriano Di Florio

*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TTree.h"
#include "TVectorD.h"
#include <vector>
#include <sstream>


//
// class declaration
//

class TrackAnayzerAOD : public edm::EDAnalyzer {
   public:
      explicit TrackAnayzerAOD(const edm::ParameterSet&);
      ~TrackAnayzerAOD() override;

      bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      void beginRun(edm::Run const&, edm::EventSetup const&) override;
      void endRun(edm::Run const&, edm::EventSetup const&) override;
      void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::TrackCollection> TrackCollection_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;

  std::string TreeName_;

  UInt_t run, event, lumi, numPrimaryVertices, numTracks;

  TLorentzVector track_p4;

  UInt_t lowMuon_NPixelHits, lowMuon_NStripHits, lowMuon_NTrackhits, lowMuon_NBPixHits;
  UInt_t lowMuon_NPixLayers, lowMuon_NTraLayers, lowMuon_NStrLayers, lowMuon_NBPixLayers;

  Int_t theTrack_NPixelHits, theTrack_NStripHits, theTrack_NTrackhits, theTrack_NBPixHits, theTrack_NPixLayers;
  Int_t theTrack_NTraLayers, theTrack_NStrLayers, theTrack_NBPixLayers, lowTrack_NPixelHits, lowTrack_NStripHits;
  Double_t theTrackFromPVBS, theTrackFromPVDZ;
  Double_t theTrack_pt, theTrack_eta, theTrack_phi, theTrack_charge, theTrack_dz, theTrack_dxy;

  Double_t theTrack_ptErr, theTrack_etaErr, theTrack_phiErr;

  Int_t theTrack_nChi2, theTrack_covVersion, theTrack_covSchema;

  Double_t theTrack_SQopQop, theTrack_SQopLam, theTrack_SQopPhi, theTrack_SQopDxy, theTrack_SQopDsz, theTrack_SLamLam;
  Double_t theTrack_SLamPhi, theTrack_SLamDxy, theTrack_SLamDsz, theTrack_SPhiPhi, theTrack_SPhiDxy, theTrack_SPhiDsz;
  Double_t theTrack_SDxyDxy, theTrack_SDxyDsz, theTrack_SDszDsz;
  Double_t covMatrix[25];
  Double_t eigVector[5];
  TTree* track_tree;

};


//
// constructors and destructor
//
TrackAnayzerAOD::TrackAnayzerAOD(const edm::ParameterSet& iConfig):
        TrackCollection_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
        thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
        TreeName_(iConfig.getParameter<std::string>("TreeName"))
{
	      edm::Service<TFileService> fs;
        track_tree = fs->make<TTree>(TreeName_.data(),"Tree of DiMuon and Four Tracks");

        track_tree->Branch("run",                &run,                "run/I");
        track_tree->Branch("event",              &event,              "event/I");
        track_tree->Branch("lumi",              &lumi,              "lumi/I");
        track_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        track_tree->Branch("numTracks",                &numTracks,                "numTracks/I");

        track_tree->Branch("track_p4",   "TLorentzVector", &track_p4);

        track_tree->Branch("theTrack_pt",         &theTrack_pt,         "theTrack_pt/D");
        track_tree->Branch("theTrack_eta",        &theTrack_eta,        "theTrack_eta/D");
        track_tree->Branch("theTrack_phi",        &theTrack_phi,        "theTrack_phi/D");
        track_tree->Branch("theTrack_charge",     &theTrack_charge,     "theTrack_charge/D");
        track_tree->Branch("theTrack_dz",         &theTrack_dz,         "theTrack_dz/D");
        track_tree->Branch("theTrack_dxy",        &theTrack_dxy,        "theTrack_dxy/D");

        track_tree->Branch("theTrack_ptErr",        &theTrack_ptErr,        "theTrack_ptErr/D");
        track_tree->Branch("theTrack_etaErr",        &theTrack_etaErr,        "theTrack_etaErr/D");
        track_tree->Branch("theTrack_phiErr",        &theTrack_phiErr,        "theTrack_phiErr/D");

        track_tree->Branch("covMatrix",&covMatrix,"covMatrix[25]/D");
        track_tree->Branch("eigVector",&eigVector,"eigVector[5]/D");

        track_tree->Branch("theTrack_nChi2",        &theTrack_nChi2,        "theTrack_nChi2/I");
        // track_tree->Branch("theTrack_covVersion",        &theTrack_covVersion,        "theTrack_covVersion/I");
        // track_tree->Branch("theTrack_covSchema",        &theTrack_covSchema,        "theTrack_covSchema/I");

        track_tree->Branch("theTrack_SQopQop",        &theTrack_SQopQop,        "theTrack_SQopQop/D");
        track_tree->Branch("theTrack_SQopLam",        &theTrack_SQopLam,        "theTrack_SQopLam/D");
        track_tree->Branch("theTrack_SQopPhi",        &theTrack_SQopPhi,        "theTrack_SQopPhi/D");
        track_tree->Branch("theTrack_SQopDxy",        &theTrack_SQopDxy,        "theTrack_SQopDxy/D");
        track_tree->Branch("theTrack_SQopDsz",        &theTrack_SQopDsz,        "theTrack_SQopDsz/D");

        track_tree->Branch("theTrack_SLamLam",        &theTrack_SLamLam,        "theTrack_SLamLam/D");
        track_tree->Branch("theTrack_SLamPhi",        &theTrack_SLamPhi,        "theTrack_SLamPhi/D");
        track_tree->Branch("theTrack_SLamDxy",        &theTrack_SLamDxy,        "theTrack_SLamDxy/D");
        track_tree->Branch("theTrack_SLamDsz",        &theTrack_SLamDsz,        "theTrack_SLamDsz/D");
        track_tree->Branch("theTrack_SPhiPhi",        &theTrack_SPhiPhi,        "theTrack_SPhiPhi/D");

        track_tree->Branch("theTrack_SPhiDxy",        &theTrack_SPhiDxy,        "theTrack_SPhiDxy/D");
        track_tree->Branch("theTrack_SPhiDsz",        &theTrack_SPhiDsz,        "theTrack_SPhiDsz/D");
        track_tree->Branch("theTrack_SDxyDxy",        &theTrack_SDxyDxy,        "theTrack_SDxyDxy/D");
        track_tree->Branch("theTrack_SDxyDsz",        &theTrack_SDxyDsz,        "theTrack_SDxyDsz/D");
        track_tree->Branch("theTrack_SDszDsz",        &theTrack_SDszDsz,        "theTrack_SDszDsz/D");

        //Tracks Flags

        track_tree->Branch("theTrack_NPixelHits",        &theTrack_NPixelHits,        "theTrack_NPixelHits/I");
        track_tree->Branch("theTrack_NStripHits",        &theTrack_NStripHits,        "theTrack_NStripHits/I");
        track_tree->Branch("theTrack_NTrackhits",        &theTrack_NTrackhits,        "theTrack_NTrackhits/I");
        track_tree->Branch("theTrack_NBPixHits",         &theTrack_NBPixHits,        "theTrack_NBPixHits/I");

        track_tree->Branch("theTrack_NPixLayers",        &theTrack_NPixLayers,        "theTrack_NPixLayers/I");
        track_tree->Branch("theTrack_NTraLayers",        &theTrack_NTraLayers,        "theTrack_NTraLayers/I");
        track_tree->Branch("theTrack_NStrLayers",        &theTrack_NStrLayers,        "theTrack_NStrLayers/I");
        track_tree->Branch("theTrack_NBPixLayers",       &theTrack_NBPixLayers,        "theTrack_NBPixLayers/I");


}

TrackAnayzerAOD::~TrackAnayzerAOD() {}

//
// member functions
//

bool TrackAnayzerAOD::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void TrackAnayzerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<reco::TrackCollection> track;
  iEvent.getByToken(TrackCollection_,track);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(thePVs_, primaryVertices_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run   = iEvent.id().run();
  event = iEvent.id().event();
  lumi  = iEvent.id().luminosityBlock();

  track_p4.SetPtEtaPhiM(0.,0.,0.,0.);

//std::cout << "Debug  1" << std::endl;

  if (!track.isValid()) std::cout<< "No tracks information " << run << "," << event <<std::endl;

  if (track.isValid()) {


    //std::cout << "Debug  2" << std::endl;
    numTracks = (UInt_t)(track->size());

    for (unsigned int i=0; i< track->size(); i++)
    {

      // std::cout << "Event: " << run << " - "
      //                         << event << " - "
      //                         << lumi << " - "
      //           << "N Tracks: " << numTracks << std::endl;

      auto tk = track->at(i);

      if (tk.pt() < 0.6) continue;

      theTrack_NPixelHits  = -1;
      theTrack_NStripHits  = -1;
      theTrack_NTrackhits  = -1;
      theTrack_NBPixHits   = -1;
      theTrack_NPixLayers  = -1;
      theTrack_NTraLayers  = -1;
      theTrack_NStrLayers  = -1;
      theTrack_NBPixLayers = -1;

      theTrack_pt      = tk.pt();
      theTrack_eta     = tk.eta();
      theTrack_phi     = tk.phi();
      theTrack_charge  = tk.charge();


      reco::TrackBase::CovarianceMatrix cov = tk.covariance();

      TMatrixDSym COV(cov.kRows);
      for (int j = 0; j < cov.kRows; j++)
          for (int k = 0; k < cov.kRows; k++)
          {
              covMatrix[j*5+k] = cov(j,k);
              COV(j,k) = cov(j,k);
            }
      TVectorD eig(cov.kRows);
      COV.EigenVectors(eig);

      for (size_t i = 0; i < 5; i++) {
        eigVector[i] = eig[i];
      }

      theTrack_dz      = tk.dz();
      theTrack_dxy     = tk.dxy();

      theTrack_ptErr  = tk.ptError();
      theTrack_etaErr = tk.etaError();
      theTrack_phiErr = tk.phiError();

      theTrack_nChi2 = tk.normalizedChi2();;
      // theTrack_covVersion = tk.covarianceVersion();;
      // theTrack_covSchema = tk.covarianceSchema();;

      theTrack_NPixelHits  = tk.hitPattern().numberOfValidPixelHits();
      theTrack_NStripHits  = tk.hitPattern().numberOfValidStripHits();
      theTrack_NTrackhits  = tk.hitPattern().numberOfValidTrackerHits();
      theTrack_NBPixHits   = tk.hitPattern().numberOfValidStripHits();
      theTrack_NPixLayers  = tk.hitPattern().pixelLayersWithMeasurement();
      theTrack_NTraLayers  = tk.hitPattern().trackerLayersWithMeasurement();
      theTrack_NStrLayers  = tk.hitPattern().stripLayersWithMeasurement();
      theTrack_NBPixLayers = tk.hitPattern().pixelBarrelLayersWithMeasurement();


      reco::TrackBase::CovarianceMatrix cm = tk.covariance();


      theTrack_SQopQop = cm( reco::TrackBase::i_qoverp, reco::TrackBase::i_qoverp );
      theTrack_SQopLam = cm( reco::TrackBase::i_qoverp, reco::TrackBase::i_lambda );
      theTrack_SQopPhi = cm( reco::TrackBase::i_qoverp, reco::TrackBase::i_phi    );
      theTrack_SQopDxy = cm( reco::TrackBase::i_qoverp, reco::TrackBase::i_dxy    );
      theTrack_SQopDsz = cm( reco::TrackBase::i_qoverp, reco::TrackBase::i_dsz    );
      theTrack_SLamLam = cm( reco::TrackBase::i_lambda, reco::TrackBase::i_lambda );
      theTrack_SLamPhi = cm( reco::TrackBase::i_lambda, reco::TrackBase::i_phi    );
      theTrack_SLamDxy = cm( reco::TrackBase::i_lambda, reco::TrackBase::i_dxy    );
      theTrack_SLamDsz = cm( reco::TrackBase::i_lambda, reco::TrackBase::i_dsz    );
      theTrack_SPhiPhi = cm( reco::TrackBase::i_phi   , reco::TrackBase::i_phi    );
      theTrack_SPhiDxy = cm( reco::TrackBase::i_phi   , reco::TrackBase::i_dxy    );
      theTrack_SPhiDsz = cm( reco::TrackBase::i_phi   , reco::TrackBase::i_dsz    );
      theTrack_SDxyDxy = cm( reco::TrackBase::i_dxy   , reco::TrackBase::i_dxy    );
      theTrack_SDxyDsz = cm( reco::TrackBase::i_dxy   , reco::TrackBase::i_dsz    );
      theTrack_SDszDsz = cm( reco::TrackBase::i_dsz   , reco::TrackBase::i_dsz    );



      track_tree->Fill();

      // cand.unpackCovariance();




      } //IsMC || onlyGen




        // dimuontt candidates are sorted by vProb
    }

}

// ------------ method called once each job just before starting event loop  ------------
void TrackAnayzerAOD::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void TrackAnayzerAOD::endJob() {}

// ------------ method called when starting to processes a run  ------------
void TrackAnayzerAOD::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void TrackAnayzerAOD::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void TrackAnayzerAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void TrackAnayzerAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackAnayzerAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnayzerAOD);
