#include "WZEdmAnalyzer.h"
#include "EdmAnalysisTools.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/IPTools/interface/IPTools.h"




/**************************************************************************
 *
 * fill necessry event information for this analysis
 * work in progress
 *
 **************************************************************************/
using namespace muon;


void 
WZEdmAnalyzer::copyVertexInfo(const reco::VertexCollection::const_iterator vertex,  const reco::MuonCollection *muons, int leadingMuonIndex, _vertex_ *myVertex) 
{


  const reco::Muon *muon =0;
  if (leadingMuonIndex>=0) muon = &(recoMuons->at(leadingMuonIndex));


  // select the leading muon

  myVertex->isValid            = vertex->isValid();
  myVertex->isFake             = vertex->isFake();

  myVertex->x                  = vertex->x();
  myVertex->y                  = vertex->y();
  myVertex->z                  = vertex->z();

  myVertex->xError             = vertex->xError();
  myVertex->yError             = vertex->yError();
  myVertex->zError             = vertex->zError();


  myVertex->chi2               = vertex->chi2();
  myVertex->normalizedChi2     = vertex->normalizedChi2();
  myVertex->ndof               = vertex->ndof();

  myVertex->tracksSize         = vertex->tracksSize();
  myVertex->muonIndex          = leadingMuonIndex;

  int num;
  float sumpt, sumpt2;

  float thrust, sphericity, planarity, aplanarity;

  // 0 threshold
  vertexActivites(vertex, muon, 0., num, sumpt, sumpt2,  thrust, sphericity, planarity, aplanarity);
  myVertex->validTracksSize = num;
  myVertex->sumPtMin = sumpt;
  myVertex->sumPt2Min = sumpt2;

  myVertex->thrust = thrust;
  myVertex->sphericity = sphericity;
  myVertex->planarity = planarity;
  myVertex->aplanarity = aplanarity;

  //  std::cout << "muon = " << leadingMuonIndex << "; t, s, p, a = "
  //    <<  thrust <<", "
  //    <<  sphericity << ", "
  //    <<  planarity << ", "
  //    <<  aplanarity
  //    << std::endl;



  //  std::cout << " num = " << num << " sumpt = " << sumpt << " sumpt2 " << sumpt2 << std::endl;

  vertexActivites(vertex, muon, 0.2, num, sumpt, sumpt2,  thrust, sphericity, planarity, aplanarity);
  myVertex->validTracksSize200MeV = num;
  myVertex->sumPtMin200MeV  = sumpt;
  myVertex->sumPt2Min200MeV  = sumpt2;
  //  std::cout << " num = " << num << " sumpt = " << sumpt << " sumpt2 " << sumpt2 << std::endl;

  vertexActivites(vertex, muon, 0.5, num, sumpt, sumpt2,  thrust, sphericity, planarity, aplanarity);
  myVertex->validTracksSize500MeV  = num;
  myVertex->sumPtMin500MeV  = sumpt;
  myVertex->sumPt2Min500MeV  = sumpt2;
  //  std::cout << " num = " << num << " sumpt = " << sumpt << " sumpt2 " << sumpt2 << std::endl;

  vertexActivites(vertex, muon, 1.0, num, sumpt, sumpt2,  thrust, sphericity, planarity, aplanarity);
  myVertex->validTracksSize1GeV  = num;
  myVertex->sumPtMin1GeV  = sumpt;
  myVertex->sumPt2Min1GeV  = sumpt2;

  //  std::cout << " num = " << num << " sumpt = " << sumpt << " sumpt2 " << sumpt2 << std::endl;


}


void 
WZEdmAnalyzer::copyTrackInfo(const reco::Track *fwTrack, _track_ *myTrack) 
{
  
  //  MC information if it's available
  //  only doing track simutrack matching with a given minimum pt
  if (!(_is_data)  
      && fwTrack->pt() > trackMinPtWithMCTruth && _reco) { 


    const SimTrack *aSimTrack  = closestSimTrack(*fwTrack, recoSimTracks);
    if (aSimTrack) {
      myTrack->mcType            = aSimTrack->type();
      myTrack->mcPt              = aSimTrack->momentum().Pt();
      myTrack->mcEta             = aSimTrack->momentum().Eta();
      myTrack->mcPhi             = aSimTrack->momentum().Phi();
    }
  }
  

  if (_is_debug) {
    std::cout << "(px, py, pz) = (" << fwTrack->px() << ", "
	      <<fwTrack->py() << ", "
	      <<fwTrack->pz() << ")"
	      <<std::endl;
  }


  //copy track information


  myTrack->chi2                = fwTrack->chi2();
  myTrack->normalizedChi2      = fwTrack->normalizedChi2();
  myTrack->charge              = (Int_t) fwTrack->charge();
  myTrack->ndof                = fwTrack->ndof();

  myTrack->pt                  = fwTrack->pt();
  myTrack->ptError             = fwTrack->ptError();
  myTrack->eta                 = fwTrack->eta();
  myTrack->etaError            = fwTrack->etaError();
  myTrack->phi                 = fwTrack->phi();
  myTrack->phiError            = fwTrack->phiError();

  myTrack->qoverp              = fwTrack->qoverp();
  myTrack->qoverpError         = fwTrack->qoverpError();
  myTrack->lambda              = fwTrack->lambda();
  myTrack->lambdaError         = fwTrack->lambdaError();

  myTrack->dxy0                = fwTrack->dxy();
  myTrack->dxy0Error           = fwTrack->dxyError();
  myTrack->dz0                 = fwTrack->dz();
  myTrack->dz0Error            = fwTrack->dzError();
  myTrack->dsz0                = fwTrack->dsz();
  myTrack->dsz0Error           = fwTrack->dszError();


  // dxy with respecting to the beam spot
  if (vertexBeamSpot) {
    myTrack->dxy                 = fwTrack->dxy(vertexBeamSpot->position());
    myTrack->dz                  = fwTrack->dz(vertexBeamSpot->position());
    myTrack->dsz                 = fwTrack->dsz(vertexBeamSpot->position());
    
  }

  if (recVtxs->size() > 0) {
    reco::VertexRef vtx(recVtxs, 0);    
    myTrack->dxyVtx = fwTrack->dxy(vtx->position());
    myTrack->dzVtx = fwTrack->dz(vtx->position());
  } else {

    myTrack->dxyVtx = fwTrack->dxy();
    myTrack->dzVtx = fwTrack->dz();
  }



  
  if (_reco)  {
    myTrack->recHitsSize         = fwTrack->recHitsSize();
  
    myTrack->outerOk             = fwTrack->outerOk();
  }

  myTrack->numberOfValidHits   = fwTrack->numberOfValidHits();




  //  myTrack->outerRadius         = fwTrack->outerRadius();
  // myTrack->innerRadius         = sqrt( pow(fwTrack->innerPosition().X(), 2) + pow(fwTrack->innerPosition().Y(), 2) );


  if (_is_debug)   std::cout << "track dz error --- " << setw(10) << fwTrack->dzError() << std::endl;



  myTrack->nStereoHits = 0;
  if (_reco) {
    for (trackingRecHit_iterator recHit = fwTrack->recHitsBegin(); 
	 recHit != fwTrack->recHitsEnd(); 
	 recHit ++) {
      
      if ( (*recHit)->type() != 0) continue;
      
	   
      //bool 	  isStereoHit = false;
      if ((( SiStripDetId )(*recHit)->geographicalId()).stereo()) {
	
	myTrack->nStereoHits ++;
	//isStereoHit = true;
      }
        

      if (_is_debug)  std::cout << (( SiStripDetId )(*recHit)->geographicalId()).stereo() << std::endl;
    }
  }

  // some debug information
  if (_is_debug) {
    if (myTrack->dz0Error < 4 && myTrack->numberOfValidHits > 8 && myTrack->dz0Error > 0.5 && myTrack->nStereoHits >=4) {

      std::cout << "run Number " << runNum << "; event Number " << evtNum << std::endl;
    }
    
    std::cout << endl;   
  }     



  myTrack->numberOfValidPixelHits           = fwTrack->hitPattern().numberOfValidPixelHits();
  myTrack->numberOfValidStripHits           = fwTrack->hitPattern().numberOfValidStripHits();
  myTrack->numberOfValidPixelBarrelHits     = fwTrack->hitPattern().numberOfValidPixelBarrelHits();
  
  myTrack->pixelLayersWithMeasurement       = fwTrack->hitPattern().pixelLayersWithMeasurement();
  myTrack->pixelBarrelLayersWithMeasurement = fwTrack->hitPattern().pixelBarrelLayersWithMeasurement() ;
  myTrack->pixelEndcapLayersWithMeasurement = fwTrack->hitPattern().pixelEndcapLayersWithMeasurement();

  myTrack->numberOfValidTrackerLayers       = fwTrack->hitPattern().trackerLayersWithMeasurement();
  myTrack->trackerLayersWithMeasurement     = fwTrack->hitPattern().trackerLayersWithMeasurement();


  // hits pattern for conversion rejection
  myTrack->hasValidHitInFirstPixelBarrel    = fwTrack->hitPattern().hasValidHitInFirstPixelBarrel();
 
  // NOTE74: this is not available in 7_4_x release
  // myTrack->trackerExpectedHitsInner_numberOfHits    = fwTrack->trackerExpectedHitsInner().numberOfHits();



  if (_is_debug) {std::cout << "finishing track information copy " << std::endl;}

 
}


void  
WZEdmAnalyzer::copyMuonInfo(  reco::MuonCollection::const_iterator fwMuon,  _muon_  *myMuon) 
{
  
  // MC information if it's available
  if (!_is_data 
      && (fwMuon->isGlobalMuon() || fwMuon->isTrackerMuon()) && _reco) {

    const SimTrack *mc_trk        = closestSimTrack(*fwMuon->track(), recoSimTracks);
    if (mc_trk) {
      myMuon->mcType              = mc_trk->type();
      myMuon->mcPt                = mc_trk->momentum().Pt();
      myMuon->mcEta               = mc_trk->momentum().Eta();
      myMuon->mcPhi               = mc_trk->momentum().Phi();

    }
  }


  myMuon->muonType                = fwMuon->type();


  //copy inner track information 
  if (fwMuon->isGlobalMuon() || fwMuon->isTrackerMuon()) {

    const reco::Track * aTrack = &(*fwMuon->innerTrack());
    this->copyTrackInfo(aTrack, (_track_ *)myMuon); 

    // access the closest primary vertex
    int vertexIndex               = whichVertex((const reco::TrackBaseRef)fwMuon->innerTrack(), recVtxs, myMuon->closestVertexDist);
    myMuon->vertexIndex           = vertexIndex;
    myMuon->associatedVertex = trackVertexAssociation( (const reco::TrackBaseRef)fwMuon->innerTrack(), recVtxs );


    //    std::cout << "accociated vertex " <<  myMuon->associatedVertex 
    //      << " ; " << myMuon->vertexIndex  << std::endl;
    // default dxy with respect to the closest primary vertex.
    //   if (vertexIndex >=0) {
    //
    //
    // myMuon->dxyVtx               = (*(fwMuon->track())).dxy( (*recVtxs.product())[vertexIndex].position());
    // myMuon->dzVtx                = (*(fwMuon->track())).dz( (*recVtxs.product())[vertexIndex].position());
    // myMuon->dszVtx               = (*(fwMuon->track())).dsz( (*recVtxs.product())[vertexIndex].position());
    // }



   if (fwMuon->isGlobalMuon()) {

     // TeV Muon ReFit
     reco::TrackRef pmcTrack;
     //pmcTrack = (reco::TrackRef)muon::tevOptimized(*fwMuon, tevMap1, tevMap2, tevMap3);//FOR CMSSW VERSION < 5_X_X
     // pmcTrack = muon::tevOptimized(*fwMuon, tevMap1, tevMap2, tevMap3).first;//FOR VERSIN GREATER THAN 5_X_X

     // NOTE74: for 74x release
     //     pmcTrack = muon::tevOptimized(*fwMuon).first;
     // tevMap1, tevMap2, tevMap3).first;//FOR VERSIN GREATER THAN 5_X_X
     // add in Oct. 27, 2015
     pmcTrack = fwMuon->tunePMuonBestTrack();

     myMuon->refitCharge     = (Int_t) pmcTrack->charge();
     myMuon->refitPt         = pmcTrack->pt();
     myMuon->refitEta        = pmcTrack->eta();
     myMuon->refitPhi        = pmcTrack->phi();

     myMuon->refitPtError         = pmcTrack->ptError();
     myMuon->refitEtaError        = pmcTrack->etaError();
     myMuon->refitPhiError        = pmcTrack->phiError();


     myMuon->isHighPtMuon         =
       muon::isHighPtMuon( *fwMuon,  (*recVtxs.product())[0]);
    
     myMuon->refitDxy             =
       pmcTrack->dxy( (*recVtxs.product())[0].position() );
     myMuon->refitDz              = 
       pmcTrack->dz( (*recVtxs.product())[0].position() );



    myMuon->globalCharge          = fwMuon->charge();
    myMuon->globalPt              = fwMuon->pt();
    myMuon->globalEta             = fwMuon->eta();
    myMuon->globalPhi             = fwMuon->phi();
    //    myMuon->globalPtError         = fwMuon->ptError();
    //   myMuon->globalEtaError        = fwMuon->etaError();
    //  myMuon->globalPhiError        = fwMuon->phiError();


    //    myMuon->globalCharge          = (*(fwMuon->globalTrack())).charge();
    //  myMuon->globalPt              = (*(fwMuon->globalTrack())).pt();
    // myMuon->globalEta             = (*(fwMuon->globalTrack())).eta();
    // myMuon->globalPhi             = (*(fwMuon->globalTrack())).phi();
    //  myMuon->globalPtError         = (*(fwMuon->globalTrack())).ptError();
    // myMuon->globalEtaError        = (*(fwMuon->globalTrack())).etaError();
    //  myMuon->globalPhiError        = (*(fwMuon->globalTrack())).phiError();

    myMuon->globalRecHitsSize     = (*(fwMuon->globalTrack())).recHitsSize();
    myMuon->globalNumberOfValidHits= (*(fwMuon->globalTrack())).numberOfValidHits();
    myMuon->globalNumberOfValidMuonHits = (*(fwMuon->globalTrack())).hitPattern().numberOfValidMuonHits();
    myMuon->globalNormalizedChi2  = (*(fwMuon->globalTrack())).normalizedChi2();


    
    myMuon->globalDxy0            = (*(fwMuon->globalTrack())).dxy();
    myMuon->globalDsz0            = (*(fwMuon->globalTrack())).dsz();
    myMuon->globalDz0             = (*(fwMuon->globalTrack())).dz();
    myMuon->globalDxy0Error       = (*(fwMuon->globalTrack())).dxyError();
    myMuon->globalDsz0Error       = (*(fwMuon->globalTrack())).dszError();
    myMuon->globalDz0Error        = (*(fwMuon->globalTrack())).dzError();




    myMuon->globalDxy             = (*(fwMuon->globalTrack())).dxy(vertexBeamSpot->position());
    myMuon->globalDsz             = (*(fwMuon->globalTrack())).dsz(vertexBeamSpot->position());
    myMuon->globalDz              = (*(fwMuon->globalTrack())).dz(vertexBeamSpot->position());


    myMuon->globalDxyVtx          = (*(fwMuon->globalTrack())).dxy( (*recVtxs.product())[0].position());
    myMuon->globalDzVtx           = (*(fwMuon->globalTrack())).dz( (*recVtxs.product())[0].position());
    myMuon->globalDszVtx          = (*(fwMuon->globalTrack())).dsz( (*recVtxs.product())[0].position());
   }
  }
      

  // standalone muon information
  if (fwMuon->isGlobalMuon() || fwMuon->isStandAloneMuon() ) {
      
    myMuon->outerRecHitsSize      = (*(fwMuon->standAloneMuon())).recHitsSize(); 
    myMuon->outerNumberOfValidHits= (*(fwMuon->standAloneMuon())).numberOfValidHits(); 
    myMuon->outerInnerOk          = (*(fwMuon->standAloneMuon())).innerOk();      
    myMuon->outerPt               = (*(fwMuon->standAloneMuon())).pt();      
    myMuon->outerEta              = (*(fwMuon->standAloneMuon())).eta();      
    myMuon->outerPhi              = (*(fwMuon->standAloneMuon())).phi();      

    myMuon->outerPtError          = (*(fwMuon->standAloneMuon())).ptError();
    myMuon->outerEtaError         = (*(fwMuon->standAloneMuon())).etaError();
    myMuon->outerPhiError         = (*(fwMuon->standAloneMuon())).phiError();
    myMuon->outerCharge           = (*(fwMuon->standAloneMuon())).charge();      
  }



  // access the muon recHits information
  if (fwMuon->isGlobalMuon()) {

    if (_is_debug) {
	  
      std::cout << "** " << std::endl;
      for (int ii = 0; ii < fwMuon->standAloneMuon()->found(); ii ++) {
	
	std::cout << "(X, Y) = (" 
		  << fwMuon->standAloneMuon()->residualX(ii)
		  <<", "
		  << fwMuon->standAloneMuon()->residualY(ii)
		  << ")" << std::endl;
      }
    }
  }
  


  myMuon->numberOfChambers       = fwMuon->numberOfChambers();
  myMuon->numberOfMatches        = fwMuon->numberOfMatches();
  myMuon->numberOfMatchedStations= fwMuon->numberOfMatchedStations();
  myMuon->stationMask            = fwMuon->stationMask();
  myMuon->stationGapMaskPull     = fwMuon->stationGapMaskPull();



  // muon timing information in muon chamber
  myMuon->isTimeValid            = fwMuon->isTimeValid();      
  myMuon->timeAtIPInOut          = fwMuon->time().timeAtIpInOut;
  myMuon->timeAtIPInOutErr       = fwMuon->time().timeAtIpInOutErr;
  myMuon->timeAtIPOutIn          = fwMuon->time().timeAtIpOutIn;
  myMuon->timeAtIPOutInErr       = fwMuon->time().timeAtIpOutInErr;
  myMuon->timenDof               = fwMuon->time().nDof;


  // compatibility
  myMuon->isCaloCompatibilityValid= fwMuon->isCaloCompatibilityValid();
  myMuon->caloCompatibility       = fwMuon->caloCompatibility();
      


  // isolation energy
  myMuon->isIsolationValid       = fwMuon->isIsolationValid();
  myMuon->isolationR03sumPt      = fwMuon->isolationR03().sumPt;
  myMuon->isolationR03emEt       = fwMuon->isolationR03().emEt;
  myMuon->isolationR03hadEt      = fwMuon->isolationR03().hadEt;
  myMuon->isolationR03hoEt       = fwMuon->isolationR03().hoEt;
  myMuon->isolationR03nTracks    = fwMuon->isolationR03().nTracks;
  myMuon->isolationR03nJets      = fwMuon->isolationR03().nJets;
  myMuon->isolationR03emVetoEt   = fwMuon->isolationR03().emVetoEt;
  myMuon->isolationR03hadVetoEt  = fwMuon->isolationR03().hadVetoEt;
  myMuon->isolationR03hoVetoEt   = fwMuon->isolationR03().hoVetoEt;


  myMuon->isolationR05sumPt      = fwMuon->isolationR05().sumPt;
  myMuon->isolationR05emEt       = fwMuon->isolationR05().emEt;
  myMuon->isolationR05hadEt      = fwMuon->isolationR05().hadEt;
  myMuon->isolationR05hoEt       = fwMuon->isolationR05().hoEt;
  myMuon->isolationR05nTracks    = fwMuon->isolationR05().nTracks;
  myMuon->isolationR05nJets      = fwMuon->isolationR05().nJets;
  myMuon->isolationR05emVetoEt   = fwMuon->isolationR05().emVetoEt;
  myMuon->isolationR05hadVetoEt  = fwMuon->isolationR05().hadVetoEt;
  myMuon->isolationR05hoVetoEt   = fwMuon->isolationR05().hoVetoEt;


  //
  // calorimeter energy deposit
  //
  myMuon->isEnergyValid          = fwMuon->isEnergyValid();
  myMuon->calEnergyem            = fwMuon->calEnergy().em;
  myMuon->calEnergyhad           = fwMuon->calEnergy().had;
  myMuon->calEnergyho            = fwMuon->calEnergy().ho;
  myMuon->calEnergy              = (fwMuon->calEnergy().em + fwMuon->calEnergy().had);

  myMuon->calEnergyecal_time     = fwMuon->calEnergy().ecal_time;
  myMuon->calEnergyecal_timeError= fwMuon->calEnergy().ecal_timeError;
  myMuon->calEnergyhcal_time     = fwMuon->calEnergy().hcal_time;
  myMuon->calEnergyhcal_timeError= fwMuon->calEnergy().hcal_timeError; 


  // additional muon selector information: implemented Nov. 25, 2013
  myMuon->improvedMuonBestTrackPt= 0;
  myMuon->improvedMuonBestTrackPtError =0;
  myMuon->improvedMuonBestTrackEta =0;
  myMuon->improvedMuonBestTrackPhi =0;
  myMuon->improvedMuonBestTrackCharge =0;
  // myMuon->isHighPtMuon = isHighPtMuon(*fwMuon, *(*recVtxs).begin());


  myMuon->isPFMuon = fwMuon->isPFMuon();
  myMuon->isSoftMuon = isSoftMuon(*fwMuon, *(*recVtxs).begin());
  myMuon->isGlobalMuon = fwMuon->isGlobalMuon();
  myMuon->isTightMuon =  isTightMuon( *fwMuon, *(*recVtxs).begin() );;
  myMuon->bestMuonBestTrackPt = (*(fwMuon->muonBestTrack())).pt();
  myMuon->bestMuonBestTrackPtError =(*(fwMuon->muonBestTrack())).ptError();
  myMuon->bestMuonBestTrackEta =(*(fwMuon->muonBestTrack())).eta();
  myMuon->bestMuonBestTrackPhi =(*(fwMuon->muonBestTrack())).phi();
  myMuon->bestMuonBestTrackCharge =(*(fwMuon->muonBestTrack())).charge();
  
  myMuon->bestMuonBestTrackDxy  = (*(fwMuon->muonBestTrack())).dxy( (*recVtxs.product())[0].position() );
 myMuon->bestMuonBestTrackDz  = (*(fwMuon->muonBestTrack())).dz( (*recVtxs.product())[0].position() );

  const MuonPFIsolation pfIso = fwMuon->pfIsolationR04();
  
  Float_t pfIso_beta = pfIso.sumChargedHadronPt;
  pfIso_beta += std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt);
  
  myMuon->pfIsolationR04beta = pfIso_beta;
  myMuon->pfIsolationR04     = pfIso.sumChargedHadronPt + pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt;

}



void
WZEdmAnalyzer::copyCaloJetInfo(    const edm::Event& iEvent,
				   const edm::EventSetup& iSetup, 
				   reco::CaloJetCollection::const_iterator jet,
				   edm::RefToBase<reco::Jet> &jetRef, 
				   JetCorrectionUncertainty *jetCorUnc, 
				   double scale, 
				   edm::Handle<reco::JetFlavourMatchingCollection> & theRecoTag, 
				   //		   edm::Handle<reco::JetTagCollection> & jetTags, 
				   _jet_ *myJet) {

  myJet->L2L3pt = scale * jet->pt();
  myJet->L2L3scale = scale;
  myJet->energy = scale * jet->energy();
  myJet->mass  = jet->mass();
  myJet->pt  = jet->pt();
  myJet->eta = jet->eta();
  myJet->phi = jet->phi();

  if (jetCorUnc ) {
    jetCorUnc->setJetEta(jet->eta());
    jetCorUnc->setJetPt(jet->pt() * scale); // here you must use the CORRECTED jet pt
    myJet->scaleUnc = jetCorUnc->getUncertainty(true);
  } else {
    myJet->scaleUnc = 0.0;

  }

  myJet->n60 = jet->n60();
  myJet->n90 = jet->n90();
  myJet->emEnergyInEE                = jet->emEnergyInEE();
  myJet->emEnergyInHF                = jet->emEnergyInHF();
  myJet->emEnergyInEB                = jet->emEnergyInEB();
  myJet->hadEnergyInHB               = jet->hadEnergyInHB();
  myJet->hadEnergyInHO               = jet->hadEnergyInHO();
  myJet->hadEnergyInHF               = jet->hadEnergyInHF();
  myJet->hadEnergyInHE               = jet->hadEnergyInHE();
  myJet->emEnergyFraction            = jet->emEnergyFraction();
  myJet->energyFractionHadronic      = jet->energyFractionHadronic();
  myJet->maxEInEmTowers              = jet->maxEInEmTowers();
  myJet->maxEInHadTowers             = jet->maxEInHadTowers();
  myJet->nConstituent                = jet->getCaloConstituents().size();



  // access via jet ID variables in the main jet loop.
  myJet->fHPD = -9999;
  myJet->fRBX = -9999;

  // access btagging for calo jet
  myJet->flavor[0]                   = btaggingAssociation( (reco::Jet)(*jet), jetTags.product());
  myJet->flavor[1]                   = btaggingAssociation( (reco::Jet)(*jet), jetTagsCSV.product());

  myJet->flavor[2]                   = btaggingAssociation( jetRef, myCaloJetTagsJP);
  myJet->flavor[3]                   = btaggingAssociation( jetRef, myCaloJetTagsTCHP);
  myJet->flavor[4]                   = btaggingAssociation( jetRef, myCaloJetTagsCSV);
  
  if (!_is_data) {

    myJet->mc_flavor = mcflavorAssociation( jetRef, theRecoTag);
  }
}




void
WZEdmAnalyzer::copyJPTJetInfo(      const edm::Event& iEvent,
				    const edm::EventSetup& iSetup, 
				    reco::JPTJetCollection::const_iterator jet,
				    edm::RefToBase<reco::Jet> &jetRef, 
				    JetCorrectionUncertainty *jetCorUnc,  
				    double scale, 
				    edm::Handle<reco::JetFlavourMatchingCollection> & theRecoTag, 
				    //  edm::Handle<reco::JetTagCollection> & jetTags, 
				    _jet_ *myJet) {

  myJet->L2L3pt    = scale * jet->pt();
  myJet->L2L3scale = scale;
  myJet->energy    = scale * jet->energy();
  myJet->mass      = jet->mass();
  myJet->pt        = jet->pt();
  myJet->eta       = jet->eta();
  myJet->phi       = jet->phi();


  if (jetCorUnc ) {
    jetCorUnc->setJetEta(jet->eta());
    jetCorUnc->setJetPt(jet->pt() * scale); // here you must use the CORRECTED jet pt
    myJet->scaleUnc = jetCorUnc->getUncertainty(true);
  } else {
    myJet->scaleUnc = 0.0;

  }


  myJet->chargedHadronEnergy         = jet->chargedHadronEnergy();
  myJet->neutralHadronEnergy         = jet->neutralHadronEnergy();
  myJet->chargedHadronEnergyFraction = jet->chargedHadronEnergyFraction();
  myJet->neutralHadronEnergyFraction = jet->neutralHadronEnergyFraction();
  myJet->chargedEmEnergy             = jet->chargedEmEnergy();
  myJet->chargedEmEnergyFraction     = jet->chargedEmEnergyFraction();
  myJet->neutralEmEnergyFraction     = jet->neutralEmEnergyFraction();
  myJet->neutralEmEnergy             = jet->neutralEmEnergy();
  myJet->neutralEmEnergyFraction     = jet->neutralEmEnergyFraction();
  myJet->muonMultiplicity            = jet->muonMultiplicity();
  myJet->electronMultiplicity        = jet->elecMultiplicity();
  myJet->chargedMultiplicity         = jet->chargedMultiplicity();
  myJet->zspCor                      = jet->getZSPCor();
  myJet->nConstituent                = jet->nConstituents();

 // access btagging for calo jet
  myJet->flavor[0]                   = btaggingAssociation( (reco::Jet)(*jet), jetTags.product());
  myJet->flavor[1]                   = btaggingAssociation( (reco::Jet)(*jet), jetTagsCSV.product());

  myJet->flavor[2]                   = btaggingAssociation( jetRef, myJPTJetTagsJP);
  myJet->flavor[3]                   = btaggingAssociation( jetRef, myJPTJetTagsTCHP);
  myJet->flavor[4]                   = btaggingAssociation( jetRef, myJPTJetTagsCSV);
  
  if (!_is_data) {

    myJet->mc_flavor = mcflavorAssociation( jetRef, theRecoTag);
  }

}



void
WZEdmAnalyzer::copyPFCHSJetInfo(   const edm::Event& iEvent,
				   const edm::EventSetup& iSetup, 
				   reco::PFJetCollection::const_iterator jet,
				   edm::RefToBase<reco::Jet> &jetRef, 
				   JetCorrectionUncertainty *jetCorUnc, 
				   double scale, 
				   edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
				   // edm::Handle<reco::JetTagCollection> & jetTags, 
				   _jet_ *myJet) {


  if (_check_jecref) {
    std::cout  << setw(15) << "PFCHS jet: ";
  }
  copyPFJetInfoCommon(iEvent, iSetup, jet, jetRef, jetCorUnc, scale, theRecoTag, myJet);

  // access btagging for PF CHS jet
  myJet->flavor[0]                   = btaggingAssociation( (reco::Jet)(*jet), jetTags.product());
  myJet->flavor[1]                   = btaggingAssociation( (reco::Jet)(*jet), jetTagsCSV.product());

  myJet->flavor[2]                   = btaggingAssociation( jetRef, myPFCHSJetTagsJP) + btaggingAssociation( jetRef, myPFCHSJetTagsTCHP);
  myJet->flavor[3]                   = btaggingAssociation( jetRef, myPFCHSJetTagsTCHP);
  myJet->flavor[4]                   = btaggingAssociation( jetRef, myPFCHSJetTagsCSV);
  myJet->flavor[5]                   = btaggingAssociation( jetRef, myPFCHSJetTagsCMVA);


}

void
WZEdmAnalyzer::copyPFJetInfo(   const edm::Event& iEvent,
				const edm::EventSetup& iSetup, 
				reco::PFJetCollection::const_iterator jet,
				edm::RefToBase<reco::Jet> &jetRef, 
				JetCorrectionUncertainty *jetCorUnc, 
				double scale, 
				edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
				//edm::Handle<reco::JetTagCollection> & jetTags, 
				_jet_ *myJet) {


  if (_check_jecref) {
    std::cout << setw(15) << "PF jet: ";
  }
  copyPFJetInfoCommon(iEvent, iSetup, jet, jetRef, jetCorUnc, scale, theRecoTag, myJet);

  // access btagging for PF jet
  myJet->flavor[0]                   = btaggingAssociation( (reco::Jet)(*jet), jetTags.product());
  myJet->flavor[1]                   = btaggingAssociation( (reco::Jet)(*jet), jetTagsCSV.product());

  myJet->flavor[2]                   = btaggingAssociation( jetRef, myPFJetTagsJP, _is_debug) +  btaggingAssociation( jetRef, myPFJetTagsTCHP, _is_debug);
  myJet->flavor[3]                   = btaggingAssociation( jetRef, myPFJetTagsTCHP, _is_debug);
  myJet->flavor[4]                   = btaggingAssociation( jetRef, myPFJetTagsCSV, _is_debug);
     myJet->flavor[5]                = btaggingAssociation( jetRef, myPFJetTagsCMVA, _is_debug);
  
}




void
WZEdmAnalyzer::copyAK5PFJetInfo(   const edm::Event& iEvent,
				const edm::EventSetup& iSetup, 
				reco::PFJetCollection::const_iterator jet,
				edm::RefToBase<reco::Jet> &jetRef, 
				JetCorrectionUncertainty *jetCorUnc, 
				double scale, 
				edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
				//edm::Handle<reco::JetTagCollection> & jetTags, 
				_jet_ *myJet) {


  if (_check_jecref) {
    std::cout << setw(15) << "PF jet: ";
  }
  copyPFJetInfoCommon(iEvent, iSetup, jet, jetRef, jetCorUnc, scale, theRecoTag, myJet);

  // access btagging for PF jet
  myJet->flavor[0]                   = btaggingAssociation( (reco::Jet)(*jet), jetTags.product());
  myJet->flavor[1]                   = btaggingAssociation( (reco::Jet)(*jet), jetTagsCSV.product());

  myJet->flavor[2]                   = btaggingAssociation( jetRef, myAK5PFJetTagsJP, _is_debug) +  btaggingAssociation( jetRef, myAK5PFJetTagsTCHP, _is_debug) ;
  myJet->flavor[3]                   = btaggingAssociation( jetRef, myAK5PFJetTagsTCHP, _is_debug);
  myJet->flavor[4]                   = btaggingAssociation( jetRef, myAK5PFJetTagsCSV, _is_debug);
  myJet->flavor[5]                   = btaggingAssociation( jetRef, myAK5PFJetTagsCMVA, _is_debug);
  
}


void
WZEdmAnalyzer::copyPFJetInfoCommon(   const edm::Event& iEvent,
				      const edm::EventSetup& iSetup, 
				      reco::PFJetCollection::const_iterator jet,
				      edm::RefToBase<reco::Jet> &jetRef, 
				      JetCorrectionUncertainty *jetCorUnc, 
				      double scale, 
				      edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag,
				      //				edm::Handle<reco::JetTagCollection> & jetTags, 
				      _jet_ *myJet) {

  myJet->L2L3pt    = scale * jet->pt();
  myJet->L2L3scale = scale;
  myJet->energy    = scale*jet->energy();
  myJet->mass      = jet->mass();
  myJet->pt        = jet->pt();
  myJet->eta       = jet->eta();
  myJet->phi       = jet->phi();

  if (jetCorUnc) {
    jetCorUnc->setJetEta((float)jet->eta());
    jetCorUnc->setJetPt((float)( jet->pt() * scale)); // here you must use the CORRECTED jet pt
  // std::cout << pfJetUnc->getUncertainty(true) << std::endl;
    myJet->scaleUnc = jetCorUnc->getUncertainty(true);
  } else {
    myJet->scaleUnc = 0.0;
  }

  if (_check_jecref) {

    std::cout << setw(10) << iEvent.id().run()
	      << setw(10) << iEvent.luminosityBlock()
	      << setw(15) << iEvent.id().event() 
	      << setw(15) <<  jet->pt()
	      << setw(15) <<  jet->eta()
	      << setw(15) <<  jet->phi()
	      << setw(15) <<  scale * jet->pt()
	      << setw(15) <<  100*myJet->scaleUnc
	      << std::endl;

  }


  myJet->chargedHadronEnergy         = jet->chargedHadronEnergy();
  myJet->neutralHadronEnergy         = jet->neutralHadronEnergy();
  myJet->chargedHadronEnergyFraction = jet->chargedHadronEnergyFraction();
  myJet->neutralHadronEnergyFraction = jet->neutralHadronEnergyFraction();
  myJet->chargedEmEnergy             = jet->chargedEmEnergy();
  myJet->chargedEmEnergyFraction     = jet->chargedEmEnergyFraction();
  myJet->neutralEmEnergyFraction     = jet->neutralEmEnergyFraction();
  myJet->neutralEmEnergy             = jet->neutralEmEnergy();
  myJet->neutralEmEnergyFraction     = jet->neutralEmEnergyFraction();
  myJet->muonMultiplicity            = jet->muonMultiplicity();
  myJet->electronMultiplicity        = jet->electronMultiplicity();
  myJet->chargedMultiplicity         = jet->chargedMultiplicity();


  // specific ones for PF jet
  myJet->photonEnergy                = jet->photonEnergy();
  myJet->photonEnergyFraction        = jet->photonEnergyFraction();
  myJet->electronEnergy              = jet->electronEnergy();
  myJet->electronEnergyFraction      = jet->electronEnergyFraction();
  myJet->muonEnergy                  = jet->muonEnergy();
  myJet->muonEnergyFraction          = jet->muonEnergyFraction();
  myJet->HFHadronEnergy              = jet->HFHadronEnergy();
  myJet->HFHadronEnergyFraction      = jet->HFHadronEnergyFraction();
  myJet->HFEMEnergy                  = jet->HFEMEnergy();
  myJet->HFEMEnergyFraction          = jet->HFEMEnergyFraction();
  myJet->chargedHadronMultiplicity   = jet->chargedHadronMultiplicity();
  myJet->neutralHadronMultiplicity   = jet->neutralHadronMultiplicity();
  myJet->photonMultiplicity          = jet->photonMultiplicity();
  myJet->HFHadronMultiplicity        = jet->HFHadronMultiplicity();
  myJet->HFEMMultiplicity            = jet->HFEMMultiplicity();
  myJet->chargedMuEnergy             = jet->chargedMuEnergy();
  myJet->chargedMuEnergyFraction     = jet->chargedMuEnergyFraction();
  myJet->neutralMultiplicity         = jet->neutralMultiplicity();

  myJet->nConstituent                = jet->nConstituents();

  myJet->pileup                      = jet->pileup();
  myJet->area                        = jet->jetArea();

  myJet->L                           = -9999;
  if (jet->jetArea()>0 && (*rhoChHandle)) myJet->L  = (jet->jetArea() * (*rhoChHandle)) * log(  scale * jet->pt()/ (jet->jetArea() * (*rhoChHandle)) );


  //  std::cout << "a jet " << std::endl;
  //
  // std::cout << " pile up = " << jet->pileup() 
  //    << " =? " <<  jet->jetArea() * (*rhoHandle)  << " ; L " << myJet->L 
  //    << std::endl;


  // quark-gluon separator
  /*
  double ptD = 0, pfqt2 = 0, pfqt = 0;
  int pfConstituents = jet->getPFConstituents().size();

  for (int ii = 0; ii < pfConstituents; ii ++) {


    double apt = 0;
    const reco::PFCandidatePtr myptr =    (jet->getPFConstituent(ii));
    apt = myptr->pt();
    pfqt2 += pow(apt, 2);
    pfqt += apt;
  }
  ptD = sqrt( pfqt2/pow(pfqt, 2) );

  myJet->quarkgluonLikelihood=  qgLikelihoodCal->computeQGLikelihoodPU(scale * jet->pt(), *rhoIsoHandle, jet->chargedHadronMultiplicity(), jet->neutralHadronMultiplicity() +jet->photonMultiplicity(), ptD );
  */

  //  if (jet->pt()>20)  std::cout<<  myJet->flavor[1] << "  : " << myJet->flavor[4]  << std::endl;

  //  double mean_dz, mean_dzerr;

  //  int counter =0;
  //  std::cout << jet->getTrackRefs().size() << std::endl;
  //  for (reco::TrackRefVector::const_iterator track = jet->getTrackRefs().begin(); track != jet->getTrackRefs().end(); track ++) {

  TH1F h1("h1", "", 100, 0, 100); h1.Sumw2();
  TH1F h2("h2", "", 100, 0, 100); h2.Sumw2();
  TH1F h3("h3", "", 100, 0, 100); h3.Sumw2();
  double sumN=0, sumS=0, sum =0;
  double totsumN=0, totsumS=0, totsum =0;

  double dissum5  =0, dissum5S  =0, dissum5N  =0;
  double dissum10 =0, dissum10S =0, dissum10N =0;



  for (unsigned int ii = 0; ii < jet->getTrackRefs().size() ; ii ++) {


    reco::TrackRef track = jet->getTrackRefs().at(ii);
    //    if ((*track).pt()<0.05) continue;

    int vertex_index = trackVertexAssociation( (const reco::TrackBaseRef)(track), recVtxs );

    if (vertex_index>=0) {
      h1.Fill(vertex_index, pow((*track).pt(), 2) );
      h2.Fill(vertex_index, (*track).pt() );
      h3.Fill(vertex_index);

      sum  += pow((*track).pt(), 2);
      sumS += (*track).pt();
      sumN += 1;

 
    }  else {

      if (recVtxs->size() > 0) {
	reco::VertexRef vtx(recVtxs, 0);    


	// always overestimate this error
	double vtxDxyCovxy = fabs( vtx->covariance(0, 1) ); // 
	double vtxDxyErr = sqrt( pow( vtx->xError(), 2) + pow( vtx->yError(), 2) + 2 * vtxDxyCovxy );
	
	// debug purpose
	//	std::cout << "noncor = " <<  sqrt( pow( vtx->xError(), 2) + pow( vtx->yError(), 2) ) << " -- " << vtxDxyErr << std::endl;
	

	// track dxy and dxy error
	double mydxy = track->dxy( vtx->position() );
	double mydxyError = track->dxyError();

	double mydz = track->dz( vtx->position() );
	double mydzError = track->dzError();
	
	
	// significance. 
	double dxySig = fabs(mydxy)/sqrt(pow(mydxyError, 2) + pow(vtxDxyErr, 2) );
	double dzSig =  fabs(mydz) /sqrt(pow( mydzError, 2) + pow(vtx->zError(), 2) );

	// matching again
	if ( dxySig<5) {

	  if (dzSig >5) {

	    dissum5  += pow((*track).pt(), 2);
	    dissum5S += (*track).pt();
	    dissum5N += 1;
	  }

	  if (dzSig>10) {

	    dissum10  += pow((*track).pt(), 2);
	    dissum10S += (*track).pt();
	    dissum10N += 1;


	  }
	}  // only consider so-called prompt tracks

      }
    }// end of matching

    totsum  += pow((*track).pt(), 2);
    totsumS += (*track).pt();
    totsumN += 1;
  }

  myJet->dissum5    = dissum5;
  myJet->dissum5S   = dissum5S;
  myJet->dissum5N   = dissum5N;

  myJet->dissum10   = dissum10;
  myJet->dissum10S  = dissum10S;
  myJet->dissum10N  = dissum10N;



  myJet->totsum = totsum; // all tracks above threshold included. 
  myJet->sum    = sum;    // only tracks are matched to one vertex included. 

  myJet->totsumS = totsumS; // all tracks above threshold included. 
  myJet->sumS    = sumS;    // only tracks are matched to one vertex included. 

  myJet->totsumN = totsumN; // all tracks above threshold included. 
  myJet->sumN    = sumN;    // only tracks are matched to one vertex included. 


  if (h1.Integral()>0) {
    myJet->jetVertexIndex =  h1.GetMaximumBin()-1 ;
    myJet->jetVertexFrac     =  1.0*h1.GetMaximum()/h1.Integral();



    int nlIndex = -1;
    float nlVal = -1; float maxVal = h1.GetMaximum();
    for (int kk = 1; kk <= h1.GetNbinsX(); kk ++) {
      
      if ( (h1.GetBinContent(kk) > nlVal)  
	   && ( h1.GetBinContent(kk) < maxVal ) ) {
	
	nlIndex = kk-1;
	nlVal = h1.GetBinContent(kk);
      }
    }

    myJet->jetVertexnlIndex =  nlIndex ;
    myJet->jetVertexnlFrac  =  1.0*nlVal  /h1.Integral();


    //    std::cout << "most active vertex = " << h1.GetMaximumBin()-1 
    //      <<  " ; fraction of tracks = " 
    //      << 1.0*h1.GetMaximum()/h1.Integral() 
    //      << " out of " << h1.Integral() << std::endl;
    
  } else {
    
    myJet->jetVertexnlIndex = -1;
    myJet->jetVertexnlFrac = -1;

    myJet->jetVertexIndex = -1;
    myJet->jetVertexFrac = -1;



  }

  // counting tracks with pt>0.1 GeV
  if (h2.Integral() >0) {

    myJet->jetVertexIndexSum =  h2.GetMaximumBin()-1 ;
    myJet->jetVertexFracSum     =  1.0*h2.GetMaximum()/h2.Integral();



    int nlIndex = -1;
    float nlVal = -1; float maxVal = h2.GetMaximum();
    for (int kk = 1; kk <= h2.GetNbinsX(); kk ++) {
      
      if ( (h2.GetBinContent(kk) > nlVal)  
	   && ( h2.GetBinContent(kk) < maxVal ) ) {
	
	nlIndex = kk-1;
	nlVal = h2.GetBinContent(kk);
      }
    }

    myJet->jetVertexnlIndexSum =  nlIndex ;
    myJet->jetVertexnlFracSum  =  1.0*nlVal  /h2.Integral();

   
    //   std::cout << "most active vertex = " << h2.GetMaximumBin()-1 
    //      <<  " ; fraction of tracks = " 
    //      << 1.0*h2.GetMaximum()/h2.Integral() 
    //      << " out of " << h2.Integral() << std::endl;

  } else {


    myJet->jetVertexnlIndexSum = -1;
    myJet->jetVertexnlFracSum = -1;
    
    myJet->jetVertexIndexSum = -1;
    myJet->jetVertexFracSum = -1;
  }




 // counting tracks with pt>0.1 GeV
  if (h3.Integral() >0) {

    myJet->jetVertexIndexNum =  h3.GetMaximumBin()-1 ;
    myJet->jetVertexFracNum     =  1.0*h3.GetMaximum()/h3.Integral();



    int nlIndex = -1;
    float nlVal = -1; float maxVal = h3.GetMaximum();
    for (int kk = 1; kk <= h3.GetNbinsX(); kk ++) {
      
      if ( (h3.GetBinContent(kk) > nlVal)  
	   && ( h3.GetBinContent(kk) < maxVal ) ) {
	
	nlIndex = kk-1;
	nlVal = h3.GetBinContent(kk);
      }
    }

    myJet->jetVertexnlIndexNum =  nlIndex ;
    myJet->jetVertexnlFracNum  =  1.0*nlVal  /h3.Integral();

   
    //   std::cout << "most active vertex = " << h3.GetMaximumBin()-1 
    //      <<  " ; fraction of tracks = " 
    //      << 1.0*h3.GetMaximum()/h3.Integral() 
    //      << " out of " << h3.Integral() << std::endl;

  } else {


    myJet->jetVertexnlIndexNum = -1;
    myJet->jetVertexnlFracNum = -1;
    
    myJet->jetVertexIndexNum = -1;
    myJet->jetVertexFracNum = -1;
  }


  if (!_is_data) {

    myJet->mc_flavor = mcflavorAssociation( jetRef, theRecoTag, myJet->mc_partonFlavor, _is_debug);
  }
}



void 
WZEdmAnalyzer::copySuperclusterInfo( reco::SuperClusterCollection::const_iterator supercluster, _supercluster_ *mySupercluster) {


  mySupercluster->et       = supercluster->energy() /TMath::CosH(supercluster->eta()); ;
  mySupercluster->eta      = supercluster->eta();
  mySupercluster->phi      = supercluster->phi();


  mySupercluster->scRawEnergy      = supercluster->rawEnergy();
  mySupercluster->scPreshowerEnergy      = supercluster->preshowerEnergy();
  mySupercluster->scPhiWidth      = supercluster->phiWidth();
  mySupercluster->scEtaWidth      = supercluster->etaWidth();
  mySupercluster->scClustersSize      = supercluster->clustersSize();


  mySupercluster->energy      = supercluster->energy();
  //  mySupercluster->calEnergy      = supercluster->calEnergy();
  mySupercluster->centroidX      = supercluster->x();
  mySupercluster->centroidY      = supercluster->y();
  mySupercluster->centroidZ      = supercluster->z();
  mySupercluster->size      = supercluster->size();
}



void 
WZEdmAnalyzer::copyPhotonInfo( reco::PhotonCollection::const_iterator photon, _photon_ *myPhoton) 
{



  myPhoton->et       = photon->et();
  myPhoton->eta      = photon->eta();
  myPhoton->phi      = photon->phi();


 myPhoton->fiducialFlags = photon->isEB();
  myPhoton->fiducialFlags +=  (photon->isEE())<<1;
  //  myPhoton->fiducialFlags +=  (photon->isGap())<<2;
  myPhoton->fiducialFlags +=  (photon->isEBEEGap())<<3;
  myPhoton->fiducialFlags +=  (photon->isEBGap())<<4;
  myPhoton->fiducialFlags +=  (photon->isEBEtaGap())<<5;
  myPhoton->fiducialFlags +=  (photon->isEBPhiGap())<<6;
  myPhoton->fiducialFlags +=  (photon->isEEGap())<<7;
  myPhoton->fiducialFlags +=  (photon->isEEDeeGap())<<8;
  myPhoton->fiducialFlags +=  (photon->isEERingGap())<<9;


  myPhoton->hasPixelSeed      = photon->hasPixelSeed();
  myPhoton->hasConversionTracks      = photon->hasConversionTracks();
  myPhoton->sigmaEtaEta      = photon->sigmaEtaEta();
  myPhoton->sigmaIetaIeta      = photon->sigmaIetaIeta();
  //  myPhoton->sigmaIetaIeta      = photon->sigmaIetaIeta();
  myPhoton->e1x5      = photon->e1x5();
  myPhoton->e2x5      = photon->e2x5();
  myPhoton->e3x3      = photon->e3x3();
  myPhoton->e5x5      = photon->e5x5();
  myPhoton->maxEnergyXtal      = photon->maxEnergyXtal();
  myPhoton->hadronicDepth1OverEm      = photon->hadronicDepth1OverEm();
  myPhoton->hadronicDepth2OverEm      = photon->hadronicDepth2OverEm();
  myPhoton->hadronicOverEm      = photon->hadronicOverEm();
  myPhoton->r1x5      = photon->r1x5();
  myPhoton->r2x5      = photon->r2x5();
  myPhoton->r9      = photon->r9();
  myPhoton->ecalRecHitSumEtConeDR04      = photon->ecalRecHitSumEtConeDR04();
  myPhoton->hcalTowerSumEtConeDR04      = photon->hcalTowerSumEtConeDR04();
  myPhoton->hcalDepth1TowerSumEtConeDR04      = photon->hcalDepth1TowerSumEtConeDR04();
  myPhoton->hcalDepth2TowerSumEtConeDR04      = photon->hcalDepth2TowerSumEtConeDR04();
  myPhoton->trkSumPtSolidConeDR04      = photon->trkSumPtSolidConeDR04();
  myPhoton->trkSumPtHollowConeDR04      = photon->trkSumPtHollowConeDR04();
  myPhoton->nTrkSolidConeDR04      = photon->nTrkSolidConeDR04();
  myPhoton->nTrkHollowConeDR04      = photon->nTrkHollowConeDR04();



  myPhoton->ecalRecHitSumEtConeDR03      = photon->ecalRecHitSumEtConeDR03();
  myPhoton->hcalTowerSumEtConeDR03      = photon->hcalTowerSumEtConeDR03();
  myPhoton->hcalDepth1TowerSumEtConeDR03      = photon->hcalDepth1TowerSumEtConeDR03();
  myPhoton->hcalDepth2TowerSumEtConeDR03      = photon->hcalDepth2TowerSumEtConeDR03();
  myPhoton->trkSumPtSolidConeDR03      = photon->trkSumPtSolidConeDR03();
  myPhoton->trkSumPtHollowConeDR03      = photon->trkSumPtHollowConeDR03();
  myPhoton->nTrkSolidConeDR03      = photon->nTrkSolidConeDR03();
  myPhoton->nTrkHollowConeDR03      = photon->nTrkHollowConeDR03();



  myPhoton->scRawEnergy      = photon->superCluster()->rawEnergy();
  myPhoton->scPreshowerEnergy      = photon->superCluster()->preshowerEnergy();
  myPhoton->scPhiWidth      = photon->superCluster()->phiWidth();
  myPhoton->scEtaWidth      = photon->superCluster()->etaWidth();
  myPhoton->scClustersSize      = photon->superCluster()->clustersSize();


  myPhoton->energy      = photon->superCluster()->energy();
  myPhoton->centroidX      = photon->superCluster()->x();
  myPhoton->centroidY      = photon->superCluster()->y();
  myPhoton->centroidZ      = photon->superCluster()->z();
  myPhoton->size      = photon->superCluster()->size();


  // compute the second moment
  // uncomment in 3_*_x release
  if (_reco) {

    CaloClusterPtr SCseed = photon->superCluster()->seed();
    if (photon->isEE()) {
      const EERecHitCollection* rhits =ecalEERecHitHandle.product();
      Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments( *SCseed, *rhits);
      

  //std::cout << moments.alpha << std::endl;
  //  myPhoton->cluster2ndMoments_alpha = moments.alpha;


      myPhoton->cluster2ndMoments_sMaj  = moments.sMaj;
      myPhoton->cluster2ndMoments_sMin  = moments.sMin;
  
    } else if ( photon->isEB() ) {


      const EERecHitCollection* rhits =ecalEBRecHitHandle.product();
      Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments( *SCseed, *rhits);


      //std::cout << moments.alpha << std::endl;
      //  myPhoton->cluster2ndMoments_alpha = moments.alpha;
      

      myPhoton->cluster2ndMoments_sMaj  = moments.sMaj;
      myPhoton->cluster2ndMoments_sMin  = moments.sMin;
      
    } else {

      myPhoton->cluster2ndMoments_sMaj  = -9999;
      myPhoton->cluster2ndMoments_sMin  = -9999;
      myPhoton->cluster2ndMoments_alpha = -9999;
    }
  }
}


void
WZEdmAnalyzer::copyElectronInfo( reco::GsfElectronCollection::const_iterator electron, 
				 _electron_ *myElectron, 
				 reco::GsfElectronRef electronRef, 
				 int ele_index) 
{

  // copy gsf electron information
  // reasign GSF electron track pt, eta, phi to GSF electron track pt, eta, phi 

  this->copyTrackInfo( ( const reco::Track * )( &(*electron->gsfTrack()) ), (_track_ *) myElectron);
  myElectron->pt    = electron->pt();
  myElectron->eta   = electron->eta();
  myElectron->phi   = electron->phi();
  myElectron->charge = electron->charge();  



    
  // photon information
  myElectron->classification      = electron->classification();

  myElectron->scPixCharge = electron->scPixCharge();
  myElectron->isGsfCtfScPixChargeConsistent = electron->isGsfCtfScPixChargeConsistent();
  myElectron->isGsfScPixChargeConsistent = electron->isGsfScPixChargeConsistent();
  myElectron->isGsfCtfChargeConsistent = electron->isGsfCtfChargeConsistent();
  myElectron->ecalDrivenSeed = electron->ecalDrivenSeed();
  myElectron->trackerDrivenSeed = electron->trackerDrivenSeed();



  // conversion rejection
  // old implementation 
  ConversionFinder convFinder;
  myElectron->isElFromConversion = 0;
  //convFinder.isFromConversion();
  // only work for a later release
  ConversionInfo convInfo =
    convFinder.getConversionInfo((*electron), tracks, bField);
  
  myElectron->dist = convInfo.dist();
  myElectron->dcot = convInfo.dcot();
  myElectron->radiusOfConversion = convInfo.radiusOfConversion();
  myElectron->deltaMissingHits= convInfo.deltaMissingHits();


  // implementation from 2012 data
  // const edm::Handle<reco::ConversionCollection> convCol; 
  myElectron->vtxFitConversion = ConversionTools::hasMatchedConversion(*electron, convCol, vertexBeamSpot->position());
  
  //NOTE74: not in 74x?
  myElectron->mHits = electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

  //trackerExpectedHitsInner().numberOfHits(); 
  //myElectron->numberOfLostHits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(); 




  // number of expected inner hit for electron ID
  //NOTE74: 
  /*
  const reco::Track *el_track = (const reco::Track*)(electron->gsfTrack().get());  
  const reco::HitPattern& p_inner = el_track->trackerExpectedHitsInner(); 
  myElectron->numberOfExpectedInnerHits = p_inner.numberOfHits();
  */

  //  std::cout << electron->ecalEnergy() 
  //    << setw(20) << electron->superCluster()->energy()
  //    << std::endl;

  myElectron->scEt                   = electron->superCluster()->energy()/TMath::CosH(electron->superCluster()->eta());
  myElectron->scEta                  = electron->superCluster()->eta();
  myElectron->scPhi                  = electron->superCluster()->phi();
  myElectron->energy                 = electron->superCluster()->energy();
  myElectron->centroidX              = electron->superCluster()->x();
  myElectron->centroidY              = electron->superCluster()->y();
  myElectron->centroidZ              = electron->superCluster()->z();
  myElectron->size                   = electron->superCluster()->size();
  myElectron->scRawEnergy            = electron->superCluster()->rawEnergy();
  myElectron->scPreshowerEnergy      = electron->superCluster()->preshowerEnergy();
  myElectron->scPhiWidth             = electron->superCluster()->phiWidth();
  myElectron->scEtaWidth             = electron->superCluster()->etaWidth();
  myElectron->scClustersSize         = electron->superCluster()->clustersSize();




  myElectron->chargeMode= (*(electron->gsfTrack())).chargeMode(); 
  myElectron->qoverpMode= (*(electron->gsfTrack())).qoverpMode(); 
  myElectron->lambdaMode= (*(electron->gsfTrack())).lambdaMode(); 
  myElectron->ptMode= (*(electron->gsfTrack())).ptMode(); 
  myElectron->phiMode= (*(electron->gsfTrack())).phiMode(); 
  myElectron->etaMode= (*(electron->gsfTrack())).etaMode(); 
  myElectron->qoverpModeError= (*(electron->gsfTrack())).qoverpModeError(); 
  myElectron->lambdaModeError= (*(electron->gsfTrack())).lambdaModeError(); 
  myElectron->ptModeError= (*(electron->gsfTrack())).ptModeError(); 
  myElectron->phiModeError= (*(electron->gsfTrack())).phiModeError(); 
  myElectron->etaModeError= (*(electron->gsfTrack())).etaModeError(); 




  
  // track cluster matching 
  if (_is_debug) {std::cout << "electron-track matching " << std::endl;}



  myElectron->eSuperClusterOverP        = electron->eSuperClusterOverP();      
  myElectron->eSeedClusterOverP         = electron->eSeedClusterOverP();         
  myElectron->eSeedClusterOverPout      = electron->eSeedClusterOverPout() ;      
  myElectron->eEleClusterOverPout       = electron->eEleClusterOverPout();       

  myElectron->deltaEtaSuperClusterTrackAtVtx = electron->deltaEtaSuperClusterTrackAtVtx(); 
  myElectron->deltaEtaSeedClusterTrackAtCalo = electron->deltaEtaSeedClusterTrackAtCalo(); 
  myElectron->deltaEtaEleClusterTrackAtCalo  = electron->deltaEtaEleClusterTrackAtCalo();  

  myElectron->deltaPhiSuperClusterTrackAtVtx = electron->deltaPhiSuperClusterTrackAtVtx(); 
  myElectron->deltaPhiEleClusterTrackAtCalo  = electron->deltaPhiEleClusterTrackAtCalo();  
  myElectron->deltaPhiSeedClusterTrackAtCalo = electron->deltaPhiSeedClusterTrackAtCalo(); 


  myElectron->fiducialFlags = electron->isEB() ? true:false;
  //  myElectron->fiducialFlags +=  (electron->isEE())<<1;
  // myElectron->fiducialFlags +=  (electron->isGap())<<2;
  // myElectron->fiducialFlags +=  (electron->isEBEEGap())<<3;
  // myElectron->fiducialFlags +=  (electron->isEBGap())<<4;
  // myElectron->fiducialFlags +=  (electron->isEBEtaGap())<<5;
  // myElectron->fiducialFlags +=  (electron->isEBPhiGap())<<6;
  // myElectron->fiducialFlags +=  (electron->isEEGap())<<7;
  // myElectron->fiducialFlags +=  (electron->isEEDeeGap())<<8;
  // myElectron->fiducialFlags +=  (electron->isEERingGap())<<9;



  //  }
  // Electron Isolation
  if (_is_debug) {std::cout << "electron isolation " << std::endl;}
  myElectron->dr03TkSumPt              = electron->dr03TkSumPt();
  myElectron->dr03EcalRecHitSumEt      = electron->dr03EcalRecHitSumEt();
  myElectron->dr03HcalDepth1TowerSumEt = electron->dr03HcalDepth1TowerSumEt();
  myElectron->dr03HcalDepth2TowerSumEt = electron->dr03HcalDepth2TowerSumEt();
  myElectron->dr03HcalTowerSumEt       = electron->dr03HcalTowerSumEt();
    
  myElectron->dr04TkSumPt              = electron->dr04TkSumPt();
  myElectron->dr04EcalRecHitSumEt      = electron->dr04EcalRecHitSumEt();
  myElectron->dr04HcalDepth1TowerSumEt = electron->dr04HcalDepth1TowerSumEt();
  myElectron->dr04HcalDepth2TowerSumEt = electron->dr04HcalDepth2TowerSumEt();
  myElectron->dr04HcalTowerSumEt       = electron->dr04HcalTowerSumEt();
 

  // pf isolation
  GsfElectron::PflowIsolationVariables pfIso = electron -> pfIsolationVariables();
  Double_t isoChargedHadrons_ = pfIso.sumChargedHadronPt ;
  Double_t isoNeutralHadrons_ = pfIso.sumNeutralHadronEt ;
  Double_t isoPhotons_        = pfIso.sumPhotonEt ;
  Double_t isoChargedFromPU_  = pfIso.sumPUPt ;
    
  double area =  _effectiveAreas.getEffectiveArea( std::abs( electron-> superCluster()->eta() ) );

  // std::cout << setw(20) <<  std::abs( electron-> superCluster()->eta() )
  //	    << setw(20) << area
  //	    << std::endl;

  float relIsoWithEA_ = ( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - (*fixGridRhoHandle) * area ) )/electron->pt();
    
  // Compute isolation with delta beta correction for PU
  float absiso        = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * isoChargedFromPU_ );
  float relIsoWithDBeta_ = absiso/electron->pt();
    

  myElectron->pfRelIsoR03    = relIsoWithDBeta_;
  myElectron->pfIsoCh        = isoChargedHadrons_;
  myElectron->pfIsoPhoton    = isoPhotons_;
  myElectron->pfIsoNeutral   = isoNeutralHadrons_;
  myElectron->pfIsoSumPUPt   = isoChargedFromPU_;
  myElectron->effArea        = area;

  myElectron->pfRelIsoR03EA  = relIsoWithEA_;

  // ECAL & HCAL isolation for preselection

  float ecalPFClusterIso  = (*electronEcalPFClusterIsolation)[electronRef];
  float hcalPFClusterIso  = (*electronHcalPFClusterIsolation)[electronRef];
  myElectron->ecalPFClusterIso = ecalPFClusterIso;
  myElectron->hcalPFClusterIso = hcalPFClusterIso;

  bool pass_preselction = electron->pt()>15 
    &&( ( std::abs(electron->superCluster()->eta()) < 1.4442 
	  && electron->full5x5_sigmaIetaIeta() < 0.012 
	  && electron->hcalOverEcal() < 0.09 
	  && (ecalPFClusterIso / electron->pt()) < 0.37 
	  && (hcalPFClusterIso / electron->pt()) < 0.25 
	  && (electron->dr03TkSumPt() / electron->pt()) < 0.18 
	  && std::abs(electron->deltaEtaSuperClusterTrackAtVtx()) < 0.0095 
	  && std::abs(electron->deltaPhiSuperClusterTrackAtVtx()) < 0.065 ) 
	|| ( std::abs(electron->superCluster()->eta()) > 1.5660 
	     && electron->full5x5_sigmaIetaIeta() < 0.033 
	     && electron->hcalOverEcal() <0.09 
	     && (ecalPFClusterIso / electron->pt()) < 0.45 
	     && (hcalPFClusterIso / electron->pt()) < 0.28 
	     && (electron->dr03TkSumPt() / electron->pt()) < 0.18 ) 
	);


 // Shower Shape
  if (_is_debug) {std::cout << "electron shower variables " << std::endl;}
  myElectron->sigmaEtaEta               = electron->sigmaEtaEta();
  myElectron->sigmaIetaIeta             = electron->sigmaIetaIeta();
  myElectron->full5x5_sigmaIetaIeta     = electron->full5x5_sigmaIetaIeta();
  myElectron->e1x5                      = electron->e1x5();
  myElectron->e2x5Max                   = electron->e2x5Max();
  myElectron->e5x5                      = electron->e5x5();
  myElectron->hcalDepth1OverEcal        = electron->hcalDepth1OverEcal();
  myElectron->hcalDepth2OverEcal        = electron->hcalDepth2OverEcal();
  myElectron->hcalOverEcal              = electron->hcalOverEcal();
  myElectron->hadronicOverEm            = electron->hadronicOverEm();




    
  myElectron->ecalDriven                 = electron->ecalDriven();
  myElectron->passingCutBasedPreselection= electron->passingCutBasedPreselection();
  myElectron->passingMvaPreselection     = electron->passingMvaPreselection();
  //NOTE74:
  //  myElectron->mva                        = electron->mva();
  myElectron->fbrem                      = electron->fbrem();
  myElectron->numberOfBrems              = electron->numberOfBrems();
  //  myElectron->isMomentumCorrected = electron->isMomentumCorrected();
  myElectron->ecalEnergy                 = electron->ecalEnergy();
  myElectron->correctedEcalEnergy        = electron->correctedEcalEnergy();
  myElectron->ecalEnergyError            = electron->ecalEnergyError();
  myElectron->trackMomentumError         = electron->trackMomentumError();
  myElectron->trackMomentumAtVtx         = electron->trackMomentumAtVtx().R();

  myElectron ->ooemoop                   = (1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy());


 
  // Oct. 22, 2015
  // new cut-based electron implementation
  // Dec. 30, 2015
  // add new MVA electron ID flag and MVA values
  unsigned int pass_loose         = (*loose_id_decisions)[electronRef];
  unsigned int pass_medium        = (*medium_id_decisions)[electronRef];
  unsigned int pass_tight         = (*tight_id_decisions)[electronRef];


  unsigned int trgmva_medium      = (*trigMvaMedium_id_decisions)[electronRef];
  unsigned int trgmva_tight       = (*trigMvaTight_id_decisions)[electronRef];

  unsigned int nontrgmva_medium   = (*nonTrigMvaMedium_id_decisions)[electronRef];
  unsigned int nontrgmva_tight    = (*nonTrigMvaTight_id_decisions)[electronRef];

  //  unsigned int trigwp70           = 0;
  // unsigned int trigtight          = 0;
  // unsigned int pass_eoverpcuts    = 0;
  myElectron -> idBitMap          = (trgmva_tight<< 6)
    + (trgmva_medium << 5)
    + (nontrgmva_tight << 4) 
    + (nontrgmva_medium << 3)
    + (pass_tight << 2) 
    + (pass_medium << 1) 
    + (pass_loose << 0);


  myElectron->mvaTrigV0         = (*trigMvaValues)[electronRef];
  myElectron->mvaTrigV0Cat      = (*trigMvaCategories)[electronRef];
  myElectron->mvaNonTrigV0 	= (*nonTrigMvaValues)[electronRef]; 
  myElectron->mvaNonTrigV0Cat   = (*nonTrigMvaCategories)[electronRef];

  myElectron -> mvaIdBitMap     = (pass_preselction<<4)
    + (trgmva_tight<< 3)
    + (trgmva_medium << 2)
    + (nontrgmva_tight << 1) 
    + (nontrgmva_medium << 0);


  myElectron->dEtaInSeed = (electron->superCluster().isNonnull() && electron->superCluster()->seed().isNonnull()) ? ( electron->deltaEtaSuperClusterTrackAtVtx() - electron->superCluster()->eta() + electron->superCluster()->seed()->eta()) : std::numeric_limits<float>::max();
  myElectron->dPhiIn     = electron->deltaPhiSuperClusterTrackAtVtx();


  //const edm::ValueMap<double> ele_regEne = (*regEne_handle.product());
  
  //    edm::Handle<edm::ValueMap<double>> regErr_handle;
  //  iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);
  //  const edm::ValueMap<double> ele_regErr = (*regErr_handle.product());


  if (_is_debug) {std::cout << "finish copying electron iformation " << std::endl;}
}



void WZEdmAnalyzer::copyTrgBits(const edm::Event& iEvent, _trg_bits_ * trgBits ) {


  if(_is_debug) std::cout << "copy L1/HLT trigger bits " << std::endl;

  ULong_t l1bits[MAX_TRIG_WORD];
  ULong_t l1tbits[MAX_TRIG_WORD];
  ULong_t hltaccbits[MAX_TRIG_WORD];
  ULong_t hltrunbits[MAX_TRIG_WORD];
  ULong_t hlterrbits[MAX_TRIG_WORD];
  for (int ii = 0; ii < MAX_TRIG_WORD; ii ++) {
    l1bits[ii] = 0;
    l1tbits[ii] = 0;
    hltaccbits[ii] = 0;
    hltrunbits[ii] = 0;
    hlterrbits[ii] = 0;

  }


  if (l1TriggerReadoutRec.isValid()) {

    //    int itrig(0);
    for (unsigned int ii = 0; ii < l1TriggerReadoutRec->decisionWord().size(); ii++) {

      if (ii >=MAX_TRIG_WORD*32) break;


      //      if (totalProcessedEvts==1) {

      //	const AlgorithmMap & algorithmMap = l1GTTriggerMenu->gtAlgorithmMap();

      //	for (CItAlgo itAlgo = algorithmMap.begin(); itAlgo != algorithmMap.end(); itAlgo ++) {

      //	  std::cout << itAlgo->first << std::endl;

      //	}

      //}
      int l1flag  = l1TriggerReadoutRec->decisionWord()[ii];
      int l1tflag = l1TriggerReadoutRec->technicalTriggerWord()[ii];

      
      if (_is_debug) {


	std::cout << "L1  trigger bits: "<< l1flag << std::endl;
	std::cout << "L1T trigger bits: "<< l1tflag << std::endl;

      }
      int index = ii/32;
      int shift = ii%32;

      if (l1flag)    l1bits[index]  = l1bits[index]  | (0x1 << shift);
      if (l1tflag)   l1tbits[index] = l1tbits[index] | (0x1 << shift);


    }

    trgBits->setL1TrgBits( l1bits);
    trgBits->setL1TTrgBits(l1tbits);
  }



  if (_is_debug) std::cout << "L1/L1T bits ... done" << std::endl;
  if (triggerResults.isValid()) {

  edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerResults);  
  //    edm::TriggerNames triggerNames(*triggerResults);  
    unsigned int npath = triggerResults->size();

    for (unsigned int ii = 0; ii < npath; ++ii) {

      if (ii >=MAX_TRIG_WORD*32) break;
      int hltacc = triggerResults->accept(ii);
      int hltrun = triggerResults->wasrun(ii);
      int hlterr = triggerResults->error(ii);

      int index = ii/32;
      int shift = ii%32;

      if (_is_debug) std::cout << "HLT: acc " << hltacc 
			       << " run " << hltrun
			       << " err " << hlterr
			       << std::endl;

      if (hltacc)   hltaccbits[index] = hltaccbits[index] | (0x1 << shift);
      if (hltrun)   hltrunbits[index] = hltrunbits[index] | (0x1 << shift);
      if (hlterr)   hlterrbits[index] = hlterrbits[index] | (0x1 << shift);
    }

    trgBits->setHLTAccTrgBits( hltaccbits);
    trgBits->setHLTRunTrgBits( hltrunbits);
    trgBits->setHLTErrTrgBits( hlterrbits);
  }
  if (_is_debug) std::cout << "HLT bits ... done" << std::endl;

  return;
}


void 
WZEdmAnalyzer::copyHLTInfo(const edm::Event& iEvent, 
			   edm::Handle<TriggerResults>         &triggerResults, 
			   std::string &triggerProcessName, 
			   edm::Handle<trigger::TriggerEvent>  &triggerObj, 
			   HLTConfigProvider &hltConfig,
			   bool &trgBit, const char *hltname,  
			   _event_ *aEvent,  
			   _vec4_ *(*ptr2AddHlt)(_event_ *), 
			   _vec4_ *(*ptr2AddHltLeg1)(_event_ *),
			   _vec4_ *(*ptr2AddHltLeg2)(_event_ *)) 
{


  if (_is_debug) std::cout << "copy HLT information " << hltname << std::endl;

  trgBit = false;

  if (!triggerResults.isValid() || !triggerObj.isValid()) {
    edm::LogInfo("FourVectorHLTOffline") << "TriggerResults not found, "
					 << "skipping event"; 
    return;
  }
  //edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerResults);  




  unsigned int triggerIndex = hltConfig.triggerIndex( std::string(hltname) ); 

  if (_is_debug) std::cout << "trigger index = " << triggerIndex << "; total paths = " << hltConfig.size() << std::endl;

  //  edm::TriggerNames triggerNames(*triggerResults);  
  //  int npath = triggerResults->size();
  //  std::vector< std::string > triggerNames = hltConfig.triggerNames();
  // int npath = triggerNames.size();

  // HLT path is not existing in current menu
  if (hltConfig.size()<=triggerIndex) return;




 
  if (triggerResults->accept(triggerIndex) ) { 
      trgBit = true;
      if (_is_debug)   std::cout << "event has been triggered " << triggerProcessName << "::" 
				 << hltname << std::endl;
  }


  // std::cout<< hltConfig.tableName() << std::endl;


    //  if (!trgBit) return;


  //  std::cout<< hltConfig.tableName() << std::endl;
  std::vector< std::string > trigModules = hltConfig.moduleLabels( std::string(hltname) );
  

  // locate relavent L3 filters
  // this only matters for double lepton and mixed e mu triggers
  std::vector< std::string > filtersL3;
  filtersL3.resize(0);


  // ee: 0
  // me: 1
  // mm: 2
  int channels =-1;
  for (unsigned int ii = 0; ii < trigModules.size() ; ii++) {
  
 
    if (   hltConfig.saveTags( std::string(trigModules[ii]) ) ) {

      if (_is_debug) {
	std::cout << "module: " << trigModules[ii]
		  << " is L3 filter? " 
		  <<      hltConfig.saveTags( std::string(trigModules[ii]) ) << std::endl;
      }



      // double muon 
      if (  ((std::string( hltname) ).find("Mu") != string::npos) 
	    && ( (std::string( hltname) ).find("Mu") !=  (std::string( hltname) ).rfind("Mu")) 
	    ) {
 
	channels = 2;
	if ( ( std::string( trigModules[ii] ).find("Filtered8") != string::npos )
	     ||  ( std::string( trigModules[ii] ).find("L3Filtered17") != string::npos )
	     ||  ( std::string( trigModules[ii] ).find("DzFiltered") != string::npos )
	     ) {


	  if (  ( std::string( trigModules[ii] ).find("L3Filtered17") != string::npos ) && filtersL3.size()>0) {

	    filtersL3.insert( filtersL3.begin(), std::string( trigModules[ii] ) );

	  } else {

	    filtersL3.push_back( std::string( trigModules[ii] ) );
	  }
	}
	
      } else if ( ((std::string( hltname) ).find("Ele") != string::npos) 
		  && ( (std::string( hltname) ).find("Ele") !=  (std::string( hltname) ).rfind("Ele")))  { // double electron trigger
	
	channels = 0;
	if ( ( std::string( trigModules[ii] ).find("TrackIsoLeg1Filter") != string::npos )
	     || ( std::string( trigModules[ii] ).find("TrackIsoLeg2Filter") != string::npos )
	     || ( std::string( trigModules[ii] ).find("DZFilter") != string::npos )
	    
	    //(( std::string( trigModules[ii] ).find("Ele17") != string::npos )
	     // &&  ( std::string( trigModules[ii] ).find("Ele8") != string::npos ) )  
	    //	     ||  
	    //(( std::string( trigModules[ii] ).find("Ele17") != string::npos )
	    // &&  ( std::string( trigModules[ii] ).find("Ele12") != string::npos ) ) 
	     )
	  filtersL3.push_back( std::string( trigModules[ii] ) );
	//filtersL3.insert(filtersL3.begin(),  std::string( trigModules[ii] ) );

      } else if (  ((std::string( hltname) ).find("Mu") != string::npos)
		   &&  ((std::string( hltname) ).find("Ele") != string::npos ) ) { // muon electron trigger
	

	channels = 1;
	if ( ( (std::string(trigModules[ii])).find("MuonlegL3IsoFiltered") != string::npos )
	     || ( (std::string(trigModules[ii])).find("ElectronlegTrackIsoFilter") != string::npos )
	     //	     || (  ( (std::string(trigModules[ii])).find("Mu") != string::npos ) && ( (std::string(trigModules[ii])).find("Ele") != string::npos )   ) 
	     ) 	  {


	  filtersL3.push_back( std::string( trigModules[ii] ) );
	}

      } else { // single lepton trigger

	channels = -1;
	if ( trigModules.size()>=2 && (ii == trigModules.size()-2) ) 
	  filtersL3.push_back( std::string( trigModules[ii] ) );
      }

    }
  }




  std::string longname;  



  // testing module
  // copy HLT four vectors. 
  const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());

  for (std::vector< std::string>::iterator it = filtersL3.begin(); it != filtersL3.end(); it ++) {
    
    
    if (_is_debug)    std::cout << *it << std::endl;
    
    longname.resize(0);
    longname.append( *it);
    longname.append("::");
    longname.append( triggerProcessName); 

    const int index = triggerObj->filterIndex(edm::InputTag(longname));

    
    if ( index >= triggerObj->sizeFilters() || index <0) { if (_is_debug) std::cout << "no trigger object associated with this path " << std::endl; continue;}

    
    const trigger::Keys & k = triggerObj->filterKeys(index);
    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {


      _vec4_ *hltobj = 0;
      _vec4_ *hltobj_dmleg = 0;


      // ee case
      if ( (channels ==0) || (channels == 2)  ) {

	if (it == filtersL3.begin()) {

	  if (ptr2AddHltLeg1) hltobj =       ptr2AddHltLeg1(aEvent);

	} else if ( it == (++filtersL3.begin()) ) {

	  if (ptr2AddHltLeg2) hltobj =       ptr2AddHltLeg2(aEvent);

	} else if ( it == (++ (++filtersL3.begin()) ) ) {
	  
	  hltobj =       ptr2AddHlt(aEvent);
	} else {
	  
	}

	//      } else if (channels == 2) { // mm case
	//
	//
	//
	//	hltobj =       ptr2AddHlt(aEvent);
	//
	//	if ( ptr2AddHltLeg1 
	//   && ptr2AddHltLeg2) {
	//
	//if (it == filtersL3.begin()) {
	//
	//  hltobj_dmleg = ptr2AddHltLeg2(aEvent);
	//} else {
	//
	//  hltobj_dmleg = ptr2AddHltLeg1(aEvent);
	//}
	//	}

      } else { // others
	
	hltobj =       ptr2AddHlt(aEvent);
      }

      if (hltobj) {
	hltobj->id = toc[*ki].id();
	hltobj->pt = toc[*ki].pt();
	hltobj->eta =  toc[*ki].eta();
	hltobj->phi =  toc[*ki].phi();
      

	if (_is_debug) {
	std::cout << "HLT trig object = (" 
		  << toc[*ki].id() << ", "
		  << toc[*ki].pt() << ", "
		  << toc[*ki].eta() << ", "
		  << toc[*ki].phi() << ")"
		  << std::endl;
	}      

      }


      if (hltobj_dmleg) {

	hltobj_dmleg->id = toc[*ki].id();
	hltobj_dmleg->pt = toc[*ki].pt();
	hltobj_dmleg->eta =  toc[*ki].eta();
	hltobj_dmleg->phi =  toc[*ki].phi();
      

	if(_is_debug) {
	std::cout << "HLT trig object = (" 
		  << toc[*ki].id() << ", "
		  << toc[*ki].pt() << ", "
		  << toc[*ki].eta() << ", "
		  << toc[*ki].phi() << ")"
		  << std::endl;
	}

      }
     
    }
  }
  

  return;


  longname.resize(0);
  if (trigModules.size()>=2) 
  longname.append( trigModules[ trigModules.size()-2 ]  );
  longname.append("::");
  longname.append( triggerProcessName); 
  if (_is_debug)  {
    std::cout <<  longname << std::endl;
    for (unsigned int ii = 0; ii < trigModules.size() ; ii++) {
  
      std::cout << "module: " << trigModules[ii]
		<< " is L3 filter? " 
		<<      hltConfig.saveTags( std::string(trigModules[ii]) ) << std::endl;

    }
  }



  // test a new way to copy HLT trigger objects
  /*
  const trigger::TriggerObjectCollection & mytoc(triggerObj->getObjects());
  for (unsigned int ii = 0; ii < trigModules.size() ; ii++) {

    std::string mylongname;
    mylongname.resize(0);

    mylongname.append( trigModules[ii] );
    mylongname.append("::");
    mylongname.append( triggerProcessName); 

    // skip level 1 and level 2s
    if (mylongname.find("hltL1")  != string::npos ) continue;
    if (mylongname.find("hltL2")  != string::npos ) continue;


    if (_is_debug)  std::cout << mylongname << std::endl;

    const int myindex = triggerObj->filterIndex(edm::InputTag(mylongname));
    if ( myindex >= triggerObj->sizeFilters() || myindex <0) { continue; }
    //if ( 

    const trigger::Keys & myk = triggerObj->filterKeys(myindex);
    if (_is_debug) {
      
      std::cout << mylongname << "  of " << myk.size() << std::endl;
    }
    
    for (trigger::Keys::const_iterator ki = myk.begin(); ki !=myk.end(); ++ki ) {
      
      _vec4_ *hltobj = ptr2AddHlt(aEvent);
      hltobj->id = mytoc[*ki].id();
      hltobj->pt = mytoc[*ki].pt();
      hltobj->eta =  mytoc[*ki].eta();
      hltobj->phi =  mytoc[*ki].phi();
      
      if (_is_debug) {


	std::cout << "HLT trig object = (" 
		  << mytoc[*ki].pt() << ", "
		  << mytoc[*ki].eta() << ", "
		  << mytoc[*ki].phi() << ")"
		  << std::endl;
	
      }

    }
  }


  return;

  */

  // copy HLT four vectors. 
  //  const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
  const int index = triggerObj->filterIndex(edm::InputTag(longname));

  if (_is_debug)  std::cout << "index " << index << "  " <<  triggerObj->sizeFilters()<< std::endl;

  if ( index >= triggerObj->sizeFilters() || index <0) { return; }


  const trigger::Keys & k = triggerObj->filterKeys(index);
  for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {

    _vec4_ *hltobj = ptr2AddHlt(aEvent);
    hltobj->id = toc[*ki].id();
    hltobj->pt = toc[*ki].pt();
    hltobj->eta =  toc[*ki].eta();
    hltobj->phi =  toc[*ki].phi();

    if (_is_debug) {


      std::cout << "HLT trig object = (" 
		<< toc[*ki].pt() << ", "
		<< toc[*ki].eta() << ", "
		<< toc[*ki].phi() << ")"
		<< std::endl;

	
    }

  }

  if (_is_debug)  std::cout << " " << std::endl;

  // special case for double muons with different threshold 
  if (  ((std::string( hltname) ).find("Mu") != string::npos) 
	&& ( (std::string( hltname) ).find("Mu") !=  (std::string( hltname) ).rfind("Mu")) 
	) {
    
    longname.resize(0);

    if (trigModules.size()>=3) 
      longname.append( trigModules[ trigModules.size()-3 ]  );
    longname.append("::");
    longname.append( triggerProcessName); 


    const int index_dmuons = triggerObj->filterIndex(edm::InputTag(longname));

    if ( index_dmuons >= triggerObj->sizeFilters() || index_dmuons <0) { return; }
    
    
    const trigger::Keys & k_dmuons = triggerObj->filterKeys(index_dmuons);
    for (trigger::Keys::const_iterator ki = k_dmuons.begin(); ki !=k_dmuons.end(); ++ki ) {
      
      _vec4_ *hltobj = ptr2AddHlt(aEvent);
      hltobj->id = toc[*ki].id();
      hltobj->pt = toc[*ki].pt();
      hltobj->eta =  toc[*ki].eta();
      hltobj->phi =  toc[*ki].phi();

      if (_is_debug) {


	std::cout << "HLT trig object = (" 
		  << toc[*ki].pt() << ", "
		  << toc[*ki].eta() << ", "
		  << toc[*ki].phi() << ")"
		  << std::endl;
	
      }
    }
    if (_is_debug)    std::cout <<"low leg trigger objects " << std::endl;
  }



  // for cross object triggers
  std::vector< std::string > myL3Filters;
  if (  (std::string( hltname) ).find("Mu") != string::npos
	&&  (std::string( hltname) ).find("Ele") != string::npos
	) {


    //!!!!! temperory solution for cross object trigger,
    // only work for HLT_Mu17* blah blah
    for (unsigned int ii = 0; ii < trigModules.size() ; ii++) {
  
      //      if ( hltConfig.saveTags( std::string(trigModules[ii]) ) ) myL3Filters.push_back( std::string( trigModules[ii] ) );

      if ( (std::string(trigModules[ii])).find("Filtered17") != string::npos )  myL3Filters.push_back( std::string( trigModules[ii] ) );
    }




    longname.resize(0);
    //    if( myL3Filters.size()>=2)  longname.append( myL3Filters.at( myL3Filters.size() -2) );
    if( myL3Filters.size()>=1)  longname.append( myL3Filters.at( myL3Filters.size() -1) );
    longname.append("::");
    longname.append( triggerProcessName); 
    const int index_xobj = triggerObj->filterIndex(edm::InputTag(longname));

    if (_is_debug)  std::cout << "additon trigger objects for (cross object triggers) with " 
			      << longname
			      << std::endl;
    if (_is_debug)  std::cout << "index " << index_xobj << "  " <<  triggerObj->sizeFilters()<< std::endl;
    
    if ( index_xobj >= triggerObj->sizeFilters() || index_xobj <0) { return; }
    
    
    const trigger::Keys & k_xobj = triggerObj->filterKeys(index_xobj);
    for (trigger::Keys::const_iterator ki = k_xobj.begin(); ki !=k_xobj.end(); ++ki ) {
      
      _vec4_ *hltobj = ptr2AddHlt(aEvent);
      hltobj->id = toc[*ki].id();
      hltobj->pt = toc[*ki].pt();
      hltobj->eta =  toc[*ki].eta();
      hltobj->phi =  toc[*ki].phi();
    }

  }

  return;
}




void 
WZEdmAnalyzer::copyHLTInfo(const edm::Event& iEvent, _hlt_info_ * hltInfo) 
{

  const char *hltname = "HLT_Mu15";
  hltInfo->HLT_Muon = false;

  if (!triggerResults.isValid() || !triggerObj.isValid()) {
    edm::LogInfo("FourVectorHLTOffline") << "TriggerResults not found, "
					 << "skipping event"; 
    return;
  }


 
  edm::TriggerNames triggerNames = iEvent.triggerNames(*triggerResults);  

  int npath = triggerResults->size();

  for (int i = 0; i < npath; ++i) {
    if (triggerNames.triggerName(i).find(hltname) != std::string::npos && triggerResults->accept(i)) { 
      hltInfo->HLT_Muon = true;
      break;
    }
  }

  // if no HLT_Mu15 fired, return
  if (!hltInfo->HLT_Muon) return;

  // copy HLT four vectors. 
  const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
  if (_is_debug) {

    for (unsigned int ii = 0; ii < triggerObj->sizeFilters(); ii ++) {
      std::cout << triggerObj->filterTag(ii).encode() << std::endl;
    }

  }
  const int index = triggerObj->filterIndex(edm::InputTag("hltSingleMu15L3PreFiltered15::HLT"));

  //  std::cout << "index " << index << "  " <<  triggerObj->sizeFilters()<< std::endl;
  if ( index >= triggerObj->sizeFilters()) {

    return;
  }


  const trigger::Keys & k = triggerObj->filterKeys(index);
  for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {

    _vec4_ *hltmuon = myEvent->addHLTMuon();
    hltmuon->id = toc[*ki].id();
    hltmuon->pt = toc[*ki].pt();
    hltmuon->eta =  toc[*ki].eta();
    hltmuon->phi =  toc[*ki].phi();
  }
}



void WZEdmAnalyzer::copyMET(const edm::Event& iEvent, const char *mettype, _met_ &amet) {

  std::string mymet(mettype); 


  if (!mymet.compare("tcMet")) { // tcmet
    
    iEvent.getByToken(tcMETToken_,             tcMEThandle);
    

    amet.pt  = (tcMEThandle->front() ).et();
    amet.phi = (tcMEThandle->front() ).phi();
    amet.sumEt = (tcMEThandle->front() ).sumEt();
    amet.e_longitudinal = (tcMEThandle->front() ).e_longitudinal();
    amet.HadronicEtFraction = 0;
    //(tcMEThandle->front() ).etFractionHadronic();
    amet.EMEtFraction = 0;
    //  (tcMEThandle->front() ).emEtFraction();
    amet.MuonEtFraction = 0;
    amet.InvisibleEtFraction = 0;

  } else if (!mymet.compare("caloMet") ) {
    
    iEvent.getByToken(caloMETToken_,          muCorrMEThandle);


    amet.pt  = (muCorrMEThandle->front() ).et();
    amet.phi = (muCorrMEThandle->front() ).phi();
    amet.sumEt = (muCorrMEThandle->front() ).sumEt();
    amet.e_longitudinal = (muCorrMEThandle->front() ).e_longitudinal();
    amet.HadronicEtFraction = (muCorrMEThandle->front() ).etFractionHadronic();
    amet.EMEtFraction = (muCorrMEThandle->front() ).emEtFraction();
    amet.MuonEtFraction = 0;
    amet.InvisibleEtFraction = 0;

  } else {



    //    iEvent.getByLabe(mettype,                          pfMEThandle);
    iEvent.getByToken(pfMETToken_,             pfMEThandle);

    amet.pt  = (pfMEThandle->front() ).et();
    amet.phi = (pfMEThandle->front() ).phi();
    amet.sumEt = (pfMEThandle->front() ).sumEt();
    amet.e_longitudinal = (pfMEThandle->front() ).e_longitudinal();
    amet.HadronicEtFraction = (pfMEThandle->front() ).ChargedHadEtFraction() + (pfMEThandle->front() ).NeutralHadEtFraction();
    amet.EMEtFraction = (pfMEThandle->front() ).ChargedEMEtFraction() + (pfMEThandle->front() ).NeutralEMEtFraction();
    
    amet.MuonEtFraction = (pfMEThandle->front() ).MuonEtFraction();
    amet.InvisibleEtFraction = 0;
    
  }


}

/*
void 
WZEdmAnalyzer::copyMETs(_mets_ *mymet) 
{


  
  mymet->pfMET.pt  = (pfMEThandle->front() ).et();
  mymet->pfMET.phi = (pfMEThandle->front() ).phi();
  mymet->pfMET.sumEt = (pfMEThandle->front() ).sumEt();
  mymet->pfMET.e_longitudinal = (pfMEThandle->front() ).e_longitudinal();
  mymet->pfMET.HadronicEtFraction = (pfMEThandle->front() ).ChargedHadEtFraction() + (pfMEThandle->front() ).NeutralHadEtFraction();
  mymet->pfMET.EMEtFraction = (pfMEThandle->front() ).ChargedEMEtFraction() + (pfMEThandle->front() ).NeutralEMEtFraction();

  mymet->pfMET.MuonEtFraction = (pfMEThandle->front() ).MuonEtFraction();
  mymet->pfMET.InvisibleEtFraction = 0;



  mymet->rawMET.pt  = (rawMEThandle->front() ).et();
  mymet->rawMET.phi = (rawMEThandle->front() ).phi();
  mymet->rawMET.sumEt = (rawMEThandle->front() ).sumEt();
  mymet->rawMET.e_longitudinal = (rawMEThandle->front() ).e_longitudinal();
  mymet->rawMET.HadronicEtFraction = (rawMEThandle->front() ).etFractionHadronic();
  mymet->rawMET.EMEtFraction = (rawMEThandle->front() ).emEtFraction();
  mymet->rawMET.MuonEtFraction = 0;
  mymet->rawMET.InvisibleEtFraction = 0;


  mymet->muonJESCorSC5MET.pt  = (muJESCorrMEThandle->front() ).et();
  mymet->muonJESCorSC5MET.phi = (muJESCorrMEThandle->front() ).phi();
  mymet->muonJESCorSC5MET.sumEt = (muJESCorrMEThandle->front() ).sumEt();
  mymet->muonJESCorSC5MET.e_longitudinal = (muJESCorrMEThandle->front() ).e_longitudinal();
  mymet->muonJESCorSC5MET.HadronicEtFraction = (muJESCorrMEThandle->front() ).etFractionHadronic();
  mymet->muonJESCorSC5MET.EMEtFraction = (muJESCorrMEThandle->front() ).emEtFraction();
  mymet->muonJESCorSC5MET.MuonEtFraction = 0;
  mymet->muonJESCorSC5MET.InvisibleEtFraction = 0;


}
*/

void WZEdmAnalyzer::copyBeamSpotInfo(const reco::BeamSpot *aBeamSpot, _beam_spot_ * myBeamSpot) {


  myBeamSpot->x0 = vertexBeamSpot->x0();
  myBeamSpot->y0 = vertexBeamSpot->y0();
  myBeamSpot->z0 = vertexBeamSpot->z0();
  myBeamSpot->x0Error = vertexBeamSpot->x0Error();
  myBeamSpot->y0Error = vertexBeamSpot->y0Error();
  myBeamSpot->z0Error = vertexBeamSpot->z0Error();
  myBeamSpot->sigmaZ = vertexBeamSpot->sigmaZ();
  myBeamSpot->sigmaZError = vertexBeamSpot->sigmaZ0Error();

  myBeamSpot->BeamWidthXError = vertexBeamSpot->BeamWidthXError();
  myBeamSpot->BeamWidthYError = vertexBeamSpot->BeamWidthYError();


}


void 
WZEdmAnalyzer::fillEventInfo(const edm::Event& iEvent,  const edm::EventSetup& iSetup) 
{
  

  /************************************************************************
   *
   * vertex infomation
   *
   ************************************************************************/
  if (vertexBeamSpot)  this->copyBeamSpotInfo(vertexBeamSpot, myEvent->getBeamSpot());


  if (_is_debug) std::cout << "copy vertex information ..." << std::endl;
  int   leadingMuonIndex = -1; 
  int   counter = -1;
  float leadingMuonPt    = -999;
  for (reco::MuonCollection::const_iterator muon = (*recoMuons).begin(); 
       muon != (*recoMuons).end(); 
       muon ++) {

    counter ++;

    //  
    if (!muon->isGlobalMuon()  ) continue;
    if (muon->innerTrack()->numberOfValidHits()<6) continue;
    if (muon->innerTrack()->normalizedChi2() > 10) continue;

    if (muon->pt() > leadingMuonPt) {

      leadingMuonPt = muon->pt(); leadingMuonIndex = counter;
    }
  }


  for(reco::VertexCollection::const_iterator v=recVtxs->begin(); 
      v!=recVtxs->end(); ++v){
    
    if (_is_debug) std::cout << "copying each vertex ... " << std::endl;
    _vertex_ *myVertex = myEvent->addVertex();
    if (myVertex)     this->copyVertexInfo(v, recoMuons, leadingMuonIndex, myVertex);
    else  std::cout << "WARNING: cannot append more vertexs in the output ntuple" << std::endl;

  }


  /***********************************************************************
   *
   *  copy HLT trigger objects
   *  
   *
   ***********************************************************************/
  //void 
  //WZEdmAnalyzer::copyHLTInfo(const edm::Event& iEvent, 
  //		   edm::Handle<TriggerResults>         &triggerResults, 
  //		   std::string &triggerProcessName, 
  //		   edm::Handle<trigger::TriggerEvent>  &triggerObj, 
  //		   HLTConfigProvider &hltConfig,
  //		   bool &trgBit, const char *hltname,  
  //		   _event_ *aEvent,  
  //		   _vec4_ *(*ptr2AddHlt)(_event_ *) ) 
  if (_is_debug) std::cout << "copy HLT information ..." << std::endl;


  if ( HLTTriggerMuons_.size()>=3) {

    if ( HLTTriggerMuons_[3].size() >0)  {
      this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
			myEvent->getHLTInfo()->HLT_MuonH,    
			HLTTriggerMuons_[3].data(),
			//		    HLTTriggerMuonL_.data(),
			//		    longname.data(), 
			myEvent, 
			&addHLTMuonH, 
			&addHLTMuonHleg1, 
			&addHLTMuonHleg2);
    }

  }

  if ( HLTTriggerMuons_[2].size() >0)  {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_MuonL,    
		      HLTTriggerMuons_[2].data(),
		      //		    HLTTriggerMuonL_.data(),
		      //		    longname.data(), 
		      myEvent, 
		      &addHLTMuonL, 
		      &addHLTMuonLleg1, 
		      &addHLTMuonLleg2);
  }


  if (HLTTriggerMuons_[1].size()>0) {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_Muon,    
		      HLTTriggerMuons_[1].data(),
		      //HLTTriggerMuon_.data(),
		    //longname.data(), 
		      myEvent, &addHLTMuon);
  }


  if (HLTTriggerMuons_[0].size()>0) {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_NonIsoMuon,    
		      HLTTriggerMuons_[0].data(),
		      //HLTTriggerMuon_.data(),
		      //longname.data(), 
		      myEvent, &addHLTNonIsoMuon);
  }



//  longname.clear();
//  longname.append(HLTTriggerElectronLLongName_); longname.append("::"); longname.append(TriggerProcess_);
//  this->copyHLTInfo(iEvent, 


  if ( HLTTriggerMuonElectrons_[0].size() >0) {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_MuonElectron,    
		      HLTTriggerMuonElectrons_[0].data(),
		      //		    longname.data(), 
		      myEvent, &addHLTMuonElectron);
  }


  if ( HLTTriggerMuonElectrons_[1].size() >0) {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_ElectronMuon,    
		      HLTTriggerMuonElectrons_[1].data(),
		      //		    longname.data(), 
		      myEvent, &addHLTElectronMuon);
  }






  if (HLTTriggerElectrons_[2].size()>0) {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_ElectronL,    
		      HLTTriggerElectrons_[2].data(),
		      //		    longname.data(), 
		      myEvent, &addHLTElectronL, 
		      &addHLTElectronLleg1, 
		      &addHLTElectronLleg2);
  }


  if (HLTTriggerElectrons_[1].size() >0)  {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_Electron,    
		      HLTTriggerElectrons_[1].data(),
		      //HLTTriggerElectron_.data(),
		      //		    longname.data(), 
		      myEvent, &addHLTElectron);
  }    
 
  
  if (HLTTriggerElectrons_[0].size() >0)  {
    this->copyHLTInfo(iEvent, triggerResults, TriggerProcess_, triggerObj, hltConfig,  
		      myEvent->getHLTInfo()->HLT_NonIsoElectron,    
		      HLTTriggerElectrons_[0].data(),
		      //HLTTriggerElectron_.data(),
		      //		    longname.data(), 
		      myEvent, &addHLTNonIsoElectron);
  }    
 



  this->copyMET(iEvent, "caloMet", myEvent->getMETs()->caloMET);
  // this->copyMET(iEvent, "corMetGlobalMuons", myEvent->getMETs()->caloMET);
  //  this->copyMET(iEvent, "tcMet", myEvent->getMETs()->tcMET);
  this->copyMET(iEvent, "pfMet", myEvent->getMETs()->pfMET);
  /*
  this->copyMET(iEvent, "pfType1CorrectedMet", myEvent->getMETs()->pfType1CorrectedMet);
  this->copyMET(iEvent, "pfType1p2CorrectedMet", myEvent->getMETs()->pfType1p2CorrectedMet);
  this->copyMET(iEvent, "pfMetT0pc", myEvent->getMETs()->pfMetT0pc);
  this->copyMET(iEvent, "pfMetT0pcT1", myEvent->getMETs()->pfMetT0pcT1);
  this->copyMET(iEvent, "pfMetT0pcTxy", myEvent->getMETs()->pfMetT0pcTxy);
  this->copyMET(iEvent, "pfMetT0pcT1Txy", myEvent->getMETs()->pfMetT0pcT1Txy);
  this->copyMET(iEvent, "pfMetT1", myEvent->getMETs()->pfMetT1);
  this->copyMET(iEvent, "pfMetT1Txy", myEvent->getMETs()->pfMetT1Txy);
  */

  //  this->copyMETs(myEvent->getMETs());
  //  this->copyDiLeadingJets(myEvent->getDiLeadingJets());
  //  this->copyTrgBits(iEvent, myEvent->getTrgBits());
  //  std::cout << "seg point 1" << std::endl;


  /************************************************************************
   *
   * PF CHS jet collection
   *
   ************************************************************************/
  if (_is_debug) std::cout << "copy PF CHS jet information ..." << std::endl;
  iEvent.getByToken(puJetIdFlagCHSToken_,puJetIdFlagCHS);
  iEvent.getByToken(puJetIdMvaCHSToken_, puJetIdMvaCHS);

  for (reco::PFJetCollection::const_iterator jet = ( *(pfchsJets.product()) ).begin(); jet != ( *(pfchsJets.product())  ).end(); jet ++) {

    int index = jet - pfchsJets->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfchsJets, index));
      

    double jec = pfchsJetCorr->correction(*jet);
    //    double jec = pfchsJetCorr->correction(*jet, iEvent, iSetup);
    //std::cout << "jet pt = " << jet->pt() << "afer " << std::endl;
  
    if ( jec * jet->pt() < jetMinPt) continue;    
    _jet_ *myjet = myEvent->addPFJet();
    this->copyPFCHSJetInfo(iEvent, iSetup, jet, jetRef, pfchsJetUnc, jec, thePFCHSJetFlavourInfos, myjet);


    int    idflag = (*puJetIdFlagCHS)[ jetRef ];
    myjet->puIDFlag   = idflag;
    myjet->puIDMva    =  (*puJetIdMvaCHS)[ jetRef ];
    myjet->puIDLoose  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose  ) ;
    myjet->puIDMedium = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) ;
    myjet->puIDTight  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight  );



    // have identical copy for the old calo jet collection
    //    _jet_ *mycalojet = myEvent->addCaloJet();
    // this->copyPFCHSJetInfo(iEvent, iSetup, jet, jetRef, pfchsJetUnc, jec, thePFCHSJetFlavourInfos, mycalojet);



    //    mycalojet->puIDFlag   = idflag;
    // mycalojet->puIDMva    =  (*puJetIdMvaCHS)[ jetRef ];
    // mycalojet->puIDLoose  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose  ) ;
    // mycalojet->puIDMedium = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) ;
    // mycalojet->puIDTight  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight  );
  }



  if (_is_debug) std::cout << "copy ak5 PF (no CHS) jet information ..." << std::endl;

  // pileup jet ID
  iEvent.getByToken(puJetIdFlagAK5Token_,puJetIdFlagAK5);
  iEvent.getByToken(puJetIdMvaAK5Token_, puJetIdMvaAK5);
  for (reco::PFJetCollection::const_iterator jet = ( *(ak5pfJets.product()) ).begin(); jet != (  *(ak5pfJets.product()) ).end(); jet ++) {


    // Jun 10, 2014
    int index = jet - ak5pfJets->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(ak5pfJets, index));

      
    double jec = pfJetCorr->correction(*jet);
    //std::cout << "jet pt = " << jet->pt() << "afer " << std::endl;
  
    if ( jec * jet->pt() < jetMinPt) continue;    
    _jet_ *myjet = myEvent->addCaloJet();

    this->copyAK5PFJetInfo(iEvent, iSetup, jet, jetRef, pfJetUnc, jec, theAK5PFJetFlavourInfos, myjet);

    int    idflag = (*puJetIdFlagAK5)[ jetRef ];
    myjet->puIDFlag   = idflag;
    myjet->puIDMva    =  (*puJetIdMvaAK5)[ jetRef ];
    myjet->puIDLoose  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose  ) ;
    myjet->puIDMedium = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) ;
    myjet->puIDTight  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight  );

  }




  /*************************************************************************
   *
   *
   * calo jet collections
   *
   *************************************************************************/
  /*
  if (_is_debug) std::cout << "copy Calo jet information ..." << std::endl;
  for (reco::CaloJetCollection::const_iterator jet = ( *(caloJets.product()) ).begin(); jet != (  *(caloJets.product()) ).end(); jet ++) {

    int index = jet - caloJets->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<CaloJetCollection>(caloJets, index));
      

    double jec = caloJetCorr->correction(*jet);
    
    if ( jec * jet->pt() < jetMinPt) continue;    
    _jet_ *myjet = myEvent->addCaloJet();
    this->copyCaloJetInfo(iEvent, iSetup, jet, jetRef, caloJetUnc, jec, theRecoCaloTag,  myjet);
    
    // additional jet ID information
    //    reco::CaloJetRef calojetref(caloJets, jet - caloJets->begin() );
    reco::JetID jetID = (*jetID_ValueMap_Handle)[jetRef];
    
    // myjet->n60  = jetID.n60Hits;
    myjet->n90  = jetID.n90Hits;
    myjet->fHPD = jetID.fHPD;
    myjet->fRBX = jetID.fRBX;
  }
  */

  /*************************************************************************
   *
   *
   * jpt jet collections
   * for this jet collection we check if it exists in the sample
   *
   *************************************************************************/
  /*
  if (jptJets.isValid() ) {
    if (_is_debug) std::cout << "copy JPT jet information ..." << std::endl;
 
   for (reco::JPTJetCollection::const_iterator jet = ( *(jptJets.product()) ).begin(); jet != (  *(jptJets.product()) ).end(); jet ++) {

      int index = jet - jptJets->begin();
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<JPTJetCollection>(jptJets, index));

    
      double jec = jptJetCorr->correction(*jet);

      
      if ( jec * jet->pt() < jetMinPt) continue;    
      _jet_ *myjet = myEvent->addJPTJet();
      this->copyJPTJetInfo(iEvent, iSetup, jet, jetRef, jptJetUnc, jec, theRecoJPTTag, myjet);


    }
  }
  */

  /*************************************************************************
   *
   *
   * pf jet collections
   *
   *************************************************************************/
  if (_is_debug) std::cout << "copy PF (no CHS) jet information ..." << std::endl;

  // pileup jet ID
  iEvent.getByToken(puJetIdFlagToken_,puJetIdFlag);
  iEvent.getByToken(puJetIdMvaToken_, puJetIdMva);
  for (reco::PFJetCollection::const_iterator jet = ( *(pfJets.product()) ).begin(); jet != (  *(pfJets.product()) ).end(); jet ++) {


    // Jun 10, 2014
    int index = jet - pfJets->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfJets, index));

      
    double jec = pfJetCorr->correction(*jet);
    //std::cout << "jet pt = " << jet->pt() << "afer " << std::endl;
  
    if ( jec * jet->pt() < jetMinPt) continue;    
    _jet_ *myjet = myEvent->addJet();

    this->copyPFJetInfo(iEvent, iSetup, jet, jetRef, pfJetUnc, jec, thePFJetFlavourInfos, myjet);

    if (QGTagsHandleMLP.isValid()) myjet->quarkGluonMLP =   (*QGTagsHandleMLP)[jetRef];
    if (QGTagsHandleLikelihood.isValid())  myjet->quarkGluonLikelihood = (*QGTagsHandleLikelihood)[jetRef];

    // NOTE74
    // pujet ID
    // jet ID
    //    continue;
    int    idflag = (*puJetIdFlag)[ jetRef ];
    myjet->puIDFlag   = idflag;
    myjet->puIDMva    =  (*puJetIdMva)[ jetRef ];
    myjet->puIDLoose  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose  ) ;
    myjet->puIDMedium = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) ;
    myjet->puIDTight  = PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight  );



    //    myjet->mva = (*puJetIdMva)[ jetRef ];
    // std::cout << "mva = " << mva << std::endl;
    // if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose  ) ) std::cout << "pass loose" << std::endl; // loose id etc
    // if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium ) )  std::cout <<"pass medium" << std::endl;
    //if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight  ) ) std::cout <<"pass tight" << std::endl;
  }






  /*************************************************************************
   *
   *
   * muon collection
   *
   *************************************************************************/
  if (_is_debug) std::cout <<"copy muon list ... "<<std::endl;  
  for (reco::MuonCollection::const_iterator muon = (*recoMuons).begin(); 
       muon != (*recoMuons).end(); 
       muon ++) {
    
    myMuon = myEvent->addMuon();
    this->copyMuonInfo(muon, myMuon);
    
  }
  if (!_reco_selection.compare("BTAG") ) return;


  /************************************************************************
   *
   * keep additional objects for analysis purpose
   *
   ************************************************************************/
  if (_is_debug) std::cout << "copy photon information ..." << std::endl;
  for (reco::PhotonCollection::const_iterator photon = (*recoPhotons).begin(); photon != (*recoPhotons).end(); photon ++) {
    
    myPhoton = myEvent->addPhoton(); 
    this->copyPhotonInfo(photon, myPhoton);      
  }



  if (_is_debug) std::cout << "copy supercluster information ..." << std::endl;
  for (reco::SuperClusterCollection::const_iterator supercluster = (*superclusters.product()).begin(); supercluster != (*superclusters.product()).end(); supercluster ++) {
    

    // keep 
    if ( (supercluster->energy() * sin( supercluster->position().theta() ) )<8) continue;
    
    _supercluster_ *mySupercluster = myEvent->addSupercluster(); 
    this->copySuperclusterInfo(supercluster, mySupercluster);      
  }
  



  if (_is_debug) std::cout << "copy gsf electron information ..." << std::endl;
  int index = 0;
  for (reco::GsfElectronCollection::const_iterator electron = (*recoElectrons).begin(); electron != (*recoElectrons).end(); electron ++, index ++) {

    //NOTE74  
    myElectron = myEvent->addElectron(); 
    reco::GsfElectronRef electronRef(electrons, index);
    this->copyElectronInfo(electron, myElectron, electronRef, electron - (*recoElectrons).begin() );      


    // find corresponding calibrated electron; the list could be different 
    int caliIndex = 0;
    bool caliEle = false;
    for (reco::GsfElectronCollection::const_iterator caliElec = calibratedElectrons->begin(); caliElec != calibratedElectrons->end(); caliElec++, caliIndex++) {

      double dist = ROOT::Math::VectorUtil::DeltaR( electron->momentum(), 
						    caliElec->momentum());

      if (dist<0.01) {


	myElectron->calibratedSCEt        = caliElec->superCluster()->energy()/TMath::CosH( caliElec->superCluster()->eta()  );
	myElectron->calibratedSCEta       = caliElec->superCluster()->eta();
	myElectron->calibratedSCPhi       = caliElec->superCluster()->phi();
	myElectron->calibratedEnergy      = caliElec->superCluster()->energy();


	myElectron->calibratedPt          =  caliElec->pt();
	myElectron->calibratedEta         =  caliElec->eta();
	myElectron->calibratedPhi         =  caliElec->phi();
	caliEle = true;
	break;
      }
    }
    if (!caliEle) {

      myElectron->calibratedPt    = 0;
      myElectron->calibratedSCEt  = 0;

      std::cout << "can't find the calibrated electron " 
		<< setw(20) << electron->pt()
		<< setw(20) << electron->eta()
		<< std::endl;
    }
    //    std::cout << setw(20) << (*regEne_handle.product()).get(index) * TMath::Sin( 2* TMath::ATan( TMath::Exp( - ((*calibratedElectrons.product())[index]).eta()  ) ) ) 
    //      << setw(20)  << ((*calibratedElectrons.product())[index]).superCluster()->energy()
    //      << setw(20)  << (*regEne_handle.product()).get(index) 
    //      << std::endl;
    /*
    std::cout << setw(20) << electron->pt() 
	      << setw(20)  << electron->eta() 
	      << setw(20) << electron->phi() 
	      << std::endl;
   std::cout << setw(20) <<  ((*calibratedElectrons.product())[index]).pt()
	      << setw(20)  <<  ((*calibratedElectrons.product())[index]).eta()
	      << setw(20) <<  ((*calibratedElectrons.product())[index]).phi()
	      << std::endl;
    */

    /* NOTE74
    myElectron->calibratedSCEt        = ((*calibratedElectrons.product())[index]).superCluster()->energy()/TMath::CosH(((*calibratedElectrons.product())[index]).superCluster()->eta());
    myElectron->calibratedSCEta       = ((*calibratedElectrons.product())[index]).superCluster()->eta();
    myElectron->calibratedSCPhi       = ((*calibratedElectrons.product())[index]).superCluster()->phi();
    myElectron->calibratedEnergy      = ((*calibratedElectrons.product())[index]).superCluster()->energy();
 

    myElectron->calibratedPt          = ((*calibratedElectrons.product())[index]).pt();
    myElectron->calibratedEta         = ((*calibratedElectrons.product())[index]).eta();
    myElectron->calibratedPhi         = ((*calibratedElectrons.product())[index]).phi();
    myElectron->regressionEnergy      = (*regEne_handle.product()).get(index);
    myElectron->regressionEnergyError = (*regErr_handle.product()).get(index);
    myElectron->mvaTrigV0             = (*mvaTrigV0_handle.product()).get(index);
    myElectron->mvaNonTrigV0          = (*mvaNonTrigV0_handle.product()).get(index);
    */


    /*********** MVA input variables **********************/
    // NOTE74
    /*
    InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
    InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));

    EcalClusterLazyTools myEcalCluster(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);


    bool validKF= false; 
    reco::TrackRef myTrackRef = electron->closestCtfTrackRef();
    validKF = (myTrackRef.isAvailable());
    validKF = (myTrackRef.isNonnull());  
    
    myElectron->validKF = validKF;
    

    // Pure tracking variables
    myElectron->fMVAVar_fbrem           =  electron->fbrem();
    myElectron->fMVAVar_kfchi2          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
    myElectron->fMVAVar_kfhits          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; 
    myElectron->fMVAVar_kfhitsall       =  (validKF) ? myTrackRef->numberOfValidHits() : -1. ;   //  save also this in your ntuple as possible alternative
    myElectron->fMVAVar_gsfchi2         =  electron->gsfTrack()->normalizedChi2();  

  
    // Geometrical matchings
    myElectron->fMVAVar_deta            =  electron->deltaEtaSuperClusterTrackAtVtx();
    myElectron->fMVAVar_dphi            =  electron->deltaPhiSuperClusterTrackAtVtx();
    myElectron->fMVAVar_detacalo        =  electron->deltaEtaSeedClusterTrackAtCalo();
    myElectron->fMVAVar_dphicalo        =  electron->deltaPhiSeedClusterTrackAtCalo();   //  save also this in your ntuple 


    // Pure ECAL -> shower shapes
    myElectron->fMVAVar_see             =  electron->sigmaIetaIeta();    //EleSigmaIEtaIEta
    std::vector<float> vCov = myEcalCluster.localCovariances(*(electron->superCluster()->seed())) ;
    if (!isnan(vCov[2])) myElectron->fMVAVar_spp = sqrt (vCov[2]);   //EleSigmaIPhiIPhi
    else myElectron->fMVAVar_spp = 0.;    
    myElectron->fMVAVar_sigmaIEtaIPhi = vCov[1];  //  save also this in your ntuple 
    
    myElectron->fMVAVar_etawidth        =  electron->superCluster()->etaWidth();
    myElectron->fMVAVar_phiwidth        =  electron->superCluster()->phiWidth();
    myElectron->fMVAVar_e1x5e5x5        =  (electron->e5x5()) !=0. ? 1.-(electron->e1x5()/electron->e5x5()) : -1. ;
    myElectron->fMVAVar_R9              =  myEcalCluster.e3x3(*(electron->superCluster()->seed())) / electron->superCluster()->rawEnergy();
    myElectron->fMVAVar_nbrems          =  fabs(electron->numberOfBrems());    //  save also this in your ntuple 
    
    // Energy matching
    myElectron->fMVAVar_HoE             =  electron->hadronicOverEm();
    myElectron->fMVAVar_EoP             =  electron->eSuperClusterOverP();
    //    myElectron->fMVAVar_IoEmIoP         =  (1.0/(electron->superCluster()->energy())) - (1.0 / electron->p());  // in the future to be changed with electron->momentumAtVtx().R()
    // myElectron->fMVAVar_IoEmIoP         =  (1.0/electron->ecalEnergy()) - (1.0 / electron->p());  // in the future to be changed with electron->gsfTrack().p()   // 24/04/2012 changed to correctly access the corrected supercluster energy from CMSSW_52X
    myElectron->fMVAVar_IoEmIoP         =  (1.0/electron->ecalEnergy()) - (1.0 / electron->p());  // in the future to be changed with electron->momentumAtVtx().R()   // 05/06/2012 changed to correctly access the electron track momentum (mode)
    myElectron->fMVAVar_eleEoPout       =  electron->eEleClusterOverPout();
    myElectron->fMVAVar_PreShowerOverRaw=  electron->superCluster()->preshowerEnergy() / electron->superCluster()->rawEnergy();
    myElectron->fMVAVar_EoPout          =  electron->eSeedClusterOverPout();     //  save also this in your ntuple 


    // Spectators
    myElectron->fMVAVar_eta             =  electron->superCluster()->eta();         
    myElectron->fMVAVar_pt              =  electron->pt();                          
    // for triggering electrons get the impact parameteres
    //d0
    if (electron->gsfTrack().isNonnull()) {
      myElectron->fMVAVar_d0 = (-1.0)*electron->gsfTrack()->dxy((*recVtxs.product())[0].position()); 
    } else if (electron->closestCtfTrackRef().isNonnull()) {
      myElectron->fMVAVar_d0 = (-1.0)*electron->closestCtfTrackRef()->dxy((*recVtxs.product())[0].position()); 
    } else {
      myElectron->fMVAVar_d0 = -9999.0;
    }
    
    //default values for IP3D
    myElectron->fMVAVar_ip3d = -999.0; 
    myElectron->fMVAVar_ip3dSig = 0.0;
    if (electron->gsfTrack().isNonnull()) {
      const double gsfsign   = ( (-electron->gsfTrack()->dxy((*recVtxs.product())[0].position()))   >=0 ) ? 1. : -1.;
      
      const reco::TransientTrack &tt = (*transientTrackBuilder).build(electron->gsfTrack()); 
      const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,(*recVtxs.product())[0]);
      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	double ip3derr = ip3dpv.second.error();  
	myElectron->fMVAVar_ip3d = ip3d; 
	myElectron->fMVAVar_ip3dSig = ip3d/ip3derr;
      }
    }
    */
    /************ end of MVA variables ***************************/
  }

  //  std::cout << std::endl;


  if (_is_debug) std::cout <<"copy track list (tracks with pt>3 GeV) ... "<< std::endl;
  for (reco::TrackCollection::const_iterator track = (*recoTracks).begin(); track != (*recoTracks).end(); track ++) {

    if (track->pt()<3) continue;

    myTrack = myEvent->addTrack();
    this->copyTrackInfo( &(*track), myTrack);

  } 
}






void WZEdmAnalyzer::fillRunInfo(void) {

  // is automatically calculated at the end of each RUN --  units in nb
  // is the precalculated one written in the cfg file -- units is pb
  // myRunInfo = myEvent->getRunInfo();

  // myRunInfo->autoXSec = runInfo->cross_section(); 
  myEvent->getRunInfo()->extXsec  = avgInstLumi; 
  // myRunInfo->filterEff= runInfo->filter_efficiency();
  
  if (_is_debug) {

    LogDebug("Analyzer::RunInfo") << " xsec (pb-1) -- "          << myRunInfo->autoXSec;
    LogDebug("Analyzer::RunInfo") << " external xsec (pb-1) -- " << myRunInfo->extXsec;
    LogDebug("Analyzer::RunInfo") << " filter efficiency -- "    << myRunInfo->filterEff;
  }
}

