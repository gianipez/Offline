// New patch to include noise and DAQ jitter effects in time resolution
// 1 batch in 2024-2025: 
// 200 keV electronics, 405 keV RIN+Dark, 200 ps DAQ JITTER
// set aCorr=6.285 bCorr=0.202
// 2 batches in 2025: 
// 200 keV electronics noise, 582 keV RIN+Dark noise, 200 ps DAQ JITTER
// set aCorr=7.911 bCorr=0.208
// Noise levels are from doc-db 35519
// DAQ jitter is from doc-db 35392 assuming ~7 DTCs in chain to include also the tracker (no DAQ jitter in the tracker at the moment)
// Crystal time resolution as function of noise levels is from Bertrand (doc-db in preparation)
// NOTE: by default the correction is switched OFF (aCorr=bCorr=0)

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "SeedService/inc/SeedService.hh"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>


namespace mu2e {


  class CaloCrystalHitFromHit : public art::EDProducer {
  public:

    explicit CaloCrystalHitFromHit(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloDigisToken_{consumes<CaloRecoDigiCollection>(pset.get<std::string>("caloDigisModuleLabel"))},
      time4Merge_        (pset.get<double>     ("time4Merge")),
      aCorr_             (pset.get<double>     ("aCorr")),
      bCorr_             (pset.get<double>     ("bCorr")),
      diagLevel_         (pset.get<int>        ("diagLevel")),
      engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
      randGauss_         (engine_)
    {
      produces<CaloCrystalHitCollection>();
    }

    void beginJob() override;
    void produce(art::Event& e) override;

  private:

    typedef art::Ptr<CaloRecoDigi> CaloRecoDigiPtr;

    art::ProductToken<CaloRecoDigiCollection> const caloDigisToken_;
    double      time4Merge_;
    double      aCorr_;
    double      bCorr_;
    int         diagLevel_;
    CLHEP::HepRandomEngine& engine_;
    CLHEP::RandGaussQ       randGauss_;

    std::vector<std::vector<const CaloRecoDigi*>> hitMap_;  //cached, expensive to create

    TH1F*  hEdep_;
    TH1F*  hTime_;
    TH1F*  hNRo_;
    TH1F*  hEdep_Cry_;
    TH1F*  hDelta_;
    TH2F*  hNRo2_;
    TH1F*  hEdep1_;
    TH1F*  hEdep2_;

    void makeCaloHits(CaloCrystalHitCollection& CaloHits,
                      const art::ValidHandle<CaloRecoDigiCollection>& recoCaloDigisHandle);

    void fillBuffer(int crystalId, int nRoid, double time, double timeErr, double eDep, double eDepErr,
                    std::vector<CaloRecoDigiPtr>& buffer, CaloCrystalHitCollection& caloHits);
  };


  //--------------------------------------------
  void CaloCrystalHitFromHit::beginJob()
  {
    if (diagLevel_ > 2) {
      art::ServiceHandle<art::TFileService> tfs;
      hEdep_     = tfs->make<TH1F>("hEdep",   "Hit energy deposition",        200,   0.,  500);
      hTime_     = tfs->make<TH1F>("hTime",   "Hit time ",                  12000,   0., 2000);
      hNRo_      = tfs->make<TH1F>("hNRo",    "Number RO ",                    10,   0.,   10);
      hEdep_Cry_ = tfs->make<TH1F>("hEdepCry","Energy deposited per crystal",2000,   0., 2000);
      hDelta_    = tfs->make<TH1F>("hDelta",  "Hit time difference",          200, -1,    1);
      hNRo2_     = tfs->make<TH2F>("hNRo2",   "Number RO ",                    5,    0., 5, 50, 0, 50);
      hEdep1_    = tfs->make<TH1F>("hEdep1",  "Hit energy deposition",        200,   0.,  100);
      hEdep2_    = tfs->make<TH1F>("hEdep2",  "Hit energy deposition",        200,   0.,  100);
    }
  }


  //------------------------------------------------------------
  void CaloCrystalHitFromHit::produce(art::Event& event)
  {
    auto const& recoCaloDigisHandle = event.getValidHandle(caloDigisToken_);
    auto caloHits = std::make_unique<CaloCrystalHitCollection>();
    makeCaloHits(*caloHits, recoCaloDigisHandle);

    event.put(std::move(caloHits));
  }


  //--------------------------------------------------------------------------------------------------------------
  void CaloCrystalHitFromHit::makeCaloHits(CaloCrystalHitCollection& caloHits,
                                           art::ValidHandle<CaloRecoDigiCollection> const& recoCaloDigisHandle)
  {
    Calorimeter const &cal = *(GeomHandle<Calorimeter>());
    auto const& recoCaloDigis = *recoCaloDigisHandle;
    CaloRecoDigi const* base = &recoCaloDigis.front(); // What if recoCaloDigis is empty?


    //extend hitMap_ if needed and clear it
    if (cal.nRO() > int(hitMap_.size()))
      for (int i = hitMap_.size(); i<= cal.nRO(); ++i) hitMap_.push_back(std::vector<const CaloRecoDigi*>());
    for (size_t i=0; i<hitMap_.size(); ++i) hitMap_[i].clear();


    // fill the map that associate for each crystal the corresponding CaloRecoDigi indexes
    for (unsigned int i=0; i< recoCaloDigis.size(); ++i)
      {
        int crystalId = cal.caloInfo().crystalByRO(recoCaloDigis[i].ROid());
        hitMap_[crystalId].push_back(&recoCaloDigis[i]);
      }


    for (unsigned int crystalId=0;crystalId<hitMap_.size();++crystalId)
      {
        std::vector<const CaloRecoDigi*> &hits = hitMap_[crystalId];
        if (hits.empty()) continue;

        std::sort(hits.begin(),hits.end(),[](const auto a, const auto b){return a->time() < b->time();});

        auto startHit = hits.begin();
        auto endHit   = hits.begin();

        std::vector<CaloRecoDigiPtr> buffer;
        double timeW(0);//,timeWtot(0);
        double eDepTot(0),eDepTotErr(0);
        int nRoid(0);

        while (endHit != hits.end())
          {
            double deltaTime = (*endHit)->time()-(*startHit)->time();
            if (diagLevel_ > 2) hDelta_->Fill(deltaTime);
		
            if (deltaTime > time4Merge_)
              {
                //double time = timeW/timeWtot;
                //double timeErr = 1.0/sqrt(timeWtot);
                double time = timeW/nRoid;
		// PATCH TO SIMULATE NOISE CONTRIBUTION AND DAQ JITTER (for 7 DTCs considering also the tracker)
		double sigmatcorr=sqrt(aCorr_*aCorr_/eDepTot/eDepTot*nRoid*nRoid+bCorr_*bCorr_);
		time += randGauss_.fire(0.0,sigmatcorr);
                double timeErr = 0;

                fillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, caloHits);

                buffer.clear();
                timeW      = 0.0;
                //timeWtot   = 0.0;
                eDepTot    = 0.0;
                eDepTotErr = 0.0;
                nRoid      = 0;
                startHit   = endHit;
              }
            else
              {
                //double wt  = 1.0/(*endHit)->timeErr()/(*endHit)->timeErr();
                //timeWtot   += wt;
                //timeW      += wt*(*endHit)->time();

                timeW      += (*endHit)->time();

                eDepTot    += (*endHit)->energyDep();
                eDepTotErr += (*endHit)->energyDepErr() * (*endHit)->energyDepErr();

                ++nRoid;

                size_t index = *endHit - base;
                buffer.push_back(art::Ptr<CaloRecoDigi>(recoCaloDigisHandle, index));

                ++endHit;
              }

          }

        //flush last buffer

        //double time = timeW/timeWtot;
        //double timeErr = 1.0/sqrt(timeWtot);
        double time = timeW/nRoid;
        double timeErr = 0;
	// PATCH TO SIMULATE NOISE CONTRIBUTION AND DAQ JITTER (for 7 DTCs considering also the tracker)
	double sigmatcorr=sqrt(aCorr_*aCorr_/eDepTot/eDepTot*nRoid*nRoid+bCorr_*bCorr_);
	time += randGauss_.fire(0.0,sigmatcorr);

        fillBuffer(crystalId, nRoid, time, timeErr, eDepTot/nRoid, eDepTotErr/nRoid, buffer, caloHits);
      }


    if ( diagLevel_ > 0 )
      {
        printf("[CaloCrystalHitFromHit::produce] produced RecoCrystalHits ");
        printf(": caloHits.size()  = %i \n", int(caloHits.size()));
      }

  }

  //--------------------------------------------------------------------------------------------------------------
  void CaloCrystalHitFromHit::fillBuffer(int const crystalId,
                                         int const nRoid,
                                         double const time,
                                         double const timeErr,
                                         double const eDep,
                                         double const eDepErr,
                                         std::vector<CaloRecoDigiPtr>& buffer,
                                         CaloCrystalHitCollection& caloHits)
  {
    caloHits.emplace_back(CaloCrystalHit(crystalId, nRoid, time, timeErr, eDep, eDepErr, buffer));

    if (diagLevel_ > 1)
      {
        std::cout<<"[CaloCrystalHitFromHit] created hit in crystal id="<<crystalId<<"\t with time="<<time<<"\t eDep="<<eDep<<"\t  from "<<nRoid<<" RO"<<std::endl;

        if (diagLevel_ > 2)
          {
            hTime_->Fill(time);
            hEdep_->Fill(eDep);
            hNRo_->Fill(nRoid);
            hEdep_Cry_->Fill(crystalId,eDep);
            hNRo2_->Fill(nRoid,eDep);
            if (nRoid==1) hEdep1_->Fill(eDep);
            if (nRoid==2) hEdep2_->Fill(eDep);
          }
      }
  }




}

DEFINE_ART_MODULE(mu2e::CaloCrystalHitFromHit);
