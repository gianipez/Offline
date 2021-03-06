
#ifndef ValCaloShowerStep_HH_
#define ValCaloShowerStep_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCaloShowerStep {

  public:
    ValCaloShowerStep(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CaloShowerStepCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _ht;
    TH1D* _hE;
    TH1D* _hposx;
    TH1D* _hposy;
    TH1D* _hposz;
  };
}


#endif
