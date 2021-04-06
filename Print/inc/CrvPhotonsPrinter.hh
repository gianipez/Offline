//
//  Utility class to print CrvPhotons
// 
#ifndef Print_inc_CrvPhotonsPrinter_hh
#define Print_inc_CrvPhotonsPrinter_hh

#include <cstring>
#include <iostream>

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/CrvPhotonsCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class CrvPhotonsPrinter : public ProductPrinter {
  public:

    CrvPhotonsPrinter() { }
    CrvPhotonsPrinter(const ProductPrinter::Config& conf):
      ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<CrvPhotonsCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<CrvPhotonsCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const CrvPhotonsCollection& coll, 
	       std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  };

}
#endif
