
#include "Print/inc/HelixSeedPrinter.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>

void 
mu2e::HelixSeedPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<HelixSeedCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<HelixSeedCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::HelixSeedPrinter::Print(const art::Handle<HelixSeedCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::HelixSeedPrinter::Print(const art::ValidHandle<HelixSeedCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::HelixSeedPrinter::Print(const HelixSeedCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "HelixSeedCollection has " << coll.size() << " helices\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::HelixSeedPrinter::Print(const art::Ptr<HelixSeed>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::HelixSeedPrinter::Print(const mu2e::HelixSeed& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0 && verbose()==1) os << std::setw(4) << ind;


  //  KalSegment seg;  // this will be filled with 0's and -1's
  // use the middle segment, which should be at the center 
  // of the tracker by default, later there will be proper selectors
  // if( obj.segments().size()>0 ) {
  //   std::size_t i = obj.segments().size()/2 +1;
  //   seg = obj.segments()[i];
  // }

  if(verbose()==1) {

    const mu2e::RobustHelix*robustHel = &obj.helix();

    float chi2xy   = robustHel->chi2dXY();
    float chi2zphi = robustHel->chi2dZPhi();
    
    int    nhits  = obj.hits().size();
    int    flag   = *((int*) &obj.status()); 
    
    double t0     = obj.t0()._t0;
    double t0err  = obj.t0()._t0err;

    double fz0    = robustHel->fz0();
    //    double phi0   = robustHel->phi0();
    double radius = robustHel->radius();
    double d0     = robustHel->rcent() - radius;

    double lambda = robustHel->lambda();
    //      double tandip = lambda/radius;
    
    double mm2MeV = 3/10.;
    double mom    = robustHel->momentum()*mm2MeV;
    double pt     = radius*mm2MeV;

    double x0     = robustHel->centerx();
    double y0     = robustHel->centery();
    
    double clusterEnergy(-1);
    const mu2e::CaloCluster*cluster = obj.caloCluster().get();
    if (cluster != 0) clusterEnergy = cluster->energyDep();

    os << std::dec 
       << " " << std::setw(11) << &obj
       << " " << std::setw( 3) << nhits 
       << " " << std::setw( 8) << std::setprecision(3) << mom
       << " " << std::setw( 8) << std::setprecision(3) << pt 
       << " " << std::setw( 7) << std::setprecision(3) << t0
       << " " << std::setw( 7) << std::setprecision(3) << t0err
       << " " << std::setw( 8) << std::setprecision(3) << d0
       << " " << std::setw( 8) << std::setprecision(3) << fz0
       << " " << std::setw( 8) << std::setprecision(3) << x0
       << " " << std::setw( 8) << std::setprecision(3) << y0
       << " " << std::setw( 8) << std::setprecision(3) << lambda
       << " " << std::setw( 8) << std::setprecision(3) << radius
       << " " << std::setw( 8) << std::setprecision(3) << clusterEnergy
       << " " << std::setw( 8) << std::setprecision(3) << chi2xy
       << " " << std::setw( 8) << std::setprecision(3) << chi2zphi
       << " " << "0x" 
       << std::setw( 8) << std::hex  << std::setfill('0') << flag
       << std::endl;


  } 

}

void 
mu2e::HelixSeedPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::HelixSeedPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind     Address   N     P        pT      T0     T0err  "
     << "   D0       FZ0       X0       Y0    Lambda   radius    ECal    chi2XY  chi2ZPhi    flag\n";
}

