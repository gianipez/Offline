//
// Plugin to read StepPoints in PS Vacuum and create ntuples
//
//  $Id: ReadAntiProtonSteps_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
//  $Author: genser $
//  $Date: 2013/10/21 20:44:04 $
//
// Original author Robert Bernstein
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "GeometryService/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadAntiProtonSteps : public art::EDAnalyzer {
  public:

    typedef vector<int> Vint;
    typedef SimParticleCollection::key_type key_type;
    double pbarMass;
    double pbarMass2;

    explicit ReadAntiProtonSteps(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _psVacuumStepPoints(pset.get<string>("psVacuumStepPoints","AntiProtonSteps")),
      _nAnalyzed(0),
      _maxPrint(pset.get<int>("maxPrint",0)),
      _ntAntiProtonSteps(0),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _stage(pset.get<int>("stage",0))
    {

      Vint const & pdg_ids = pset.get<Vint>("savePDG", Vint());
      if( pdg_ids.size()>0 ) {
        cout << "ReadAntiProtonSteps: save following particle types in the ntuple: ";
        for( size_t i=0; i<pdg_ids.size(); ++i ) {
          pdg_save.insert(pdg_ids[i]);
          cout << pdg_ids[i] << ", ";
        }
        cout << endl;
      }

      nt = new float[1000];
      pbarMass = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::anti_proton).ref().mass().value();
      pbarMass2 = pbarMass*pbarMass;
    }

    virtual ~ReadAntiProtonSteps() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  private:

    // Name of the VD and TVD StepPoint collections
    std::string  _psVacuumStepPoints;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;

    TNtuple* _ntAntiProtonSteps;

    float *nt; // Need this buffer to fill TTree ntpsVacuum

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    int _diagLevel;

    int _stage;
  };

  void ReadAntiProtonSteps::beginJob(){

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _ntAntiProtonSteps = tfs->make<TNtuple>( "ntpbars", "AntiProtonSteps ntuple",
					     "run:evt:trk:pdg:time:ptime:protonEndPx:protonEndPy:protonEndPz:initialX:initialY:initialZ:initialPx:initialPy:initialPz:xsecInitialCosTheta:mu2eInitialCosTheta:endX_s1:endY_s1:endZ_s1:endPx_s1:endPy_s1:endPz_s1:endX_s2:endY_s2:endZ_s2:endPx_s2:endPy_s2:endPz_s2:endX_s3:endY_s3:endZ_s3:endPx_s3:endPy_s3:endPz_s3:endX_s4:endY_s4:endZ_s4:endPx_s4:endPy_s4:endPz_s4:posZ_VD");
  }

  void ReadAntiProtonSteps::beginRun(art::Run const& run){

  }

  void ReadAntiProtonSteps::analyze(const art::Event& event) {

    ++_nAnalyzed;

    if (_diagLevel > 1)
      {
	std::cout << " \n \n \n hi, nAnalyzed = " << _nAnalyzed << std::endl;
      }


    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_psVacuumStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    const SimParticleCollection& simParticlesColl(*simParticles);
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());
    std::string tag = simParticles.provenance()->productDescription().branchName();


    art::Handle<GenParticleCollection> genParticleHandle;
    event.getByLabel(_generatorModuleLabel, genParticleHandle);
    bool haveGenPart = genParticleHandle.isValid();
    GenParticleCollection const& genParticles(*genParticleHandle);
    if ( haveGenPart ) haveGenPart = !(genParticles.empty());

    CLHEP::HepLorentzVector initialProtonFourMomentum(0.,0.,0.,0.);
    //
    // look at the genParticles and see that the geantino with the initial proton is still there.
    if (_diagLevel > 1){
      std::cout << "found gen particle: " << haveGenPart << std::endl;
    }

    // in new style, we don't use the gen particle!
    if (haveGenPart)
      {
	for (auto iGen : genParticles )
	  {
	    if (_diagLevel > 1)
	      {
		std::cout << " particle id " << iGen.pdgId() << " particle momentum " << iGen.momentum() << " position " << iGen.position() << std::endl;
	      }
	    if (iGen.pdgId() == PDGCode::proton)
	      {
		initialProtonFourMomentum = iGen.momentum();
	      }
	  }
      } 

 
    if (_diagLevel > 1){
      std::cout << "about to print size of hits" << std::endl;
      std::cout << " size of hits " << hits.isValid() << " " << hits->size()  << std::endl;
    }

    // Loop over all hits.
    if( hits.isValid() )
      //    if (hits.isValid() && hits->size() ==1) 
      {
	for ( size_t i=0; i< hits->size(); ++i ){ 
	  // Alias, used for readability.
	  const StepPointMC& hit = (*hits)[i];

	  // Get the hit information.
	  const CLHEP::Hep3Vector& pos = hit.position();
	  //const CLHEP::Hep3Vector& mom = hit.momentum();
	  CLHEP::HepLorentzVector startPbarFourMomentum(0.,0.,0.,0.);
	  CLHEP::HepLorentzVector endProtonFourMomentum(0.,0.,0.,0.);

	  // Get track info
	  key_type trackId = hit.trackId();
	  int pdgId = 0;
	  bool goodPDG = false;


	  //CLHEP::HepLorentzVector startingFourMomentum(0.,0.,0.,0.);
	  //CLHEP::Hep3Vector startingPosition(0.,0.,0.);

	  //double currentKE(0.);



	  //print of all information of the SimParticleCollection
	  if (_diagLevel > 1)
	    {
	      std::cout<< "Next collection from "<< tag << std::endl;
	      std::cout << "SimParticleCollection has " << simParticlesColl.size() << " particles\n" << std::endl;  //works
	    }
   


	  if ( haveSimPart ){
	    if( !simParticles->has(trackId) ) {
	      pdgId = 0;
	    } else {
	      SimParticle const& sim = simParticles->at(trackId);
	      //SimParticle prot =  simParticles->at(trackId); 
	      //SimParticle* prot = &prot_temp;
	      pdgId = sim.pdgId();
	      for (auto iPDG : pdg_save)        //look only at antiprotons
		{
		  if (pdgId == iPDG)
		    {
		      goodPDG = true;
		    }
		  if (goodPDG)
		    {
		      //
		      // what was the starting momentum and direction of the track?  
		      // just use angle wrt z as a close-enough estimate
		      //auto originalParticle = sim.originParticle();
		      //startingFourMomentum = originalParticle.startMomentum();
		      //startingPosition     = originalParticle.startPosition();
		      //currentKE = sqrt( mom.mag()*mom.mag() + pbarMass2) - pbarMass;
		      //while (prot.hasParent()) {prot = prot.parent();}
		      
		      if (_diagLevel > 1)
			{
			  std::cout << "_nAnalyzed, pdg, propertime, time, volume " << _nAnalyzed << " " << pdgId << " " << hit.properTime() << " " << hit.time() << " " << hit.volumeId() <<"\n"
				    << std::setprecision(15) << " position" <<   " " << pos.x() << " " << pos.y() << " " << pos.z() 
				    << std::endl;
			}
		      
		  
		    }
		}
	    }
	  }

	  if (goodPDG)//hack to get rid of genParticles
	    {

	      nt[0]  = event.id().run();
	      nt[1]  = event.id().event();
	      nt[2]  = trackId.asInt();
	      nt[3]  = pdgId;
	      nt[4]  = hit.time();
	      nt[5]  = hit.properTime();

	      //to read various VD
	      nt[41] = pos.z();

	      

	      int flag{0};
	      for (const auto& iSim: simParticlesColl){
		//parent of the simParticle
		art::Ptr<SimParticle> const& pptr = iSim.second.parent();
		int pkey = -1;
		if (pptr) pkey = int(pptr.key());
		int mykey = iSim.first.asUint();

		if (pkey == -1 && iSim.second.stoppingCode().name() == "protonInelastic" && iSim.second.pdgId() == 2212){
		  endProtonFourMomentum = iSim.second.endMomentum();
		  nt[6] = endProtonFourMomentum.vect().x();   //components for proton
		  nt[7] = endProtonFourMomentum.vect().y();
		  nt[8] = endProtonFourMomentum.vect().z();
		}

		if (pkey == 1 && iSim.second.stoppingCode().name() == "mu2eProtonInelastic"  && iSim.second.pdgId() == -2212){
		  startPbarFourMomentum = iSim.second.startMomentum();
		  nt[9]  = iSim.second.startPosition().x();   //position of generated pbar
		  nt[10] = iSim.second.startPosition().y();
		  nt[11] = iSim.second.startPosition().z();
		  nt[12] = startPbarFourMomentum.vect().x();   //momentum components for generated pbar
		  nt[13] = startPbarFourMomentum.vect().y();
		  nt[14] = startPbarFourMomentum.vect().z();
		  //nt[12] = startingFourMomentum.vect().x();   
		  //nt[13] = startingFourMomentum.vect().y();
		  //nt[14] = startingFourMomentum.vect().z();
		  //wrt proton direction at the end
		  nt[15] = (endProtonFourMomentum.vect().dot(startPbarFourMomentum.vect()))	/   (endProtonFourMomentum.vect().mag()*startPbarFourMomentum.vect().mag());
		  //wrt mu2e z
		  nt[16] = startPbarFourMomentum.cosTheta();
		}
		//end s1
		if ((mykey>100000 && mykey<200000) && _stage>=1 && iSim.second.stoppingCode().name() == "mu2eKillerVolume"  && iSim.second.pdgId() == pdgId && iSim.second.endVolumeIndex() == 14846){
		  nt[17] = iSim.second.endPosition().x();
		  nt[18] = iSim.second.endPosition().y();
		  nt[19] = iSim.second.endPosition().z();
		  nt[20] = iSim.second.endMomentum().vect().x();
		  nt[21] = iSim.second.endMomentum().vect().y();
		  nt[22] = iSim.second.endMomentum().vect().z();
		}
		//end s2
		if ((mykey>200000 && mykey<300000) && _stage>=2 && iSim.second.stoppingCode().name() == "mu2eKillerVolume"  && iSim.second.pdgId() == pdgId){
		  nt[23] = iSim.second.endPosition().x();
		  nt[24] = iSim.second.endPosition().y();
		  nt[25] = iSim.second.endPosition().z();
		  nt[26] = iSim.second.endMomentum().vect().x();
		  nt[27] = iSim.second.endMomentum().vect().y();
		  nt[28] = iSim.second.endMomentum().vect().z();
		}
		//end s3
		if ((mykey>300000 && mykey<400000) && _stage>=3 && iSim.second.stoppingCode().name() == "mu2eKillerVolume"  && iSim.second.pdgId() == pdgId){
		  nt[29] = iSim.second.endPosition().x();
		  nt[30] = iSim.second.endPosition().y();
		  nt[31] = iSim.second.endPosition().z();
		  nt[32] = iSim.second.endMomentum().vect().x();
		  nt[33] = iSim.second.endMomentum().vect().y();
		  nt[34] = iSim.second.endMomentum().vect().z();
		}
		//end s4
		if ((mykey>400000 && mykey<500000) && _stage>=4 && iSim.second.stoppingCode().name() == "hFritiofCaptureAtRest"  && iSim.second.pdgId() == pdgId){
		  nt[35] = iSim.second.endPosition().x();
		  nt[36] = iSim.second.endPosition().y();
		  nt[37] = iSim.second.endPosition().z();
		  nt[38] = iSim.second.endMomentum().vect().x();
		  nt[39] = iSim.second.endMomentum().vect().y();
		  nt[40] = iSim.second.endMomentum().vect().z();
		}

		if ((mykey>400000 && mykey<500000) && _stage>=4  && iSim.second.pdgId() == pdgId && pos.z()>5000.){
		  nt[38]=hit.momentum().x();
		  nt[39]=hit.momentum().y();
		  nt[40]=hit.momentum().z();
		}

		if(_diagLevel > 0 && (startPbarFourMomentum.cosTheta()<=-0.8 && startPbarFourMomentum.vect().mag()>=3000. && flag ==0)) {
		  std::cout << "key      parent    pdgId       Start  Position                  Start P               CosTheta              EndPosition                End P                  vol   process\n" << std::endl;
		  flag =1;
		}

		if(_diagLevel > 0 && (startPbarFourMomentum.cosTheta()<=-0.8 && startPbarFourMomentum.vect().mag()>=3000.)) {
		  std::cout 
		    << " " << std::setw(7) << mykey
		    << " " << std::setw(7) << pkey
		    << " " << std::setw(8) << iSim.second.pdgId()
		    << " | " << std::setw(8)  << std::setprecision(5) << iSim.second.startPosition().x()
		    << " " << std::setw(8)  << std::setprecision(5) << iSim.second.startPosition().y()
		    << " " << std::setw(8)  << std::setprecision(5) << iSim.second.startPosition().z()
		    << " | " << std::setw(9) << std::setprecision(5) << iSim.second.startMomentum().vect().x()
		    << " " << std::setw(9) << std::setprecision(5) << iSim.second.startMomentum().vect().y()
		    << " " << std::setw(9) << std::setprecision(5) << iSim.second.startMomentum().vect().z()
		    << " | " << std::setw(9) << std::setprecision(5) << iSim.second.startMomentum().vect().cosTheta()
		    << "   "
		    << " | " << std::setw(8)  << std::setprecision(5) << iSim.second.endPosition().x()
		    << " " << std::setw(8)  << std::setprecision(5) << iSim.second.endPosition().y()
		    << " " << std::setw(8)  << std::setprecision(5) << iSim.second.endPosition().z()
		    << " | " << std::setw(9) << std::setprecision(5) << iSim.second.endMomentum().vect().x()
		    << " " << std::setw(9) << std::setprecision(5) << iSim.second.endMomentum().vect().y()
		    << " " << std::setw(9) << std::setprecision(5) << iSim.second.endMomentum().vect().z()
		    << " | " << std::setw(6) << iSim.second.endVolumeIndex()
		    << "  "
		    << " " << std::setiosflags(std::ios::left) << iSim.second.stoppingCode().name() 
		    << std::endl;
		}

	      }//end loop over simParticles

	      if (_stage<4 || pos.z()>5000) _ntAntiProtonSteps->Fill(nt);
	    }//end goodPDG

	} // end loop over hits.
      }
  }

}  // end namespace mu2e

using mu2e::ReadAntiProtonSteps;
DEFINE_ART_MODULE(ReadAntiProtonSteps);
