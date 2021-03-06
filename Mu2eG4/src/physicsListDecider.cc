//
// Decide which physics list to use.
//
//
// Original author Rob Kutschke
//
// Notes:
// 1) Extract the name of the requested physics list from the
//    config file, instantiate the requested physics list and
//    return a bare pointer to it.
//
// 2) The caller receives the pointer and immediately passes it
//    to G4, which takes ownership of the physics list object.
//    The G4 interface requires a bare pointer.
//
// 3) There are same special names:
//     Minimal - the original Mu2e minimal physics list
//     *_MU2E* Mu2e custom lists
//
// 4) All other names are presumed to be valid names for physics lists that
//    can be created by the PhysListFactory.  At this writing ( April 2010),
//    the PhysListFactory has a default if it does not recognize name as
//    one of its known lists.  The default is QGSP_BERT 3.3.
//
// C++ includes
#include <string>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/DecayMuonsWithSpin.hh"
#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/MinDEDXPhysicsList.hh"
#include "Mu2eG4/inc/StepLimiterPhysConstructor.hh"
#include "Mu2eG4/inc/Mu2eG4CustomizationPhysicsConstructor.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// G4 includes
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4ErrorPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"

#if G4VERSION>4103
#include "G4EmParameters.hh"
#endif

using namespace std;

namespace mu2e{

  G4VUserPhysicsList* physicsListDecider(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug) {

    G4VModularPhysicsList* tmpPL(nullptr);

    const string name = phys.physicsListName();

    debug.diagLevel()>-1 && G4cout << __func__ << " invoked with " << name << G4endl;

    // special cases
    if ( name  == "Minimal" ) {
      return new MinimalPhysicsList();
    }

    else if ( name  == "MinDEDX" ) {
      return new MinDEDXPhysicsList(); // limited EM Processes
    }

    else if ( name  == "ErrorPhysicsList" ) {
      // rather special case of G4VUserPhysicsList for Track Error
      // Propagation, with special Energy Loss implementation
      // (see User's Guide: For Application Developers)
      return new G4ErrorPhysicsList();
    }

    // General case
    else {

      G4PhysListFactory physListFactory;
      physListFactory.SetVerbose(debug.diagLevel());
      tmpPL = physListFactory.GetReferencePhysList(name);

    }

    if ( tmpPL==nullptr ) {
      throw cet::exception("G4CONTROL")
        << "Unable to load physics list named: "
        << name
        << "\n";
    }

    // The modular physics list takes ownership of the StepLimiterPhysConstructor.
    tmpPL->RegisterPhysics( new StepLimiterPhysConstructor() );

    // Mu2e Customizations
    tmpPL->RegisterPhysics( new Mu2eG4CustomizationPhysicsConstructor(&phys, &debug));

    if (phys.turnOffRadioactiveDecay()) {
      tmpPL->RemovePhysics("G4RadioactiveDecay");
    }

    if ( phys.turnOffRadioactiveDecay() && phys.turnOnRadioactiveDecay() ) {
      mf::LogError("Config") << "Inconsistent config";
      G4cout << "Error: turnOnRadioactiveDecay & turnOffRadioactiveDecay on" << G4endl;
      throw cet::exception("BADINPUT")<<" decide on turnOn/OffRadioactiveDecay\n";
    }

    if (phys.turnOnRadioactiveDecay()) {
      tmpPL->RegisterPhysics(new G4RadioactiveDecayPhysics(debug.diagLevel()));
    }

#if G4VERSION>4104

    // Changing MSC model transition energy if requested
    { double mscModelTransitionEnergy(std::numeric_limits<double>::max());  // initializing
      // to something distinct before fetching the requested value if present and only using it then
      if (phys.mscModelTransitionEnergy(mscModelTransitionEnergy)) {
        if (debug.diagLevel()>0) {
          G4cout << __func__
                 << " Changing MscEnergyLimit to "
                 << mscModelTransitionEnergy << " MeV" << G4endl;
        }
        G4EmParameters* emParams = G4EmParameters::Instance();
        emParams->SetMscEnergyLimit(mscModelTransitionEnergy*CLHEP::MeV);
      }
    }

    if ( phys.useEmOption4InTracker() && (name.find("_EMZ") == std::string::npos) ) {
      // assign EmStandard_opt4 to the tracker
      if (debug.diagLevel()>0) {
        G4cout << __func__
               << " Assigning EmStandardPhysics_option4 to the TrackerMother" << G4endl;
      }
      G4EmParameters* emParams = G4EmParameters::Instance();
      emParams->AddPhysics("TrackerMother", "G4EmStandard_opt4");
    }

#endif

    // Muon Spin and Radiative decays plus pion muons with spin
    if ( phys.decayMuonsWithSpin() ) {

      // requires spin tracking: G4ClassicalRK4WSpin
      if ( phys.stepper() != "G4ClassicalRK4WSpin" &&
           phys.stepper() != "G4DormandPrince745WSpin" ) {
        mf::LogError("Config") << "Inconsistent config";
        G4cout << "Error: DecayMuonsWithSpin requires enabling spin tracking" << G4endl;
        throw cet::exception("BADINPUT")<<" DecayMuonsWithSpin requires enabling spin tracking\n";
      }

      tmpPL->RegisterPhysics( new DecayMuonsWithSpin(debug.diagLevel()));
    }

    G4double productionCut = phys.minRangeCut();
    G4double protonProductionCut = phys.protonProductionCut();
    mf::LogInfo("GEOM_MINRANGECUT")
      << "Setting production cut to " << productionCut
      << ", protonProductionCut to " << protonProductionCut << " mm";
    if (debug.diagLevel() > 0) {
      G4cout << __func__ << " Setting gamma, e- and e+ production cut"
             << " to " << productionCut << " mm and for proton to "
             << protonProductionCut << " mm" << G4endl;
    }
    //setCutCmd equivalent:
    tmpPL->SetDefaultCutValue(productionCut);
    tmpPL->SetCutValue(protonProductionCut, "proton");
    // regional cuts (if any) are set during the geometry construction

    if (debug.diagLevel() > 0) tmpPL->DumpCutValuesTable();

    return dynamic_cast<G4VUserPhysicsList*>(tmpPL);

  }

} // end namespace mu2e
