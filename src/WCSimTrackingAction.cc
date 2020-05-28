#include "WCSimTrackingAction.hh"
#include "WCSimTrajectory.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "WCSimTrackInformation.hh"
#include "WCSimTrackingMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


WCSimTrackingAction::WCSimTrackingAction()
{

  ProcessList.insert("Decay") ;                         // Michel e- from pi+ and mu+
  ProcessList.insert("conv") ;                         // Products of gamma conversion

  //ProcessList.insert("muMinusCaptureAtRest") ;          // Includes Muon decay from K-shell: for Michel e- from mu0. This dominates/replaces the mu- decay (noticed when switching off this process in PhysicsList)                                                   // TF: IMPORTANT: ONLY USE FROM G4.9.6 onwards because buggy/double counting before.
  ////////// ToDo: switch ON the above when NuPRISM uses G4 >= 4.9.6
  ProcessList.insert("nCapture");

//   ProcessList.insert("conv");
  ParticleList.insert(111); // pi0
  ParticleList.insert(211); // pion+
  ParticleList.insert(-211);
  ParticleList.insert(321);
  ParticleList.insert(-321); // kaon-
  ParticleList.insert(311); // kaon0
  ParticleList.insert(-311); // kaon0 bar
  // don't put gammas there or there'll be too many

  //TF: add protons and neutrons
  ParticleList.insert(2212);
  ParticleList.insert(2112);

  percentageOfCherenkovPhotonsToDraw = 0.0;

  messenger = new WCSimTrackingMessenger(this);
}

WCSimTrackingAction::~WCSimTrackingAction(){;}

void WCSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //TF: userdefined now
  //G4float percentageOfCherenkovPhotonsToDraw = 100.0;
  // if larger than zero, will keep trajectories of many secondaries as well
  // and store them in output file. Difficult to control them all, so best only
  // use for visualization, not for storing in ROOT.



  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()
       || G4UniformRand() < percentageOfCherenkovPhotonsToDraw/100. )
    {
      WCSimTrajectory* thisTrajectory = new WCSimTrajectory(aTrack);
      fpTrackingManager->SetTrajectory(thisTrajectory);
      fpTrackingManager->SetStoreTrajectory(true);
    }
  else 
    fpTrackingManager->SetStoreTrajectory(false);
}

void WCSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // added by M Fechner
  const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
  //  if ( creatorProcess )


  WCSimTrackInformation* anInfo;
  if (aTrack->GetUserInformation())
    anInfo = (WCSimTrackInformation*)(aTrack->GetUserInformation());   //eg. propagated to all secondaries blelow.
  else anInfo = new WCSimTrackInformation();
  /** TF's particle list (ToDo: discuss/converge)
  
   // is it a primary ?
  // is the process in the set ? eg. Michel e-, but only keep e-
  // is the particle in the set ? eg. pi0, pi+-, K, p, n
  // is it a gamma above 50 MeV ? OR gamma from nCapture? ToDo: Oxygen de-excitation
  // is it a muon that can still produce Cherenkov light?
  // due to lazy evaluation of the 'or' in C++ the order is important  
  if( aTrack->GetParentID()==0 ||                                                                            // Primary
      ((creatorProcess!=0) && ProcessList.count(creatorProcess->GetProcessName())                            
       && aTrack->GetMomentum().mag() > .5*MeV && abs(aTrack->GetDefinition()->GetPDGEncoding()) == 11) ||   
      (ParticleList.count(aTrack->GetDefinition()->GetPDGEncoding()) )                                       
      || (aTrack->GetDefinition()->GetPDGEncoding()==22 && aTrack->GetVertexKineticEnergy() > 4.*MeV)       // Bugfixed: need vertex kinetic energy for gamma, rest is zero. ToDo: check whether pi0 gamma's are saved!
      || (aTrack->GetDefinition()->GetPDGEncoding()==22 && creatorProcess->GetProcessName() == "nCapture")   
      || (abs(aTrack->GetDefinition()->GetPDGEncoding()== 13) && aTrack->GetMomentum().mag() > 110.0*MeV) //mu+- above Cherenkov Threshold in water (119 MeV/c)

  **/
    // TF: Currently use the nuPRISM one

    // is it a primary ?
    // is the process in the set ? 
    // is the particle in the set ?
    // is it a gamma 
    // due to lazy evaluation of the 'or' in C++ the order is important
    if( aTrack->GetParentID()==0 
	|| ((creatorProcess!=0) && ProcessList.count(creatorProcess->GetProcessName()))
	|| (ParticleList.count(aTrack->GetDefinition()->GetPDGEncoding()))
	|| (aTrack->GetDefinition()->GetPDGEncoding()==22 && aTrack->GetTotalEnergy() > 1.0*MeV)
      || (creatorProcess->GetProcessName() == "muMinusCaptureAtRest" && aTrack->GetTotalEnergy() > 1.0*MeV)
      )
    {
    // if so the track is worth saving
    anInfo->WillBeSaved(true);
    
    //      G4cout << "track # " << aTrack->GetTrackID() << " is worth saving\n";
   
	 G4cout << "It is a " <<aTrack->GetDefinition()->GetParticleName() << G4endl;


    //For Ch hits: use Parent ID of actual mother process:
    // Decay: keep both decaying particle and Michel e-, Hit Parent ID should be Track ID from Michel e-
    // Pi0 : keep both pi0 and gamma's, Hit Parent ID should be Track ID from gamma
    // nCapture: keep both n and gamma, Hit Parent ID should be Track ID from gamma.
    // Hits from LE electrons from muIonization, etc. : Hit Parent ID should be from Mother particle, not secondary track ID.
    
    if (aTrack->GetDefinition()->GetPDGEncoding()==111)
      pi0List.insert(aTrack->GetTrackID()); // list of all pi0-s 
    
    // Be careful with gamma's. I want the ones closest to the actual mother process, not all secondaries.
    if (aTrack->GetDefinition()->GetPDGEncoding() == 22){
      // also use lazy evaluation of "or" here:
      if( aTrack->GetParentID() == 0  || // then this gamma has no creator process (eg. nRooTracker particles)
	  pi0List.count(aTrack->GetParentID()) ||
	  (creatorProcess->GetProcessName() == "nCapture") ||
	  (creatorProcess->GetProcessName() == "NeutronInelastic")	  
	  )
    anInfo->SetPrimaryParentID(aTrack->GetTrackID());  
    }
    //TF: crucial bugfix: I want this for all tracks that I save to match Ch hits with tracks that can
    // produce Cherenkov light.
    else
      anInfo->SetPrimaryParentID(aTrack->GetTrackID());
      
	}
  else {
    anInfo->WillBeSaved(false);
	}


  G4Track* theTrack = (G4Track*)aTrack;
  theTrack->SetUserInformation(anInfo);




 //changed from here
 if( aTrack -> GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
/*	G4ThreeVector PrePosition = aTrack -> GetStep() -> GetPreStepPoint() -> GetPosition();
	std::cout << "PrePosition_x " << PrePosition(0)/cm << std::endl;
	std::cout << "PrePosition_z " << PrePosition(2)/cm << std::endl;
	G4double preposition_r = sqrt(pow(PrePosition(0)/cm,2)+pow(PrePosition(2)/cm,2));
        std::cout << "radius of Pre step position  = " << preposition_r/cm << std::endl;

	G4String materialName = aTrack ->GetMaterial()-> GetName();
	std::cout << "PrePositionMaterial - " << materialName << std::endl;
	

	G4ThreeVector PostPosition = aTrack -> GetStep() -> GetPostStepPoint() -> GetPosition();
	std::cout << "PostPosition_x " << PostPosition(0)/cm << std::endl;
	std::cout << "PostPosition_z " << PostPosition(2)/cm << std::endl;
	G4double postposition_r = sqrt(pow(PostPosition(0)/cm,2)+pow(PostPosition(2)/cm,2));
        std::cout << "radius of Post step position  = " << postposition_r/cm << std::endl;

	G4String PostPositionName = aTrack -> GetMaterial() -> GetName();
	std::cout << "PostPositionMaterial - " << PostPositionName << std::endl;


	G4ThreeVector PrePosition_y = aTrack -> GetStep() -> GetPreStepPoint() -> GetPosition(1);
	G4ThreeVector PrePosition_z = aTrack -> GetStep() -> GetPreStepPoint() -> GetPosition(2);

	G4double preposition_r = sqrt(pow(PrePosition_x,2)+pow(PrePosition_z,2));
	std::cout << "radius of Pre step position  = " << preposition_r/cm << std::endl;
			
	G4ThreeVector PostPosition_x = aTrack -> GetStep() -> GetPostStepPoint() -> GetPosition(0);
	G4ThreeVector PostPosition_y = aTrack -> GetStep() -> GetPostStepPoint() -> GetPosition(1);
	G4ThreeVector PostPosition_z = aTrack -> GetStep() -> GetPostStepPoint() -> GetPosition(2);

	G4double postposition_r = sqrt(pow(PostPosition_x,2)+pow(PostPosition_z,2));
	std::cout << "radius of Post step position = " << postposition_r/cm << std::endl;
		
	G4double stepLength = aTrack -> GetStep() -> GetStepLength();
	std::cout << "Step length = " << stepLength << std::endl;
*/
	
	//G4String processName = aTrack ->GetStep() -> Get 
	}
 
 //till here





 // pass primary parent ID to children
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      { 
	WCSimTrackInformation* infoSec = new WCSimTrackInformation(anInfo);
	if(anInfo->isSaved()){ // Parent is primary, so we want start pos & time of this secondary
        infoSec->SetPhotonStartTime((*secondaries)[i]->GetGlobalTime());
        infoSec->SetPhotonStartPos((*secondaries)[i]->GetPosition());
	}
	infoSec->WillBeSaved(false); // ADDED BY MFECHNER, temporary, 30/8/06
	(*secondaries)[i]->SetUserInformation(infoSec);
      }
    } 
  }
  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() )
    //   if (aTrack->GetDefinition()->GetPDGCharge() == 0) 
  {
    WCSimTrajectory *currentTrajectory = 
      (WCSimTrajectory*)fpTrackingManager->GimmeTrajectory();

    G4ThreeVector currentPosition      = aTrack->GetPosition();
    G4VPhysicalVolume* currentVolume   = aTrack->GetVolume();
	
    currentTrajectory->SetStoppingPoint(currentPosition);
    currentTrajectory->SetStoppingVolume(currentVolume);

    if (anInfo->isSaved())
      currentTrajectory->SetSaveFlag(true);// mark it for WCSimEventAction ;
    else currentTrajectory->SetSaveFlag(false);// mark it for WCSimEventAction ;
	std::cout << "hey again" << std::endl;
  }
}





