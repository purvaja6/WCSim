#include <stdlib.h>
#include <stdio.h>

#include "WCSimSteppingAction.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "WCSimWCSD.hh"
#include "G4UnitsTable.hh"

G4int WCSimSteppingAction::n_photons_through_mPMTLV = 0;
G4int WCSimSteppingAction::n_photons_through_acrylic = 0;
G4int WCSimSteppingAction::n_photons_through_gel = 0;
G4int WCSimSteppingAction::n_photons_on_blacksheet = 0;
G4int WCSimSteppingAction::n_photons_on_smallPMT = 0;


void WCSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	//DISTORTION must be used ONLY if INNERTUBE or INNERTUBEBIG has been defined in BidoneDetectorConstruction.cc
	WCSimWCSD *pSD = WCSimWCSD::aSDPointer;// me:for logicreflector
	const G4Track* track       = aStep->GetTrack();
	//const G4VProcess* creatorProcess = track->GetCreatorProcess();

	// Not used:
	//const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
	//G4VPhysicalVolume* volume  = track->GetVolume();
	//G4SDManager* SDman   = G4SDManager::GetSDMpointer();
	//G4HCofThisEvent* HCE = evt->GetHCofThisEvent();


	// Debug for photon tracking
	G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
	//me:changed from here
	G4ThreeVector PrePosition = aStep->GetPreStepPoint()->GetPosition();
	// std::cout <<"PrePosition_x " << PrePosition(0) << std::endl;
	// std::cout << "PrePosition_z " << PrePosition(2) << std::endl;

	G4TouchableHandle theTouchable = thePrePoint->GetTouchableHandle(); //me:will tell you where in the geometry you are
	G4String VolumeName = theTouchable->GetVolume()->GetLogicalVolume()->GetName();
	G4ParticleDefinition *particleDefinition =    aStep->GetTrack()->GetDefinition();


	if(VolumeName == "reflectorCone" && particleDefinition == G4OpticalPhoton::OpticalPhotonDefinition())
	{
		bool tagflag=true;

		if(tagflag)
		{
			pSD->reflectortag.push_back(track->GetTrackID());
		}

		for(unsigned ij=0; ij < pSD->reflectortag.size(); ij++) //me: this will check the trackID of all the photons in an event (of any particle) 
		{
			if(track->GetTrackID() == pSD->reflectortag[ij]) 
			{
				tagflag=false;
				break;
			}		

		}
	}	

	G4double PrePosition_x = PrePosition(0) + 314.159;
	G4double PrePosition_z = PrePosition(2) + 314.159;
	G4double preposition_r = sqrt(pow(PrePosition_x,2) + pow(PrePosition_z,2));

	G4ThreeVector PreMomentum = aStep -> GetPreStepPoint() -> GetMomentum();
	G4ThreeVector PreMomentumDirection = aStep -> GetPreStepPoint() -> GetMomentumDirection();
	G4double X = PreMomentumDirection(0);
	G4double Y = PreMomentumDirection(1);
	G4double Z = PreMomentumDirection(2);
	//std::cout << "radius of Pre step position  = " << preposition_r << std::endl;
	//till here 

	G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
	//me: changed from here
	//G4Material* PreMaterial = thePrePoint->GetMaterial();
	//std::cout << "PreStepPointMaterial " << thePrePoint->GetMaterial()->GetName() << std::endl;
	//std::cout << "PreStepPointPosition " << thePrePoint<< std::endl;
	//me:till here

	G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
	//me:changed from here
	G4ThreeVector PostPosition = aStep->GetPostStepPoint()->GetPosition();  
	//std::cout << "PostPosition_x " << PostPosition(0) << std::endl;
	//std::cout << "PostPosition_z " << PostPosition(2) << std::endl;
	G4double PostPosition_x = PostPosition(0) + 314.159;
	G4double PostPosition_z = PostPosition(2) + 314.159;
	G4double postposition_r = sqrt(pow(PostPosition_x,2)+pow(PostPosition_z,2));
	G4ThreeVector PostMomentum = aStep -> GetPostStepPoint() -> GetMomentum();
	G4ThreeVector PostMomentumDirection = aStep -> GetPostStepPoint() -> GetMomentumDirection();
	G4double x = PostMomentumDirection(0);
	G4double y = PostMomentumDirection(1);
	G4double z = PostMomentumDirection(2);

	//std::cout << "radius of Post step position  = " << postposition_r << std::endl;

	//std::cout << "PostStepPointMaterial " << thePostPoint->GetMaterial()->GetName() << std::endl;
	//  std::cout << "PostStepPointPosition " << thePostPoint << std::endl;
	//if(creatorProcess != NULL)
	//{
	//G4String  PreProcessName = creatorProcess->GetProcessName();

	//	G4ProcessType  PreProcessType = creatorProcess->GetProcessType();
	//	std::cout << "PreStepProcessName " << PreProcessName << std::endl;
	//	std::cout << "PreStepProcessType " << PreProcessType << std::endl;
	//	G4String  PostProcessName = creatorProcess->GetProcessType();
	//	G4ProcessType  PostProcessType = creatorProcess->GetProcessType();
	//	std::cout << "PostStepProcessName " << PostProcessName << std::endl;
	//	std::cout << "PostStepProcessType " << PostProcessType << std::endl;
	//}
	/*		if(preposition_r >= 49.5983 && preposition_r <=50.4975)
			{
			if(thePrePoint->GetMaterial()->GetName()=="SilGel" && thePostPoint->GetMaterial()->GetName()=="SilGel" )
			{
			*/
	//	if((preposition_r >= 49.5983 && preposition_r <=50.85 )||( postposition_r >= 49.5983 && postposition_r <= 50.85))
	//	{
	G4ThreeVector vertex_position = aStep -> GetTrack()->GetVertexPosition();
	G4ThreeVector vertex_direction = aStep -> GetTrack() -> GetVertexMomentumDirection();
	G4double vertex_x = vertex_position(0) + 314.159*CLHEP::mm;
	G4double vertex_y = vertex_position(1);
	G4double vertex_z = vertex_position(2) + 314.159*CLHEP::mm;
	G4double vertex_r = sqrt(pow(vertex_x,2)+pow(vertex_z,2));
	//aStep -> GetTrack()->  SetVertexPosition(vertex_x, vertex_y, vertex_z);

	//		std::cout << "vertex_position_x = "  << vertex_x << std::endl;
	//		std::cout << "vertex_position_y = "  << vertex_y << std::endl;
	//		std::cout << "vertex_position_z = "  << vertex_z << std::endl;
	//		std::cout << "vertex_r = " << vertex_r << std::endl;
	//		std::cout << "vertex_direction = " << vertex_direction << std::endl;
	/*	if(vertex_r >= 49.5983 && vertex_r <= 50.4975)
		{
		std::cout << "PreStep Radius  = " << preposition_r << std::endl;
		std::cout << "PreStep Material " << thePrePoint->GetMaterial()->GetName() << std::endl;
		std::cout << "PreStep Momentum = " << PreMomentum << std::endl;
		std::cout << "PreStep MomentumDirection = " << PreMomentumDirection << std::endl;


		std::cout << "PostStep Radius  = " << postposition_r << std::endl;
		std::cout << "PostStep Material " << thePostPoint->GetMaterial()->GetName() << std::endl;
		std::cout << "PostStep Momentum = " << PostMomentum << std::endl;
		std::cout << "PostStep MomentumDirection = " << PostMomentumDirection << std::endl;
		G4double angle = acos((x*X + y*Y + z*Z)/((sqrt(x*x + y*y + z*z)) * (sqrt(X*X + Y*Y + Z*Z))));
		std::cout << " angle_calculated = " << angle << std::endl;
		G4double Angle = PostMomentumDirection.angle(PreMomentumDirection);
		std::cout << " angle_function = " << Angle << std::endl;
	//		aStep ->GetTrack() ->  SetTrackStatus(fStopAndKill);		

	}
	//	}
	*/	//till here  

	//G4OpBoundaryProcessStatus boundaryStatus=Undefined;
	//static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;  //doesn't work and needs #include tls.hh from Geant4.9.6 and beyond
	G4OpBoundaryProcess* boundary=NULL;
	/*if((preposition_r >= 49.5983 && preposition_r <=50.85 )||( postposition_r >= 49.5983 && postposition_r <= 50.85))
	  {
	  G4ThreeVector vertex_position = aStep -> GetTrack()->GetVertexPosition();
	  G4ThreeVector vertex_direction = aStep -> GetTrack() -> GetVertexMomentumDirection();
	  std::cout << "vertex_position = "  << vertex_position << std::endl;
	  std::cout << "vertex_direction = " << vertex_direction << std::endl; 
	  */
	if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition()){ //&& energyDeposition == 0.0) {

	} else {
		//G4Track* aTrack = aStep->GetTrack();
		if(aStep -> GetTrack()->GetParentID()==0) {
			//std::cout<<"-x-x-x-x-x-x-x-x-x-"<<std::endl;
			//std::cout<<"StepNum = "<<track->GetCurrentStepNumber()<<" volume "<<VolumeName<<" momDir "<<track->GetMomentumDirection() << " KE "<< G4BestUnit(track->GetKineticEnergy(),"Energy") << std::endl;
			const  G4VProcess* proc1 = thePostPoint->GetProcessDefinedStep();
			if(proc1) {
				//	std::cout<<"Process Defined step "<<proc1->GetProcessName()<<std::endl;
			}
		} // if(aTrack->GetParentID()==0) {
	}

	if(!boundary)
	{

		G4ProcessManager* pm
			= aStep->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		G4int i;
		for( i=0;i<nprocesses;i++)
		{

			if( (*pv)[i]->GetProcessType()==fOptical && (*pv)[i]->GetProcessSubType() == 32 )
			{
				G4cout<<" ProcessName:" << (*pv)[i]->GetProcessName() <<G4endl;				
				boundary = (G4OpBoundaryProcess*)(*pv)[i];
				G4cout<<"  G4OpBoundaryProcessStatus:" << boundary->GetStatus() <<G4endl;
				//		if(boundary-> GetStatus() ==4)
				//		{			
				//		if(vertex_r >= 49.5983 && vertex_r <= 50.4975)
				//			if(thePrePoint->GetMaterial()->GetName()=="SilGel" && thePostPoint->GetMaterial()->GetName()=="Air")

				//			{
				std::cout << "vertex_x = " <<  vertex_x << std::endl;
				std::cout << "vertex_y = " <<  vertex_y << std::endl;
				std::cout << "vertex_z = " <<  vertex_z << std::endl;

				//	if(vertex_x>=-340 && vertex_x<=-290 && vertex_z>=-345 && vertex_z<=-315)
				//	{							
				std::cout << "preposition_x = " << PrePosition_x << std::endl;
				std::cout << "preposition_y = " << PrePosition(1) + 2838.02*CLHEP::mm << std::endl;
				std::cout << "preposition_z = " << PrePosition_z << std::endl;
				std::cout << "PreStep Radius  = " << preposition_r << std::endl;
				std::cout << "PreStep Material " << thePrePoint->GetMaterial()->GetName() << std::endl;
				std::cout << "PreStep Momentum = " << PreMomentum << std::endl;
				std::cout << "PreStep MomentumDirection = " << PreMomentumDirection << std::endl;

				std::cout << "postposition_x = " << PostPosition_x << std::endl;
				std::cout << "postposition_y = " << PostPosition(1) + 2838.02*CLHEP::mm << std::endl;
				std::cout << "postposition_z = " << PostPosition_z << std::endl;
				std::cout << "PostStep Radius  = " << postposition_r << std::endl;
				std::cout << "PostStep Material " << thePostPoint->GetMaterial()->GetName() << std::endl;
				std::cout << "PostStep Momentum = " << PostMomentum << std::endl;
				std::cout << "PostStep MomentumDirection = " << PostMomentumDirection << std::endl;
				//G4double angle = acos((x*X + y*Y + z*Z)/((sqrt(x*x + y*y + z*z)) * (sqrt(X*X + Y*Y + Z*Z))));
				//std::cout << " angle_calculated = " << angle << std::endl;
				G4double Angle = PostMomentumDirection.angle(PreMomentumDirection);
				std::cout << " angle_function = " << Angle << std::endl;
				G4double steplength = aStep -> GetStepLength();
				std::cout << "steplength = " << steplength << std::endl;
				//					aStep ->GetTrack() ->  SetTrackStatus(fKillTrackAndSecondaries);
				//			}

				break;

			}

			}
			//}
			//}	
			for(G4int j=0;j<nprocesses;j++)
			{

				if( (*pv)[j]->GetProcessType()==fOptical && (*pv)[j]->GetProcessSubType() == 32 )
				{
					boundary = (G4OpBoundaryProcess*)(*pv)[i];
					std::cout << "boundary" << boundary << std::endl;
					if(boundary-> GetStatus() ==4)
					{			
						if(thePrePoint->GetMaterial()->GetName()=="SilGel" && thePostPoint->GetMaterial()->GetName()=="Air")

						{
							aStep ->GetTrack() ->  SetTrackStatus(fKillTrackAndSecondaries);
						}
					}
				}
			}			
	}
	//changed from here
	/*	if( (*pv)[i]->GetProcessType()==fOptical && (*pv)[i]->GetProcessSubType() == 32 )
		{
	//		if((preposition_r >= 49.5983 && preposition_r <=50.4975 )||( postposition_r >= 49.5983 && postposition_r <= 50.4975))
	//		{
	G4cout<<" ProcessName:" << (*pv)[i]->GetProcessName() <<G4endl;
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	G4cout<<"  G4OpBoundaryProcessStatus:" << boundary->GetStatus() <<G4endl;
	std::cout << "PreStep Radius  = " << preposition_r << std::endl;
	std::cout << "PreStep Material " << thePrePoint->GetMaterial()->GetName() << std::endl;
	std::cout << "PreStep Momentum = " << PreMomentum << std::endl;
	std::cout << "PreStep MomentumDirection = " << PreMomentumDirection << std::endl;


	std::cout << "PostStep Radius  = " << postposition_r << std::endl;
	std::cout << "PostStep Material " << thePostPoint->GetMaterial()->GetName() << std::endl;
	std::cout << "PostStep Momentum = " << PostMomentum << std::endl;
	std::cout << "PostStep MomentumDirection = " << PostMomentumDirection << std::endl;
	G4double angle = acos((x*X + y*Y + z*Z)/((sqrt(x*x + y*y + z*z)) * (sqrt(X*X + Y*Y + Z*Z))));
	std::cout << " angle_calculated = " << angle << std::endl;
	G4double Angle = PostMomentumDirection.angle(PreMomentumDirection);
	std::cout << " angle_function = " << Angle << std::endl;
	//till here*/




	G4ParticleDefinition *particleType = track->GetDefinition();
	if(particleType == G4OpticalPhoton::OpticalPhotonDefinition())
	{
		if( (thePrePV->GetName().find("MultiPMT") != std::string::npos) &&
				(thePostPV->GetName().find("vessel") != std::string::npos) )
			n_photons_through_mPMTLV++;

		if( (thePostPV->GetName().find("container") != std::string::npos) &&
				(thePrePV->GetName().find("vessel") != std::string::npos))
			n_photons_through_acrylic++;

		if( (thePrePV->GetName().find("container") != std::string::npos) )
			n_photons_through_gel++;

		if( (thePrePV->GetName().find("container") != std::string::npos) &&
				(thePostPV->GetName().find("inner") != std::string::npos) )
			n_photons_on_blacksheet++;

		if( (thePrePV->GetName().find("container") != std::string::npos) &&
				(thePostPV->GetName().find("pmt") != std::string::npos) )
			n_photons_on_smallPMT++;


		/*
		   if( (thePrePV->GetName().find("pmt") != std::string::npos)){
		   std::cout << "Photon between " << thePrePV->GetName() <<
		   " and " << thePostPV->GetName() << " because " << 
		   thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
		   " and boundary status: " <<  boundary->GetStatus() << " with track status " << track->GetTrackStatus() << std::endl;

		   }*/



		/*			if(track->GetTrackStatus() == fStopAndKill)
					{
					if(boundary->GetStatus() == NoRINDEX)
					{
					std::cout << "Optical photon is killed because of missing refractive index in either " << thePrePoint->GetMaterial()->GetName() << " or " << thePostPoint->GetMaterial()->GetName() << 
					" : could also be caused by Overlaps with volumes with logicalBoundaries." << std::endl;

					}*/
		/* Debug :  
		   if( (thePrePV->GetName().find("PMT") != std::string::npos) ||
		   (thePrePV->GetName().find("pmt") != std::string::npos)){

		   if(boundary->GetStatus() != StepTooSmall){
		//	if(thePostPoint->GetProcessDefinedStep()->GetProcessName() != "Transportation")
		std::cout << "Killed photon between " << thePrePV->GetName() <<
		" and " << thePostPV->GetName() << " because " << 
		thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
		" and boundary status: " <<  boundary->GetStatus() <<
		std::endl;
		}
		}	
		*/

		//}



		//debugging 
		//  G4Track* theTrack = aStep->GetTrack();
		//  const G4DynamicParticle* aParticle = theTrack->GetDynamicParticle();
		//  G4ThreeVector aMomentum = aParticle->GetMomentumDirection();
		//  G4double vx = aMomentum.x();
		//  G4int ix = std::isnan(vx);
		//  if(ix != 0){
		//    G4cout << " PROBLEM! " << theTrack->GetCreatorProcess()->GetProcessName() <<
		//  std::flush << G4endl;
		//  }



	}
	}
	G4int WCSimSteppingAction::G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
			G4ThreeVector lArPos,
			G4ThreeVector start,
			G4int i)
	{
		G4double x0 = start[0]-lArPos[0];
		G4double y0 = start[1]-lArPos[1];
		//   G4double y0 = 2121.3;//mm
		G4double z0 = start[2]-lArPos[2];

		G4double dt=0.8;//mm
		//   G4double midt = 2651.625;
		G4double pitch = 3;//mm
		//   G4double midwir = 1207.10;
		G4double c45 = 0.707106781;
		G4double s45 = 0.707106781;

		G4double w1;
		G4double w2;
		G4double t;

		//   G4double xField(0.);
		//   G4double yField(0.);
		//   G4double zField(0.);

		//   if(detector->getElectricFieldDistortion())
		//     {
		//       Distortion(pos3d->getX(),
		// 		 pos3d->getY());

		//       xField = ret[0];
		//       yField = ret[1];
		//       zField = pos3d->getZ();

		//       w1 = (int)(((zField+z0)*c45 + (x0-xField)*s45)/pitch); 
		//       w2 = (int)(((zField+z0)*c45 + (x0+xField)*s45)/pitch); 
		//       t = (int)(yField+1);

		//       //G4cout<<" x orig "<<pos3d->getX()<<" y orig "<<(pos3d->getY()+y0)/dt<<G4endl;
		//       //G4cout<<" x new "<<xField<<" y new "<<yField<<G4endl;
		//     }
		//   else 
		//     {

		w1 = (int) (((pos3d->getZ()+z0)*c45 + (x0-pos3d->getX())*s45)/pitch); 
		w2 = (int)(((pos3d->getZ()+z0)*c45 + (x0+pos3d->getX())*s45)/pitch); 
		t  = (int)((pos3d->getY()+y0)/dt +1);
		//     }

		if (i==0)
			return (int)w1;
		else if (i==1)
			return (int)w2;
		else if (i==2)
			return (int)t;
		else return 0;
	} 


	void WCSimSteppingAction::Distortion(G4double /*x*/,G4double /*y*/)
	{

		//   G4double theta,steps,yy,y0,EvGx,EvGy,EField,velocity,tSample,dt;
		//   y0=2121.3;//mm
		//   steps=0;//1 mm steps
		//   tSample=0.4; //micros
		//   dt=0.8;//mm
		//   LiquidArgonMedium medium;  
		//   yy=y;
		//   while(y<y0 && y>-y0 )
		//     {
		//       EvGx=FieldLines(x,y,1);
		//       EvGy=FieldLines(x,y,2);
		//       theta=atan(EvGx/EvGy);
		//       if(EvGy>0)
		// 	{
		// 	  x+=sin(theta);
		// 	  y+=cos(theta);
		// 	}
		//       else
		// 	{
		// 	  y-=cos(theta);
		// 	  x-=sin(theta);
		// 	}
		//       EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
		//       velocity=medium.DriftVelocity(EField*10);// mm/microsec
		//       steps+=1/(tSample*velocity);

		//       //G4cout<<" step "<<steps<<" x "<<x<<" y "<<y<<" theta "<<theta<<" Gx "<<eventaction->Gx->Eval(x,y)<<" Gy "<<eventaction->Gy->Eval(x,y)<<" EField "<<EField<<" velocity "<<velocity<<G4endl;
		//     }

		//   //numbers
		//   //EvGx=FieldLines(0,1000,1);
		//   //EvGy=FieldLines(0,1000,2);
		//   //EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
		//   //velocity=medium.DriftVelocity(EField*10);// mm/microsec
		//   //G4double quenching;
		//   //quenching=medium.QuenchingFactor(2.1,EField*10);
		//   //G4cout<<" Gx "<<EvGx<<" Gy "<<EvGy<<" EField "<<EField<<" velocity "<<velocity<<" quenching "<<quenching<<G4endl;


		//ret[0]=5;
		//   if(yy>0)
		//     ret[1]=2*y0/dt -steps;
		//   else
		//     ret[1]=steps; 
	}


	double WCSimSteppingAction::FieldLines(G4double /*x*/,G4double /*y*/,G4int /*coord*/)
	{ //0.1 kV/mm = field
		//G4double Radius=302;//mm
		//   G4double Radius=602;//mm
		//   if(coord==1) //x coordinate
		//     return  (0.1*(2*Radius*Radius*x*abs(y)/((x*x+y*y)*(x*x+y*y))));
		//   else //y coordinate
		//     return 0.1*((abs(y)/y)*(1-Radius*Radius/((x*x+y*y)*(x*x+y*y))) + abs(y)*(2*Radius*Radius*y/((x*x+y*y)*(x*x+y*y))));
		return 0;
	}

