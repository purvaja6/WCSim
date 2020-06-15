#include "WCSimDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SDManager.hh"
#include "WCSimWCSD.hh"
#include "WCSimPMTObject.hh"

#include "G4SystemOfUnits.hh"

#include "G4LogicalSkinSurface.hh"

//PMT logical volume construction.
//WCSimDetectorConstruction::PMTMap_t WCSimDetectorConstruction::PMTLogicalVolumes;

//ToDo: use enum instead of string for PMTtype, and PMTs should be objects of one simple class.
G4LogicalVolume* WCSimDetectorConstruction::ConstructPMT(G4String PMTName, G4String CollectionName)
{
	PMTKey_t key(PMTName,CollectionName);

	// Return pre-created PMT Logical Volume if it already exists.
	PMTMap_t::iterator it = PMTLogicalVolumes.find(key);
	if (it != PMTLogicalVolumes.end()) {
		//G4cout << "Restore PMT" << G4endl;
		return it->second;
	}

	//G4cout << "Create PMT" << G4endl;


	if (Vis_Choice == "RayTracer"){
		// Blue wireframe visual style
		// Used in the RayTracer visualizer
		G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
		WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
		WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
	}


	G4double expose = 0.;
	G4double radius = 0.;
	G4double glassThickness = 0.;

	WCSimPMTObject *PMT = GetPMTPointer(CollectionName);
	//expose = PMT->GetExposeHeight();
	expose = 15.31*CLHEP::mm;
	std::cout << "expose " << expose << std::endl;
	//radius = PMT->GetRadius();                            //r at height = expose
	radius = 36*CLHEP::mm;
	std::cout << "radius " << radius << std::endl;
	glassThickness = PMT->GetPMTGlassThickness();
	std::cout <<"glass thickness " << glassThickness << std::endl;
	//sphereRadius R: radius of curvature, based on spherical approx near exposeHeight
	//if radius of spherical cap at exposeHeight = PMTradius r (eg. 20"/2)
	//then from r = sqrt(R*R - (R-expose)*(R-expose)), we get:
	G4double sphereRadius = (expose*expose + radius*radius)/(2*expose);
	std::cout << "SphereRadius " << sphereRadius << std::endl;
	G4double PMTOffset =  sphereRadius - expose;
	std::cout << "PMToffset " << PMTOffset << std::endl; 

	//Optional reflectorCone:
	G4double reflectorRadius = radius + id_reflector_height * tan(id_reflector_angle); // PMT radius+ r = h * tan (theta)
	std::cout << "reflector height" << id_reflector_height << std::endl;
	//G4double reflectorThickness = 0.5*CLHEP::mm;
	G4double reflectorThickness = 0.328*CLHEP::mm; //me: the actual reflector thickness is 0.5 mm but due to solid works design, I am taking the horizontal component to calculate the reflector radius. Refer my CAD drawing.
	std::cout << "reflector thickness" << reflectorThickness << std::endl;
	if((reflectorRadius - radius) < 1.*CLHEP::mm)
		reflectorThickness = 0.*CLHEP::mm;


	G4cout << "========================================================" << G4endl;
	G4cout << "Expose height for PMT " << CollectionName << " : " << expose << " mm" << G4endl;
	G4cout << "Radius at expose height : " << radius << " mm" << G4endl;
	G4cout << "Radius of curvature : " << sphereRadius << " mm" << G4endl;
	G4cout << "Reflector radius : " << reflectorRadius << " mm" << G4endl;
	G4cout << "Reflector thickness : " << reflectorThickness << " mm" << G4endl;
	G4cout << "========================================================" << G4endl;


	G4VSolid *solidWCPMT;
	G4double position_z_offset = 0.;  // for positioning PMTs
	G4double wcpmt_z_offset = 0.;     // for positioning single PMT support (mPMT)
	G4bool addPMTBase = false; 
	//G4double pmtModuleHeight = 54.*CLHEP::mm; //includes puck and single PMT support, not PMT base. //default wcsim design
	G4double pmtModuleHeight = 59.62 *CLHEP::mm; //the height of pmt module for solid works design

	if(nID_PMTs == 1){

		//All components of the PMT are now contained in a single logical volume logicWCPMT.
		//Origin is on the blacksheet, faces positive z-direction.

		G4double PMTHolderZ[2] = {0, std::max(expose,id_reflector_height+id_reflector_z_offset)};
		G4double PMTHolderR[2] = {std::max(radius,reflectorRadius) + reflectorThickness, 
			std::max(radius,reflectorRadius) + reflectorThickness};
		G4double PMTHolderr[2] = {0,0};
		//  std::cout << "pmtholderr" << 
		// std::cout << "PmtholderR" << 

		solidWCPMT = 
			new G4Polycone("WCPMT",                    
					0.0*deg,
					360.0*deg,
					2,
					PMTHolderZ,
					PMTHolderr, // R Inner
					PMTHolderR);// R Outer

	} else{

		// NEW: reduce overlaps by using curved surfaces, but only for mPMTs because
		// the single PMTs will rest on a curved PMT support (in WCSim as the polygon is
		// quite non trivial to model and not so important to get right initially, later
		// we'll directly use the SolidWorks exported STL to CADMesh as Tesselated surface
		// but that's not urgent now)

		//////////////////////////////
		/// Optional: Add PMT base ///
		//////////////////////////////

		// ToDo: Base is PMT property! Should not hard coded here.
		G4double baseHeight = 0.;
		G4double baseRadius = 0.;
		if(addPMTBase){
			G4double basePinLength = 28.*CLHEP::mm; //TF measurement: pins + base + transformer.
			baseHeight = 97.*CLHEP::mm - expose + basePinLength; //97mm includes the PMT top as well.
			baseRadius = 26.*CLHEP::mm; //for R121990-02
		} else {
			// version without a base but with optional reflectorCone
			baseHeight = expose;
			baseRadius = std::max(radius,reflectorRadius) + reflectorThickness; //radius + reflectorThickness; too tight
		}

		// ToDo: extend the shell 
		/////////////////////////////


		// for single PMTs no overlap issues and will rest on flat black sheets
		// Need spherical shell for curved top to better match (and reduce overlaps) with
		// gel layer ("container") and vessel.
		solidWCPMT =
			new G4Sphere("WCPMT",
					vessel_radius_curv - mPMT_outer_material_d - pmtModuleHeight, //me:rMin = 342 - 10 - 54 = 278 mm, 54mm is position of support structure
					vessel_radius_curv - mPMT_outer_material_d,          //me: rMax = 332 mm
					0.0*deg,                                             //phiStart
					360.0*deg,                                           //Deltaphi
					0.0*deg,                                             //thetaStart
					mPMT_pmt_openingAngle);	//me: 8.7 deg                                            //Deltatheta
		// Once we have a full PMT support in (not urgent),
		// replacing the current "inner" black sphere support by using
		// G4Tesselated from CADMesh'ing an STL,
		// Cutoff the bottom to make it flat.



		//z=0 for spherical shell needs to be shifted
		wcpmt_z_offset = vessel_radius_curv - mPMT_outer_material_d - pmtModuleHeight;
		position_z_offset = vessel_radius_curv - mPMT_outer_material_d - expose - dist_pmt_vessel;
		std::cout << "wcpmt_z_offset " << wcpmt_z_offset << std::endl;
		std::cout << "position_z_offset " << position_z_offset << std::endl;
	}

	/*
	   G4Material *material_around_pmt = G4Material::GetMaterial("Water");
	   if(id_reflector_height > 0.1*CLHEP::mm && 
	   (vessel_radius > 1.*CLHEP::mm || vessel_cyl_height > 1.*CLHEP::mm)) //or make this a user option? 
	   material_around_pmt = ; //"SilGel"
	   */

	//ToDo: define defaults for single PMTs and rename variable (optional gel for single 20" PMT)
	G4LogicalVolume* logicWCPMT =
		new G4LogicalVolume(    solidWCPMT,
				G4Material::GetMaterial(mPMT_material_pmtAssembly),        // Default: SilGel for n > 1, Water for n = 1 
				"WCPMT",
				0,0,0);



	//////////////////////////
	// Visualize logicWCPMT //
	//////////////////////////

	// if (Vis_Choice == "RayTracer"){
	// Makes the volume containing the PMT visible, solid, and forces the auxiliary edges to be viewed.
	G4VisAttributes* WCPMTVisAtt5 = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	WCPMTVisAtt5->SetForceSolid(true); // force the object to be visualized with a surface
	WCPMTVisAtt5->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

	logicWCPMT->SetVisAttributes(WCPMTVisAtt5);
	/* } else{
	   G4VisAttributes* WCPMTVisAtt2 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
	   WCPMTVisAtt2->SetForceSolid(true);
	   if(addPMTBase) 
	   logicWCPMT->SetVisAttributes(WCPMTVisAtt2);
	   else{
	// Makes the volume containg the PMT invisible for normal visualization
	G4VisAttributes* WCPMTVisAttGrey = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
	WCPMTVisAttGrey->SetForceSolid(true);
	// logicWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
	logicWCPMT->SetVisAttributes(WCPMTVisAttGrey);                   //TF debug overlaps with this volume
	}
	}*/

	//Need a volume to cut away excess behind blacksheet
	G4Box* solidCutOffTubs =
		new G4Box(    "cutOffTubs",
				sphereRadius+1.*cm,
				sphereRadius+1.*cm,
				PMTOffset);

	////////////////////////////
	///  Create PMT Interior ///
	////////////////////////////

	G4Sphere* tmpSolidInteriorWCPMT =
		new G4Sphere(    "tmpInteriorWCPMT",
				0.0*m,(sphereRadius-glassThickness), //me: 50 - 1 = 49 mm
				0.0*deg,360.0*deg,
				0.0*deg,90.0*deg);

	G4SubtractionSolid* solidInteriorWCPMT =
		new G4SubtractionSolid(    "InteriorWCPMT",
				tmpSolidInteriorWCPMT,
				solidCutOffTubs);

	// "Air" here is not true air, but a modified material
	// with n = 1 and a very short absorption length
	G4LogicalVolume* logicInteriorWCPMT =
		new G4LogicalVolume(    solidInteriorWCPMT,
				G4Material::GetMaterial("Air"), //me: vacuum
				"InteriorWCPMT",
				0,0,0);

	G4VPhysicalVolume* physiInteriorWCPMT =
		new G4PVPlacement(0,
				G4ThreeVector(0, 0, -1.0*PMTOffset + position_z_offset),
				logicInteriorWCPMT,
				"InteriorWCPMT",
				logicWCPMT,
				false,
				0,
				checkOverlaps);


	//////////////////////////////////
	// Visualize logicInteriorWCPMT //
	/////////////////////////////////

	//  if (Vis_Choice == "RayTracer"){
	// Adding color and forcing the inner portion of the PMT's to be viewed
	G4VisAttributes* WCPMTVisAtt6 = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	WCPMTVisAtt6->SetForceSolid(true); // force the object to be visualized with a surface
	WCPMTVisAtt6->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 

	logicInteriorWCPMT->SetVisAttributes(WCPMTVisAtt6);
	/* } else {
	// Making the inner portion of the detector invisible for OGLSX visualization
	logicInteriorWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
	}*/

	/////////////////////////////
	/// Create PMT Glass Face ///
	/////////////////////////////

	G4Sphere* tmpGlassFaceWCPMT =
		new G4Sphere(    "tmpGlassFaceWCPMT",
				(sphereRadius-glassThickness),
				sphereRadius,
				0.0*deg,360.0*deg,
				0.0*deg,90.0*deg);

	G4SubtractionSolid* solidGlassFaceWCPMT =
		new G4SubtractionSolid(    CollectionName,
				tmpGlassFaceWCPMT,
				solidCutOffTubs); 

	G4LogicalVolume *logicGlassFaceWCPMT =
		new G4LogicalVolume(    solidGlassFaceWCPMT,
				G4Material::GetMaterial("Glass"),
				CollectionName,
				0,0,0);

	G4VPhysicalVolume* physiGlassFaceWCPMT =
		new G4PVPlacement(0,
				G4ThreeVector(0, 0, -1.0*PMTOffset+position_z_offset),
				logicGlassFaceWCPMT,
				CollectionName,
				logicWCPMT,
				false,
				0,
				checkOverlaps);

	// For either visualization type, logicGlassFaceWCPMT will either be visible 
	// or invisible depending on which line is commented at the end of the respective if statements.

	if (Vis_Choice == "OGLSX")
	{ // Gray wireframe visual style
		// used in OGLSX visualizer
		//G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
		G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
		WCPMTVisAtt->SetForceWireframe(true);
		//logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
		logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);}

	if (Vis_Choice == "RayTracer"){
		// Blue wireframe visual style
		// Used in the RayTracer visualizer
		G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
		WCPMTVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
		WCPMTVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown 
		//logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);

		logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);
	} else { 

		// Gray wireframe visual style
		// used in OGLSX visualizer
		//G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));  //better for seeing hits
		//WCPMTVisAtt->SetForceWireframe(true);
		G4VisAttributes* WCPMTVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //better for seeing geometry
		WCPMTVisAtt->SetForceSolid(true);
		//logicGlassFaceWCPMT->SetVisAttributes(G4VisAttributes::Invisible);
		logicGlassFaceWCPMT->SetVisAttributes(WCPMTVisAtt);

	}

	//Add Logical Border Surface
	new G4LogicalBorderSurface("GlassCathodeSurface",
			physiGlassFaceWCPMT,
			physiInteriorWCPMT,
			OpGlassCathodeSurface);



	//////////////////////////
	/// Optional Reflector ///
	//////////////////////////

	//me: reflector configuration for solid works design

	/*Some details:
	 *                  *               1.1mm is the gap between reflector and PMT for mechanical construction 
	 *                                   *               (KM3NeT support matrix value)
	 *                                                    * three degrees of freedom: height, z position and opening angle
	 *                                                                     */

	if(id_reflector_height > 0.1*CLHEP::mm
			&& (reflectorRadius-radius) > -5*CLHEP::mm){

		G4Cons* reflectorCone =
			new G4Cons("WCPMT_reflect",
					radius - 0.72*CLHEP::mm,                               //rmin
					radius - 0.72*CLHEP::mm + reflectorThickness,          //rmax
					reflectorRadius - 0.72*CLHEP::mm,                      //Rmin
					reflectorRadius - 0.72*CLHEP::mm + reflectorThickness, //Rmax
					id_reflector_height/2,                                //z/2
					0, 2*CLHEP::pi);
		std::cout << "reflectorradiusmin " << reflectorRadius - 0.72*CLHEP::mm << std::endl;
		std::cout << "reflectorradiusmax " << reflectorRadius - 0.72*CLHEP::mm + reflectorThickness << std::endl;

		std::cout << "radiusmin" << radius - 0.72*CLHEP::mm << std::endl;
		std::cout << "radiusmax" << radius - 0.72*CLHEP::mm + reflectorThickness << std::endl;
		G4LogicalVolume* logicReflector =
			new G4LogicalVolume(    reflectorCone,
					G4Material::GetMaterial("Aluminum"), //It actually is Al+ Ag evaporation
					"reflectorCone",
					0,0,0);
		G4VisAttributes* WCPMTVisAtt3 = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
		WCPMTVisAtt3->SetForceSolid(true);
		logicReflector->SetVisAttributes(WCPMTVisAtt3);

		new G4LogicalSkinSurface("ReflectorLogSkinSurface",logicReflector,ReflectorSkinSurface);
		new G4PVPlacement(0,
				G4ThreeVector(0, 0, id_reflector_z_offset+id_reflector_height/2+position_z_offset),
				logicReflector, // me: current logical
				"reflectorWCPMT",
				logicWCPMT, // me: mother logical
				false,
				0,
				checkOverlaps);
	}
	//me: till here I have made changes in the reflector geometry based on solid works design

	// me:I have commented out from here

	if(nID_PMTs > 1){
		////////////////////
		/// 1-PMT support //
		////////////////////

		//me:changed the lines below to adapt solid works design

		//support

		if(id_reflector_height > 0.1*CLHEP::mm
				&& (reflectorRadius-radius) > -5*CLHEP::mm){

			G4double ReflectorHolderZ[4] = {0, 23.18, 41.64, 46.36};
			G4double ReflectorHolderR[4] = {39.735, 43.19, 45.92, 46.64};

			G4double ReflectorHolderr[4] = {26.75, 40.52, 40.52, 45.96};

			G4Polycone * solidWCPMTsupport2 =
				new G4Polycone("WCPMTsupport2",
						0.0*deg,
						360.0*deg,
						4,
						ReflectorHolderZ,
						ReflectorHolderr, // R Inner
						ReflectorHolderR);// R Outer

			G4LogicalVolume* logicWCPMTsupport2 =
				new G4LogicalVolume(solidWCPMTsupport2,
						G4Material::GetMaterial("Blacksheet"),
						"WCPMTsupport",
						0,0,0);

			// Visualize
			G4VisAttributes* WCPMTVisAtt_sup = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
			WCPMTVisAtt_sup->SetForceSolid(true);
			//                         logicWCPMTsupport->SetVisAttributes(WCPMTVisAtt_sup);
			new G4LogicalSkinSurface("FoamLogSkinSurface2",logicWCPMTsupport2,OpGelFoamSurface);
			logicWCPMTsupport2->SetVisAttributes(WCPMTVisAtt_sup);

			new G4PVPlacement(0,
					G4ThreeVector(0, 0, 272.36), // value of the z=0 plane with respect to the origin
					logicWCPMTsupport2,
					"WCPMTsupport2",
					logicWCPMT,
					false,
					0,
					checkOverlaps);


		}
	}
	// me: till here code has been modified to adapt solid works design



	/////////////////////////
	// Sensitive detector ///
	/////////////////////////

	// Instantiate a new sensitive detector 
	// and register this sensitive detector volume with the SD Manager. 
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDName = "/WCSim/";
	SDName += CollectionName;

	// If there is no such sensitive detector with that SDName yet,
	// make a new one
	if( ! SDman->FindSensitiveDetector(SDName, false) ) {

		aWCPMT = new WCSimWCSD(CollectionName,SDName,this );
		SDman->AddNewDetector( aWCPMT );
	}

	logicGlassFaceWCPMT->SetSensitiveDetector( aWCPMT );


	PMTLogicalVolumes[key] = logicWCPMT;

	return logicWCPMT;
}

