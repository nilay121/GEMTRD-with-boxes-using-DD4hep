//==========================================================================
//  dRICh: Dual Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: Christopher Dilks (Duke University)
//
// - Design Adapted from Standalone Fun4all and GEMC implementations
//   [ Evaristo Cisbani, Cristiano Fanelli, Alessio Del Dotto, et al. ]
//
//==========================================================================

#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "GeometryHelpers.h"
#include "Math/Point2D.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <iostream>

using namespace dd4hep;
using namespace dd4hep::rec;

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();

  DetElement det(detName, detID);
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // attributes -----------------------------------------------------------
  // - vessel
  double  vesselZmin       =  dims.attr<double>(_Unicode(zmin));
  double  vesselLength     =  dims.attr<double>(_Unicode(length));
  double  vesselRmin0      =  dims.attr<double>(_Unicode(rmin0));
  double  vesselRmin1      =  dims.attr<double>(_Unicode(rmin1));
  double  vesselRmax0      =  dims.attr<double>(_Unicode(rmax0));
  double  vesselRmax1      =  dims.attr<double>(_Unicode(rmax1));
  double  vesselRmax2      =  dims.attr<double>(_Unicode(rmax2));
  double  snoutLength      =  dims.attr<double>(_Unicode(snout_length));
  int     nSectors         =  dims.attr<int>(_Unicode(nsectors));
  double  wallThickness    =  dims.attr<double>(_Unicode(wall_thickness));
  double  windowThickness  =  dims.attr<double>(_Unicode(window_thickness));
  auto    vesselMat        =  desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto    gasvolMat        =  desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto    vesselVis        =  desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto    gasvolVis        =  desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));
  // - radiator (applies to aerogel and filter)
  auto    radiatorElem        =  detElem.child(_Unicode(radiator));
  double  radiatorRmin        =  radiatorElem.attr<double>(_Unicode(rmin));
  double  radiatorRmax        =  radiatorElem.attr<double>(_Unicode(rmax));
  double  radiatorPhiw        =  radiatorElem.attr<double>(_Unicode(phiw));
  double  radiatorPitch       =  radiatorElem.attr<double>(_Unicode(pitch));
  double  radiatorFrontplane  =  radiatorElem.attr<double>(_Unicode(frontplane));
  // - aerogel
  auto    aerogelElem       =  radiatorElem.child(_Unicode(aerogel));
  auto    aerogelMat        =  desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  auto    aerogelVis        =  desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  double  aerogelThickness  =  aerogelElem.attr<double>(_Unicode(thickness));
  // - filter
  auto    filterElem       =  radiatorElem.child(_Unicode(filter));
  auto    filterMat        =  desc.material(filterElem.attr<std::string>(_Unicode(material)));
  auto    filterVis        =  desc.visAttributes(filterElem.attr<std::string>(_Unicode(vis)));
  double  filterThickness  =  filterElem.attr<double>(_Unicode(thickness));
  // - mirror
  auto    mirrorElem       =  detElem.child(_Unicode(mirror));
  auto    mirrorMat        =  desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto    mirrorVis        =  desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto    mirrorSurf       =  surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  double  mirrorBackplane  =  mirrorElem.attr<double>(_Unicode(backplane));
  double  mirrorThickness  =  mirrorElem.attr<double>(_Unicode(thickness));
  double  mirrorRmin       =  mirrorElem.attr<double>(_Unicode(rmin));
  double  mirrorRmax       =  mirrorElem.attr<double>(_Unicode(rmax));
  double  mirrorPhiw       =  mirrorElem.attr<double>(_Unicode(phiw));
  double  focusTuneZ       =  mirrorElem.attr<double>(_Unicode(focus_tune_z));
  double  focusTuneX       =  mirrorElem.attr<double>(_Unicode(focus_tune_x));
  // - sensor module
  auto    sensorElem       =  detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto    sensorMat        =  desc.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto    sensorVis        =  desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto    sensorSurf       =  surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  double  sensorSide       =  sensorElem.attr<double>(_Unicode(side));
  double  sensorGap        =  sensorElem.attr<double>(_Unicode(gap));
  double  sensorThickness  =  sensorElem.attr<double>(_Unicode(thickness));
  // - sensor sphere
  auto    sensorSphElem     =  detElem.child(_Unicode(sensors)).child(_Unicode(sphere));
  double  sensorSphRadius   =  sensorSphElem.attr<double>(_Unicode(radius));
  double  sensorSphCenterX  =  sensorSphElem.attr<double>(_Unicode(centerx));
  double  sensorSphCenterZ  =  sensorSphElem.attr<double>(_Unicode(centerz));
  // - sensor sphere patch cuts
  auto    sensorSphPatchElem  =  detElem.child(_Unicode(sensors)).child(_Unicode(sphericalpatch));
  double  sensorSphPatchPhiw  =  sensorSphPatchElem.attr<double>(_Unicode(phiw));
  double  sensorSphPatchRmin  =  sensorSphPatchElem.attr<double>(_Unicode(rmin));
  double  sensorSphPatchRmax  =  sensorSphPatchElem.attr<double>(_Unicode(rmax));
  double  sensorSphPatchZmin  =  sensorSphPatchElem.attr<double>(_Unicode(zmin));
  // - debugging switches
  int   debug_optics_mode  =  detElem.attr<int>(_Unicode(debug_optics));
  bool  debug_mirror       =  mirrorElem.attr<bool>(_Unicode(debug));
  bool  debug_sensors      =  sensorSphElem.attr<bool>(_Unicode(debug));

  // if debugging optics, override some settings
  bool debug_optics = debug_optics_mode > 0;
  if(debug_optics) {
    printout(WARNING,"DRich_geo","DEBUGGING DRICH OPTICS");
    switch(debug_optics_mode) {
      case 1: vesselMat = aerogelMat = filterMat = sensorMat = gasvolMat = desc.material("VacuumOptical"); break;
      case 2: vesselMat = aerogelMat = filterMat = sensorMat = desc.material("VacuumOptical"); break;
      default: printout(FATAL,"DRich_geo","UNKNOWN debug_optics_mode"); return det;
    };
    aerogelVis = sensorVis = mirrorVis;
    gasvolVis = vesselVis = desc.invisible();
  };


  // BUILD VESSEL ====================================================================
  /* - `vessel`: aluminum enclosure, the mother volume of the dRICh
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   * - the dRICh vessel geometry has two regions: the snout refers to the conic region
   *   in the front, housing the aerogel, while the tank refers to the cylindrical
   *   region, housing the rest of the detector components
   */

  // derived attributes
  double tankLength = vesselLength - snoutLength;
  double vesselZmax = vesselZmin + vesselLength;

  // snout solids
  double boreDelta = vesselRmin1 - vesselRmin0;
  double snoutDelta = vesselRmax1 - vesselRmax0;
  Cone vesselSnout(
      snoutLength/2.0,
      vesselRmin0,
      vesselRmax0,
      vesselRmin0 + boreDelta * snoutLength / vesselLength,
      vesselRmax1
      );
  Cone gasvolSnout(
      /* note: `gasvolSnout` extends a bit into the tank, so it touches `gasvolTank`
       * - the extension distance is equal to the tank `windowThickness`, so the
       *   length of `gasvolSnout` == length of `vesselSnout`
       * - the extension backplane radius is calculated using similar triangles
       */
      snoutLength/2.0,
      vesselRmin0 + wallThickness,
      vesselRmax0 - wallThickness,
      vesselRmin0 + boreDelta * (snoutLength-windowThickness) / vesselLength + wallThickness,
      vesselRmax1 - wallThickness + windowThickness * (vesselRmax1 - vesselRmax0) / snoutLength
      );

  // tank solids
  Cone vesselTank(
      tankLength/2.0,
      vesselSnout.rMin2(),
      vesselRmax2,
      vesselRmin1,
      vesselRmax2
      );
  Cone gasvolTank(
      tankLength/2.0 - windowThickness,
      gasvolSnout.rMin2(),
      vesselRmax2 - wallThickness,
      vesselRmin1 + wallThickness,
      vesselRmax2 - wallThickness
      );

  // snout + tank solids
  UnionSolid vesselUnion(
      vesselTank,
      vesselSnout,
      Position(0., 0., -vesselLength/2.)
      );
  UnionSolid gasvolUnion(
      gasvolTank,
      gasvolSnout,
      Position(0., 0., -vesselLength/2. + windowThickness)
      );

  //  extra solids for `debug_optics` only
  Box vesselBox(1001,1001,1001);
  Box gasvolBox(1000,1000,1000);

  // choose vessel and gasvol solids (depending on `debug_optics_mode` (0=disabled))
  Solid vesselSolid, gasvolSolid;
  switch(debug_optics_mode) {
    case 0:  vesselSolid=vesselUnion;  gasvolSolid=gasvolUnion;  break; // `!debug_optics`
    case 1:  vesselSolid=vesselBox;    gasvolSolid=gasvolBox;    break;
    case 2:  vesselSolid=vesselBox;    gasvolSolid=gasvolUnion;  break;
  };

  // volumes
  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName+"_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  // reference positions
  // - the vessel is created such that the center of the cylindrical tank volume
  //   coincides with the origin; this is called the "origin position" of the vessel
  // - when the vessel (and its children volumes) is placed, it is translated in
  //   the z-direction to be in the proper ATHENA-integration location
  // - these reference positions are for the frontplane and backplane of the vessel,
  //   with respect to the vessel origin position
  auto originFront = Position(0., 0., -tankLength/2.0 - snoutLength );
  auto originBack =  Position(0., 0., tankLength/2.0 );

  // initialize sensor centroids (used for mirror parameterization below); this is
  // the average (x,y,z) of the placed sensors, w.r.t. originFront
  double sensorCentroidX = 0;
  double sensorCentroidZ = 0;
  int sensorCount = 0;


  // sensitive detector type
  sens.setType("tracker");
  
  
  //Radiator material+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // BUILD RADIATOR ====================================================================

  // attributes
  double airGap = 0.01*mm; // air gap between aerogel and filter (FIXME? actually it's currently a gas gap)

  // solid and volume: create aerogel and filter
   
  Box radiator(100*cm,100*cm,100*cm);
   
  //radiator material
  Material radiatorMat = desc.material("radiatorMatCustom");  
  Volume aerogelVol( detName+"_aerogel", radiator, radiatorMat );
  //Volume filterVol(  detName+"_filter",  filterSolid,  filterMat );
  aerogelVol.setVisAttributes(aerogelVis);
  //filterVol.setVisAttributes(filterVis);

  // aerogel placement and surface properties
  // TODO [low-priority]: define skin properties for aerogel and filter
  
  //position of the radiator --> coordinates
  auto aerogelPV = gasvolVol.placeVolume(aerogelVol,Position(0*cm, 0*cm, 0*cm));
  DetElement aerogelDE(det, "aerogel_de", 0);
  aerogelDE.setPlacement(aerogelPV);
  // place gas volume
  //PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol,Position(0, 0, 0));
  //DetElement gasvolDE(det, "gasvol_de", 0);
  //gasvolDE.setPlacement(gasvolPV);

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(gasvolVol,
      Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
};

// clang-format off
DECLARE_DETELEMENT(athena_DRICH, createDetector)

