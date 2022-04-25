/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
 *
 * \ingroup trackers
 *
 * @{
 */
#include <map>
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"
#include "XML/Layering.h"

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Endcap Trapezoidal Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     vacuum   = description.vacuum();
  Material                     radiatorMat = description.material("radiatorMatCustom");
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  bool                         reflect  = x_det.reflect(false);
  DetElement                   sdet(det_name, det_id);
  Assembly                     assembly(det_name);

  Material  air  = description.material("Air");
  // Volume      assembly    (det_name,Box(10000,10000,10000),vacuum);
  Volume                  motherVol = description.pickMotherVolume(sdet);
  int                     m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume>     modules;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>>        volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume            pv;


  // ACTS extension
 /* {
    Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();

    for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
      xml_comp_t x_boundary_material = bmat;
      Acts::xmlToProtoSurfaceMaterial(x_boundary_material, *detWorldExt, "boundary_material");
    }
    sdet.addExtension<Acts::ActsExtension>(detWorldExt);
  }

  assembly.setVisAttributes(description.invisible());
//  sens.setType("tracker");*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) 
  {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    xml_comp_t trd   = x_mod.trd();

    double     posY;
    double     x1 = trd.x1();
    double     x2 = trd.x2();
    double     z  = trd.z();
    double     total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    double     thickness_so_far = 0.0;
    double     thickness_sum    = -total_thickness / 2.0;
    double     y1               = total_thickness / 2;
    double     y2               = total_thickness / 2;
    Trapezoid m_solid(x1, x2, y1, y2, z);
    Volume m_volume(m_nam, m_solid, radiatorMat);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));


    
    modules[m_nam] = m_volume;
  }
///////////////////////////////////////////////////////////////////////
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int        l_id    = x_layer.id();
    int        mod_num = 1;

    xml_comp_t l_env      = x_layer.child(_U(envelope));
    string     layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis    = l_env.attr<std::string>(_Unicode(vis));
    double      layer_rmin   = l_env.attr<double>(_Unicode(rmin));
    double      layer_rmax   = l_env.attr<double>(_Unicode(rmax));
    double      layer_length = l_env.attr<double>(_Unicode(length));
    double      layer_zstart = l_env.attr<double>(_Unicode(zstart));
    double      layer_center_z =  layer_zstart + layer_length/2.0;

    Tube       layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    Volume     layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv =
          assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_N";
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_P";
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);
    
 //   Acts::ActsExtension* layerExtension = new Acts::ActsExtension();

    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t  x_ring   = ri;
      double      r        = x_ring.r();
      double      phi0     = x_ring.phi0(0);
      double      zstart   = x_ring.zstart();
      double      dz       = x_ring.dz(0);
      int         nmodules = x_ring.nmodules();
      string      m_nam    = x_ring.moduleStr();
      Volume      m_vol    = modules[m_nam];
      double      iphi     = 2 * M_PI / nmodules;
      double      phi      = phi0;
      Placements& sensVols = sensitives[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        string     m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        double     x      = -r * std::cos(phi);
        double     y      = -r * std::sin(phi);

        if (!reflect) {
          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(
              m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, zstart + dz)));
          pv.addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
            comp_elt.setPlacement(sens_pv);
        //std::cout << " adding ACTS extension" << "\n";
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension("XZY");
            comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        } else {
          pv = layer_vol.placeVolume(
              m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2), Position(x, y, -zstart - dz)));
          pv.addPhysVolID("module", mod_num);
          DetElement r_module(layer_element, m_base + "_neg", det_id);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement   comp_elt(r_module, sens_pv.volume().name(), mod_num);
            comp_elt.setPlacement(sens_pv);
        //std::cout << " adding ACTS extension" << "\n";
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension("XZY");
            comp_elt.addExtension<Acts::ActsExtension>(moduleExtension);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        }
        dz = -dz;
        phi += iphi;
        ++mod_num;
      }
    }
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  pv = motherVol.placeVolume(assembly,Position(0,0,0) );
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
//DECLARE_DETELEMENT(athena_TrapEndcapTracker, create_detector)
DECLARE_DETELEMENT(athena_radiator, create_detector)
