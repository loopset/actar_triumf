#ifndef triumf_cxx
#define triumf_cxx
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include <string>
using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;


void do_simu(const std::string& beam, const std::string& target, const std::string& light, double Tbeam, double Ex,
             bool inspect)
{
    // Set number of iterations
    auto niter {static_cast<int>(1e5)};

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons.conf");
    // We have to centre the silicons with the beam input
    double zOfBeam {256. / 2};
    for(auto& [name, layer] : sils->GetLayers())
        layer.MoveZTo(zOfBeam, {5, 6, 7, 8});
    // This means: make the Z of silicons {5,6,...} be that zOfBeam.
    // shift the others accordingly
    sils->DrawGeo();

    for(int it = 0; it < niter; it++)
    {
        // Sample vertex
        XYZPoint vertex {};
        // Slow with straggling beam
        // Sample Ex/Ecm/thetaCM...
        // Generate lab kinematics for protons
        XYZVector direction {};
        // How to check whether tracks would read the silicons with new class:
        auto [silIndex0, silPoint0] {sils->FindSPInLayer("f0", vertex, direction)};
        // "f0": key name of layer to check for SP
        // silIndex == -1 if NO SP
        // else, returns the silicon index
        // silPoint0: SP in millimetres. in standard ACTAR frame (same as analysis)
    }
}
#endif
