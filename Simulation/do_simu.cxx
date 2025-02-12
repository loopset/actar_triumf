#ifndef triumf_cxx
#define triumf_cxx
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include <string>
void do_simu(const std::string& beam, const std::string& target, const std::string& light, double Tbeam, double Ex,
             bool inspect)
{
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    ActPhysics::SilSpecs sils;
    sils.ReadFile("../configs/silicons.conf");
    sils.DrawGeo();
}
#endif
