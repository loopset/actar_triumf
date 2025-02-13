#ifndef Histos_h
#define Histos_h
#include "ROOT/RDF/HistoModels.hxx"
namespace Histos
{
// Ex
const ROOT::RDF::TH1DModel Ex {"hEx", "Ex;Ex [MeV];Counts", 300, -5, 10};

// Ecm
const ROOT::RDF::TH1DModel Ecm {"hEcm", "Ecm;E_{CM} [MeV];Counts", 300, -5, 10};

// SP
const ROOT::RDF::TH2DModel SP {"hSP", "SP;X or Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300};

// RP
const ROOT::RDF::TH2DModel RP {"hRP", "RP;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

// Kin
const ROOT::RDF::TH2DModel Kin {"hKin", "Kinematics;#theta_{light, Lab} [#circ];E_{light} [MeV]", 350, 0, 95, 350, 0,
                                40};

// Efficiency
const ROOT::RDF::TH1DModel ThetaCM {"hThetaCM", "ThetaCM;#theta_{CM} [#circ]", 600, 0, 180};
} // namespace Histos

#endif // !Histos_h
