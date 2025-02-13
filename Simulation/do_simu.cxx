#ifndef triumf_cxx
#define triumf_cxx
#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"

#include <iostream>
#include <string>

#include "../Histos.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

XYZPoint SampleVertex(ActRoot::TPCParameters* tpc)
{
    // Define sigmas along Y and Z
    double sigmaY {4};
    double sigmaZ {4};
    auto y {gRandom->Gaus(tpc->Y() / 2, sigmaY)};
    auto z {gRandom->Gaus(tpc->Z() / 2, sigmaZ)};
    auto x {gRandom->Uniform() * tpc->X()};
    return {x, y, z};
}

std::pair<double, double> SampleCM()
{
    auto theta {TMath::ACos(gRandom->Uniform(-1, 1))};
    auto phi {gRandom->Uniform(0, TMath::TwoPi())};
    return {theta, phi};
}

void do_simu(const std::string& beam, const std::string& target, const std::string& light, double Tbeam, double Ex,
             bool inspect)
{
    // Set number of iterations
    auto niter {static_cast<int>(1e6)};

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons.conf");
    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    for(auto& [name, layer] : sils->GetLayers())
        layer.MoveZTo(tpc.Z() / 2, {5, 6, 7, 8, 9});
    std::cout << "Sils Z centred at : " << tpc.Z() / 2 << " mm" << '\n';
    // This means: make the Z of silicons {5,6,...} be that zOfBeam.
    // shift the others accordingly
    // sils->DrawGeo();

    // SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", TString::Format("../SRIM files/%s_Butane_160Torr.txt", beam.c_str()).Data());
    srim->ReadTable("light", TString::Format("../SRIM files/%s_Butane_160Torr.txt", light.c_str()).Data());
    srim->ReadTable("lightInSil", TString::Format("../SRIM files/%s_silicon.txt", light.c_str()).Data());

    // Kinematics
    auto* kin {new ActPhysics::Kinematics {beam, target, light, Tbeam, Ex}};

    // Declare histograms
    auto hKin {Histos::Kin.GetHistogram()};
    auto hSP {Histos::SP.GetHistogram()};
    auto hRP {Histos::RP.GetHistogram()};
    auto hThetaCMAll {Histos::ThetaCM.GetHistogram()};
    auto hThetaCM {Histos::ThetaCM.GetHistogram()};

    for(int it = 0; it < niter; it++)
    {
        // Sample vertex
        auto vertex {SampleVertex(&tpc)};
        // Slow beam with straggling
        auto TbeamCorr {srim->SlowWithStraggling("beam", Tbeam, vertex.X())};
        kin->SetBeamEnergy(TbeamCorr);
        // Sample Ex/Ecm/thetaCM...
        // So far traditional approach. We have to change for Ecm sampling == sampling vertex.X() i think
        // thetaCM would be fixed in that case
        auto [thetaCM, phiCM] {SampleCM()};
        // Generate lab kinematics for protons
        kin->ComputeRecoilKinematics(thetaCM, phiCM);
        // Fill thetaCMall
        hThetaCMAll->Fill(thetaCM * TMath::RadToDeg());
        // Extract direction
        auto T3Lab {kin->GetT3Lab()};
        auto theta3Lab {kin->GetTheta3Lab()};
        auto phi3Lab {kin->GetPhi3Lab()};
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        // How to check whether tracks would read the silicons with new class:
        auto [silIndex0, silPoint0] {sils->FindSPInLayer("f0", vertex, direction)};
        // "f0": key name of layer to check for SP
        // silIndex == -1 if NO SP
        // else, returns the silicon index
        // silPoint0: SP in millimetres. in standard ACTAR frame (same as analysis)
        if(silIndex0 == -1)
            continue;

        // Test
        hKin->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
        hRP->Fill(vertex.X(), vertex.Y());
        hSP->Fill(silPoint0.Y(), silPoint0.Z());
        hThetaCM->Fill(thetaCM * TMath::RadToDeg()); // only thetaCm that enter our cuts
    }

    // Draw if not running for multiple Exs
    if(inspect)
    {
        auto* c0 {new TCanvas {"c0", "Sim inspect 0"}};
        c0->DivideSquare(6);
        c0->cd(1);
        hKin->DrawClone("colz");
        // Draw theo kin
        kin->SetBeamEnergyAndEx(Tbeam, Ex);
        auto* gtheo {kin->GetKinematicLine3()};
        gtheo->Draw("l");
        c0->cd(2);
        hSP->DrawClone("colz");
        sils->GetLayer("f0").GetSilMatrix()->Draw();
        c0->cd(3);
        hRP->DrawClone("colz");
        c0->cd(4);
        hThetaCMAll->SetTitle("All #theta_{CM}");
        hThetaCMAll->DrawClone();
        c0->cd(5);
        hThetaCM->SetTitle("#theta_{CM} in cuts");
        hThetaCM->DrawClone();
    }
}
#endif
