#ifndef triumf_cxx
#define triumf_cxx
#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>

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

void ApplyNaN(double& e, double t = 0, const std::string& comment = "stopped")
{
    if(e <= t)
        e = std::nan(comment.c_str());
}

void ApplyThetaRes(double& theta)
{
    double sigma {0.95 / 2.355}; // FWHM to sigma
    theta = gRandom->Gaus(theta, sigma * TMath::DegToRad());
}

void do_simu(const std::string& beam, const std::string& target, const std::string& light, double Tbeam, double Ex,
             const std::unordered_map<std::string, double>& opts, bool inspect)
{
    // Set number of iterations
    auto niter {static_cast<int>(1e6)};

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons.conf");
    const double sigmaSil {0.060 / 2.355}; // Si resolution
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return sigmaSil * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
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
    hKin->SetTitle("Sampled kinematics");
    auto hKinRec {Histos::Kin.GetHistogram()};
    hKinRec->SetTitle("Reconstructed kinematics");
    auto hSP {Histos::SP.GetHistogram()};
    auto hRP {Histos::RP.GetHistogram()};
    auto hThetaCMAll {Histos::ThetaCM.GetHistogram()};
    auto hThetaCM {Histos::ThetaCM.GetHistogram()};
    auto hEx {Histos::Ex.GetHistogram()};

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
        // Save without resolution
        auto theta3LabSampled {theta3Lab};
        // Apply angle resolution
        ApplyThetaRes(theta3Lab);
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

        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // CHeck if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
            continue;
        // Slow down in silicon
        auto normal {sils->GetLayer("f0").GetNormal()};
        auto angleWithNormal {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        auto T3AfterSil0 {srim->SlowWithStraggling("lightInSil", T3AtSil, sils->GetLayer("f0").GetUnit().GetThickness(),
                                                   angleWithNormal)};
        auto T3AfterSil0Res {gRandom->Gaus(T3AfterSil0, silRes->Eval(T3AfterSil0))};
        auto eLoss0 {T3AtSil - T3AfterSil0Res};
        ApplyNaN(eLoss0, sils->GetLayer("f0").GetThresholds().at(silIndex0));
        if(std::isnan(eLoss0))
            continue;

        // Reconstruct!
        bool isOk {T3AfterSil0 == 0}; // no punchthrouhg
        if(isOk)
        {
            // Assuming no punchthrough!
            auto T3Rec {srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R())};
            auto ExRec {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};

            // Fill
            hKin->Fill(theta3LabSampled * TMath::RadToDeg(), T3Lab); // before uncertainties implemented
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec);     // after reconstruction
            hEx->Fill(ExRec);
            hRP->Fill(vertex.X(), vertex.Y());
            hSP->Fill(silPoint0.Y(), silPoint0.Z());
            hThetaCM->Fill(thetaCM * TMath::RadToDeg()); // only thetaCm that enter our cuts
        }
    }

    // Compute efficiency
    auto* eff {new TEfficiency {*hThetaCM, *hThetaCMAll}};
    eff->SetNameTitle("eff", "#theta_{CM} efficiency;#epsilon;#theta_{CM} [#circ]");

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
        hKinRec->DrawClone("colz");
        gtheo->Draw("l");
        c0->cd(3);
        hSP->DrawClone("colz");
        sils->GetLayer("f0").GetSilMatrix()->Draw();
        c0->cd(4);
        hRP->DrawClone("colz");
        c0->cd(5);
        hThetaCMAll->SetTitle("All #theta_{CM}");
        hThetaCMAll->DrawClone();
        c0->cd(6);
        hThetaCM->SetTitle("#theta_{CM} in cuts");
        hThetaCM->DrawClone();

        auto* c1 {new TCanvas {"c1", "Sim inspect 1"}};
        c1->DivideSquare(6);
        c1->cd(1);
        eff->Draw("apl");
        c1->cd(2);
        hEx->DrawClone();
    }
}
#endif
