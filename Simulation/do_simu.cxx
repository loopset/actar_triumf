#ifndef triumf_cxx
#define triumf_cxx
#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

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

TF1* GetEcmSampler(std::string filename)
{
    auto graphEcm {new TGraphErrors("../Inputs/Mg20p_150_A_tot30keV.dat", "%lg %lg")};

    int nPoints = graphEcm->GetN();
    double lastX, lastY, firstX, firstY;
    graphEcm->GetPoint(0, firstX, firstY);
    graphEcm->GetPoint(nPoints - 1, lastX, lastY);

    std::cout<<nPoints<<std::endl;
    std::cout<<lastX<<std::endl;

    double step = 0.01;
    // Expand to the right if necessary
    if (lastX < 6.4)
    {
        for (double x = lastX + step; x <= 6.3; x += step)
        {
            graphEcm->AddPoint(x, lastY);
        }
    }
    // Expand to the left if necessary
    if (firstX > 0.1)
    {
        for (double x = firstX - step; x >= 0; x -= step)
        {
            graphEcm->AddPoint(x, firstY);
        }
    }

    int nPoints1 = graphEcm->GetN();
    double lastX1, lastY1;
    graphEcm->GetPoint(nPoints1 - 1, lastX1, lastY1);

    std::cout<<nPoints1<<std::endl;
    std::cout<<lastX1<<std::endl;

    auto functionEcm {
        new TF1("fEcm", [=](double* x, double* p) { return graphEcm->Eval(x[0], nullptr, "s"); }, 0, 6.2, 0)};
    return functionEcm;
}

TF1* GetEcmRPRelation(ActPhysics::Kinematics* kin, ActPhysics::SRIM* srim, double Tbeam, double Ex)
{
    auto* gRPxEcm {new TGraphErrors};
    gRPxEcm->SetTitle("#font[12]{Resonant} E_{CM};RP.X() [mm];E_{CM} [MeV]");
    gRPxEcm->SetLineWidth(2);
    auto* gEcmRPx {new TGraphErrors};

    double x0 {0};
    double x1 {256};
    double step {1};
    kin->GetParticle(1).Print();
    kin->GetParticle(2).Print();
    kin->GetParticle(3).Print();
    for(double x = x0; x <= x1; x += step)
    {
        auto TbeamCorr {srim->Slow("beam", Tbeam, x)};
        auto beamThreshold {ActPhysics::Kinematics(kin->GetParticle(1), kin->GetParticle(2), kin->GetParticle(3), -1, Ex).GetT1Thresh()};
        // std::cout<<"TbeamCorr: "<<TbeamCorr<< " Beam Threshold: "<< beamThreshold << " x: "<< x <<std::endl;
        if(std::isnan(TbeamCorr) || TbeamCorr < beamThreshold){
            continue;}
        kin->SetBeamEnergy(TbeamCorr);
        gRPxEcm->AddPoint(x, kin->GetResonantECM());
        if(TbeamCorr > 0)
            gEcmRPx->AddPoint(kin->GetResonantECM(), x);
    }

    auto funcEcmRPx {
        new TF1("funcEcmRPx", [=](double* x, double* p) { return gEcmRPx->Eval(x[0], nullptr, "S"); }, 0, 20, 0)};
    return funcEcmRPx;
}

void do_simu(const std::string& beam, const std::string& target, const std::string& light, double Tbeam, double Ex,
             const std::unordered_map<std::string, double>& opts, bool inspect)
{
    gRandom->SetSeed(0);
    // Set number of iterations
    auto niter {static_cast<int>(1e5)};

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
        layer.MoveZTo(tpc.Z() / 2, {4, 5, 6, 7});
    sils->DrawGeo();
    std::cout << "Sils Z centred at : " << tpc.Z() / 2 << " mm" << '\n';
    // This means: make the Z of silicons {5,6,...} be that zOfBeam.
    // shift the others accordingly
    // sils->DrawGeo();

    // SRIM
    double pressure {opts.at("pressure")};
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", TString::Format("../SRIM files/%s_Butane_%.0fmbar.txt", beam.c_str(), pressure).Data());
    srim->ReadTable("light", TString::Format("../SRIM files/%s_Butane_%.0fmbar.txt", light.c_str(), pressure).Data());
    srim->ReadTable("lightInSil", TString::Format("../SRIM files/%s_silicon.txt", light.c_str()).Data());

    // Kinematics
    auto* kin {new ActPhysics::Kinematics {beam, target, light, Tbeam, Ex}};

    // Ecm sampler and Ecm-rp.X relation
    auto* ecmSampler {GetEcmSampler("../Inputs/Mg20p_150_A_tot30keV.dat")};
    auto* eCMtoRP {GetEcmRPRelation(kin, srim, Tbeam, Ex)};

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
    auto hThetaCMThetaLab {Histos::ThetaCMThetaLab.GetHistogram()};
    auto hRPxEfirstSil {Histos::RP_E.GetHistogram()};
    hRPxEfirstSil->SetTitle("RPvsE if just 1 Sil");
    auto hRPxEbothSil {Histos::RP_E.GetHistogram()};
    hRPxEbothSil->SetTitle("RPvsE if 2 Sil");
    auto hRPxEStoppedGas {Histos::RP_E.GetHistogram()};
    auto hRPeffAll {Histos::RP_eff.GetHistogram()};
    auto hRPeffIn {Histos::RP_eff.GetHistogram()};

    // Output file
    TString fileName {
        TString::Format("../Outputs/%.1fAMeV/transfer_TRIUMF_%s_Eex_%.3f.root", Tbeam / 20, beam.c_str(), Ex)};
    auto* outFile {new TFile(fileName, "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Eex_tree {};
    outTree->Branch("Eex", &Eex_tree);
    double RP_tree {};
    outTree->Branch("RP", &RP_tree);
    double T3Rec_tree {};
    outTree->Branch("T3Rec", &T3Rec_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double phi3CM_tree {};
    outTree->Branch("phi3CM", &phi3CM_tree);
    double T4_tree {};
    outTree->Branch("T4Lab", &T4_tree);
    double range4_tree {};
    outTree->Branch("Range4", &range4_tree);

    for(int it = 0; it < niter; it++)
    {
        // Sample vertex and Ecm
        auto Ecm {ecmSampler->GetRandom()};
        // std::cout<<"Searching for: "<<Ecm<<std::endl;
        auto vertex {SampleVertex(&tpc)};
        vertex.SetX(eCMtoRP->Eval(Ecm));
        hRPeffAll->Fill(vertex.X());
        // std::cout<<vertex.X()<<std::endl;
        
        // Slow beam with straggling
        auto TbeamCorr {srim->SlowWithStraggling("beam", Tbeam, vertex.X())};
        auto beamThreshold {ActPhysics::Kinematics(beam, target, light, -1, Ex).GetT1Thresh()};
        if(std::isnan(TbeamCorr) || TbeamCorr < beamThreshold){
            continue;}
        kin->SetBeamEnergy(TbeamCorr);
        // Sample Ex/Ecm/thetaCM...
        // So far traditional approach. We have to change for Ecm sampling == sampling vertex.X() i think
        // thetaCM would be fixed in that case
        auto [thetaCM, phiCM] {SampleCM()};
        thetaCM = 150 * TMath::DegToRad(); // Fixed in 150 for Bea's data
        // Generate lab kinematics for protons
        kin->ComputeRecoilKinematics(thetaCM, phiCM);
        // Fill thetaCMall
        hThetaCMAll->Fill(kin->GetThetaCM() * TMath::RadToDeg());
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
        // Heavy particle information
        auto T4Lab {kin->GetT4Lab()};
        auto range4 {srim->EvalRange("beam", T4Lab)};

        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // Check if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
        {
            hRPxEStoppedGas->Fill(vertex.X(), T3Lab);
            continue;
        }
        // Slow down in silicon
        auto normal {sils->GetLayer("f0").GetNormal()};
        auto angleWithNormal {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        auto T3AfterSil0 {srim->SlowWithStraggling("lightInSil", T3AtSil, sils->GetLayer("f0").GetUnit().GetThickness(),
                                                   angleWithNormal)};
        auto eLoss0preSilRes {T3AtSil - T3AfterSil0};
        auto eLoss0 {gRandom->Gaus(eLoss0preSilRes, silRes->Eval(eLoss0preSilRes))}; // after silicon resolution
        ApplyNaN(eLoss0, sils->GetLayer("f0").GetThresholds().at(silIndex0));
        if(std::isnan(eLoss0))
            continue;


        // Apply 2nd layer of silicons
        double T3AfterInterGas {};
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double eLoss1 {};
        if(T3AfterSil0 > 0.) 
        {
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f1", vertex, direction);
            if(silIndex1 == -1)
            {} // If a silicon is not reached, don't continue with punchthough calculation
            else
            {
                T3AfterInterGas = {srim->SlowWithStraggling("light", T3AfterSil0, (silPoint0 - silPoint1).R())};
                if(T3AfterInterGas == 0)
                {} // If slow in gas don't continue with calculation
                else
                {
                    auto T3AfterSil1 {srim->SlowWithStraggling("lightInSil", T3AfterInterGas, sils->GetLayer("f1").GetUnit().GetThickness(),
                                                   angleWithNormal)};
                    auto eLoss1preSilRes {T3AfterInterGas - T3AfterSil1};
                    eLoss1 = gRandom->Gaus(eLoss1preSilRes, silRes->Eval(eLoss1preSilRes)); // after silicon resolution
                    ApplyNaN(eLoss1, sils->GetLayer("f1").GetThresholds().at(silIndex1));
                    if(std::isnan(eLoss1))
                        eLoss1 = 0;
                }
            }
        }
        // Reconstruct!
        bool isOk {T3AfterSil0 == 0}; // no punchthrouhg
        if(true)
        {
            // Assuming no punchthrough!
            double T3Rec {};
            if(eLoss1 == 0)
            {
                T3Rec = srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R());
                hRPxEfirstSil->Fill(vertex.X(), T3Rec);
            }
            else
            {
                auto T3Rec0 {srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R())};
                hRPxEfirstSil->Fill(vertex.X(), T3Rec0); // This is to show what would be obtained if there was just one sil

                // Reconstruction ob T3 with 2 silicon layers
                auto T3Rec1 {srim->EvalInitialEnergy("light", eLoss1, (silPoint1 - silPoint0).R())};
                T3Rec = srim->EvalInitialEnergy("light", eLoss0+T3Rec1, (silPoint0 - vertex).R());
            }

            auto ExRec {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};

            // Fill
            hKin->Fill(theta3LabSampled * TMath::RadToDeg(), T3Lab); // before uncertainties implemented
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec);     // after reconstruction
            hEx->Fill(ExRec);
            hRP->Fill(vertex.X(), vertex.Y());
            hSP->Fill(silPoint0.Y(), silPoint0.Z());
            hThetaCM->Fill(thetaCM * TMath::RadToDeg()); // only thetaCm that enter our cuts
            hThetaCMThetaLab->Fill(thetaCM * TMath::RadToDeg(), theta3Lab * TMath::RadToDeg());
            hRPxEbothSil->Fill(vertex.X(), T3Rec);
            hRPeffIn->Fill(vertex.X());

            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = thetaCM * TMath::RadToDeg();
            RP_tree = vertex.X();
            T3Rec_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phiCM;
            T4_tree = T4Lab;
            range4_tree = range4;
            outTree->Fill();
        }
    }

    // Compute efficiency
    auto* eff {new TEfficiency {*hThetaCM, *hThetaCMAll}};
    eff->SetNameTitle("eff", "#theta_{CM} efficiency;#epsilon;#theta_{CM} [#circ]");
    // Now efficiency in intervals
    auto* effIntervals {new TEfficiency {*hRPeffIn, *hRPeffAll}};
    effIntervals->SetNameTitle("effIntervals", "#RPx efficiency;#epsilon;#RPx [#mm]");

    // SAVING
    outFile->cd();
    outTree->Write();
    outFile->Close();
    delete outFile;
    outFile = nullptr;

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
        c1->cd(3);
        hThetaCMThetaLab->DrawClone("colz");
        auto* gCMLab {kin->GetThetaLabvsThetaCMLine()};
        gCMLab->Draw("l");
        c1->cd(4);
        hRPxEfirstSil->DrawClone();
        c1->cd(5);
        hRPxEbothSil->DrawClone();
        c1->cd(6);
        hRPxEStoppedGas->DrawClone();

        auto* cEff {new TCanvas {"cEff", "Eff in RPx intervals"}};
        effIntervals->Draw("apl");

        auto cEcm = new TCanvas("cEcm", "Resonant E_{CM} vs RP.X()");
        cEcm->DivideSquare(2);
        cEcm->cd(1);
        eCMtoRP->SetLineColor(kBlue);
        eCMtoRP->SetLineWidth(2);
        eCMtoRP->Draw();
        cEcm->cd(2);
        ecmSampler->Draw();
    }
}
#endif
