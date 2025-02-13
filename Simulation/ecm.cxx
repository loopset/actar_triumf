#include "ActKinematics.h"
#include "ActSRIM.h"

#include "TCanvas.h"
#include "TGraphErrors.h"

void ecm()
{
    // SRIM tables
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", "../SRIM files/20Mg_Butane_160Torr.txt");

    // Kinematics
    auto* kin {new ActPhysics::Kinematics {"20Mg(p,p)@130"}};
    double Tbeam {130};
    // Sum of masses
    auto massSum {kin->GetMass(1) + kin->GetMass(2)};

    auto* gcm {new TGraphErrors};
    gcm->SetTitle("Plain E_{CM};RP.X() [mm];E_{CM} [MeV]");
    auto* gcorr {new TGraphErrors};
    gcorr->SetTitle("#font[12]{Resonant} E_{CM};RP.X() [mm];E_{CM} [MeV]");
    for(auto g : {gcm, gcorr})
        g->SetLineWidth(2);

    double x0 {0};
    double x1 {256};
    double step {1};
    for(double x = x0; x <= x1; x += step)
    {
        auto TbeamCorr {srim->Slow("beam", Tbeam, x)};
        kin->SetBeamEnergy(TbeamCorr);
        auto cm {kin->GetECM()};
        gcm->AddPoint(x, cm);
        gcorr->AddPoint(x, cm - massSum);
        // gcorr->AddPoint(x, kin->GetResonantECM());
        // This is equivalent to calling now kin->GetResonantECM()
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "ECM tests"}};
    c0->DivideSquare(2);
    c0->cd(1);
    gcm->Draw("apl");
    c0->cd(2);
    gcorr->Draw("apl");
}
