#include "ActKinematics.h"
#include "ActSRIM.h"
#include "TF1.h"

#include "TCanvas.h"
#include "TGraphErrors.h"

void ecm()
{
    // SRIM tables
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", "../SRIM files/20Mg_Butane_160mbar.txt");

    // Kinematics
    auto* kin {new ActPhysics::Kinematics {"20Mg(p,p)@130"}};
    double Tbeam {130};
    // Sum of masses
    auto massSum {kin->GetMass(1) + kin->GetMass(2)};

    auto* gcm {new TGraphErrors};
    gcm->SetTitle("Plain E_{CM};RP.X() [mm];E_{CM} [MeV]");
    auto* gcorr {new TGraphErrors};
    gcorr->SetTitle("#font[12]{Resonant} E_{CM};E_{CM} [MeV]; RP.X() [mm]");
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
        gcorr->AddPoint(cm - massSum, x);
        // gcorr->AddPoint(x, kin->GetResonantECM());
        // This is equivalent to calling now kin->GetResonantECM()
    }

    auto fEcmToRP {new TF1("fEcmToRP", [=](double* x, double* p){return gcorr->Eval(x[0], nullptr, "s");}, 0, 6, 0)};
    std::cout<<fEcmToRP->Eval(5)<<std::endl;

    // Draw
    auto* c0 {new TCanvas {"c0", "ECM tests"}};
    c0->DivideSquare(2);
    c0->cd(1);
    gcm->Draw("apl");
    c0->cd(2);
    gcorr->Draw("apl");
    fEcmToRP->Draw("same");
}
