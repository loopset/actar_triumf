#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <ROOT/RDataFrame.hxx>

#include "../Histos.h"

void separationNaMg()
{
    TString fileNa ("../Outputs/6.5AMeV/transfer_TRIUMF_20Na_Eex_0.000.root");
    TString fileMg ("../Outputs/6.5AMeV/transfer_TRIUMF_20Mg_Eex_0.000.root");
    TString fileMg2plus ("../Outputs/6.5AMeV/transfer_TRIUMF_20Mg_Eex_1.598.root");

    ROOT::RDataFrame dfNa ("SimulationTTree", fileNa);
    ROOT::RDataFrame dfMg ("SimulationTTree", fileMg);
    ROOT::RDataFrame dfMg2plus ("SimulationTTree", fileMg2plus);

    auto h2Na {Histos::RP_E.GetHistogram()};
    h2Na->SetTitle("20Na Data");
    h2Na->SetXTitle("Range of heavy [mm]");
    auto h2Mg {Histos::RP_E.GetHistogram()};
    h2Mg->SetTitle("20Mg Data");
    h2Mg->SetXTitle("Range of heavy [mm]");
    auto h2Mg2plus {Histos::RP_E.GetHistogram()};
    h2Mg2plus->SetTitle("20Mg 2+ Data");
    h2Mg2plus->SetXTitle("Range of heavy [mm]");

    dfNa.Foreach([&](double x, double y) { h2Na->Fill(x, y); }, {"Range4", "T3Rec"});
    dfMg.Foreach([&](double x, double y) { h2Mg->Fill(x, y); }, {"Range4", "T3Rec"});
    dfMg2plus.Foreach([&](double x, double y) { h2Mg2plus->Fill(x, y); }, {"Range4", "T3Rec"});

    TCanvas *c = new TCanvas("c", "Separation Na vs Mg", 800, 600);
    c->cd();
    
    h2Na->SetMarkerColor(kRed);
    h2Na->DrawClone("COLZ");

    h2Mg->SetMarkerColor(kBlue);
    h2Mg->DrawClone("SAME COLZ");

    h2Mg2plus->SetMarkerColor(kGreen);
    h2Mg2plus->DrawClone("SAME COLZ");
}