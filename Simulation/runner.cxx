#include "TString.h"

#include <string>
#include <unordered_map>
#include <vector>

#include "./do_simu.cxx"

void runner(TString what = "plot", bool inspect = true)
{
    std::string beam {"20Mg"};
    std::string target {"1H"};
    std::string light {"1H"};
    // Beam energy
    double Tbeam {130}; // MeV
    // Vector with Exs
    std::vector<double> Exs {0}; // 20Mg gs only
    // How to pass other options to simulation?
    // I propose this
    std::unordered_map<std::string, double> opts {{"pressure", 160}};

    // Run simu or plot
    if(what.Contains("simu"))
    {
        for(const auto& ex : Exs)
        {

            do_simu(beam, target, light, Tbeam, ex, opts, inspect);
            // auto str {TString::Format("root -l -b -x -q \'triumf.cxx(\"%s\",\"%s\",\"%s\",%f,%f,%d)\'", beam.c_str(),
            //                           target.c_str(), light.c_str(), Tbeam, ex, inspect)};
            if(inspect)
                break; // inspect: to debug simulation
        }
    }
    else
    {
        std::cout << "No plot method implemented yet" << '\n';
    }
}
