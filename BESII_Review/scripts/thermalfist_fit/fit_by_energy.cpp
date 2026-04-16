#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "HRGBase.h"
#include "HRGFit.h"

namespace {

struct EnergyInput {
  double energyGeV;
  std::string fitFile;
  int nPoints;
};

std::vector<EnergyInput> readManifest(const std::string &manifestPath) {
  std::vector<EnergyInput> rows;
  std::ifstream fin(manifestPath);
  if (!fin.is_open()) {
    throw std::runtime_error("Cannot open manifest: " + manifestPath);
  }

  std::string line;
  std::getline(fin, line); // header

  while (std::getline(fin, line)) {
    if (line.empty()) {
      continue;
    }
    std::stringstream ss(line);
    std::string energyStr, fileStr, pointsStr;
    if (!std::getline(ss, energyStr, ',')) {
      continue;
    }
    if (!std::getline(ss, fileStr, ',')) {
      continue;
    }
    if (!std::getline(ss, pointsStr, ',')) {
      continue;
    }

    EnergyInput row;
    row.energyGeV = std::stod(energyStr);
    row.fitFile = fileStr;
    row.nPoints = std::stoi(pointsStr);
    rows.push_back(row);
  }

  return rows;
}

std::vector<thermalfist::FittedQuantity> readQuantities(const std::string &fitFile) {
  std::vector<thermalfist::FittedQuantity> quantities;
  std::ifstream fin(fitFile);
  if (!fin.is_open()) {
    return quantities;
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    int fitFlag = 1;
    long long pdg1 = 0;
    long long pdg2 = 0;
    int feed1 = 3;
    int feed2 = 0;
    double value = 0.0;
    double error = 0.0;
    if (!(iss >> fitFlag >> pdg1 >> pdg2 >> feed1 >> feed2 >> value >> error)) {
      continue;
    }
    if (error <= 0.0) {
      continue;
    }

    if (pdg2 == 0) {
      thermalfist::FittedQuantity q(thermalfist::ExperimentMultiplicity(
          pdg1, value, error, static_cast<thermalfist::Feeddown::Type>(feed1)));
      q.toFit = (fitFlag != 0);
      quantities.push_back(q);
    } else {
      thermalfist::FittedQuantity q(thermalfist::ExperimentRatio(
          pdg1, pdg2, value, error,
          static_cast<thermalfist::Feeddown::Type>(feed1),
          static_cast<thermalfist::Feeddown::Type>(feed2)));
      q.toFit = (fitFlag != 0);
      quantities.push_back(q);
    }
  }
  return quantities;
}

void printUsage(const char *prog) {
  std::cerr << "Usage: " << prog << " <manifest.csv> <particle_list.dat> <output.csv> [points_outdir]\n";
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 4) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string manifestPath = argv[1];
  const std::string particleListPath = argv[2];
  const std::string outputPath = argv[3];
  const std::string pointsOutDir = (argc >= 5) ? argv[4] : "";

  std::vector<EnergyInput> energies;
  try {
    energies = readManifest(manifestPath);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 2;
  }

  if (energies.empty()) {
    std::cerr << "No energies found in manifest: " << manifestPath << "\n";
    return 3;
  }

  thermalfist::ThermalParticleSystem tps(particleListPath);
  thermalfist::ThermalModelIdeal model(&tps);

  model.SetUseWidth(thermalfist::ThermalParticle::BWTwoGamma);
  model.SetStatistics(true);
  model.ConstrainMuS(false);
  model.ConstrainMuQ(true);
  model.SetQoverB(0.4);

  thermalfist::ThermalModelFit fitter(&model);

  // Fit T, muB, muS, and R (R -> dV/dy conversion in output),
  // while keeping gammaS fixed to 1.
  fitter.SetParameterFitFlag("T", true);
  fitter.SetParameterFitFlag("muB", true);
  fitter.SetParameterFitFlag("muS", true);
  fitter.SetParameterFitFlag("gammaS", false);
  fitter.SetParameterFitFlag("R", true);

  fitter.SetParameterFitFlag("muQ", false);
  fitter.SetParameterFitFlag("muC", false);
  fitter.SetParameterFitFlag("Rc", false);
  fitter.SetParameterFitFlag("gammaq", false);
  fitter.SetParameterFitFlag("gammaC", false);
  fitter.SetParameterFitFlag("Tkin", false);

  std::ofstream fout(outputPath);
  if (!fout.is_open()) {
    std::cerr << "Cannot open output file: " << outputPath << "\n";
    return 4;
  }

  fout << "energy_GeV,n_points,T_MeV,T_err_MeV,muB_MeV,muB_err_MeV,muS_MeV,muS_err_MeV,gammaS,gammaS_err,R_fm,R_err_fm,dVdy_fm3,dVdy_err_fm3,chi2,ndf,chi2_ndf,input_file\n";

  constexpr double kPi = 3.14159265358979323846;
  if (!pointsOutDir.empty()) {
    std::filesystem::create_directories(pointsOutDir);
  }

  for (const auto &entry : energies) {
    if (std::abs(entry.energyGeV - 130.0) < 1e-6) {
      std::cout << "Skipping sqrt(sNN)=130 GeV by request.\n";
      continue;
    }

    const std::vector<thermalfist::FittedQuantity> quantities = readQuantities(entry.fitFile);
    if (quantities.empty()) {
      std::cerr << "Skipping " << entry.energyGeV << " GeV: no quantities loaded from " << entry.fitFile << "\n";
      continue;
    }

    fitter.SetQuantities(quantities);

    const double muBInit = 1.308 / (1.0 + 0.273 * entry.energyGeV);
    fitter.SetParameter("T", 0.150, 0.020, 0.090, 0.190);
    fitter.SetParameter("muB", muBInit, 0.050, 0.0, 0.900);
    fitter.SetParameter("muS", 0.0, 0.020, -0.300, 0.300);
    fitter.SetParameter("gammaS", 1.0, 0.0, 1.0, 1.0);
    fitter.SetParameter("R", 8.0, 1.0, 0.5, 20.0);

    const thermalfist::ThermalModelFitParameters result = fitter.PerformFit(false);
    const auto &fitted = fitter.FittedQuantities();

    const double R = result.R.value;
    const double Rerr = result.R.error;
    const double dVdy = (4.0 / 3.0) * kPi * R * R * R;
    const double dVdyErr = 4.0 * kPi * R * R * Rerr;

    fout << std::fixed << std::setprecision(6)
         << entry.energyGeV << ","
         << entry.nPoints << ","
         << result.T.value * 1000.0 << ","
         << result.T.error * 1000.0 << ","
         << result.muB.value * 1000.0 << ","
         << result.muB.error * 1000.0 << ","
         << result.muS.value * 1000.0 << ","
         << result.muS.error * 1000.0 << ","
         << result.gammaS.value << ","
         << result.gammaS.error << ","
         << result.R.value << ","
         << result.R.error << ","
         << dVdy << ","
         << dVdyErr << ","
         << result.chi2 << ","
         << result.ndf << ","
         << result.chi2ndf << ","
         << entry.fitFile
         << "\n";

    std::cout << "fit sqrt(sNN)=" << entry.energyGeV << " GeV: "
              << "T=" << result.T.value * 1000.0 << " MeV, "
              << "muB=" << result.muB.value * 1000.0 << " MeV, "
              << "muS=" << result.muS.value * 1000.0 << " MeV, "
              << "gammaS=" << result.gammaS.value << ", "
              << "dV/dy=" << dVdy << " fm^3, "
              << "chi2/ndf=" << result.chi2ndf << "\n";

    if (!pointsOutDir.empty()) {
      std::ostringstream eoss;
      eoss << std::fixed << std::setprecision(1) << entry.energyGeV;
      std::string etag = eoss.str();
      for (char &c : etag) {
        if (c == '.') c = 'p';
      }
      std::string pointsFile = pointsOutDir + "/sqrts_" + etag + "GeV_points.csv";
      std::ofstream fpt(pointsFile);
      if (fpt.is_open()) {
        fpt << "energy_GeV,pdg_id,particle_name,data_yield,data_err,model_yield\n";
        for (int i = 0; i < fitter.ModelDataSize(); ++i) {
          if (i >= static_cast<int>(fitted.size())) {
            continue;
          }
          const auto &q = fitted[i];
          if (q.type != thermalfist::FittedQuantity::Multiplicity) {
            continue;
          }
          const long long pdg = q.mult.fPDGID;
          const std::string name = model.TPS()->GetNameFromPDG(pdg);
          fpt << std::fixed << std::setprecision(6)
              << entry.energyGeV << ","
              << pdg << ","
              << "\"" << name << "\"" << ","
              << q.mult.fValue << ","
              << q.mult.fError << ","
              << fitter.ModelData(i) << "\n";
        }
      }
    }
  }

  std::cout << "Wrote fit summary: " << outputPath << "\n";
  return 0;
}
