#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCollection.h"

#include "TTMParameter.h"
#include "TTMParameterSetBSQ.h"
#include "TTMParticle.h"
#include "TTMParticleSet.h"
#include "TTMYield.h"
#include "TThermalFitBSQ.h"

namespace {

struct EnergyInput {
  double energyGeV = 0.0;
  std::string fitFile;
  int nPoints = 0;
};

struct ObsPoint {
  int pdg = 0;
  std::string label;
  double data = 0.0;
  double error = 0.0;
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

std::vector<ObsPoint> readObsPoints(const std::string &fitFile) {
  std::vector<ObsPoint> points;
  std::ifstream fin(fitFile);
  if (!fin.is_open()) {
    return points;
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    ObsPoint p;
    if (!(iss >> p.pdg >> p.label >> p.data >> p.error)) {
      continue;
    }
    if (p.error <= 0.0) {
      continue;
    }
    points.push_back(p);
  }

  return points;
}

void printUsage(const char *prog) {
  std::cerr << "Usage: " << prog
            << " <manifest.csv> <particle_list.txt> <decays_dir> <output.csv> [points_outdir]\n";
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 5) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string manifestPath = argv[1];
  const std::string particleListPath = argv[2];
  const std::string decaysDir = argv[3];
  const std::string outputPath = argv[4];
  const std::string pointsOutDir = (argc >= 6) ? argv[5] : "";

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

  std::ofstream fout(outputPath);
  if (!fout.is_open()) {
    std::cerr << "Cannot open output file: " << outputPath << "\n";
    return 4;
  }

  if (!pointsOutDir.empty()) {
    std::filesystem::create_directories(pointsOutDir);
  }

  fout << "energy_GeV,n_points,T_MeV,T_err_MeV,muB_MeV,muB_err_MeV,muS_MeV,muS_err_MeV,gammaS,gammaS_err,R_fm,R_err_fm,dVdy_fm3,dVdy_err_fm3,chi2,ndf,chi2_ndf,input_file\n";

  constexpr double kPi = 3.14159265358979323846;
  for (const auto &entry : energies) {
    if (std::abs(entry.energyGeV - 130.0) < 1e-6) {
      std::cout << "Skipping sqrt(sNN)=130 GeV by request.\n";
      continue;
    }

    std::vector<ObsPoint> points = readObsPoints(entry.fitFile);
    if (points.empty()) {
      std::cerr << "Skipping " << entry.energyGeV << " GeV: no points from " << entry.fitFile << "\n";
      continue;
    }

    TTMParticleSet set(particleListPath.c_str(), true);
    set.InputDecays(decaysDir.c_str(), true);

    const double t0 = 0.150;
    const double muB0 = 1.308 / (1.0 + 0.273 * entry.energyGeV);
    const double muS0 = 0.0;
    const double muQ0 = 0.0;
    const double gammaS0 = 1.0;
    const double r0 = 8.0;

    TTMParameterSetBSQ par(t0, muB0, muS0, muQ0, gammaS0, r0);
    par.FitT(t0, 0.090, 0.190, 0.001);
    par.FitMuB(muB0, 0.0, 0.900, 0.001);
    par.FitMuS(muS0, 0.0, 0.300, 0.001);
    par.FixMuQ(0.0, 0.0);
    par.FixGammas(1.0, 0.0);
    par.FitRadius(r0, 0.5, 20.0, 0.05);
    par.FixMuC(0.0, 0.0);
    par.FixGammac(1.0, 0.0);
    par.FixMub(0.0, 0.0);
    par.FixGammab(1.0, 0.0);

    TTMThermalFitBSQ fit(&set, &par, entry.fitFile.c_str());
    fit.SetQStats(kFALSE);
    fit.SetWidth(kFALSE);
    fit.FitData(0);
    fit.GenerateYields();

    const TTMParameter *pT = par.GetTPar();
    const TTMParameter *pMuB = par.GetMuBPar();
    const TTMParameter *pMuS = par.GetMuSPar();
    const TTMParameter *pGammaS = par.GetGammasPar();
    const TTMParameter *pR = par.GetRadiusPar();

    const double dVdy = (4.0 / 3.0) * kPi * std::pow(pR->GetValue(), 3);
    const double dVdyErr = 4.0 * kPi * std::pow(pR->GetValue(), 2) * pR->GetError();

    int nFitPoints = 0;
    TIter next(fit.GetYields());
    while (TTMYield *y = static_cast<TTMYield *>(next())) {
      if (y->GetFit()) {
        ++nFitPoints;
      }
    }
    const int nFitParams = 4;
    const int ndf = nFitPoints - nFitParams;
    const double chi2 = fit.GetChiSquare();
    const double chi2ndf = (ndf > 0) ? (chi2 / static_cast<double>(ndf)) : 0.0;

    fout << std::fixed << std::setprecision(6) << entry.energyGeV << "," << nFitPoints << ","
         << pT->GetValue() * 1000.0 << "," << pT->GetError() * 1000.0 << ","
         << pMuB->GetValue() * 1000.0 << "," << pMuB->GetError() * 1000.0 << ","
         << pMuS->GetValue() * 1000.0 << "," << pMuS->GetError() * 1000.0 << ","
         << pGammaS->GetValue() << "," << pGammaS->GetError() << ","
         << pR->GetValue() << "," << pR->GetError() << ","
         << dVdy << "," << dVdyErr << ","
         << chi2 << "," << ndf << "," << chi2ndf << "," << entry.fitFile << "\n";

    std::cout << "fit sqrt(sNN)=" << entry.energyGeV << " GeV: "
              << "T=" << pT->GetValue() * 1000.0 << " MeV, "
              << "muB=" << pMuB->GetValue() * 1000.0 << " MeV, "
              << "muS=" << pMuS->GetValue() * 1000.0 << " MeV, "
              << "gammaS=" << pGammaS->GetValue() << ", "
              << "dV/dy=" << dVdy << " fm^3, "
              << "chi2/ndf=" << chi2ndf << "\n";

    if (!pointsOutDir.empty()) {
      std::ostringstream eoss;
      eoss << std::fixed << std::setprecision(1) << entry.energyGeV;
      std::string etag = eoss.str();
      for (char &c : etag) {
        if (c == '.') {
          c = 'p';
        }
      }
      std::string pointsFile = pointsOutDir + "/sqrts_" + etag + "GeV_points.csv";
      std::ofstream fpt(pointsFile);
      if (!fpt.is_open()) {
        continue;
      }
      fpt << "energy_GeV,pdg_id,particle_name,data_yield,data_err,model_yield\n";

      for (const auto &p : points) {
        TTMYield *y = fit.GetYield(p.pdg, 0, p.label.c_str());
        double modelValue = std::nan("");
        if (y) {
          modelValue = y->GetModelValue();
        }

        fpt << std::fixed << std::setprecision(6) << entry.energyGeV << ","
            << p.pdg << ",\"" << p.label << "\","
            << p.data << "," << p.error << "," << modelValue << "\n";
      }
    }
  }

  std::cout << "Wrote fit summary: " << outputPath << "\n";
  return 0;
}
