#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "TCollection.h"

#include "TTMParameter.h"
#include "TTMParameterSetBQ.h"
#include "TTMParticleSet.h"
#include "TTMYield.h"
#include "TThermalFitBQ.h"

namespace {

struct ObsPoint {
  int pdg = 0;
  std::string label;
  double data = 0.0;
  double error = 0.0;
};

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
            << " <fit_input.txt> <particle_list.txt> <decays_dir> <summary.csv> <points.csv>\n";
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 6) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string fitFile = argv[1];
  const std::string particleListPath = argv[2];
  const std::string decaysDir = argv[3];
  const std::string summaryPath = argv[4];
  const std::string pointsPath = argv[5];

  std::vector<ObsPoint> points = readObsPoints(fitFile);
  if (points.empty()) {
    std::cerr << "No valid points found in " << fitFile << "\n";
    return 2;
  }

  TTMParticleSet set(particleListPath.c_str(), true);
  set.InputDecays(decaysDir.c_str(), true);

  constexpr double kSqrtS = 3.0;
  constexpr double kB2Q = 197.0 / (2.0 * 79.0);
  constexpr double kPi = 3.14159265358979323846;

  const double t0 = 0.100;
  const double muB0 = 1.308 / (1.0 + 0.273 * kSqrtS);
  const double muQ0 = 0.0;
  const double gammaS0 = 1.0;
  const double rc0 = 3.0;
  const double r0 = 5.0;

  TTMParameterSetBQ par(t0, muB0, muQ0, gammaS0, rc0, r0);
  par.FitT(t0, 0.040, 0.190, 0.001);
  par.FitMuB(muB0, 0.050, 1.200, 0.001);
  par.ConstrainMuQ(kB2Q);
  par.FixGammas(1.0, 0.0);
  par.FitCanRadius(rc0, 0.1, 15.0, 0.05);
  par.FitRadius(r0, 0.1, 20.0, 0.05);

  TTMThermalFitBQ fit(&set, &par, const_cast<char *>(fitFile.c_str()));
  fit.SetNonStrangeQStats(kFALSE);
  fit.SetWidth(kFALSE);
  fit.FitData(0);
  fit.GenerateYields();

  const TTMParameter *pT = par.GetTPar();
  const TTMParameter *pMuB = par.GetMuBPar();
  const TTMParameter *pMuQ = par.GetMuQPar();
  const TTMParameter *pGammaS = par.GetGammasPar();
  const TTMParameter *pRc = par.GetCanRadiusPar();
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
  const int nFitParams = 4; // T, muB, Rc, R ; muQ constrained, gammaS fixed
  const int ndf = nFitPoints - nFitParams;
  const double chi2 = fit.GetChiSquare();
  const double chi2ndf = (ndf > 0) ? (chi2 / static_cast<double>(ndf)) : 0.0;

  {
    std::ofstream fout(summaryPath);
    if (!fout.is_open()) {
      std::cerr << "Cannot open summary output: " << summaryPath << "\n";
      return 3;
    }
    fout << "energy_GeV,n_points,T_MeV,T_err_MeV,muB_MeV,muB_err_MeV,muQ_MeV,muQ_err_MeV,gammaS,gammaS_err,"
            "Rc_fm,Rc_err_fm,R_fm,R_err_fm,dVdy_fm3,dVdy_err_fm3,chi2,ndf,chi2_ndf,input_file,ensemble,muQ_constraint_B2Q\n";
    fout << std::fixed << std::setprecision(6)
         << kSqrtS << "," << nFitPoints << ","
         << pT->GetValue() * 1000.0 << "," << pT->GetError() * 1000.0 << ","
         << pMuB->GetValue() * 1000.0 << "," << pMuB->GetError() * 1000.0 << ","
         << pMuQ->GetValue() * 1000.0 << "," << pMuQ->GetError() * 1000.0 << ","
         << pGammaS->GetValue() << "," << pGammaS->GetError() << ","
         << pRc->GetValue() << "," << pRc->GetError() << ","
         << pR->GetValue() << "," << pR->GetError() << ","
         << dVdy << "," << dVdyErr << ","
         << chi2 << "," << ndf << "," << chi2ndf << ","
         << fitFile << ",SCanonical," << kB2Q << "\n";
  }

  {
    std::filesystem::create_directories(std::filesystem::path(pointsPath).parent_path());
    std::ofstream fpt(pointsPath);
    if (!fpt.is_open()) {
      std::cerr << "Cannot open points output: " << pointsPath << "\n";
      return 4;
    }
    fpt << "energy_GeV,pdg_id,particle_name,data_yield,data_err,model_yield\n";
    for (const auto &p : points) {
      TTMYield *y = fit.GetYield(p.pdg, 0, const_cast<char *>(p.label.c_str()));
      double modelValue = std::numeric_limits<double>::quiet_NaN();
      if (y) {
        modelValue = y->GetModelValue();
      }
      fpt << std::fixed << std::setprecision(6)
          << kSqrtS << "," << p.pdg << ",\"" << p.label << "\","
          << p.data << "," << p.error << "," << modelValue << "\n";
    }
  }

  std::cout << std::fixed << std::setprecision(3)
            << "SCE fit sqrt(sNN)=" << kSqrtS << " GeV: "
            << "T=" << pT->GetValue() * 1000.0 << " MeV, "
            << "muB=" << pMuB->GetValue() * 1000.0 << " MeV, "
            << "muQ=" << pMuQ->GetValue() * 1000.0 << " MeV, "
            << "Rc=" << pRc->GetValue() << " fm, "
            << "R=" << pR->GetValue() << " fm, "
            << "chi2/ndf=" << chi2ndf << "\n";
  std::cout << "Wrote fit summary: " << summaryPath << "\n";
  std::cout << "Wrote points table: " << pointsPath << "\n";
  return 0;
}
