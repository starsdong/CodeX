#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "TTMParameter.h"
#include "TTMParameterSetBQ.h"
#include "TTMParticleSet.h"
#include "TTMYield.h"
#include "TThermalFitBQ.h"

namespace {

struct PredictionInput {
  double energyGeV = 0.0;
  std::string predictionFile;
  double tMeV = 0.0;
  double muBMeV = 0.0;
  double rcFm = 0.0;
  double rFm = 0.0;
  double gammaS = 1.0;
  double b2q = 0.0;
  std::string interpolationNote;
};

struct ObsPoint {
  int pdg = 0;
  std::string label;
};

std::vector<PredictionInput> readManifest(const std::string &manifestPath) {
  std::vector<PredictionInput> rows;
  std::ifstream fin(manifestPath);
  if (!fin.is_open()) {
    throw std::runtime_error("Cannot open manifest: " + manifestPath);
  }

  std::string line;
  std::getline(fin, line);
  while (std::getline(fin, line)) {
    if (line.empty()) {
      continue;
    }
    std::stringstream ss(line);
    std::string field;
    PredictionInput row;
    std::getline(ss, field, ',');
    row.energyGeV = std::stod(field);
    std::getline(ss, row.predictionFile, ',');
    std::getline(ss, field, ',');
    row.tMeV = std::stod(field);
    std::getline(ss, field, ',');
    row.muBMeV = std::stod(field);
    std::getline(ss, field, ',');
    row.rcFm = std::stod(field);
    std::getline(ss, field, ',');
    row.rFm = std::stod(field);
    std::getline(ss, field, ',');
    row.gammaS = std::stod(field);
    std::getline(ss, field, ',');
    row.b2q = std::stod(field);
    std::getline(ss, row.interpolationNote);
    while (!row.interpolationNote.empty() &&
           (row.interpolationNote.back() == '\r' || row.interpolationNote.back() == '\n')) {
      row.interpolationNote.pop_back();
    }
    rows.push_back(row);
  }
  return rows;
}

std::vector<ObsPoint> readPredictionPoints(const std::string &predictionFile) {
  std::vector<ObsPoint> points;
  std::ifstream fin(predictionFile);
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
    double dummyValue = 0.0;
    double dummyError = 0.0;
    if (!(iss >> p.pdg >> p.label >> dummyValue >> dummyError)) {
      continue;
    }
    points.push_back(p);
  }
  return points;
}

std::string energyTag(double energyGeV) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << energyGeV;
  std::string tag = oss.str();
  for (char &c : tag) {
    if (c == '.') {
      c = 'p';
    }
  }
  return tag;
}

void printUsage(const char *prog) {
  std::cerr << "Usage: " << prog
            << " <manifest.csv> <particle_list.txt> <decays_dir> <summary.csv> <points_outdir>\n";
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 6) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string manifestPath = argv[1];
  const std::string particleListPath = argv[2];
  const std::string decaysDir = argv[3];
  const std::string summaryPath = argv[4];
  const std::string pointsOutDir = argv[5];

  std::vector<PredictionInput> energies;
  try {
    energies = readManifest(manifestPath);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 2;
  }
  if (energies.empty()) {
    std::cerr << "No prediction energies found in manifest: " << manifestPath << "\n";
    return 3;
  }

  std::ofstream fout(summaryPath);
  if (!fout.is_open()) {
    std::cerr << "Cannot open output file: " << summaryPath << "\n";
    return 4;
  }

  std::filesystem::create_directories(pointsOutDir);
  fout << "energy_GeV,T_MeV,muB_MeV,muQ_MeV,gammaS,Rc_fm,R_fm,dVdy_fm3,"
          "prediction_file,ensemble,muQ_constraint_B2Q,interpolation_note\n";

  constexpr double kPi = 3.14159265358979323846;
  for (const auto &entry : energies) {
    std::vector<ObsPoint> points = readPredictionPoints(entry.predictionFile);
    if (points.empty()) {
      std::cerr << "Skipping " << entry.energyGeV
                << " GeV: no prediction species in " << entry.predictionFile << "\n";
      continue;
    }

    TTMParticleSet set(particleListPath.c_str(), true);
    set.InputDecays(decaysDir.c_str(), true);

    TTMParameterSetBQ par(
        entry.tMeV / 1000.0,
        entry.muBMeV / 1000.0,
        0.0,
        entry.gammaS,
        entry.rcFm,
        entry.rFm);
    par.FixT(entry.tMeV / 1000.0, 0.0);
    par.FixMuB(entry.muBMeV / 1000.0, 0.0);
    par.ConstrainMuQ(entry.b2q);
    par.FixGammas(entry.gammaS, 0.0);
    par.FixCanRadius(entry.rcFm, 0.0);
    par.FixRadius(entry.rFm, 0.0);

    TTMThermalFitBQ fit(&set, &par, const_cast<char *>(entry.predictionFile.c_str()));
    fit.SetNonStrangeQStats(kFALSE);
    fit.SetWidth(kFALSE);
    fit.GenerateYields();

    const TTMParameter *pMuQ = par.GetMuQPar();
    const double dVdy = (4.0 / 3.0) * kPi * std::pow(entry.rFm, 3);

    fout << std::fixed << std::setprecision(6)
         << entry.energyGeV << ","
         << entry.tMeV << ","
         << entry.muBMeV << ","
         << pMuQ->GetValue() * 1000.0 << ","
         << entry.gammaS << ","
         << entry.rcFm << ","
         << entry.rFm << ","
         << dVdy << ","
         << entry.predictionFile << ",SCanonical,"
         << entry.b2q << ",\"" << entry.interpolationNote << "\"\n";

    std::string pointsFile = pointsOutDir + "/sqrts_" + energyTag(entry.energyGeV) + "GeV_points.csv";
    std::ofstream fpt(pointsFile);
    if (!fpt.is_open()) {
      std::cerr << "Cannot open points output: " << pointsFile << "\n";
      continue;
    }
    fpt << "energy_GeV,pdg_id,particle_name,model_yield\n";
    for (const auto &p : points) {
      TTMYield *y = fit.GetYield(p.pdg, 0, const_cast<char *>(p.label.c_str()));
      double modelValue = std::numeric_limits<double>::quiet_NaN();
      if (y) {
        modelValue = y->GetModelValue();
      }
      fpt << std::fixed << std::setprecision(6)
          << entry.energyGeV << "," << p.pdg << ",\"" << p.label << "\","
          << modelValue << "\n";
    }
  }

  std::cout << "Wrote SCE prediction summary: " << summaryPath << "\n";
  std::cout << "Wrote SCE prediction points: " << pointsOutDir << "\n";
  return 0;
}
