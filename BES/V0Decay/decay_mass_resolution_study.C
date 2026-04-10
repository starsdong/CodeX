#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TVector3.h"

struct DecayConfig {
  Long64_t nEvents = 500000;
  double motherMass = 5.279;
  double daughterMass1 = 0.13957;
  double daughterMass2 = 0.13957;

  double etaMin = -1.0;
  double etaMax = 1.0;
  double phiMin = -TMath::Pi();
  double phiMax = TMath::Pi();

  double ptMin = 0.2;
  double ptMax = 5.0;
  double ptSpectrumScale = 1.0;

  double Cpt = 0.005;      // sigma(pT)/pT = Cpt * pT
  double Ctheta = 0.0008;  // sigma(theta) constant [rad]
  double Cphi = 0.0008;    // sigma(phi) constant [rad]

  int nPtBins = 20;
  int nMassBins = 200;     // for truth/reco mass histograms
  double dMMin = -0.01;    // -10 MeV
  double dMMax = 0.01;     // +10 MeV
  bool autoDeltaMRange = false;
  double deltaMBinWidthTarget = 1.0e-5; // 10 keV
  int nDeltaMBinsMin = 120;
  int nDeltaMBinsMax = 6000;
  Long64_t nPilotForRange = 100000;

  int randomSeed = 12345;
};

struct DecayChannel {
  std::string name;
  double motherMass;
  double daughterMass1;
  double daughterMass2;
};

TLorentzVector SmearTrack(const TLorentzVector &truth, double mass, TRandom3 &rng,
                          double Cpt, double Ctheta, double Cphi) {
  const TVector3 p3 = truth.Vect();
  const double pt = p3.Pt();
  const double theta = p3.Theta();
  const double phi = p3.Phi();

  const double relSigmaPt = std::max(1e-12, Cpt * pt);
  double ptReco = pt * (1.0 + rng.Gaus(0.0, relSigmaPt));
  if (ptReco < 1e-6) ptReco = 1e-6;

  double thetaReco = theta + rng.Gaus(0.0, Ctheta);
  thetaReco = std::clamp(thetaReco, 1e-6, TMath::Pi() - 1e-6);

  const double phiReco = TVector2::Phi_mpi_pi(phi + rng.Gaus(0.0, Cphi));

  const double sinTheta = std::sin(thetaReco);
  const double pReco = ptReco / std::max(1e-9, sinTheta);
  const double px = ptReco * std::cos(phiReco);
  const double py = ptReco * std::sin(phiReco);
  const double pz = pReco * std::cos(thetaReco);
  const double eReco = std::sqrt(pReco * pReco + mass * mass);

  return TLorentzVector(px, py, pz, eReco);
}

std::pair<TLorentzVector, TLorentzVector> TwoBodyDecayRestFrame(double motherMass,
                                                                 double m1,
                                                                 double m2,
                                                                 TRandom3 &rng) {
  const double term1 = motherMass * motherMass - (m1 + m2) * (m1 + m2);
  const double term2 = motherMass * motherMass - (m1 - m2) * (m1 - m2);
  const double pstar = 0.5 / motherMass * std::sqrt(std::max(0.0, term1 * term2));

  const double cosTheta = rng.Uniform(-1.0, 1.0);
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const double phi = rng.Uniform(-TMath::Pi(), TMath::Pi());

  const double px = pstar * sinTheta * std::cos(phi);
  const double py = pstar * sinTheta * std::sin(phi);
  const double pz = pstar * cosTheta;

  TLorentzVector d1(px, py, pz, std::sqrt(pstar * pstar + m1 * m1));
  TLorentzVector d2(-px, -py, -pz, std::sqrt(pstar * pstar + m2 * m2));
  return {d1, d2};
}

double SampleMotherPt(TRandom3 &rng, const DecayConfig &cfg) {
  // Requested spectrum: f(pT) ~ pT * exp(-pT/scale), sampled as
  // a Gamma(k=2, theta=scale) variate and truncated to [ptMin, ptMax].
  while (true) {
    const double u1 = std::max(1e-12, rng.Uniform());
    const double u2 = std::max(1e-12, rng.Uniform());
    const double pt = -cfg.ptSpectrumScale * std::log(u1 * u2);
    if (pt >= cfg.ptMin && pt <= cfg.ptMax) return pt;
  }
}

void GenerateEventObservables(const DecayConfig &cfg, TRandom3 &rng,
                              double &pt, double &eta, double &phi,
                              double &mTruth, double &mReco, double &dM) {
  pt = SampleMotherPt(rng, cfg);
  eta = rng.Uniform(cfg.etaMin, cfg.etaMax);
  phi = rng.Uniform(cfg.phiMin, cfg.phiMax);

  TLorentzVector pMother;
  pMother.SetPtEtaPhiM(pt, eta, phi, cfg.motherMass);

  auto [d1Rest, d2Rest] = TwoBodyDecayRestFrame(cfg.motherMass, cfg.daughterMass1,
                                                cfg.daughterMass2, rng);

  const TVector3 beta = pMother.BoostVector();
  TLorentzVector d1Truth = d1Rest;
  TLorentzVector d2Truth = d2Rest;
  d1Truth.Boost(beta);
  d2Truth.Boost(beta);

  TLorentzVector d1Reco = SmearTrack(d1Truth, cfg.daughterMass1, rng, cfg.Cpt, cfg.Ctheta, cfg.Cphi);
  TLorentzVector d2Reco = SmearTrack(d2Truth, cfg.daughterMass2, rng, cfg.Cpt, cfg.Ctheta, cfg.Cphi);

  const TLorentzVector pMotherTruth = d1Truth + d2Truth;
  const TLorentzVector pMotherReco = d1Reco + d2Reco;

  mTruth = pMotherTruth.M();
  mReco = pMotherReco.M();
  dM = mReco - mTruth;
}

void run_decay_study(const char *outFile = "decay_study.root", Long64_t nEvents = 500000,
                     double motherMass = 5.279, double daughterMass1 = 0.13957,
                     double daughterMass2 = 0.13957, double Cpt = 0.005,
                     double Ctheta = 0.0008, double Cphi = 0.0008,
                     int randomSeed = 12345) {
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  DecayConfig cfg;
  cfg.nEvents = nEvents;
  cfg.motherMass = motherMass;
  cfg.daughterMass1 = daughterMass1;
  cfg.daughterMass2 = daughterMass2;
  cfg.Cpt = Cpt;
  cfg.Ctheta = Ctheta;
  cfg.Cphi = Cphi;
  cfg.ptMin = 0.2;
  cfg.ptMax = 5.0;
  cfg.ptSpectrumScale = 0.5 * cfg.motherMass;
  cfg.randomSeed = randomSeed;

  if (cfg.motherMass < cfg.daughterMass1 + cfg.daughterMass2) {
    std::cerr << "Invalid masses: motherMass < m1+m2. Abort." << std::endl;
    return;
  }

  double dMMin = cfg.dMMin;
  double dMMax = cfg.dMMax;
  int nDeltaMBins = std::clamp(
      static_cast<int>(std::ceil((dMMax - dMMin) / cfg.deltaMBinWidthTarget)),
      cfg.nDeltaMBinsMin, cfg.nDeltaMBinsMax);

  if (cfg.autoDeltaMRange) {
    const Long64_t nPilot = std::clamp(cfg.nPilotForRange, Long64_t(20000), cfg.nEvents);
    TRandom3 rngPilot(cfg.randomSeed + 100003);

    double mean = 0.0;
    double m2 = 0.0;
    Long64_t n = 0;

    for (Long64_t i = 0; i < nPilot; ++i) {
      double pt = 0.0, eta = 0.0, phi = 0.0, mTruth = 0.0, mReco = 0.0, dM = 0.0;
      GenerateEventObservables(cfg, rngPilot, pt, eta, phi, mTruth, mReco, dM);
      ++n;
      const double delta = dM - mean;
      mean += delta / static_cast<double>(n);
      m2 += delta * (dM - mean);
    }

    const double sigma = (n > 1) ? std::sqrt(std::max(1e-16, m2 / static_cast<double>(n - 1))) : 2.5e-4;
    const double half = 4.0 * sigma;
    dMMin = mean - half;
    dMMax = mean + half;

    const double minWidth = 2.0 * cfg.deltaMBinWidthTarget * cfg.nDeltaMBinsMin;
    if ((dMMax - dMMin) < minWidth) {
      const double mid = 0.5 * (dMMin + dMMax);
      dMMin = mid - 0.5 * minWidth;
      dMMax = mid + 0.5 * minWidth;
    }

    nDeltaMBins = std::clamp(
        static_cast<int>(std::ceil((dMMax - dMMin) / cfg.deltaMBinWidthTarget)),
        cfg.nDeltaMBinsMin, cfg.nDeltaMBinsMax);
  }

  TRandom3 rng(cfg.randomSeed);

  const double massWindowHalf = 0.02;

  auto *hPtMotherTruth = new TH1D("hPtMotherTruth", "Mother truth p_{T};p_{T}^{truth} [GeV];Events",
                                  cfg.nPtBins, cfg.ptMin, cfg.ptMax);
  auto *hEtaMotherTruth = new TH1D("hEtaMotherTruth", "Mother truth #eta;#eta^{truth};Events",
                                   120, cfg.etaMin, cfg.etaMax);
  auto *hPhiMotherTruth = new TH1D("hPhiMotherTruth", "Mother truth #phi;#phi^{truth} [rad];Events",
                                   128, -TMath::Pi(), TMath::Pi());
  auto *hMassTruth = new TH1D("hMassTruth", "Truth mother mass;M_{truth} [GeV];Events",
                              cfg.nMassBins, cfg.motherMass - massWindowHalf, cfg.motherMass + massWindowHalf);
  auto *hMassReco = new TH1D("hMassReco", "Reco mother mass;M_{reco} [GeV];Events",
                             cfg.nMassBins, cfg.motherMass - massWindowHalf, cfg.motherMass + massWindowHalf);
  auto *hDeltaMass = new TH1D("hDeltaMass", "Mass residual;#Delta M=M_{reco}-M_{truth} [GeV];Events",
                              nDeltaMBins, dMMin, dMMax);

  auto *h2DeltaMassVsPt = new TH2D("h2DeltaMassVsPt",
                                   "Mass residual vs mother p_{T}^{truth};p_{T}^{truth} [GeV];#Delta M [GeV]",
                                   cfg.nPtBins, cfg.ptMin, cfg.ptMax,
                                   nDeltaMBins, dMMin, dMMax);

  for (Long64_t i = 0; i < cfg.nEvents; ++i) {
    double pt = 0.0, eta = 0.0, phi = 0.0, mTruth = 0.0, mReco = 0.0, dM = 0.0;
    GenerateEventObservables(cfg, rng, pt, eta, phi, mTruth, mReco, dM);

    hPtMotherTruth->Fill(pt);
    hEtaMotherTruth->Fill(eta);
    hPhiMotherTruth->Fill(phi);
    hMassTruth->Fill(mTruth);
    hMassReco->Fill(mReco);
    hDeltaMass->Fill(dM);
    h2DeltaMassVsPt->Fill(pt, dM);
  }

  TFile fout(outFile, "RECREATE");
  hPtMotherTruth->Write();
  hEtaMotherTruth->Write();
  hPhiMotherTruth->Write();
  hMassTruth->Write();
  hMassReco->Write();
  hDeltaMass->Write();
  h2DeltaMassVsPt->Write();
  fout.Close();
  TH1::AddDirectory(oldAddDir);

  std::cout << "Wrote simulation output to: " << outFile << std::endl;
  std::cout << "Used pT shape ~ pT*exp(-pT/scale), scale=" << cfg.ptSpectrumScale
            << " GeV, pT range=[" << cfg.ptMin << ", " << cfg.ptMax << "] GeV" << std::endl;
  std::cout << (cfg.autoDeltaMRange ? "Auto" : "Fixed")
            << " #DeltaM range: [" << dMMin << ", " << dMMax << "] GeV with "
            << nDeltaMBins << " bins (~"
            << 1.0e6 * (dMMax - dMMin) / std::max(1, nDeltaMBins) << " keV/bin)" << std::endl;
}

void BuildOffsetGraphFromMode(TH2D *h2, TGraphErrors *g) {
  const int nBinsX = h2->GetNbinsX();
  int p = 0;
  for (int ib = 1; ib <= nBinsX; ++ib) {
    TH1D *proj = h2->ProjectionY(Form("hMode_proj_%d", ib), ib, ib, "e");
    if (proj->GetEntries() < 20) {
      delete proj;
      continue;
    }

    const int maxBin = proj->GetMaximumBin();
    const double x = h2->GetXaxis()->GetBinCenter(ib);
    const double ex = 0.5 * h2->GetXaxis()->GetBinWidth(ib);
    const double y = proj->GetXaxis()->GetBinCenter(maxBin);
    const double ey = 0.5 * proj->GetXaxis()->GetBinWidth(maxBin);

    g->SetPoint(p, x, y);
    g->SetPointError(p, ex, ey);
    ++p;
    delete proj;
  }
}

void BuildOffsetGraphFromGaussian(TH2D *h2, TGraphErrors *g, bool constrainedRange) {
  const int nBinsX = h2->GetNbinsX();
  int p = 0;

  for (int ib = 1; ib <= nBinsX; ++ib) {
    TH1D *proj = h2->ProjectionY(Form("hGaus_proj_%d_%d", constrainedRange ? 1 : 0, ib), ib, ib, "e");
    if (proj->GetEntries() < 40) {
      delete proj;
      continue;
    }

    const double yBin = proj->GetXaxis()->GetBinWidth(1);
    const double mean0 = proj->GetMean();
    const double rms0 = std::max(2.0 * yBin, proj->GetRMS());

    TF1 gfit("gfit", "gaus", mean0 - 2.5 * rms0, mean0 + 2.5 * rms0);
    gfit.SetParameters(std::max(1.0, proj->GetMaximum()), mean0, rms0);
    TFitResultPtr r0 = proj->Fit(&gfit, "QSN0");
    if (int(r0) != 0) {
      delete proj;
      continue;
    }

    double fitMean = gfit.GetParameter(1);
    double fitSigma = std::fabs(gfit.GetParameter(2));

    if (constrainedRange) {
      const double halfRange = std::max(4.0 * yBin, 1.2 * std::max(1e-12, fitSigma));
      TF1 gfit2("gfit2", "gaus", fitMean - halfRange, fitMean + halfRange);
      gfit2.SetParameters(gfit.GetParameter(0), fitMean, fitSigma);
      TFitResultPtr r1 = proj->Fit(&gfit2, "QSN0");
      if (int(r1) == 0) {
        fitMean = gfit2.GetParameter(1);
        fitSigma = std::fabs(gfit2.GetParameter(2));
      }
    }

    const double x = h2->GetXaxis()->GetBinCenter(ib);
    const double ex = 0.5 * h2->GetXaxis()->GetBinWidth(ib);
    const double ey = std::max(yBin / std::sqrt(12.0), fitSigma / std::sqrt(std::max(1.0, proj->GetEntries())));

    g->SetPoint(p, x, fitMean);
    g->SetPointError(p, ex, ey);
    ++p;
    delete proj;
  }
}

double RoundTo10keV(double xGeV) {
  const double step = 1e-5; // 10 keV in GeV
  return step * std::round(xGeV / step);
}

void PrintGraphInKeV(const TGraphErrors *g, const char *label) {
  std::cout << "\n" << label << " offsets (rounded to 10 keV):" << std::endl;
  std::cout << "  pT_center[GeV]  offset[keV]  statErr[keV]" << std::endl;
  for (int i = 0; i < g->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    g->GetPoint(i, x, y);
    const double yRounded = RoundTo10keV(y);
    const double ey = g->GetErrorY(i);
    std::cout << "  " << std::fixed << std::setprecision(3) << x
              << "          " << std::setprecision(1) << yRounded * 1.0e6
              << "        " << ey * 1.0e6 << std::endl;
  }
}

void DrawAndWriteGraph(TGraphErrors *g, const char *name, const char *title,
                       int color, const char *outPdf) {
  g->SetName(name);
  g->SetTitle(title);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);

  TCanvas c(name, name, 900, 700);
  c.SetGrid();
  g->Draw("APL");
  c.SaveAs(outPdf);
}

TGraphErrors *ConvertGraphYToMeV(const TGraphErrors *gIn, const char *name) {
  auto *gOut = new TGraphErrors();
  gOut->SetName(name);
  for (int i = 0; i < gIn->GetN(); ++i) {
    double x = 0.0, y = 0.0;
    gIn->GetPoint(i, x, y);
    gOut->SetPoint(i, x, 1000.0 * y);
    gOut->SetPointError(i, gIn->GetErrorX(i), 1000.0 * gIn->GetErrorY(i));
  }
  return gOut;
}

void analyze_decay_study(const char *inFile = "decay_study.root",
                         const char *plotPrefix = "decay_mass_offset") {
  gStyle->SetOptStat(0);

  TFile fin(inFile, "READ");
  if (fin.IsZombie()) {
    std::cerr << "Cannot open input file: " << inFile << std::endl;
    return;
  }

  auto *h2 = dynamic_cast<TH2D *>(fin.Get("h2DeltaMassVsPt"));
  if (!h2) {
    std::cerr << "Histogram h2DeltaMassVsPt is missing in " << inFile << std::endl;
    return;
  }

  auto *gMode = new TGraphErrors();
  auto *gGaus = new TGraphErrors();
  auto *gGausCR = new TGraphErrors();

  BuildOffsetGraphFromMode(h2, gMode);
  BuildOffsetGraphFromGaussian(h2, gGaus, false);
  BuildOffsetGraphFromGaussian(h2, gGausCR, true);

  auto *gModeMeV = ConvertGraphYToMeV(gMode, "gOffsetMode");
  auto *gGausMeV = ConvertGraphYToMeV(gGaus, "gOffsetGaus");
  auto *gGausCRMeV = ConvertGraphYToMeV(gGausCR, "gOffsetGausConstrained");

  const std::string p1 = std::string(plotPrefix) + "_mode.pdf";
  const std::string p2 = std::string(plotPrefix) + "_gaus.pdf";
  const std::string p3 = std::string(plotPrefix) + "_gaus_constrained.pdf";

  DrawAndWriteGraph(gModeMeV, "gOffsetMode",
                    "Mass offset vs p_{T} (peak/bin mode);p_{T}^{truth} [GeV];Offset #DeltaM [MeV]",
                    kBlue + 1, p1.c_str());

  DrawAndWriteGraph(gGausMeV, "gOffsetGaus",
                    "Mass offset vs p_{T} (Gaussian fit);p_{T}^{truth} [GeV];Offset #DeltaM [MeV]",
                    kRed + 1, p2.c_str());

  DrawAndWriteGraph(gGausCRMeV, "gOffsetGausConstrained",
                    "Mass offset vs p_{T} (constrained Gaussian fit);p_{T}^{truth} [GeV];Offset #DeltaM [MeV]",
                    kGreen + 2, p3.c_str());

  TCanvas cCmp("cOffsetCompare", "cOffsetCompare", 1000, 750);
  cCmp.SetGrid();
  gModeMeV->SetTitle("Mass offset comparison vs p_{T};p_{T}^{truth} [GeV];Offset #DeltaM [MeV]");
  gModeMeV->Draw("APL");
  gGausMeV->SetLineColor(kRed + 1);
  gGausMeV->SetMarkerColor(kRed + 1);
  gGausMeV->SetMarkerStyle(21);
  gGausMeV->Draw("PL SAME");
  gGausCRMeV->SetLineColor(kGreen + 2);
  gGausCRMeV->SetMarkerColor(kGreen + 2);
  gGausCRMeV->SetMarkerStyle(22);
  gGausCRMeV->Draw("PL SAME");

  TLegend leg(0.14, 0.72, 0.50, 0.89);
  leg.AddEntry(gModeMeV, "Peak/bin mode", "lp");
  leg.AddEntry(gGausMeV, "Gaussian fit", "lp");
  leg.AddEntry(gGausCRMeV, "Constrained Gaussian fit", "lp");
  leg.Draw();

  const std::string p4 = std::string(plotPrefix) + "_compare.pdf";
  cCmp.SaveAs(p4.c_str());

  fin.Close();

  TFile fout(inFile, "UPDATE");
  gModeMeV->Write("gOffsetMode", TObject::kOverwrite);
  gGausMeV->Write("gOffsetGaus", TObject::kOverwrite);
  gGausCRMeV->Write("gOffsetGausConstrained", TObject::kOverwrite);
  fout.Close();

  std::cout << "Analysis plots written: " << p1 << ", " << p2
            << ", " << p3 << ", " << p4 << std::endl;
  std::cout << "Graphs saved into ROOT file: " << inFile << std::endl;

  PrintGraphInKeV(gMode, "Mode");
  PrintGraphInKeV(gGaus, "Gaussian");
  PrintGraphInKeV(gGausCR, "Constrained Gaussian");
}

void plot_mass_diff_fits_per_ptbin(const char *inFile = "decay_study.root",
                                   const char *outPdf = "mass_diff_fits_per_ptbin.pdf") {
  gStyle->SetOptStat(0);

  TFile fin(inFile, "READ");
  if (fin.IsZombie()) {
    std::cerr << "Cannot open input file: " << inFile << std::endl;
    return;
  }

  auto *h2 = dynamic_cast<TH2D *>(fin.Get("h2DeltaMassVsPt"));
  if (!h2) {
    std::cerr << "Histogram h2DeltaMassVsPt is missing in " << inFile << std::endl;
    return;
  }

  const int nBinsX = h2->GetNbinsX();
  TCanvas c("cMassDiffFits", "cMassDiffFits", 900, 700);
  c.SetGrid();
  c.SaveAs((std::string(outPdf) + "(").c_str());

  for (int ib = 1; ib <= nBinsX; ++ib) {
    TH1D *projGeV = h2->ProjectionY(Form("hDM_ptbin_%d", ib), ib, ib, "e");
    if (projGeV->GetEntries() < 50) {
      delete projGeV;
      continue;
    }

    const int nY = projGeV->GetNbinsX();
    TH1D *proj = new TH1D(Form("hDM_ptbin_mev_%d", ib),
                          Form("#DeltaM in p_{T} bin %d;#DeltaM [MeV];Counts", ib),
                          nY,
                          1000.0 * projGeV->GetXaxis()->GetXmin(),
                          1000.0 * projGeV->GetXaxis()->GetXmax());
    for (int j = 1; j <= nY; ++j) {
      proj->SetBinContent(j, projGeV->GetBinContent(j));
      proj->SetBinError(j, projGeV->GetBinError(j));
    }
    delete projGeV;

    const double yBin = proj->GetXaxis()->GetBinWidth(1);
    const double mean0 = proj->GetMean();
    const double rms0 = std::max(2.0 * yBin, proj->GetRMS());
    const double mode0 = proj->GetXaxis()->GetBinCenter(proj->GetMaximumBin());

    TF1 fStd(Form("fStd_%d", ib), "gaus", mean0 - 2.5 * rms0, mean0 + 2.5 * rms0);
    fStd.SetParameters(std::max(1.0, proj->GetMaximum()), mean0, rms0);
    const int fitStatusStd = int(proj->Fit(&fStd, "QSN0"));

    double meanStd = fStd.GetParameter(1);
    double sigmaStd = std::fabs(fStd.GetParameter(2));
    if (fitStatusStd != 0) {
      meanStd = mean0;
      sigmaStd = rms0;
    }

    const double halfRange = std::max(4.0 * yBin, 1.2 * std::max(1e-12, sigmaStd));
    TF1 fCR(Form("fCR_%d", ib), "gaus", meanStd - halfRange, meanStd + halfRange);
    fCR.SetParameters(std::max(1.0, proj->GetMaximum()), meanStd, sigmaStd);
    const int fitStatusCR = int(proj->Fit(&fCR, "QSN0"));

    TF1 fMode(Form("fMode_%d", ib), "gaus", mode0 - 2.5 * rms0, mode0 + 2.5 * rms0);
    fMode.SetParameters(std::max(1.0, proj->GetMaximum()), mode0, rms0);
    const int fitStatusMode = int(proj->Fit(&fMode, "QSN0"));

    TF1 fDouble(Form("fDouble_%d", ib),
                "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
                meanStd - std::max(4.0 * yBin, 3.0 * sigmaStd),
                meanStd + std::max(4.0 * yBin, 3.0 * sigmaStd));
    fDouble.SetParameters(std::max(1.0, 0.7 * proj->GetMaximum()),
                          (fitStatusStd == 0 ? meanStd : mean0),
                          std::max(yBin, (fitStatusStd == 0 ? sigmaStd : rms0)),
                          std::max(1.0, 0.3 * proj->GetMaximum()),
                          std::max(2.0 * yBin, 2.0 * (fitStatusStd == 0 ? sigmaStd : rms0)));
    fDouble.SetParLimits(2, 0.5 * yBin, 20.0 * rms0);
    fDouble.SetParLimits(4, 0.5 * yBin, 40.0 * rms0);
    const int fitStatusDouble = int(proj->Fit(&fDouble, "QSN0"));

    c.Clear();
    proj->SetLineColor(kBlack);
    proj->SetMarkerStyle(20);
    proj->SetMarkerSize(0.8);
    proj->Draw("E");

    fStd.SetLineColor(kRed + 1);
    fStd.SetLineWidth(2);
    if (fitStatusStd == 0) fStd.Draw("SAME");

    fCR.SetLineColor(kBlue + 1);
    fCR.SetLineWidth(2);
    if (fitStatusCR == 0) fCR.Draw("SAME");

    fMode.SetLineColor(kGreen + 2);
    fMode.SetLineWidth(2);
    if (fitStatusMode == 0) fMode.Draw("SAME");

    fDouble.SetLineColor(kMagenta + 2);
    fDouble.SetLineWidth(2);
    if (fitStatusDouble == 0) fDouble.Draw("SAME");

    TLegend leg(0.14, 0.12, 0.84, 0.34);
    leg.SetBorderSize(0);
    leg.SetNColumns(2);
    leg.AddEntry(proj, Form("Data: N=%.0f", proj->GetEntries()), "lep");
    if (fitStatusStd == 0) {
      leg.AddEntry(&fStd, Form("Std Gaus: #mu=%.3f MeV, #sigma=%.3f MeV",
                               fStd.GetParameter(1), std::fabs(fStd.GetParameter(2))), "l");
    } else {
      leg.AddEntry((TObject *)nullptr, "Std Gaus: fit failed", "");
    }
    if (fitStatusCR == 0) {
      leg.AddEntry(&fCR, Form("Constrained Gaus: #mu=%.3f MeV, #sigma=%.3f MeV",
                              fCR.GetParameter(1), std::fabs(fCR.GetParameter(2))), "l");
    } else {
      leg.AddEntry((TObject *)nullptr, "Constrained Gaus: fit failed", "");
    }
    if (fitStatusMode == 0) {
      leg.AddEntry(&fMode, Form("Mode-seeded Gaus: #mu=%.3f MeV, #sigma=%.3f MeV",
                                fMode.GetParameter(1), std::fabs(fMode.GetParameter(2))), "l");
    } else {
      leg.AddEntry((TObject *)nullptr, "Mode-seeded Gaus: fit failed", "");
    }
    if (fitStatusDouble == 0) {
      leg.AddEntry(&fDouble, Form("Double Gaus: #mu=%.3f MeV", fDouble.GetParameter(1)), "l");
    } else {
      leg.AddEntry((TObject *)nullptr, "Double Gaus: fit failed", "");
    }
    leg.Draw();

    TLatex txt;
    txt.SetNDC(true);
    txt.SetTextSize(0.032);
    txt.DrawLatex(0.14, 0.93, Form("%s | p_{T} in [%.3f, %.3f] GeV",
                                   inFile,
                                   h2->GetXaxis()->GetBinLowEdge(ib),
                                   h2->GetXaxis()->GetBinUpEdge(ib)));

    c.SaveAs(outPdf);
    delete proj;
  }

  c.SaveAs((std::string(outPdf) + ")").c_str());
  fin.Close();

  std::cout << "Wrote per-pT-bin mass-difference fit pages to: " << outPdf << std::endl;
}

void run_full_decay_mass_resolution_study(const char *outFile = "decay_study.root",
                                          Long64_t nEvents = 500000,
                                          double motherMass = 5.279,
                                          double daughterMass1 = 0.13957,
                                          double daughterMass2 = 0.13957,
                                          double Cpt = 0.005,
                                          double Ctheta = 0.0008,
                                          double Cphi = 0.0008,
                                          int randomSeed = 12345,
                                          const char *plotPrefix = "decay_mass_offset") {
  run_decay_study(outFile, nEvents, motherMass, daughterMass1, daughterMass2,
                  Cpt, Ctheta, Cphi, randomSeed);
  analyze_decay_study(outFile, plotPrefix);
}

void run_requested_channels(Long64_t nEventsPerChannel = 500000,
                            double Cpt = 0.005,
                            double Ctheta = 0.0008,
                            double Cphi = 0.0008,
                            int randomSeed = 12345) {
  // Masses in GeV. 3HL/4HL masses can be updated to your preferred references.
  const std::vector<DecayChannel> channels = {
      {"Kshort_to_pip_pim", 0.497611, 0.13957039, 0.13957039},
      {"Lambda_to_p_pim", 1.115683, 0.93827208816, 0.13957039},
      {"3HL_to_He3_pim", 2.99109, 2.80940, 0.13957039},
      {"4HL_to_He4_pim", 3.92250, 3.72840, 0.13957039},
  };

  for (size_t i = 0; i < channels.size(); ++i) {
    const auto &ch = channels[i];
    const std::string outRoot = ch.name + ".root";
    const std::string plotPrefix = ch.name + "_offset";
    const int seed = randomSeed + static_cast<int>(i);

    std::cout << "\n=== Running channel: " << ch.name << " ===" << std::endl;
    run_full_decay_mass_resolution_study(outRoot.c_str(), nEventsPerChannel,
                                         ch.motherMass, ch.daughterMass1, ch.daughterMass2,
                                         Cpt, Ctheta, Cphi, seed, plotPrefix.c_str());
  }
}

void plot_mass_diff_fits_for_requested_channels() {
  const std::vector<std::string> channels = {
      "Kshort_to_pip_pim",
      "Lambda_to_p_pim",
      "3HL_to_He3_pim",
      "4HL_to_He4_pim",
  };

  for (const auto &ch : channels) {
    const std::string inRoot = ch + ".root";
    const std::string outPdf = ch + "_mass_diff_fits_per_ptbin.pdf";
    plot_mass_diff_fits_per_ptbin(inRoot.c_str(), outPdf.c_str());
  }
}
