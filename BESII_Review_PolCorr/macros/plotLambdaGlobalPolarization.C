#include "../../BESII_Review/FreezeOut/style.C"
#include "../../BESII_Review/FreezeOut/draw.C"

#include "TBox.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"

#include <cmath>
#include <vector>

struct DataSet {
  const char *name;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> stat;
  std::vector<double> sysLow;
  std::vector<double> sysHigh;
  int marker;
  double markerSize;
  int color;
  bool openMarker;
};

TGraphErrors *MakeStatGraph(const DataSet &d)
{
  const int n = d.x.size();
  std::vector<double> ex(n, 0.0);
  TGraphErrors *g = new TGraphErrors(n, d.x.data(), d.y.data(), ex.data(), d.stat.data());
  g->SetName(Form("gStat_%s", d.name));
  g->SetMarkerStyle(d.marker);
  g->SetMarkerSize(d.markerSize);
  g->SetMarkerColor(d.color);
  g->SetLineColor(d.color);
  g->SetLineWidth(2);
  if (d.openMarker) {
    g->SetMarkerColor(d.color);
    g->SetLineColor(d.color);
  }
  return g;
}

TGraphAsymmErrors *MakeSysGraph(const DataSet &d)
{
  const int n = d.x.size();
  std::vector<double> exLow(n, 0.0), exHigh(n, 0.0);
  TGraphAsymmErrors *g = new TGraphAsymmErrors(n, d.x.data(), d.y.data(), exLow.data(), exHigh.data(), d.sysLow.data(), d.sysHigh.data());
  g->SetName(Form("gSys_%s", d.name));
  return g;
}

TGraph *MakeCurve(const char *name, const std::vector<double> &x, const std::vector<double> &y, int color, int style, int width)
{
  TGraph *g = new TGraph(x.size(), x.data(), y.data());
  g->SetName(name);
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
  return g;
}

void DrawSysBoxes(TGraphAsymmErrors *g, double fracWidth, int color, double alpha)
{
  if (!g) return;
  for (int i = 0; i < g->GetN(); ++i) {
    double x, y;
    g->GetPoint(i, x, y);
    TBox *box = new TBox(x * (1.0 - fracWidth), y - g->GetErrorYlow(i), x * (1.0 + fracWidth), y + g->GetErrorYhigh(i));
    box->SetFillColorAlpha(color, alpha);
    box->SetLineColorAlpha(color, alpha);
    box->SetLineWidth(1);
    box->Draw("same");
  }
}

void plotLambdaGlobalPolarization()
{
  style();

  const Double_t EMIN = 1.9;
  const Double_t EMAX = 7000.0;
  const Double_t PMIN = -1.0;
  const Double_t PMAX = 10.0;

  // Published data points. Errors are in percent; systematic boxes are separate from stat bars.
  DataSet hades10 = {
    "HADES_10_40",
    {2.42, 2.55},
    {5.3, 4.4},
    {1.0, 0.3},
    {1.3, 0.4},
    {1.3, 0.4},
    25,
    1.6,
    kRed + 1,
    true
  };
  DataSet hades20 = {
    "HADES_20_40",
    {2.42, 2.55},
    {6.8, 6.2},
    {1.3, 0.4},
    {2.1, 0.6},
    {2.1, 0.6},
    21,
    1.45,
    kRed + 1,
    false
  };
  DataSet star = {
    "STAR",
    {3.0, 7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {4.90823, 1.80, 1.19, 1.17, 0.9152, 0.7168, 0.45, 1.20, 0.243},
    {0.81389, 0.60, 0.35, 0.40, 0.0503, 0.0553, 0.40, 1.00, 0.035},
    {0.15485, 0.18, 0.18, 0.27, 0.0172, 0.0143, 0.18, 0.0, 0.043},
    {0.15485, 0.18, 0.18, 0.27, 0.0172, 0.0143, 0.18, 0.0, 0.034},
    29,
    1.9,
    kBlue + 1,
    false
  };
  DataSet alice = {
    "ALICE",
    {2760.0, 5020.0},
    {-0.07, 0.12},
    {0.09, 0.10},
    {0.04, 0.04},
    {0.04, 0.04},
    33,
    1.8,
    kGreen + 2,
    false
  };

  TGraphErrors *gHades10 = MakeStatGraph(hades10);
  TGraphErrors *gHades20 = MakeStatGraph(hades20);
  TGraphErrors *gStar = MakeStatGraph(star);
  TGraphErrors *gAlice = MakeStatGraph(alice);
  TGraphAsymmErrors *sHades10 = MakeSysGraph(hades10);
  TGraphAsymmErrors *sHades20 = MakeSysGraph(hades20);
  TGraphAsymmErrors *sStar = MakeSysGraph(star);
  TGraphAsymmErrors *sAlice = MakeSysGraph(alice);

  TGraph *g3fdHadronic = MakeCurve(
    "g3fdHadronic",
    {2.376, 2.626, 2.876, 3.376, 3.876, 4.876, 5.876, 6.876, 8.876, 11.376},
    {2.80, 3.10, 3.20, 3.30, 3.10, 2.60, 1.80, 1.10, 0.50, 0.15},
    kGreen + 2, 1, 3);
  TGraph *g3fdCross = MakeCurve(
    "g3fdCross",
    {2.426, 2.626, 2.776, 2.976, 3.376, 3.876, 4.876, 5.876, 7.376, 8.876, 10.876},
    {3.70, 4.60, 5.10, 4.70, 3.60, 3.00, 2.10, 1.70, 1.25, 1.00, 0.80},
    kBlue - 6, 1, 3);
  TGraph *g3fd1pt = MakeCurve(
    "g3fd1pt",
    {2.426, 2.676, 3.076, 3.676, 4.376, 5.376, 6.876, 8.876, 11.876, 17.876, 26.876},
    {3.10, 3.25, 3.45, 3.60, 3.65, 3.20, 2.60, 1.60, 1.05, 0.75, 0.65},
    kOrange + 7, 1, 3);
  TGraph *gAmptLow = MakeCurve(
    "gAmptLow",
    {2.376, 2.676, 3.076, 3.676, 4.876, 6.876, 9.876, 11.876, 15.876, 21.876},
    {0.65, 0.70, 0.75, 0.82, 0.90, 1.02, 1.05, 1.00, 0.85, 0.72},
    kOrange + 1, 1, 4);
  TGraph *gUrqmd = MakeCurve(
    "gUrqmd",
    {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {1.45, 1.12, 0.95, 0.72, 0.55, 0.41, 0.30, 0.17},
    kBlue + 1, 1, 2);
  TGraph *gUrqmdFd = MakeCurve(
    "gUrqmdFd",
    {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {1.25, 1.02, 0.86, 0.63, 0.48, 0.36, 0.26, 0.14},
    kBlue + 1, 2, 2);
  TGraph *gAmpt = MakeCurve(
    "gAmpt",
    {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {3.10, 2.25, 1.75, 1.20, 0.95, 0.72, 0.55, 0.25},
    kMagenta + 1, 1, 2);
  TGraph *gAmptFd = MakeCurve(
    "gAmptFd",
    {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {2.55, 1.95, 1.50, 1.05, 0.82, 0.62, 0.47, 0.22},
    kMagenta + 1, 2, 2);
  TGraph *gChiral = MakeCurve(
    "gChiral",
    {7.7, 10.0, 14.5, 19.6, 27.0, 39.0, 62.4, 100.0, 150.0, 200.0},
    {2.45, 2.35, 2.20, 2.05, 1.85, 1.55, 1.25, 1.00, 0.82, 0.70},
    kAzure + 7, 3, 2);

  TCanvas *c1 = new TCanvas("cLambdaGlobalPolarization", "", 900, 700);
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.13);
  c1->SetTopMargin(0.04);
  c1->SetRightMargin(0.03);
  c1->SetLogx(1);
  c1->Draw();

  TH1D *frame = new TH1D("frame", "", 1, EMIN, EMAX);
  frame->SetMinimum(PMIN);
  frame->SetMaximum(PMAX);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitle("P_{#Lambda} (%)");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->Draw();

  drawLine(EMIN, 0.0, EMAX, 0.0, 2, 2, kGray + 2);

  g3fdHadronic->Draw("C SAME");
  g3fdCross->Draw("C SAME");
  g3fd1pt->Draw("C SAME");
  gAmptLow->Draw("C SAME");
  gUrqmd->Draw("C SAME");
  gUrqmdFd->Draw("C SAME");
  gAmpt->Draw("C SAME");
  gAmptFd->Draw("C SAME");
  gChiral->Draw("C SAME");

  DrawSysBoxes(sHades10, 0.045, kRed + 1, 0.16);
  DrawSysBoxes(sHades20, 0.045, kRed + 1, 0.22);
  DrawSysBoxes(sStar, 0.045, kBlue + 1, 0.16);
  DrawSysBoxes(sAlice, 0.045, kGreen + 2, 0.18);

  gHades10->Draw("P SAME");
  gHades20->Draw("P SAME");
  gStar->Draw("P SAME");
  gAlice->Draw("P SAME");

  drawText(2.1, 9.35, "Global #Lambda Polarization", 42, 0.045);
  drawText(2.1, 8.78, "stat. bars, sys. boxes; theory curves digitized from publication figures", 42, 0.026, 0, kGray + 2);

  TLegend *legData = new TLegend(0.15, 0.65, 0.47, 0.88);
  legData->SetFillStyle(4000);
  legData->SetBorderSize(0);
  legData->SetTextSize(0.028);
  legData->AddEntry(gHades20, "HADES 20-40%", "p");
  legData->AddEntry(gHades10, "HADES 10-40%", "p");
  legData->AddEntry(gStar, "STAR 20-50%", "p");
  legData->AddEntry(gAlice, "ALICE 15-50%", "p");
  legData->Draw();

  TLegend *legTheory = new TLegend(0.55, 0.52, 0.95, 0.88);
  legTheory->SetFillStyle(4000);
  legTheory->SetBorderSize(0);
  legTheory->SetTextSize(0.023);
  legTheory->AddEntry(g3fdHadronic, "3FD hadronic EoS", "l");
  legTheory->AddEntry(g3fdCross, "3FD crossover EoS", "l");
  legTheory->AddEntry(g3fd1pt, "3FD 1PT EoS", "l");
  legTheory->AddEntry(gAmptLow, "AMPT low-energy", "l");
  legTheory->AddEntry(gUrqmd, "UrQMD+vHLLE primary", "l");
  legTheory->AddEntry(gUrqmdFd, "UrQMD+vHLLE feed-down", "l");
  legTheory->AddEntry(gAmpt, "AMPT primary", "l");
  legTheory->AddEntry(gAmptFd, "AMPT feed-down", "l");
  legTheory->AddEntry(gChiral, "Chiral kinetic", "l");
  legTheory->Draw();

  drawHistBox(EMIN, EMAX, PMIN, PMAX);

  gSystem->mkdir("figures", kTRUE);
  c1->Update();
  c1->SaveAs("figures/lambda_global_polarization_vs_energy_ROOT.pdf");
  c1->SaveAs("figures/lambda_global_polarization_vs_energy_ROOT.png");
}
