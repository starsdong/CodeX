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
#include "TLatex.h"
#include "TPad.h"
#include "TSystem.h"

#include <cmath>
#include <vector>

// Hand-editable ROOT macro for the global Lambda and anti-Lambda polarization plot.
// Units:
//   x = sqrt(s_NN) in GeV
//   y/stat/sys = P_Lambda in percent
//
// Run from BESII_Review_PolCorr:
//   root -l -b -q macros/plotLambdaGlobalPolarization_withAntiLambda.C

namespace LambdaPolWithAnti {

// -------------------------
// Global plot knobs
// -------------------------
const bool kLogX = true;
const bool kCompressMiddleX = false; // one-panel mode: ALICE is drawn at sqrt(sNN)/10
const double kXMin = 1.9;
const double kXMax = 700.0;
const double kXBreakLow = 300.0;
const double kXBreakHigh = 2000.0;
const double kYMin = -0.3;
const double kYMax = 10.0;
const char *kOutputBase = "figures/lambda_global_polarization_with_antilambda_scaled_alice";
const double kAliceXScale = 0.10;
const double kScaledAliceBreakX = 230.0; // visual marker between 200 and the scaled ALICE region
const double kSinglePanelXLabelSize = 0.052;
const double kSinglePanelXLabelY = 0.138;

const double kLeftMargin = 0.14;
const double kRightMargin = 0.03;
const double kTopMargin = 0.04;
const double kBottomMargin = 0.16;
const double kCompressedInnerMargin = 0.006;
const double kMiddleGuideHalfGap = 0.008;
const double kRightPadXLabelSize = 0.16;
const double kRightPadXLabelY = 0.098;

const double kSysBoxFracWidth = 0.045; // box half-width as fraction of x for log axis
const double kSysBoxAlpha = 0.18;
const double kMarkerScale = 1.30;
const double kStarAliceMarkerScale = 1.48;
const int kThinTheoryLineWidth = 1;

// Toggle groups here.
const bool kDrawHades10to40 = true;
const bool kDrawHades20to40 = false;
const bool kDrawStar = true;
const bool kDrawStarAntiLambda = true;
const bool kDrawAlice = true;
const bool kDrawAliceAntiLambda = true;
const bool kDrawTheory = true;

// -------------------------
// Data point blocks
// -------------------------
struct DataSet {
  const char *legend;
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

DataSet Hades20to40()
{
  return {
    "HADES #Lambda",
    {2.42, 2.55},
    {6.8, 6.2},
    {1.3, 0.4},
    {2.1, 0.6},
    {2.1, 0.6},
    21, 1.45 * kMarkerScale, kBlue + 1, false
  };
}

DataSet Hades10to40()
{
  return {
    "HADES #Lambda",
    {2.42, 2.55},
    {5.3, 4.4},
    {1.0, 0.3},
    {1.3, 0.4},
    {1.3, 0.4},
    25, 1.60 * kMarkerScale, kBlue + 1, true
  };
}

DataSet Star20to50()
{
  return {
    "STAR #Lambda",
    {3.0, 7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {4.90823, 1.80, 1.19, 1.17, 0.9152, 0.7168, 0.45, 1.20, 0.243},
    {0.81389, 0.60, 0.35, 0.40, 0.0503, 0.0553, 0.40, 1.00, 0.035},
    {0.15485, 0.18, 0.18, 0.27, 0.0172, 0.0143, 0.18, 0.0, 0.043},
    {0.15485, 0.18, 0.18, 0.27, 0.0172, 0.0143, 0.18, 0.0, 0.034},
    29, 1.90 * kStarAliceMarkerScale, kRed + 1, false
  };
}

DataSet StarAntiLambda20to50()
{
  return {
    "STAR #bar{#Lambda}",
    {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
    {7.70, 1.59, 2.01, 0.8976, 0.8257, 0.83, 1.50, 0.210},
    {3.20, 1.10, 1.10, 0.1170, 0.1046, 0.50, 1.40, 0.039},
    {0.90, 0.13, 0.13, 0.0171, 0.0163, 0.13, 0.0, 0.035},
    {0.0, 0.0, 0.35, 0.0171, 0.0163, 0.0, 0.0, 0.053},
    30, 2.00 * kStarAliceMarkerScale, kRed + 1, true
  };
}

DataSet Alice15to50()
{
  return {
    "ALICE #Lambda",
    {2760.0, 5020.0},
    {-0.07, 0.12},
    {0.09, 0.10},
    {0.04, 0.04},
    {0.04, 0.04},
    33, 1.80 * kStarAliceMarkerScale, kGreen + 2, false
  };
}

DataSet AliceAntiLambda15to50()
{
  return {
    "ALICE #bar{#Lambda}",
    {2760.0, 5020.0},
    {0.05, -0.13},
    {0.09, 0.11},
    {0.03, 0.03},
    {0.03, 0.03},
    27, 1.95 * kStarAliceMarkerScale, kGreen + 2, true
  };
}

DataSet ScaleDataX(DataSet d, double scale)
{
  for (double &x : d.x) x *= scale;
  return d;
}

DataSet Alice15to50Plot()
{
  return ScaleDataX(Alice15to50(), kAliceXScale);
}

DataSet AliceAntiLambda15to50Plot()
{
  return ScaleDataX(AliceAntiLambda15to50(), kAliceXScale);
}

// -------------------------
// Theory curve blocks
// -------------------------
struct CurveSet {
  const char *legend;
  std::vector<double> x;
  std::vector<double> y;
  int color;
  int lineStyle;
  int lineWidth;
};

std::vector<CurveSet> TheoryCurves()
{
  return {
    {"3FD hadronic EoS",
     {2.376, 2.626, 2.876, 3.376, 3.876, 4.876, 5.876, 6.876, 8.876, 11.376},
     {2.80, 3.10, 3.20, 3.30, 3.10, 2.60, 1.80, 1.10, 0.50, 0.15},
     kBlack, 1, 3},
    {"3FD crossover EoS",
     {2.426, 2.626, 2.776, 2.976, 3.376, 3.876, 4.876, 5.876, 7.376, 8.876, 10.876},
     {3.70, 4.60, 5.10, 4.70, 3.60, 3.00, 2.10, 1.70, 1.25, 1.00, 0.80},
     kBlack, 2, 3},
    {"3FD 1PT EoS",
     {2.426, 2.676, 3.076, 3.676, 4.376, 5.376, 6.876, 8.876, 11.876, 17.876, 26.876},
     {3.10, 3.25, 3.45, 3.60, 3.65, 3.20, 2.60, 1.60, 1.05, 0.75, 0.65},
     kBlack, 7, 3},
    {"AMPT low-energy",
     {2.376, 2.676, 3.076, 3.676, 4.876, 6.876, 9.876, 11.876, 15.876, 21.876},
     {0.65, 0.70, 0.75, 0.82, 0.90, 1.02, 1.05, 1.00, 0.85, 0.72},
     kOrange + 1, 1, 4},
    {"AMPT",
     {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
     {2.55, 1.95, 1.50, 1.05, 0.82, 0.62, 0.47, 0.22},
     kMagenta + 1, 2, 2},
    {"UrQMD+vHLLE",
     {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0},
     {1.25, 1.02, 0.86, 0.63, 0.48, 0.36, 0.26, 0.14},
     kBlue + 1, 2, 2},
    {"Chiral kinetic",
     {7.7, 10.0, 14.5, 19.6, 27.0, 39.0, 62.4, 100.0, 150.0, 200.0},
     {2.45, 2.35, 2.20, 2.05, 1.85, 1.55, 1.25, 1.00, 0.82, 0.70},
     kAzure + 7, 3, 2}
  };
}

TGraphErrors *MakeStatGraph(const DataSet &d)
{
  const int n = d.x.size();
  std::vector<double> ex(n, 0.0);
  TGraphErrors *g = new TGraphErrors(n, d.x.data(), d.y.data(), ex.data(), d.stat.data());
  g->SetMarkerStyle(d.marker);
  g->SetMarkerSize(d.markerSize);
  g->SetMarkerColor(d.color);
  g->SetLineColor(d.color);
  g->SetLineWidth(2);
  return g;
}

TGraphAsymmErrors *MakeSysGraph(const DataSet &d)
{
  const int n = d.x.size();
  std::vector<double> exLow(n, 0.0), exHigh(n, 0.0);
  return new TGraphAsymmErrors(n, d.x.data(), d.y.data(), exLow.data(), exHigh.data(), d.sysLow.data(), d.sysHigh.data());
}

TGraph *MakeCurve(const CurveSet &c)
{
  TGraph *g = new TGraph(c.x.size(), c.x.data(), c.y.data());
  g->SetLineColor(c.color);
  g->SetLineStyle(c.lineStyle);
  g->SetLineWidth(kThinTheoryLineWidth);
  return g;
}

DataSet SelectDataRange(const DataSet &d, double xmin, double xmax)
{
  DataSet out = {d.legend, {}, {}, {}, {}, {}, d.marker, d.markerSize, d.color, d.openMarker};
  for (size_t i = 0; i < d.x.size(); ++i) {
    if (d.x[i] < xmin || d.x[i] > xmax) continue;
    out.x.push_back(d.x[i]);
    out.y.push_back(d.y[i]);
    out.stat.push_back(d.stat[i]);
    out.sysLow.push_back(d.sysLow[i]);
    out.sysHigh.push_back(d.sysHigh[i]);
  }
  return out;
}

CurveSet SelectCurveRange(const CurveSet &c, double xmin, double xmax)
{
  CurveSet out = {c.legend, {}, {}, c.color, c.lineStyle, c.lineWidth};
  for (size_t i = 0; i < c.x.size(); ++i) {
    if (c.x[i] < xmin || c.x[i] > xmax) continue;
    out.x.push_back(c.x[i]);
    out.y.push_back(c.y[i]);
  }
  return out;
}

void DrawSysBoxes(TGraphAsymmErrors *g, int color, double alpha = kSysBoxAlpha)
{
  if (!g) return;
  for (int i = 0; i < g->GetN(); ++i) {
    double x, y;
    g->GetPoint(i, x, y);
    const double eyLow = g->GetErrorYlow(i);
    const double eyHigh = g->GetErrorYhigh(i);
    if (eyLow <= 0.0 && eyHigh <= 0.0) continue;
    TBox *box = new TBox(x * (1.0 - kSysBoxFracWidth), y - eyLow,
                         x * (1.0 + kSysBoxFracWidth), y + eyHigh);
    box->SetFillColorAlpha(color, alpha);
    box->SetLineColorAlpha(color, 0.0);
    box->SetLineWidth(0);
    box->Draw("same");
  }
}

void DrawData(const DataSet &d, TLegend *leg)
{
  if (d.x.empty()) return;
  TGraphAsymmErrors *sys = MakeSysGraph(d);
  TGraphErrors *stat = MakeStatGraph(d);
  DrawSysBoxes(sys, d.color, d.openMarker ? 0.13 : kSysBoxAlpha);
  stat->Draw("P SAME");
  if (leg) leg->AddEntry(stat, d.legend, "p");
}

TH1D *DrawFrame(const char *name, double xmin, double xmax, bool showYAxis)
{
  TH1D *frame = new TH1D(name, "", 1, xmin, xmax);
  frame->SetMinimum(kYMin);
  frame->SetMaximum(kYMax);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  frame->GetXaxis()->SetTitleOffset(1.12);
  frame->GetYaxis()->SetTitle("P_{#Lambda} (%)");
  frame->GetYaxis()->SetTitleOffset(1.15);
  if (!showYAxis) {
    frame->GetYaxis()->SetLabelSize(0.0);
    frame->GetYaxis()->SetTitleSize(0.0);
  }
  frame->Draw();
  return frame;
}

void DrawPlotObjects(TLegend *legData, TLegend *legTheory, double xmin = kXMin, double xmax = kXMax)
{
  if (kDrawTheory) {
    for (const auto &curve : TheoryCurves()) {
      CurveSet selected = SelectCurveRange(curve, xmin, xmax);
      if (selected.x.size() < 2) continue;
      TGraph *g = MakeCurve(selected);
      g->Draw("C SAME");
      if (legTheory) legTheory->AddEntry(g, curve.legend, "l");
    }
  }

  if (kDrawHades20to40) DrawData(SelectDataRange(Hades20to40(), xmin, xmax), legData);
  if (kDrawHades10to40) DrawData(SelectDataRange(Hades10to40(), xmin, xmax), legData);
  if (kDrawStar) DrawData(SelectDataRange(Star20to50(), xmin, xmax), legData);
  if (kDrawStarAntiLambda) DrawData(SelectDataRange(StarAntiLambda20to50(), xmin, xmax), legData);
  if (kDrawAlice) DrawData(SelectDataRange(Alice15to50Plot(), xmin, xmax), legData);
  if (kDrawAliceAntiLambda) DrawData(SelectDataRange(AliceAntiLambda15to50Plot(), xmin, xmax), legData);
}

void DrawXAxisBreak(bool leftPad)
{
  const double y = gPad->GetBottomMargin();
  const double dy = 0.022;
  const double dx = 0.014;
  const double x = leftPad ? 0.992 : 0.008;

  TLine *l1 = new TLine(x - dx, y - dy * 0.2, x + dx, y + dy);
  TLine *l2 = new TLine(x - dx, y + dy * 0.7, x + dx, y + dy * 1.9);
  l1->SetNDC();
  l2->SetNDC();
  l1->SetLineWidth(1);
  l2->SetLineWidth(1);
  l1->SetLineColor(kGray + 1);
  l2->SetLineColor(kGray + 1);
  l1->Draw("same");
  l2->Draw("same");
}

void DrawRightPadXLabel(double xValue, const char *label)
{
  const double left = gPad->GetLeftMargin();
  const double right = gPad->GetRightMargin();
  const double usable = 1.0 - left - right;
  const double u = (std::log10(xValue) - std::log10(kXBreakHigh)) / (std::log10(kXMax) - std::log10(kXBreakHigh));
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(23);
  text.SetTextFont(42);
  text.SetTextSize(kRightPadXLabelSize);
  text.DrawLatex(left + usable * u, kRightPadXLabelY, label);
}

double LogXToPadNDC(double x, double xmin, double xmax)
{
  const double left = gPad->GetLeftMargin();
  const double right = gPad->GetRightMargin();
  const double usable = 1.0 - left - right;
  const double u = (std::log10(x) - std::log10(xmin)) / (std::log10(xmax) - std::log10(xmin));
  return left + usable * u;
}

void DrawSinglePanelXLabel(double xValue, const char *label, int align = 23, double textSize = kSinglePanelXLabelSize)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(align);
  text.SetTextFont(42);
  text.SetTextSize(textSize);
  text.DrawLatex(LogXToPadNDC(xValue, kXMin, kXMax), kSinglePanelXLabelY, label);
}

void DrawSinglePanelXLabels()
{
  DrawSinglePanelXLabel(2.0, "2");
  DrawSinglePanelXLabel(4.0, "4");
  DrawSinglePanelXLabel(10.0, "10");
  DrawSinglePanelXLabel(20.0, "20");
  DrawSinglePanelXLabel(40.0, "40");
  DrawSinglePanelXLabel(100.0, "100");
  DrawSinglePanelXLabel(500.0, "5000");
}

void DrawScaledAliceXAxisBreak()
{
  const double x = LogXToPadNDC(kScaledAliceBreakX, kXMin, kXMax);
  const double y = gPad->GetBottomMargin();
  const double offsets[] = {-0.006, 0.006};
  const double dx = 0.006;
  const double dy = 0.014;

  for (const double offset : offsets) {
    TLine *slash = new TLine(x + offset - dx, y - dy, x + offset + dx, y + dy);
    slash->SetNDC();
    slash->SetLineColor(kBlack);
    slash->SetLineWidth(3);
    slash->Draw("same");
  }
}

void DrawPartialHistBox(double x1, double x2, double y1, double y2, bool drawLeft, bool drawRight)
{
  TLine *bottom = new TLine(x1, y1, x2, y1);
  TLine *top = new TLine(x1, y2, x2, y2);
  bottom->SetLineWidth(3);
  top->SetLineWidth(3);
  bottom->Draw("same");
  top->Draw("same");

  if (drawLeft) {
    TLine *left = new TLine(x1, y1, x1, y2);
    left->SetLineWidth(3);
    left->Draw("same");
  }

  if (drawRight) {
    TLine *right = new TLine(x2, y1, x2, y2);
    right->SetLineWidth(3);
    right->Draw("same");
  }
}

void HideLeftPanelRightBoundary()
{
  const double x = 1.0 - gPad->GetRightMargin();
  TLine cover(x, 0.0, x, 1.0);
  cover.SetNDC();
  cover.SetLineColor(kWhite);
  cover.SetLineWidth(9);
  cover.Draw("same");
}

void HideLeftPanelRightTicks()
{
  const double x = 1.0 - gPad->GetRightMargin();
  TLine cover(x - 0.002, gPad->GetBottomMargin(), x + 0.002, 1.0 - gPad->GetTopMargin());
  cover.SetNDC();
  cover.SetLineColor(kWhite);
  cover.SetLineWidth(12);
  cover.Draw("same");
}

void HideRightPanelLeftBoundary()
{
  const double left = gPad->GetLeftMargin();
  TLine cover(left, gPad->GetBottomMargin(), left, 1.0 - gPad->GetTopMargin());
  cover.SetNDC();
  cover.SetLineColor(kWhite);
  cover.SetLineWidth(18);
  cover.Draw("same");
}

void HideCanvasMiddleBoundary(double split)
{
  const double xEdges[] = {
    split * (1.0 - kCompressedInnerMargin),
    split,
    split + (1.0 - split) * kCompressedInnerMargin
  };

  for (const double x : xEdges) {
    TLine *cover = new TLine(x, kBottomMargin, x, 1.0 - kTopMargin);
    cover->SetNDC();
    cover->SetLineColor(kWhite);
    cover->SetLineWidth(8);
    cover->Draw("same");
  }

  const double xGuideLines[] = {
    split - kMiddleGuideHalfGap,
    split + kMiddleGuideHalfGap
  };

  for (const double x : xGuideLines) {
    TLine *guide = new TLine(x, kBottomMargin, x, 1.0 - kTopMargin);
    guide->SetNDC();
    guide->SetLineColorAlpha(kGray + 1, 0.50);
    guide->SetLineWidth(1);
    guide->Draw("same");
  }
}

} // namespace LambdaPolWithAnti

void plotLambdaGlobalPolarization_withAntiLambda()
{
  using namespace LambdaPolWithAnti;

  style();

  TCanvas *c = new TCanvas("cLambdaGlobalPolarizationWithAntiLambda", "", 900, 700);
  c->SetLeftMargin(kLeftMargin);
  c->SetRightMargin(kRightMargin);
  c->SetTopMargin(kTopMargin);
  c->SetBottomMargin(kBottomMargin);
  if (kLogX) c->SetLogx(1);

  if (kCompressMiddleX) {
    c->Clear();
    const double split = 0.79;

    TPad *pLeft = new TPad("pLeftCompressed", "", 0.0, 0.0, split, 1.0);
    TPad *pRight = new TPad("pRightCompressed", "", split, 0.0, 1.0, 1.0);
    pLeft->SetBorderMode(0);
    pLeft->SetBorderSize(0);
    pLeft->SetFrameBorderMode(0);
    pRight->SetBorderMode(0);
    pRight->SetBorderSize(0);
    pRight->SetFrameBorderMode(0);
    pLeft->SetLeftMargin(kLeftMargin);
    pLeft->SetRightMargin(kCompressedInnerMargin);
    pLeft->SetTopMargin(kTopMargin);
    pLeft->SetBottomMargin(kBottomMargin);
    pRight->SetLeftMargin(kCompressedInnerMargin);
    pRight->SetRightMargin(kRightMargin);
    pRight->SetTopMargin(kTopMargin);
    pRight->SetBottomMargin(kBottomMargin);
    pLeft->SetTickx(1);
    pLeft->SetTicky(0);
    pRight->SetTickx(1);
    pRight->SetTicky(0);
    if (kLogX) {
      pLeft->SetLogx(1);
      pRight->SetLogx(1);
    }
    pLeft->Draw();
    pRight->Draw();

    pLeft->cd();
    TH1D *frameLeft = DrawFrame("frameLambdaGlobalPolarizationLeft", kXMin, kXBreakLow, true);
    frameLeft->GetXaxis()->SetTitle("");
    gPad->Update();
    gPad->GetFrame()->SetLineColor(kWhite);
    drawLine(kXMin, 0.0, kXBreakLow, 0.0, 2, 2, kGray + 2);

    TLegend *legTheory = new TLegend(0.55, 0.53, 0.93, 0.89);
    legTheory->SetFillStyle(4000);
    legTheory->SetBorderSize(0);
    legTheory->SetTextSize(0.022);

    TLegend *legData = new TLegend(0.16, 0.53, 0.44, 0.89);
    legData->SetFillStyle(4000);
    legData->SetBorderSize(0);
    legData->SetTextSize(0.028);

    DrawPlotObjects(nullptr, legTheory, kXMin, kXBreakLow);

//    drawText(2.1, 7.35, "Global #Lambda Polarization", 42, 0.045);
    if (kDrawTheory) legTheory->Draw();
    DrawPartialHistBox(kXMin, kXBreakLow, kYMin, kYMax, true, false);
    HideLeftPanelRightTicks();
    HideLeftPanelRightBoundary();
    DrawXAxisBreak(true);

    pRight->cd();
    TH1D *frameRight = DrawFrame("frameLambdaGlobalPolarizationRight", kXBreakHigh, kXMax, false);
    frameRight->GetXaxis()->SetTitle("");
    frameRight->GetXaxis()->SetLabelSize(0.0);
    frameRight->GetYaxis()->SetTickLength(0.0);
    gPad->Update();
    gPad->GetFrame()->SetLineColor(kWhite);
    drawLine(kXBreakHigh, 0.0, kXMax, 0.0, 2, 2, kGray + 2);
    DrawPlotObjects(nullptr, nullptr, kXBreakHigh, kXMax);
    DrawPartialHistBox(kXBreakHigh, kXMax, kYMin, kYMax, false, true);
    HideRightPanelLeftBoundary();
    DrawXAxisBreak(false);
    DrawRightPadXLabel(3000.0, "3000");
    DrawRightPadXLabel(5000.0, "5000");

    TLegend *legDataRight = new TLegend(0.04, 0.55, 0.97, 0.82);
    legDataRight->SetFillStyle(4000);
    legDataRight->SetBorderSize(0);
    legDataRight->SetTextSize(0.09);
    if (kDrawHades10to40) legDataRight->AddEntry(MakeStatGraph(Hades10to40()), Hades10to40().legend, "p");
    if (kDrawStar) legDataRight->AddEntry(MakeStatGraph(Star20to50()), Star20to50().legend, "p");
    if (kDrawStarAntiLambda) legDataRight->AddEntry(MakeStatGraph(StarAntiLambda20to50()), StarAntiLambda20to50().legend, "p");
    if (kDrawAlice) legDataRight->AddEntry(MakeStatGraph(Alice15to50Plot()), Alice15to50().legend, "p");
    if (kDrawAliceAntiLambda) legDataRight->AddEntry(MakeStatGraph(AliceAntiLambda15to50Plot()), AliceAntiLambda15to50().legend, "p");
    legDataRight->Draw();

    c->cd();
    HideCanvasMiddleBoundary(split);
    TLatex xTitle;
    xTitle.SetNDC();
    xTitle.SetTextFont(42);
    xTitle.SetTextSize(0.055);
    xTitle.DrawLatex(0.42, 0.015, "#sqrt{s_{NN}} (GeV)");
    gSystem->mkdir("figures", kTRUE);
    c->Update();
    c->SaveAs(Form("%s.pdf", kOutputBase));
    c->SaveAs(Form("%s.png", kOutputBase));
    return;
  }

  TH1D *frame = new TH1D("frameLambdaGlobalPolarizationManual", "", 1, kXMin, kXMax);
  frame->SetMinimum(kYMin);
  frame->SetMaximum(kYMax);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitle("Global P_{#Lambda,#bar{#Lambda}} (%)");
  frame->GetYaxis()->SetTitleOffset(0.90);
  frame->Draw();

  drawLine(kXMin, 0.0, kXMax, 0.0, 2, 2, kGray + 2);

  TLegend *legTheory = new TLegend(0.68, 0.50, 0.99, 0.90);
  legTheory->SetFillStyle(4000);
  legTheory->SetBorderSize(0);
  legTheory->SetTextSize(0.030);

  TLegend *legData = new TLegend(0.43, 0.58, 0.68, 0.90);
  legData->SetFillStyle(4000);
  legData->SetBorderSize(0);
  legData->SetTextSize(0.036);

  DrawPlotObjects(legData, legTheory, kXMin, kXMax);

//  drawText(2.1, 7.35, "Global #Lambda Polarization", 42, 0.045);

  legData->Draw();
  if (kDrawTheory) legTheory->Draw();
  drawHistBox(kXMin, kXMax, kYMin, kYMax);
  DrawScaledAliceXAxisBreak();
  DrawSinglePanelXLabels();

  gSystem->mkdir("figures", kTRUE);
  c->Update();
  c->SaveAs(Form("%s.pdf", kOutputBase));
  c->SaveAs(Form("%s.png", kOutputBase));
}
