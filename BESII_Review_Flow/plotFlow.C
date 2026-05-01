#include "../BESII_Review/FreezeOut/draw.C"
#include "../BESII_Review/FreezeOut/style.C"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct FlowPoint {
  std::string experiment;
  std::string particle;
  double s = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

std::vector<std::string> splitCSV(const std::string &line)
{
  std::vector<std::string> out;
  std::string field;
  bool quoted = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (c == ',' && !quoted) {
      out.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  out.push_back(field);
  return out;
}

double valueOrZero(const std::vector<std::string> &fields, const std::map<std::string, int> &index, const std::string &name)
{
  const auto it = index.find(name);
  if (it == index.end() || it->second >= (int)fields.size() || fields[it->second].empty()) return 0.0;
  return atof(fields[it->second].c_str());
}

std::string stringValue(const std::vector<std::string> &fields, const std::map<std::string, int> &index, const std::string &name)
{
  const auto it = index.find(name);
  if (it == index.end() || it->second >= (int)fields.size()) return "";
  return fields[it->second];
}

std::vector<FlowPoint> readFlowCSV(const char *path)
{
  std::vector<FlowPoint> points;
  std::ifstream in(path);
  if (!in.good()) {
    printf("Cannot open %s\n", path);
    return points;
  }

  std::string line;
  if (!std::getline(in, line)) return points;
  std::vector<std::string> header = splitCSV(line);
  std::map<std::string, int> index;
  for (int i = 0; i < (int)header.size(); ++i) index[header[i]] = i;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<std::string> fields = splitCSV(line);
    const double statP = valueOrZero(fields, index, "stat_err_plus");
    const double statM = valueOrZero(fields, index, "stat_err_minus");
    const double systP = valueOrZero(fields, index, "syst_err_plus");
    const double systM = valueOrZero(fields, index, "syst_err_minus");
    const double errP = std::sqrt(statP * statP + systP * systP);
    const double errM = std::sqrt(statM * statM + systM * systM);

    FlowPoint p;
    p.experiment = stringValue(fields, index, "experiment");
    p.particle = stringValue(fields, index, "particle");
    p.s = valueOrZero(fields, index, "sqrt_s_NN_GeV");
    p.y = valueOrZero(fields, index, "value");
    p.ey = 0.5 * (errP + errM);
    if (p.s > 0.0) points.push_back(p);
  }
  return points;
}

void appendPoints(std::vector<FlowPoint> &out, const std::vector<FlowPoint> &in)
{
  out.insert(out.end(), in.begin(), in.end());
}

bool isStarLike(const FlowPoint &p)
{
  return p.experiment == "STAR" || p.experiment == "STAR FXT";
}

std::vector<FlowPoint> selectParticle(const std::vector<FlowPoint> &all, const char *particle, bool starLike)
{
  std::vector<FlowPoint> selected;
  for (const auto &p : all) {
    if (p.particle == particle && isStarLike(p) == starLike) selected.push_back(p);
  }
  std::sort(selected.begin(), selected.end(), [](const FlowPoint &a, const FlowPoint &b) { return a.s < b.s; });
  return selected;
}

std::vector<FlowPoint> selectExperiment(const std::vector<FlowPoint> &all, const char *experiment)
{
  std::vector<FlowPoint> selected;
  for (const auto &p : all) {
    if (p.experiment == experiment) selected.push_back(p);
  }
  std::sort(selected.begin(), selected.end(), [](const FlowPoint &a, const FlowPoint &b) {
    if (std::fabs(a.s - b.s) > 1e-8) return a.s < b.s;
    return a.particle < b.particle;
  });
  return selected;
}

TGraphErrors *makeGraph(const std::vector<FlowPoint> &points, const char *name)
{
  TGraphErrors *g = new TGraphErrors(points.size());
  g->SetName(name);
  for (int i = 0; i < (int)points.size(); ++i) {
    g->SetPoint(i, points[i].s, points[i].y);
    g->SetPointError(i, 0.0, points[i].ey);
  }
  g->SetLineWidth(2);
  return g;
}

void setGraph(TGraphErrors *g, int markerStyle, int color, double markerSize = 1.25)
{
  if (!g) return;
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerColor(color);
  g->SetMarkerSize(markerSize);
  g->SetLineColor(color);
  g->SetLineWidth(2);
}

TH1D *drawFrame(const char *name, double xmin, double xmax, double ymin, double ymax, const char *xtitle, const char *ytitle)
{
  TH1D *h = new TH1D(name, "", 1, xmin, xmax);
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNoExponent(kTRUE);
  h->GetXaxis()->SetTitleSize(0.060);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.04);
  h->GetYaxis()->SetTitleSize(0.060);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.02);
  h->Draw();
  return h;
}

void drawZeroLine(double xmin, double xmax)
{
  drawLine(xmin, 0.0, xmax, 0.0, 1, 2, kGray + 2);
}

void plotFlow()
{
  style();

  std::vector<FlowPoint> dv1 = readFlowCSV("data/flow_dv1dy_vs_energy.csv");
  appendPoints(dv1, readFlowCSV("data/flow_dv1dy_low_energy_expanded.csv"));
  std::vector<FlowPoint> v2 = readFlowCSV("data/flow_v2_low_energy_expanded.csv");

  const double XMIN = 1.75;
  const double XMAX = 260.0;
  const double DV1MIN = -0.12;
  const double DV1MAX = 0.43;
  const double V2MIN = -0.095;
  const double V2MAX = 0.085;

  TCanvas *c1 = new TCanvas("c1", "", 1200, 560);
  c1->Draw();

  TPad *p1 = new TPad("p1", "", 0.00, 0.00, 0.50, 1.00);
  TPad *p2 = new TPad("p2", "", 0.50, 0.00, 1.00, 1.00);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.03);
  p1->SetTopMargin(0.05);
  p1->SetBottomMargin(0.16);
  p2->SetLeftMargin(0.13);
  p2->SetRightMargin(0.03);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.16);
  p1->SetLogx(1);
  p2->SetLogx(1);
  p1->Draw();
  p2->Draw();

  const int NDV1 = 7;
  const char *dv1Particle[NDV1] = {"proton", "net-proton", "antiproton", "pi+", "pi-", "Lambda", "K0S"};
  const char *dv1Label[NDV1] = {"p", "net-p", "#bar{p}", "#pi^{+}", "#pi^{-}", "#Lambda", "K^{0}_{S}"};
  const int dv1Color[NDV1] = {kBlue + 1, kBlack, kRed + 1, kGreen + 2, kViolet + 1, kOrange + 1, kCyan + 2};
  const int dv1Marker[NDV1] = {20, 33, 21, 22, 23, 29, 34};
  const int dv1OpenMarker[NDV1] = {24, 27, 25, 26, 32, 30, 28};

  p1->cd();
  drawFrame("h_dv1", XMIN, XMAX, DV1MIN, DV1MAX, "#sqrt{s_{NN}} (GeV)", "dv_{1}/dy |_{y=0}");
  drawZeroLine(XMIN, XMAX);

  TGraphErrors *gDv1Main[NDV1] = {0};
  TGraphErrors *gDv1Other[NDV1] = {0};
  for (int i = 0; i < NDV1; ++i) {
    std::vector<FlowPoint> starLike = selectParticle(dv1, dv1Particle[i], true);
    std::vector<FlowPoint> other = selectParticle(dv1, dv1Particle[i], false);
    gDv1Main[i] = makeGraph(starLike, Form("g_dv1_%d_main", i));
    gDv1Other[i] = makeGraph(other, Form("g_dv1_%d_other", i));
    setGraph(gDv1Main[i], dv1Marker[i], dv1Color[i], i < 5 ? 1.25 : 1.5);
    setGraph(gDv1Other[i], dv1OpenMarker[i], dv1Color[i], 1.25);
    if (gDv1Main[i]->GetN() > 1) gDv1Main[i]->Draw("LP SAME");
    else if (gDv1Main[i]->GetN() == 1) gDv1Main[i]->Draw("P SAME");
    if (gDv1Other[i]->GetN() > 0) gDv1Other[i]->Draw("P SAME");
  }
  drawText(2.02, 0.385, "(a) Directed flow", 42, 0.045);
  drawText(2.02, 0.340, "Au+Au/Pb+Pb, near midrapidity", 42, 0.030, 0, kGray + 2);

  TLegend *legDv1 = new TLegend(0.51, 0.56, 0.94, 0.92);
  legDv1->SetFillStyle(4000);
  legDv1->SetBorderSize(0);
  legDv1->SetTextSize(0.030);
  legDv1->SetNColumns(2);
  for (int i = 0; i < NDV1; ++i) {
    if (gDv1Main[i] && gDv1Main[i]->GetN() > 0) legDv1->AddEntry(gDv1Main[i], dv1Label[i], "p");
  }
  TGraphErrors *gOpen = new TGraphErrors(1);
  setGraph(gOpen, 24, kGray + 2, 1.25);
  legDv1->AddEntry(gOpen, "older/SIS data", "p");
  legDv1->Draw();
  drawHistBox(XMIN, XMAX, DV1MIN, DV1MAX);

  const int NV2 = 8;
  const char *v2Experiment[NV2] = {"FOPI", "EOS/E895/E877", "STAR FXT", "CERES", "NA49", "STAR", "PHENIX", "PHOBOS"};
  const char *v2Label[NV2] = {"FOPI", "EOS/E895/E877", "STAR FXT", "CERES", "NA49", "STAR", "PHENIX", "PHOBOS"};
  const int v2Color[NV2] = {kBlue + 1, kOrange + 7, kBlack, kOrange + 1, kCyan + 2, kGreen + 2, kRed + 1, kViolet + 1};
  const int v2Marker[NV2] = {20, 33, 25, 21, 34, 22, 23, 32};

  p2->cd();
  drawFrame("h_v2", XMIN, XMAX, V2MIN, V2MAX, "#sqrt{s_{NN}} (GeV)", "v_{2}");
  drawZeroLine(XMIN, XMAX);

  TGraphErrors *gV2[NV2] = {0};
  for (int i = 0; i < NV2; ++i) {
    std::vector<FlowPoint> selected = selectExperiment(v2, v2Experiment[i]);
    gV2[i] = makeGraph(selected, Form("g_v2_%d", i));
    setGraph(gV2[i], v2Marker[i], v2Color[i], 1.25);
    if (std::string(v2Experiment[i]) == "FOPI" || std::string(v2Experiment[i]) == "EOS/E895/E877") {
      gV2[i]->Draw("LP SAME");
    } else {
      gV2[i]->Draw("P SAME");
    }
  }
  drawText(2.02, 0.070, "(b) Elliptic flow", 42, 0.045);
  drawText(2.02, 0.056, "centrality/species differ across data sets", 42, 0.030, 0, kGray + 2);

  TLegend *legV2 = new TLegend(0.50, 0.16, 0.94, 0.48);
  legV2->SetFillStyle(4000);
  legV2->SetBorderSize(0);
  legV2->SetTextSize(0.030);
  legV2->SetNColumns(2);
  for (int i = 0; i < NV2; ++i) {
    if (gV2[i] && gV2[i]->GetN() > 0) legV2->AddEntry(gV2[i], v2Label[i], "p");
  }
  legV2->Draw();
  drawHistBox(XMIN, XMAX, V2MIN, V2MAX);

  gSystem->mkdir("fig", kTRUE);
  c1->Update();
  c1->SaveAs("fig/Flow_dv1_v2_two_panel.pdf");
  c1->SaveAs("fig/Flow_dv1_v2_two_panel.png");
}
