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

struct NCQPoint {
  double energy = 0.0;
  std::string particle;
  std::string centrality;
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

std::vector<std::string> splitNCQCSV(const std::string &line)
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

double getDouble(const std::vector<std::string> &fields, const std::map<std::string, int> &index, const char *name)
{
  const auto it = index.find(name);
  if (it == index.end() || it->second >= (int)fields.size() || fields[it->second].empty()) return 0.0;
  return atof(fields[it->second].c_str());
}

std::string getString(const std::vector<std::string> &fields, const std::map<std::string, int> &index, const char *name)
{
  const auto it = index.find(name);
  if (it == index.end() || it->second >= (int)fields.size()) return "";
  return fields[it->second];
}

std::string canonicalParticle(std::string particle)
{
  if (particle == "proton") return "p";
  if (particle == "pion") return "#pi";
  return particle;
}

std::vector<NCQPoint> readNCQCSV(const char *path)
{
  std::vector<NCQPoint> points;
  std::ifstream in(path);
  if (!in.good()) {
    printf("Cannot open %s\n", path);
    return points;
  }

  std::string line;
  if (!std::getline(in, line)) return points;
  std::vector<std::string> header = splitNCQCSV(line);
  std::map<std::string, int> index;
  for (int i = 0; i < (int)header.size(); ++i) index[header[i]] = i;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<std::string> fields = splitNCQCSV(line);
    NCQPoint p;
    p.energy = getDouble(fields, index, "energy");
    p.particle = canonicalParticle(getString(fields, index, "particle"));
    p.centrality = getString(fields, index, "centrality");
    p.x = getDouble(fields, index, "mt_minus_m0_over_nq");
    p.y = getDouble(fields, index, "v2_over_nq");
    const double stat = getDouble(fields, index, "stat_over_nq");
    const double syst = getDouble(fields, index, "syst_over_nq");
    const double glob = getDouble(fields, index, "glob_over_nq");
    p.ey = std::sqrt(stat * stat + syst * syst + glob * glob);
    if (p.energy > 0.0) points.push_back(p);
  }
  return points;
}

std::vector<NCQPoint> selectNCQ(const std::vector<NCQPoint> &all, double energy, const char *particle)
{
  std::vector<NCQPoint> selected;
  for (const auto &p : all) {
    if (std::fabs(p.energy - energy) < 1e-6 && p.particle == particle) selected.push_back(p);
  }
  std::sort(selected.begin(), selected.end(), [](const NCQPoint &a, const NCQPoint &b) { return a.x < b.x; });
  return selected;
}

bool skipPanelSpecies(double energy, const char *particle)
{
  std::string sp = particle;
  if (std::fabs(energy - 3.0) < 1e-6 && (sp == "K0S" || sp == "Lambda")) return true;
  if (std::fabs(energy - 7.7) < 1e-6 && (sp == "K-" || sp == "pbar")) return true;
  return false;
}

std::string centralityForEnergy(const std::vector<NCQPoint> &all, double energy)
{
  for (const auto &p : all) {
    if (std::fabs(p.energy - energy) < 1e-6) return p.centrality;
  }
  return "";
}

TGraphErrors *makeNCQGraph(const std::vector<NCQPoint> &points, const char *name)
{
  TGraphErrors *g = new TGraphErrors(points.size());
  g->SetName(name);
  for (int i = 0; i < (int)points.size(); ++i) {
    g->SetPoint(i, points[i].x, points[i].y);
    g->SetPointError(i, 0.0, points[i].ey);
  }
  g->SetLineWidth(2);
  return g;
}

void setNCQStyle(TGraphErrors *g, int markerStyle, int color, double markerSize = 1.15)
{
  if (!g) return;
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerColor(color);
  g->SetMarkerSize(markerSize);
  g->SetLineColor(color);
  g->SetLineWidth(2);
}

TH1D *drawNCQFrame(const char *name, double xmin, double xmax, double ymin, double ymax, bool leftColumn, bool bottomRow)
{
  TH1D *h = new TH1D(name, "", 1, xmin, xmax);
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
  h->GetYaxis()->SetTitle("v_{2}/n_{q}");
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleSize(0.060);
  h->GetXaxis()->SetLabelSize(0.050);
  h->GetXaxis()->SetTitleOffset(1.02);
  h->GetYaxis()->SetTitleSize(0.064);
  h->GetYaxis()->SetLabelSize(0.050);
  h->GetYaxis()->SetTitleOffset(1.12);
  h->Draw();
  return h;
}

void plotStarV2NCQ2x2()
{
  style();

  std::vector<NCQPoint> all = readNCQCSV("data/star_identified_v2_ncq.csv");

  const int NE = 4;
  const double energy[NE] = {3.0, 4.5, 7.7, 19.6};
  const char *energyLabel[NE] = {"3.0", "4.5", "7.7", "19.6"};
  const double xmax[NE] = {0.65, 0.95, 1.15, 1.75};
  const double ymin[NE] = {-0.040, -0.005, -0.005, -0.005};
  const double ymax[NE] = {0.008, 0.072, 0.082, 0.082};

  const int NS = 15;
  const char *species[NS] = {"pi+", "pi-", "K+", "K-", "K0S", "p", "pbar", "phi", "Lambda", "Lambdabar", "Xi-", "Xibar+", "Omega-", "Omegabar+", "#pi"};
  const char *label[NS] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "K^{0}_{S}", "p", "#bar{p}", "#phi", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#bar{#Xi}^{+}", "#Omega^{-}", "#bar{#Omega}^{+}", "#pi"};
  const int color[NS] = {kRed + 1, kRed - 4, kBlue + 1, kBlue - 7, kCyan + 2, kBlack, kGray + 2, kViolet + 1, kGreen + 2, kGreen - 6, kOrange + 7, kOrange + 1, kMagenta + 2, kMagenta - 7, kRed + 1};
  const int marker[NS] = {20, 24, 21, 25, 33, 22, 26, 34, 23, 32, 29, 30, 47, 46, 20};
  const double markerSize[NS] = {1.05, 1.05, 1.05, 1.05, 1.28, 1.05, 1.05, 1.18, 1.18, 1.18, 1.28, 1.28, 1.35, 1.35, 1.05};

  TCanvas *c1 = new TCanvas("c1", "", 1100, 980);
  c1->Draw();

  TPad *pad[NE];
  pad[0] = new TPad("p_3gev", "", 0.00, 0.54, 0.50, 0.93);
  pad[1] = new TPad("p_4p5gev", "", 0.50, 0.54, 1.00, 0.93);
  pad[2] = new TPad("p_7p7gev", "", 0.00, 0.02, 0.50, 0.50);
  pad[3] = new TPad("p_19p6gev", "", 0.50, 0.02, 1.00, 0.50);

  for (int i = 0; i < NE; ++i) {
    pad[i]->SetLeftMargin(0.16);
    pad[i]->SetRightMargin(0.03);
    pad[i]->SetTopMargin(i < 2 ? 0.06 : 0.04);
    pad[i]->SetBottomMargin(0.17);
    pad[i]->Draw();
  }

  TGraphErrors *graphs[NE][NS] = {{0}};

  for (int ie = 0; ie < NE; ++ie) {
    pad[ie]->cd();
    const bool left = (ie % 2 == 0);
    const bool bottom = (ie >= 2);
    drawNCQFrame(Form("h_ncq_%d", ie), 0.0, xmax[ie], ymin[ie], ymax[ie], left, bottom);
    drawLine(0.0, 0.0, xmax[ie], 0.0, 1, 2, kGray + 2);

    for (int is = 0; is < NS; ++is) {
      if (skipPanelSpecies(energy[ie], species[is])) continue;
      std::vector<NCQPoint> selected = selectNCQ(all, energy[ie], species[is]);
      if (selected.empty()) continue;
      graphs[ie][is] = makeNCQGraph(selected, Form("g_ncq_%d_%d", ie, is));
      setNCQStyle(graphs[ie][is], marker[is], color[is], markerSize[is]);
      graphs[ie][is]->Draw("P SAME");
    }

    const double yrange = ymax[ie] - ymin[ie];
    drawText(0.04 * xmax[ie], ymax[ie] - 0.12 * yrange, Form("(%c) #sqrt{s_{NN}} = %s GeV", 'a' + ie, energyLabel[ie]), 42, 0.060);
    drawHistBox(0.0, xmax[ie], ymin[ie], ymax[ie]);
  }

  c1->cd();
  TLegend *leg = new TLegend(0.08, 0.935, 0.98, 0.995);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.030);
  leg->SetNColumns(8);
  for (int is = 0; is < NS; ++is) {
    bool present = false;
    for (int ie = 0; ie < NE; ++ie) {
      if (graphs[ie][is] && graphs[ie][is]->GetN() > 0) present = true;
    }
    if (!present) continue;
    TGraphErrors *gdummy = new TGraphErrors(1);
    gdummy->SetPoint(0, 0.0, 0.0);
    setNCQStyle(gdummy, marker[is], color[is], markerSize[is]);
    leg->AddEntry(gdummy, label[is], "p");
  }
  leg->Draw();

  gSystem->mkdir("fig", kTRUE);
  c1->Update();
  c1->SaveAs("fig/STAR_v2_NCQ_2x2_3_4p5_7p7_19p6.pdf");
  c1->SaveAs("fig/STAR_v2_NCQ_2x2_3_4p5_7p7_19p6.png");
}
