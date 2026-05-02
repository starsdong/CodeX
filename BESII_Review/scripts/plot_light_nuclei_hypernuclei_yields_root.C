#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

struct CsvRow {
  std::map<std::string, std::string> fields;
};

struct Series {
  std::string label;
  int marker = 20;
  Color_t color = kBlack;
  double marker_size = 1.4;
  int line_width = 2;
  bool draw_line = true;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
};

struct ModelCurve {
  std::string label;
  int color = kBlack;
  int style = 1;
  int width = 2;
  std::vector<double> x;
  std::vector<double> y;
};

std::vector<std::string> splitCsvLine(const std::string &line)
{
  std::vector<std::string> fields;
  std::string field;
  bool in_quotes = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (c == ',' && !in_quotes) {
      fields.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  fields.push_back(field);
  return fields;
}

bool parseDouble(const std::string &text, double &value)
{
  if (text.empty()) return false;
  char *end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return end != text.c_str();
}

std::string field(const CsvRow &row, const std::string &name)
{
  const auto it = row.fields.find(name);
  if (it == row.fields.end()) return "";
  return it->second;
}

bool hasObservable(const CsvRow &row, const std::set<std::string> &observables)
{
  return observables.count(field(row, "observable")) > 0;
}

double rowValue(const CsvRow &row, const std::string &name, double fallback = 0.0)
{
  double value = fallback;
  parseDouble(field(row, name), value);
  return value;
}

double rowError(const CsvRow &row)
{
  double total = 0.0;
  if (parseDouble(field(row, "total_err"), total)) return total;
  const double stat = rowValue(row, "stat_err", 0.0);
  const double sys = rowValue(row, "sys_err", 0.0);
  return std::sqrt(stat * stat + sys * sys);
}

double multiplicityHigh(const CsvRow &row)
{
  const std::string note = field(row, "kinematic_note");
  const std::size_t open = note.find('[');
  const std::size_t comma = note.find(',', open == std::string::npos ? 0 : open);
  const std::size_t close = note.find(']', comma == std::string::npos ? 0 : comma);
  if (open == std::string::npos || comma == std::string::npos || close == std::string::npos) {
    return -1.0;
  }
  double value = -1.0;
  parseDouble(note.substr(comma + 1, close - comma - 1), value);
  return value;
}

bool isCentralSummaryRow(const CsvRow &row)
{
  const std::string cent = field(row, "centrality");
  return cent == "0-10%" || cent == "0-20%" || cent == "max Npart row (most central)";
}

std::string multiplicityGroupKey(const CsvRow &row)
{
  return field(row, "experiment") + "|" + field(row, "system") + "|" + field(row, "sqrt_sNN_GeV") + "|" +
         field(row, "particle") + "|" + field(row, "matter") + "|" + field(row, "observable");
}

std::vector<CsvRow> readCsv(const TString &path)
{
  std::vector<CsvRow> rows;
  std::ifstream input(path.Data());
  if (!input) {
    std::cerr << "Cannot open " << path << std::endl;
    return rows;
  }

  std::string line;
  if (!std::getline(input, line)) return rows;
  const std::vector<std::string> header = splitCsvLine(line);
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> values = splitCsvLine(line);
    CsvRow row;
    for (std::size_t i = 0; i < header.size(); ++i) {
      row.fields[header[i]] = i < values.size() ? values[i] : "";
    }
    rows.push_back(row);
  }
  return rows;
}

std::vector<CsvRow> selectSummaryRows(const std::vector<CsvRow> &rows)
{
  const std::set<std::string> rapidity_observables = {
      "dN/dy", "BR*dN/dy", "dN/dy (BR=0.23 assumed)", "dN/dy (BR=0.25)",
  };

  std::vector<CsvRow> selected;
  std::map<std::string, std::vector<CsvRow>> multiplicity_groups;

  for (const CsvRow &row : rows) {
    const std::string obs = field(row, "observable");
    if (obs == "4pi_yield") {
      selected.push_back(row);
      continue;
    }
    if (!rapidity_observables.count(obs)) continue;

    if (field(row, "centrality").empty() && multiplicityHigh(row) >= 0.0) {
      multiplicity_groups[multiplicityGroupKey(row)].push_back(row);
      continue;
    }
    if (isCentralSummaryRow(row)) selected.push_back(row);
  }

  for (const auto &entry : multiplicity_groups) {
    const auto best = std::max_element(
        entry.second.begin(),
        entry.second.end(),
        [](const CsvRow &a, const CsvRow &b) { return multiplicityHigh(a) < multiplicityHigh(b); });
    if (best != entry.second.end()) selected.push_back(*best);
  }

  return selected;
}

std::string seriesKey(const CsvRow &row)
{
  return field(row, "particle") + "|" + field(row, "matter") + "|" + field(row, "observable");
}

std::string displayLabel(const CsvRow &row)
{
  const std::string particle = field(row, "particle");
  const std::string obs = field(row, "observable");

  if (particle == "d") return "d";
  if (particle == "anti-d") return "anti-d";
  if (particle == "t" && obs == "4pi_yield") return "t 4#pi";
  if (particle == "t") return "t dN/dy (STAR)";
  if (particle == "anti-t") return "anti-t";
  if (particle == "He3") return "{}^{3}He";
  if (particle == "anti-He3") return "anti-{}^{3}He";
  if (particle == "He4") return "{}^{4}He";
  if (particle == "anti-He4") return "anti-{}^{4}He";
  if (particle == "H4Lambda") return "{}_{#Lambda}^{4}H BR#timesdN/dy";
  if (particle == "H3Lambda" && obs == "BR*dN/dy") return "{}_{#Lambda}^{3}H BR#timesdN/dy";
  if (particle == "H3Lambda") return "{}_{#Lambda}^{3}H dN/dy, BR=0.23";
  if (particle == "(H3Lambda+anti-H3Lambda)/2") return "({}_{#Lambda}^{3}H+anti)/2, BR=0.25";
  return particle;
}

void applyStyleForRow(const CsvRow &row, Series &series)
{
  const std::string particle = field(row, "particle");
  const std::string obs = field(row, "observable");
  series.label = displayLabel(row);
  series.marker_size = 1.35;
  series.draw_line = true;

  if (particle == "d") {
    series.color = kBlue + 1;
    series.marker = 20;
  } else if (particle == "anti-d") {
    series.color = kBlue + 1;
    series.marker = 24;
  } else if (particle == "t") {
    series.color = kOrange + 7;
    series.marker = 21;
    series.marker_size = 1.55;
    series.line_width = 3;
  } else if (particle == "anti-t") {
    series.color = kOrange + 7;
    series.marker = 25;
  } else if (particle == "He3") {
    series.color = kGreen + 2;
    series.marker = 22;
  } else if (particle == "anti-He3") {
    series.color = kGreen + 2;
    series.marker = 26;
  } else if (particle == "He4") {
    series.color = kRed + 1;
    series.marker = 33;
    series.marker_size = 1.55;
  } else if (particle == "anti-He4") {
    series.color = kRed + 1;
    series.marker = 27;
    series.marker_size = 1.55;
  } else if (particle == "H4Lambda") {
    series.color = kOrange + 3;
    series.marker = 34;
    series.marker_size = 1.45;
  } else if (particle == "H3Lambda" && obs == "BR*dN/dy") {
    series.color = kViolet + 5;
    series.marker = 29;
    series.marker_size = 1.65;
  } else if (particle == "H3Lambda") {
    series.color = kViolet + 5;
    series.marker = 23;
    series.marker_size = 1.45;
  } else if (particle == "(H3Lambda+anti-H3Lambda)/2") {
    series.color = kViolet + 5;
    series.marker = 30;
    series.marker_size = 1.45;
  }
}

void sortSeries(Series &series)
{
  struct Point {
    double x;
    double y;
    double ey;
  };
  std::vector<Point> points;
  for (std::size_t i = 0; i < series.x.size(); ++i) {
    points.push_back({series.x[i], series.y[i], series.ey[i]});
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    series.x[i] = points[i].x;
    series.y[i] = points[i].y;
    series.ey[i] = points[i].ey;
  }
}

std::map<std::string, Series> buildSeries(const std::vector<CsvRow> &rows, bool four_pi)
{
  const std::set<std::string> rapidity_observables = {
      "dN/dy", "BR*dN/dy", "dN/dy (BR=0.23 assumed)", "dN/dy (BR=0.25)",
  };

  std::map<std::string, Series> series_map;
  for (const CsvRow &row : rows) {
    if (field(row, "experiment") == "ALICE") continue;

    const std::string obs = field(row, "observable");
    if (four_pi) {
      if (obs != "4pi_yield") continue;
    } else if (!rapidity_observables.count(obs)) {
      continue;
    }

    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "sqrt_sNN_GeV"), x) || !parseDouble(field(row, "value"), y) || y <= 0.0) {
      continue;
    }
    if (x > 300.0) continue;

    const std::string key = seriesKey(row);
    Series &series = series_map[key];
    if (series.label.empty()) applyStyleForRow(row, series);
    series.x.push_back(x);
    series.y.push_back(y);
    series.ey.push_back(rowError(row));
  }

  for (auto &entry : series_map) sortSeries(entry.second);
  return series_map;
}

std::vector<std::string> sortedKeys(const std::map<std::string, Series> &series_map)
{
  std::vector<std::string> keys;
  for (const auto &entry : series_map) keys.push_back(entry.first);
  std::sort(keys.begin(), keys.end(), [&](const std::string &a, const std::string &b) {
    const std::string &la = series_map.at(a).label;
    const std::string &lb = series_map.at(b).label;
    return la < lb;
  });
  return keys;
}

TGraphErrors *makeGraph(const Series &series)
{
  if (series.x.empty()) return nullptr;
  std::vector<double> ex(series.x.size(), 0.0);
  TGraphErrors *graph = new TGraphErrors(
      static_cast<int>(series.x.size()),
      const_cast<double *>(series.x.data()),
      const_cast<double *>(series.y.data()),
      ex.data(),
      const_cast<double *>(series.ey.data()));
  graph->SetMarkerStyle(series.marker);
  graph->SetMarkerSize(series.marker_size);
  graph->SetMarkerColor(series.color);
  graph->SetLineColor(series.color);
  graph->SetLineWidth(series.line_width);
  return graph;
}

ModelCurve readModelCurve(const TString &path, const std::string &label, int color, int style)
{
  ModelCurve curve;
  curve.label = label;
  curve.color = color;
  curve.style = style;
  curve.width = 3;

  std::ifstream input(path.Data());
  if (!input) return curve;

  std::string a;
  std::string b;
  while (input >> a >> b) {
    double x = 0.0;
    double y = 0.0;
    if (parseDouble(a, x) && parseDouble(b, y) && y > 0.0) {
      curve.x.push_back(x);
      curve.y.push_back(y);
    }
  }
  return curve;
}

TGraph *makeModelGraph(const ModelCurve &curve)
{
  if (curve.x.empty()) return nullptr;
  TGraph *graph = new TGraph(
      static_cast<int>(curve.x.size()),
      const_cast<double *>(curve.x.data()),
      const_cast<double *>(curve.y.data()));
  graph->SetLineColor(curve.color);
  graph->SetLineStyle(curve.style);
  graph->SetLineWidth(curve.width);
  return graph;
}

void applyFreezeOutStyle()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  gStyle->SetTickLength(0.04, "X");
  gStyle->SetTickLength(0.04, "Y");
  gStyle->SetGridWidth(1);
  gStyle->SetGridStyle(2);
  gStyle->SetGridColor(kGray + 1);

  gStyle->SetFillColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);

  gStyle->SetNdivisions(108, "X");
  gStyle->SetNdivisions(108, "Y");
  gStyle->SetLabelOffset(0.005, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetLabelFont(42, "X");
  gStyle->SetLabelFont(42, "Y");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTitleSize(0.07, "X");
  gStyle->SetTitleSize(0.07, "Y");
  gStyle->SetTitleFont(42, "X");
  gStyle->SetTitleFont(42, "Y");
  gStyle->SetMarkerSize(1.8);
  gStyle->SetMarkerStyle(20);
  gStyle->SetLegendFillColor(10);
  gROOT->ForceStyle();
}

void drawFrameBox(double x1, double x2, double y1, double y2)
{
  TLine *bottom = new TLine(x1, y1, x2, y1);
  TLine *top = new TLine(x1, y2, x2, y2);
  TLine *left = new TLine(x1, y1, x1, y2);
  TLine *right = new TLine(x2, y1, x2, y2);
  for (TLine *line : {bottom, top, left, right}) {
    line->SetLineColor(kBlack);
    line->SetLineWidth(3);
    line->Draw("same");
  }
}

void drawText(double x, double y, const char *text, double size = 0.045)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetTextFont(42);
  latex->SetTextSize(size);
  latex->Draw("same");
}

void drawPanel(
    const char *frame_name,
    const std::map<std::string, Series> &series_map,
    const std::vector<ModelCurve> &models,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    const char *ylabel,
    const char *title,
    bool draw_x_title,
    double title_x,
    double title_y)
{
  TH1D *frame = new TH1D(frame_name, "", 1, xmin, xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetYaxis()->SetTitle(ylabel);
  frame->GetXaxis()->SetTitleSize(draw_x_title ? 0.055 : 0.0);
  frame->GetXaxis()->SetLabelSize(draw_x_title ? 0.045 : 0.0);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleSize(0.057);
  frame->GetYaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetLabelSize(0.043);
  frame->Draw();

  drawText(title_x, title_y, title, 0.038);

  std::vector<TGraph *> model_graphs;
  for (const ModelCurve &curve : models) {
    TGraph *graph = makeModelGraph(curve);
    if (!graph) continue;
    model_graphs.push_back(graph);
    graph->Draw("C SAME");
  }

  std::vector<TGraphErrors *> data_graphs;
  for (const std::string &key : sortedKeys(series_map)) {
    const Series &series = series_map.at(key);
    TGraphErrors *graph = makeGraph(series);
    if (!graph) continue;
    data_graphs.push_back(graph);
    graph->Draw(series.draw_line && series.x.size() > 1 ? "LP SAME" : "P SAME");
  }

  drawFrameBox(xmin, xmax, ymin, ymax);
}

void drawLegends(
    const std::map<std::string, Series> &rapidity,
    const std::vector<ModelCurve> &models)
{
  TLegend *model_leg = new TLegend(0.72, 0.74, 0.99, 0.92);
  model_leg->SetFillStyle(4000);
  model_leg->SetBorderSize(0);
  model_leg->SetTextSize(0.030);
  model_leg->SetHeader("Model curves", "C");
  for (const ModelCurve &curve : models) {
    TGraph *dummy = new TGraph();
    dummy->SetLineColor(curve.color);
    dummy->SetLineStyle(curve.style);
    dummy->SetLineWidth(curve.width);
    model_leg->AddEntry(dummy, curve.label.c_str(), "l");
  }
  model_leg->Draw();

  TLegend *data_leg = new TLegend(0.72, 0.08, 0.99, 0.72);
  data_leg->SetFillStyle(4000);
  data_leg->SetBorderSize(0);
  data_leg->SetTextSize(0.023);
  data_leg->SetNColumns(1);
  data_leg->SetHeader("Data", "C");

  std::set<std::string> seen;
  for (const std::string &key : sortedKeys(rapidity)) {
    const Series &series = rapidity.at(key);
    if (seen.count(series.label)) continue;
    seen.insert(series.label);
    TGraphErrors *dummy = new TGraphErrors();
    dummy->SetMarkerStyle(series.marker);
    dummy->SetMarkerSize(series.marker_size * 0.85);
    dummy->SetMarkerColor(series.color);
    dummy->SetLineColor(series.color);
    dummy->SetLineWidth(series.line_width);
    data_leg->AddEntry(dummy, series.label.c_str(), "p");
  }
  data_leg->Draw();
}

}  // namespace

void plot_light_nuclei_hypernuclei_yields_root()
{
  applyFreezeOutStyle();

  TString macro_path = __FILE__;
  if (!macro_path.BeginsWith("/")) {
    macro_path = TString(gSystem->WorkingDirectory()) + "/" + macro_path;
  }
  const TString script_dir = gSystem->DirName(macro_path);
  const TString repo_dir = gSystem->DirName(script_dir);
  const TString input_path = repo_dir + "/data/light_nuclei_hypernuclei_yields_vs_energy.csv";
  const TString model_dir = repo_dir + "/data/data_for_Xin/model";
  const TString out_pdf = repo_dir + "/data/light_nuclei_hypernuclei_yields_root_style.pdf";
  const TString out_png = repo_dir + "/data/light_nuclei_hypernuclei_yields_root_style.png";

  const std::vector<CsvRow> rows = readCsv(input_path);
  const std::vector<CsvRow> summary = selectSummaryRows(rows);
  if (summary.empty()) {
    std::cerr << "No summary rows loaded from " << input_path << std::endl;
    return;
  }

  std::map<std::string, Series> rapidity = buildSeries(summary, false);
  std::vector<ModelCurve> models = {
      readModelCurve(model_dir + "/Thermal_fist_H3L_dNdy.txt", "Thermal-FIST", kBlack, 1),
      readModelCurve(model_dir + "/UrQMD_coal_H3L_dNdy.txt", "UrQMD+coal.", kTeal + 3, 2),
  };

  TCanvas *canvas = new TCanvas("cLightNucleiHypernucleiYields", "", 1100, 700);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.30);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  drawPanel(
      "frameLightNucleiRapidity",
      rapidity,
      models,
      2.4,
      300.0,
      1.0e-4,
      4.0e1,
      "yield near midrapidity",
      "Light nuclei and hypernuclei rapidity-density yields",
      true,
      4.2,
      7.5);
  drawLegends(rapidity, models);

  canvas->Update();
  canvas->SaveAs(out_pdf);
  canvas->SaveAs(out_png);

  std::cout << "wrote " << out_pdf << std::endl;
  std::cout << "wrote " << out_png << std::endl;
}
