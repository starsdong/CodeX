#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace {

TH1 *GetHist(TFile &file, const char *path)
{
  TH1 *hist = dynamic_cast<TH1 *>(file.Get(path));
  if (!hist) {
    std::cerr << "Missing histogram: " << file.GetName() << ":" << path << std::endl;
  }
  return hist;
}

TGraphErrors *MakeGraph(TH1 *hist, const char *name, bool useXErrors)
{
  auto *graph = new TGraphErrors();
  graph->SetName(name);

  if (!hist) return graph;

  int point = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double y = hist->GetBinContent(bin);
    const double ey = hist->GetBinError(bin);
    if (!std::isfinite(y) || y <= 0.0) continue;

    const double x = hist->GetBinCenter(bin);
    const double ex = useXErrors ? 0.5 * hist->GetBinWidth(bin) : 0.0;
    graph->SetPoint(point, x, y);
    graph->SetPointError(point, ex, std::isfinite(ey) ? ey : 0.0);
    ++point;
  }

  return graph;
}

void StyleGraph(TGraphErrors *graph, Color_t color, Style_t marker, Style_t line, int width)
{
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.05);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineStyle(line);
  graph->SetLineWidth(width);
}

void SavePoints(std::ofstream &csv, const char *component, const char *series, TH1 *hist)
{
  if (!hist) return;

  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double y = hist->GetBinContent(bin);
    if (!std::isfinite(y) || y <= 0.0) continue;

    csv << component << ","
        << series << ","
        << hist->GetBinLowEdge(bin) << ","
        << hist->GetBinLowEdge(bin) + hist->GetBinWidth(bin) << ","
        << hist->GetBinCenter(bin) << ","
        << y << ","
        << hist->GetBinError(bin) << "\n";
  }
}

void DrawOverlay(const char *component,
                 const char *title,
                 const char *yTitle,
                 TH1 *hData,
                 TH1 *hFastSim,
                 const char *outDir,
                 double xMax,
                 double yMax,
                 TFile &outFile,
                 std::ofstream &csv)
{
  auto *gData = MakeGraph(hData, Form("g_data_%s", component), true);
  auto *gFastSim = MakeGraph(hFastSim, Form("g_fastsim_%s", component), true);

  StyleGraph(gData, kBlack, 20, 1, 2);
  StyleGraph(gFastSim, kRed + 1, 24, 1, 3);

  TCanvas canvas(Form("c_data_vs_fastsim_%s", component),
                 Form("Data vs FastSim %s", component),
                 900, 700);
  canvas.SetLeftMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.08);
  canvas.SetBottomMargin(0.12);
  canvas.SetGrid(1, 1);

  TH1F frame(Form("frame_%s", component),
             Form("%s;p or <p> (GeV/c);%s", title, yTitle),
             100, 0.0, xMax);
  frame.SetMinimum(0.0);
  frame.SetMaximum(yMax);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleOffset(1.05);
  frame.GetYaxis()->SetTitleOffset(1.15);
  frame.Draw();

  gFastSim->Draw("PZ SAME");
  gData->Draw("PZ SAME");

  TLegend legend(0.45, 0.70, 0.92, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.037);
  legend.AddEntry(gData, "Data, iterative Gaussian", "lep");
  legend.AddEntry(gFastSim, "FastSim, reconstructed DCA", "lep");
  legend.Draw();

  const TString base = TString::Format("%s/data_vs_fastsim_two_track_dca_%s", outDir, component);
  canvas.SaveAs(base + ".pdf");
  canvas.SaveAs(base + ".png");

  outFile.cd();
  canvas.Write();
  gData->Write();
  gFastSim->Write();

  SavePoints(csv, component, "data_iterative_gaussian", hData);
  SavePoints(csv, component, "fastsim_reconstructed_dca", hFastSim);
}

}  // namespace

void OverlayDataFastSimDca(const char *dataResolutionFile = "resolution_vs_x_methods_final/output_PlottingMacro_data_resolution_methods.root",
                           const char *fastSimFile = "/Users/starsdong/Work/CodeX/FastSim/outputs/alltracks_dcasp_1M/histograms.root",
                           const char *outDir = "resolution_vs_x_methods_final/data_vs_fastsim_two_track_dca",
                           double xMax = 3.0,
                           double yMaxXY = 0.015,
                           double yMaxZ = 0.015)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);

  gSystem->mkdir(outDir, kTRUE);

  TFile dataFile(dataResolutionFile, "READ");
  if (dataFile.IsZombie()) {
    std::cerr << "Cannot open data resolution file: " << dataResolutionFile << std::endl;
    return;
  }

  TFile simFile(fastSimFile, "READ");
  if (simFile.IsZombie()) {
    std::cerr << "Cannot open FastSim file: " << fastSimFile << std::endl;
    return;
  }

  TFile outFile(TString::Format("%s/data_vs_fastsim_two_track_dca.root", outDir), "RECREATE");
  std::ofstream csv(TString::Format("%s/data_vs_fastsim_two_track_dca.csv", outDir).Data());
  csv << "component,series,x_low,x_high,x_center,sigma_cm,sigma_error_cm\n";

  TH1 *hDataXY = GetHist(dataFile, "hdcaxy_vs_preco_allparticles/iterative_sigma");
  TH1 *hDataZ = GetHist(dataFile, "hdcaz_vs_preco_allparticles/iterative_sigma");
  TH1 *hFastXY = GetHist(simFile, "rec_dca_xy_vs_pavg_smeared_2");
  TH1 *hFastZ = GetHist(simFile, "rec_dca_z_vs_pavg_smeared_2");

  DrawOverlay("xy",
              "Two-track DCA_{xy} resolution: data vs FastSim",
              "#sigma(DCA_{xy}) (cm)",
              hDataXY, hFastXY, outDir, xMax, yMaxXY, outFile, csv);

  DrawOverlay("z",
              "Two-track DCA_{z} resolution: data vs FastSim",
              "#sigma(DCA_{z}) (cm)",
              hDataZ, hFastZ, outDir, xMax, yMaxZ, outFile, csv);

  csv.close();
  outFile.Close();

  std::cout << "Wrote overlay plots to " << outDir << std::endl;
}
