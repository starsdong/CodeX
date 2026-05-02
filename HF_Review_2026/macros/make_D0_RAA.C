#include <iomanip>

using namespace std;

void make_D0_RAA()
{
  gROOT->Reset();
  gSystem->mkdir("fig", kTRUE);

  const TString kHFInputDir = "/Users/starsdong/Work/work/HF";

  const Double_t global_sys = sqrt(0.081*0.081+0.052*0.052); // global uncertainty
  const Int_t N = 6;
  const Char_t *FileName[N] = {"Duke","LBT","Catania","TAMU","PHSD","Nantes"};
  const Char_t *Label[N] = {"Duke","LBT","Catania","TAMU","PHSD","SUBATECH"};
  const Int_t kModelColor[N] = {kAzure+10, kBlue+1, kOrange+7, kGreen+2, kMagenta+1, kRed+1};
  const Int_t kModelStyle[N] = {10, 4, 5, 1, 8, 7};

  TGraph *gr[N];
  for(int i=0;i<N;i++) {
    gr[i] = new TGraph(Form("/Users/starsdong/work/datapoints/HQJET_201808/%s.txt",FileName[i]),"%lg %lg");
  }

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);

   // c1->SetGridx();
   // c1->SetGridy();
   c1->SetTickx();
   c1->SetTicky();

   c1->SetLeftMargin(0.16);
   c1->SetBottomMargin(0.16);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = 0;
   double x2 = 8.0;
   double y1 = 0.0;
   double y2 = 2.0;

   TH1D *d0 = new TH1D("d0","",1,x1,x2);
   d0->SetMinimum(y1);
   d0->SetMaximum(y2);
   d0->GetXaxis()->SetNdivisions(208);
   d0->GetXaxis()->CenterTitle();
   d0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   d0->GetXaxis()->SetTitleOffset(1.1);
   d0->GetXaxis()->SetTitleSize(0.06);
   d0->GetXaxis()->SetLabelOffset(0.01);
   d0->GetXaxis()->SetLabelSize(0.05);
   d0->GetXaxis()->SetLabelFont(42);
   d0->GetXaxis()->SetTitleFont(42);
   d0->GetYaxis()->SetNdivisions(505);
   d0->GetYaxis()->CenterTitle();
   d0->GetYaxis()->SetTitle("R_{AA}");
   d0->GetYaxis()->SetTitleOffset(1.0);
   d0->GetYaxis()->SetTitleSize(0.06);
   d0->GetYaxis()->SetLabelOffset(0.015);
   d0->GetYaxis()->SetLabelSize(0.05);
   d0->GetYaxis()->SetLabelFont(42);
   d0->GetYaxis()->SetTitleFont(42);
   d0->Draw("c");

   TLine *l0 = new TLine(x1,1,x2,1);
   l0->SetLineWidth(2);
   l0->SetLineStyle(2);
   l0->Draw("same");

   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   double xg1 = x2-0.2;
   double xg2 = x2-0.02;
   double yg1 = 1 - global_sys;
   double yg2 = 1 + global_sys;
   
   TBox *box = new TBox(xg1,yg1,xg2,yg2);
   box->SetFillColor(3);
   box->Draw("same");

   for(int i=0;i<N;i++) {
     gr[i]->SetLineColor(kModelColor[i]);
     gr[i]->SetLineWidth(2);
     gr[i]->SetLineStyle(kModelStyle[i]);
     gr[i]->Draw("c");
   }


   TFile *f1 = new TFile(Form("%s/root/D0RAA_Run14HFT_QM180510.root", kHFInputDir.Data()));
   TGraphErrors *gr_HFT_0_10_stat = (TGraphErrors *)f1->Get("D0_RAA_err_0_10");
   TGraphAsymmErrors *gr_HFT_0_10_sys = (TGraphAsymmErrors *)f1->Get("D0_RAA_sys_combine_0_10");
   TFile *f2 = new TFile(Form("%s/root/D0_RAA_AsymError1_run11_0507.root", kHFInputDir.Data()));
   TGraphErrors *gr_TPC_0_10_stat = (TGraphErrors *)f2->Get("RAA_0_10c_err");
   TGraphAsymmErrors *gr_TPC_0_10_sys = (TGraphAsymmErrors *)f2->Get("RAA_0_10c_sys");

   int n_hft = gr_HFT_0_10_stat->GetN();
   int n_tpc = gr_TPC_0_10_stat->GetN();

   const int color_tpc = 4;
   if(0) {
   for(int i=0;i<n_tpc;i++) {
     double x = gr_TPC_0_10_sys->GetX()[i];
     double y = gr_TPC_0_10_sys->GetY()[i];
     double ye_up = gr_TPC_0_10_sys->GetEYhigh()[i];
     double ye_down = gr_TPC_0_10_sys->GetEYlow()[i];
     
    double x1 = x - 0.08;
    double x2 = x + 0.08;
    double y1 = y - ye_down;
    double y2 = y + ye_up;
    
    TLine *la = new TLine(x1, y1, x1, y1+0.03);
    la->SetLineColor(color_tpc);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.03);
    lb->SetLineColor(color_tpc);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.03);
    lc->SetLineColor(color_tpc);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.03);
    ld->SetLineColor(color_tpc);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineColor(color_tpc);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineColor(color_tpc);
    lf->SetLineWidth(2);
    lf->Draw("same");
     
   }
   }

   const int color_hft = 4;
   for(int i=0;i<n_hft;i++) {
     double x = gr_HFT_0_10_sys->GetX()[i];
     double y = gr_HFT_0_10_sys->GetY()[i];
     double ye_up = gr_HFT_0_10_sys->GetEYhigh()[i];
     double ye_down = gr_HFT_0_10_sys->GetEYlow()[i];
     
    double x1 = x - 0.08;
    double x2 = x + 0.08;
    double y1 = y - ye_down;
    double y2 = y + ye_up;
    
    TLine *la = new TLine(x1, y1, x1, y1+0.03);
    la->SetLineColor(color_hft);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.03);
    lb->SetLineColor(color_hft);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.03);
    lc->SetLineColor(color_hft);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.03);
    ld->SetLineColor(color_hft);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineColor(color_hft);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineColor(color_hft);
    lf->SetLineWidth(2);
    lf->Draw("same");
     
   }

   gr_HFT_0_10_stat->SetMarkerStyle(20);
   gr_HFT_0_10_stat->SetMarkerColor(color_hft);
   gr_HFT_0_10_stat->SetMarkerSize(1.8);
   gr_HFT_0_10_stat->SetLineColor(color_hft);
   gr_HFT_0_10_stat->SetLineWidth(2);
   gr_HFT_0_10_stat->Draw("p");
   //   gr_HFT_0_10_sys->Draw("p");
   // gr_TPC_0_10_stat->Draw("p");
   //   gr_TPC_0_10_sys->Draw("p");


   TLatex * tex = new TLatex(0.6, 1.75, "Au+Au 200 GeV, 0-10%");
   tex->SetTextFont(42);
   tex->SetTextSize(0.055);
   tex->Draw("same");

   TLegend *leg = new TLegend(0.72, 0.65, 0.96, 0.94);
   leg->SetFillColor(10);
   leg->SetFillStyle(10);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.034);
   leg->SetMargin(0.28);
   for(int i=0;i<N;i++) {
     leg->AddEntry(gr[i],Label[i],"l");
   }
   leg->Draw();

   c1->SaveAs("fig/D0_RAA.pdf");
   c1->SaveAs("fig/D0_RAA.png");
   


   
   // Model test with chi2
   
   double chi2_HFT[N];
   for(int i=0;i<N;i++) chi2_HFT[i] = 0.0;
   double chi2_TPC[N];
   for(int i=0;i<N;i++) chi2_TPC[i] = 0.0;

   
   for(int i=0;i<n_hft;i++) {
     double x = gr_HFT_0_10_stat->GetX()[i];
     double y = gr_HFT_0_10_stat->GetY()[i];
     double ye = gr_HFT_0_10_stat->GetEY()[i];
     double yes_g = y*global_sys;
     for(int j=0;j<N;j++) {
       double y_th = gr[j]->Eval(x);
       double yes = 1.e9;
       if(y>y_th) {
	 yes = gr_HFT_0_10_sys->GetEYlow()[i];
       } else {
	 yes = gr_HFT_0_10_sys->GetEYhigh()[i];
       }
       chi2_HFT[j] += (y-y_th)*(y-y_th)/(ye*ye+yes*yes+yes_g*yes_g);
       if(j==0) {
	 cout << "x=" << x << "\ty=" << y << "\ty_th=" << "\t" << y_th << "\tye=" << ye << "\tyes=" << yes << "\tchi2=" << chi2_HFT[j] << endl;
       }
     }     
   }
   
   for(int i=0;i<n_tpc;i++) {
     double x = gr_TPC_0_10_stat->GetX()[i];
     double y = gr_TPC_0_10_stat->GetY()[i];
     double ye = gr_TPC_0_10_stat->GetEY()[i];
     double yes_g = y*global_sys;
     for(int j=0;j<N;j++) {
       double y_th = gr[j]->Eval(x);
       double yes = 1.e9;
       if(y>y_th) {
	 yes = gr_TPC_0_10_sys->GetEYlow()[i];
       } else {
	 yes = gr_TPC_0_10_sys->GetEYhigh()[i];
       }
       chi2_TPC[j] += (y-y_th)*(y-y_th)/(ye*ye+yes*yes+yes_g*yes_g);

       if(j==0) {
	 cout << "x=" << x << "\ty=" << y << "\ty_th=" << "\t" << y_th << "\tye=" << ye << "\tyes=" << yes << "\tchi2=" << chi2_TPC[j] << endl;
       }
     }     
   }

   cout << " ======= Chi2 Test ======= " << endl;
   for(int j=0;j<N;j++) {
     cout << Label[j] << "\t HFT\t" << chi2_HFT[j] << "/" << n_hft << "\t Prob=" << setw(10) << TMath::Prob(chi2_HFT[j],n_hft) << "\t TPC\t" << chi2_TPC[j] << "/" << n_tpc << "\t Prob=" << setw(10) << TMath::Prob(chi2_TPC[j],n_tpc) << endl;
   }

   
}
