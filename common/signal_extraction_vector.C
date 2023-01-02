
using namespace RooStats;
using namespace RooFit;
#include "Riostream.h"

Double_t Pvalue(Double_t significance)
{
  return ROOT::Math::chisquared_cdf_c(pow(significance, 2), 1) / 2;
}

void signal_extraction_vector()
{
  std::vector<std::string> framePaths = {"ct12/eff0.71_pol1", "ct24/eff0.76_pol1", "ct46/eff0.80_pol1", "ct68/eff0.83_pol1",
                                         "ct810/eff0.80_pol1", "ct1014/eff0.80_pol1", "ct1418/eff0.80_pol1", "ct1823/eff0.74_pol1", "ct2335/eff0.62_pol1"};
  std::vector<TCanvas> canvasVec;
  TFile fit_file("../Results/2Body/2body_analysis_upd_signal_extraction_dscb.root");
  TFile outFile("outfile.root", "RECREATE");

  for (unsigned int i = 0; i < framePaths.size(); i++)
  {
    RooPlot *frame = static_cast<RooPlot *>(fit_file.Get(framePaths[i].c_str()));
    // frame->remove();
    frame->SetTitle("");

    constexpr double kYsize = 0.4;

    gStyle->SetOptStat(1);
    gStyle->SetOptDate(0);
    gStyle->SetOptFit(1);
    gStyle->SetLabelSize(0.04, "xyz"); // size of axis value font
    gStyle->SetTitleSize(0.05, "xyz"); // size of axis title font
    gStyle->SetTitleFont(42, "xyz");   // font option
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleOffset(1.05, "x");
    gStyle->SetTitleOffset(1.1, "y");
    // default canvas options
    gStyle->SetCanvasDefW(800);
    gStyle->SetCanvasDefH(600);
    gStyle->SetPadBottomMargin(0.12); // margins...
    gStyle->SetPadTopMargin(0.01);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadGridX(0); // grids, tickmarks
    gStyle->SetPadGridY(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPaperSize(20, 24); // US letter size
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetEndErrorSize(0.);
    gStyle->SetMarkerSize(10);
    // gStyle->SetLineWidth(4);

    auto cv = TCanvas(Form("cv%i", i), Form("cv%i", i), 1500, 1500);
    cv.cd();

    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleOffset(0.9);
    frame->GetXaxis()->SetTitleOffset(1.1);
    // frame->GetYaxis()->SetTitle("Entries / (2.35 MeV/#it{c}^{2})");
    frame->GetXaxis()->SetTitle("#it{M}(^{3}He + #pi^{-} and c.c.)   (GeV/#it{c}^{2})");
    frame->SetMinimum(0.01);
    // frame->SetMarkerSize();
    TPaveText *pinfo_alice = new TPaveText(0.14, 0.75, 0.44, 0.93, "NDC");
    pinfo_alice->SetBorderSize(0);
    pinfo_alice->SetFillStyle(0);
    pinfo_alice->SetTextFont(42);
    pinfo_alice->AddText("ALICE, Pb#minusPb");
    pinfo_alice->AddText("0#minus90%, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pinfo_alice->AddText("2 #leq #it{p}_{T} < 9 GeV/#it{c}");

    frame->setDrawOptions("data", "PE0");
    frame->Draw();

    TLegend *leg1 = new TLegend(0.14, 0.52, 0.42, 0.72);
    leg1->SetBorderSize(0);
    leg1->AddEntry("data", "{}_{#Lambda}^{3}H + {}_{#bar{#Lambda}}^{3}#bar{H}", "PE");
    leg1->AddEntry("pol1_total_pdf", "Signal + Background", "L");
    leg1->AddEntry("pol1_total_pdf_Norm[m]_Comp[bkg]_Range[fit_nll_pol1_total_pdf_data]_NormRange[fit_nll_pol1_total_pdf_data]", "Background", "L");
    leg1->SetMargin(0.15);

    if (i == 0)
    {
      leg1->Draw();
      pinfo_alice->Draw();
    }
    frame->Print();

    outFile.cd();
    cv.Write();
  }
  outFile.Close();
}