#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TF1.h"

void plotFinal2D()
{

  std::vector<float> hildX{0.0032608695652174002, 0.026086956521739174, 0.05869565217391309, 0.09130434782608698, 0.12608695652173912, 0.15869565217391307, 0.19130434782608702, 0.25652173913043486, 0.3891304347826088, 0.5217391304347827, 0.7521739130434782, 0.8521739130434783, 0.9826086956521739, 1.3456521739130434, 1.4782608695652173};
  std::vector<float> hildY{0.9949999999999999, 0.986, 0.98, 0.9769999999999999, 0.9739999999999999, 0.971, 0.9694999999999999, 0.965, 0.9574999999999999, 0.95, 0.938, 0.9319999999999999, 0.9259999999999999, 0.9079999999999999, 0.9019999999999999};
  for (auto &y : hildY)
  {
    y = 1.f / y;
  }

  std::vector<float> hildXeD{0.004347826086956552, 0.017391304347826125, 0.03043478260869567, 0.052173913043478265, 0.07608695652173919, 0.10652173913043478, 0.13478260869565223, 0.17173913043478267, 0.22391304347826085, 0.3630434782608696, 0.45869565217391306, 0.5239130434782608, 0.5956521739130435, 0.6804347826086957, 0.7391304347826088, 0.7978260869565218, 1.017391304347826, 1.2043478260869565, 1.3478260869565215, 1.4804347826086954, 1.5};
  std::vector<float> hildYeD{0.965, 0.9409999999999998, 0.9139999999999999, 0.8899999999999999, 0.8689999999999999, 0.8449999999999999, 0.8269999999999998, 0.8059999999999998, 0.7789999999999999, 0.722, 0.6889999999999998, 0.6679999999999999, 0.6469999999999999, 0.6229999999999999, 0.6079999999999998, 0.5929999999999999, 0.5419999999999998, 0.5029999999999997, 0.47599999999999976, 0.45199999999999974, 0.44899999999999984};
  for (auto &y : hildYeD)
  {
    y = 1.f / y;
  }

  std::vector<float> hildXeU{0.004347826086956552, 0.0195652173913044, 0.04782608695652177, 0.07391304347826094, 0.13478260869565223, 0.19565217391304351, 0.28695652173913044, 0.34347826086956523, 0.40652173913043477, 0.5, 0.5608695652173914, 0.6543478260869565, 0.717391304347826, 0.7804347826086957, 0.9086956521739132, 1.0282608695652176, 1.1956521739130435, 1.402173913043478, 1.5};
  std::vector<float> hildYeU{1.0219999999999998, 1.043, 1.0699999999999998, 1.0879999999999999, 1.121, 1.145, 1.175, 1.19, 1.205, 1.226, 1.238, 1.2530000000000001, 1.265, 1.274, 1.292, 1.307, 1.325, 1.343, 1.3519999999999999};
  for (auto &y : hildYeU)
  {
    y = 1.f / y;
  }

  TCanvas *cv = new TCanvas();
  cv->SetTopMargin(0.12);
  cv->SetBottomMargin(0.14);
  cv->SetRightMargin(0.025);

  auto pad = cv->DrawFrame(-0., 0.5, 0.5, 1.3, ";B_{#Lambda} (MeV);#tau/#tau_{#Lambda}");
  pad->GetXaxis()->SetTitleSize(0.07);
  pad->GetXaxis()->SetTitleOffset(0.9);

  pad->GetYaxis()->SetTitleSize(0.08);
  pad->GetYaxis()->SetTitleOffset(0.78);

  pad->GetYaxis()->SetLabelSize(0.05);
  pad->GetXaxis()->SetLabelSize(0.05);

  TGraph *gr = new TGraph(hildX.size(), hildX.data(), hildY.data());
  TGraph *grmin = new TGraph(hildXeU.size(), hildXeU.data(), hildYeU.data());
  TGraph *grmax = new TGraph(hildXeD.size(), hildXeD.data(), hildYeD.data());
  TGraph *grshade = new TGraph(2 * hildXeU.size());
  for (int i = 0; i < hildXeU.size(); i++)
  {
    grshade->SetPoint(i, hildXeU[i], hildYeU[i]);
    grshade->SetPoint(hildXeU.size() + i, hildXeU[hildXeU.size() - i - 1], grmax->Eval(hildXeU[hildXeU.size() - i - 1]));
  }
  grshade->SetFillStyle(3013);
  grshade->SetFillColor(kRed - 9);
  gr->SetLineColor(kRed);
  gr->Draw("L same");
  gr->SetFillStyle(3013);
  gr->SetFillColor(kRed - 9);
  grshade->Draw("f same");


  TF1 *dalitz = new TF1("dalitz", "1. / (1. + 0.14 * std::sqrt(x))", 0, 0.7);
  dalitz->SetLineColor(kMagenta);
  dalitz->SetLineWidth(3);
  dalitz->SetLineStyle(kDashed);
  dalitz->Draw("same");


  std::vector<float> galX{69.e-3, 135.e-3, 159.e-3, 410.e-3};
  std::vector<float> galY{234. / 263.2, 190. / 263.2, 180. / 263.2, 163. / 263.2};
  std::vector<float> galEY{27. / 263.2, 22. / 263.2, 21. / 263.2, 18. / 263.2};
  TGraphErrors *galGr = new TGraphErrors(galX.size(), galX.data(), galY.data(), 0x0, galEY.data());
  galGr->SetFillStyle(3354);
  galGr->SetFillColor(kOrange - 3);
  galGr->SetLineColor(kOrange - 3);
  galGr->SetLineWidth(3);
  galGr->Draw("lf3");


  float measX[1]{0.102};
  float measXstat[1]{0.063};
  float measXsyst[1]{0.067};
  float measY[1]{253. / 263.2};
  float measYstat[1]{11. / 263.2};
  float measYsyst[1]{6. / 263.2};

  float prd70X[1]{-0.24};
  float prd70Xstat[1]{0.22};
  float prd70Y[1]{264. / 263.2};
  float prd70Ylow[1]{52. / 263.2};
  float prd70Yup[1]{84. / 263.2};
  TGraphAsymmErrors *prd70Stat = new TGraphAsymmErrors(1, prd70X, prd70Y, prd70Xstat, prd70Xstat, prd70Ylow, prd70Yup);
  prd70Stat->SetLineColor(kBlack);
  prd70Stat->SetMarkerColor(kBlack);
  prd70Stat->SetMarkerStyle(20);

  TLatex measText;
  measText.SetTextFont(42);
  measText.SetTextColor(kBlue);
  measText.SetTextAlign(21);
  measText.SetTextSize(0.06);
  measText.DrawText(measX[0] - 0.03, 1.2, "ALICE");





  TGraphErrors *measStat = new TGraphErrors(1, measX, measY, measXstat, measYstat);
  TGraphErrors *measSyst = new TGraphErrors(1, measX, measY, measXsyst, measYsyst);
  measStat->SetLineColor(kBlue);
  measStat->SetMarkerColor(kBlue);
  measStat->SetMarkerStyle(20);
  measSyst->SetLineColor(kBlue);
  measSyst->SetMarkerColor(kBlue);
  measSyst->SetFillStyle(0);
  measSyst->SetFillColor(kBlue);
  measSyst->SetMarkerStyle(20);
  measStat->Draw("pZ");
  measSyst->Draw("p2Z");

  TLegend *leg = new TLegend(0.12, 0.88, 0.975, .98);
  leg->SetMargin(0.12);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->AddEntry(dalitz, "Nuo. Cim. A 46 (1966) 786", "lf");
  leg->AddEntry(galGr, "PLB 811 (2020) 135916", "lf");
  leg->AddEntry(gr, "PRC 102 (2020) 064002", "lf");
  leg->SetTextSize(0.05);
  leg->Draw();
  cv->SaveAs("plot.pdf");
  cv->SaveAs("plot.png");




}