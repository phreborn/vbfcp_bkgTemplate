#include "utils.h"

#ifdef __CLING__
#include "/scratchfs/atlas/chenhr/atlaswork/ATLAS_style/atlasrootstyle/AtlasLabels.C"
#include "/scratchfs/atlas/chenhr/atlaswork/ATLAS_style/atlasrootstyle/AtlasUtils.C"
#endif

void ratioPlot(){
   TString dirname = "templatePlots/";
   TString CR = "invID_invIso";
   //TString CR = "invID";
   TString name = dirname+CR+"_uncertainty";

  char *cf_cats = (char*)"cats.cfg";
  map<TString, string> catCuts;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;

   SetAtlasStyle();
   gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);

   map<TString, pair<float, float>> bins;
   bins["OOincl"] = make_pair(-999999999, 999999999);
   bins["b1"] = make_pair(-999999999, -2);
   bins["b2"] = make_pair(-2, -1);
   bins["b3"] = make_pair(-1, 0);
   bins["b4"] = make_pair(0, 1);
   bins["b5"] = make_pair(1, 2);
   bins["b6"] = make_pair(2, 99999999);

   TFile *f_in = new TFile("template.root","read");

   //for(auto bin : bins){
   for(auto cat : catCuts){

     //cout<<"========  "<<bin.first<<"  ========"<<endl;
     cout<<"========  "<<cat.first<<"  ========"<<endl;

     //TH1F *h1 = (TH1F*) f_in->Get("data_nom_"+bin.first);
     //TH1F *h2 = (TH1F*) f_in->Get("template_uncer_"+bin.first);
     //TH1F *h4 = (TH1F*) f_in->Get("yj_component_"+bin.first);
     TH1F *h1 = (TH1F*) f_in->Get("data_nom_"+cat.first);
     TH1F *h2 = (TH1F*) f_in->Get("template_uncer_"+cat.first);
     TH1F *h4 = (TH1F*) f_in->Get("yj_component_"+cat.first);

     TH1F *h2_noErr = (TH1F*) h2->Clone(Form("%s_noErr", h2->GetName()));
     for(int i=1; i<=h2_noErr->GetNbinsX(); i++){
       h2_noErr->SetBinError(i, 0.);
     }
  
     TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  
     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
     pad1->SetBottomMargin(0); // Upper and lower plot are joined
     pad1->SetGridx();         // Vertical grid
     pad1->Draw();             // Draw the upper pad: pad1
     pad1->cd();               // pad1 becomes the current pad
     h1->SetStats(0);          // No statistics on upper plot
     h1->SetMinimum(0.);
     h1->Draw("ep");               // Draw h1
     h2->Draw("same hist");         // Draw h2 on top of h1
  
     h4->SetLineWidth(0);
     h4->SetFillColor(TColor::GetColor(255,215,0));
     h4->Draw("same hist");
  
     ATLASLabel(0.22,0.90,"Internal");
     myText(0.22, 0.85, 1, "#sqrt{s}= 13 TeV, 139 fb^{-1}");
  
     h1->GetYaxis()->SetLabelSize(.04);
     TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
     axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
     axis->SetLabelSize(15);
     axis->Draw();
  
     TLegend* lg = new TLegend(0.6, 0.75, 0.95, 0.9);
     if(h1!=nullptr) lg->AddEntry(h1, "Data(2015-2018)", "lp");
     if(h2!=nullptr) lg->AddEntry(h2, "Background templ.", "f");
     if(h4!=nullptr) lg->AddEntry(h4, "Reducible(#gamma+jet)", "f");
     lg -> SetFillStyle(0);
     lg -> SetBorderSize(0);
     lg -> Draw("same");
  
     c->cd();          // Go back to the main canvas before defining pad2
     TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
     pad2->SetTopMargin(0);
     pad2->SetBottomMargin(0.2);
     pad2->SetGridx(); // vertical grid
     pad2->Draw();
     pad2->cd();       // pad2 becomes the current pad
  
     TH1F *h3 = (TH1F*)h1->Clone("h3");
     h3->SetLineColor(kBlack);
     h3->SetMinimum(0.5);  // Define Y ..
     h3->SetMaximum(1.5); // .. range
     h3->Sumw2();
     h3->SetStats(0);      // No statistics on lower plot
     h3->Divide(h2_noErr);
     h3->SetMarkerStyle(20);
     h3->Draw("ep");       // Draw the ratio plot
  
     TH1F *h5 = (TH1F*)h2->Clone("h5");
     h5->SetLineColor(kBlack);
     h5->Divide(h2_noErr);
     h5->SetMarkerStyle(20);
     h5->SetMarkerSize(0);
     h5->SetFillStyle(3001);
     h5->SetFillColor(kBlack);
     h5->Draw("same e2"); 
  
     if(h5!=nullptr) lg->AddEntry(h5, "Background unc.", "f");
  
     h1->SetLineColor(kBlack);
     h1->SetLineWidth(2);
     h1->GetYaxis()->SetTitleSize(20);
     h1->GetYaxis()->SetTitleFont(43);
     h1->GetYaxis()->SetTitleOffset(1.55);
  
     h3->SetTitle(""); // Remove the ratio title
  
     h3->GetYaxis()->SetTitle("#gamma#gamma fraction");
     h3->GetYaxis()->SetNdivisions(505);
     h3->GetYaxis()->SetTitleSize(20);
     h3->GetYaxis()->SetTitleFont(43);
     h3->GetYaxis()->SetTitleOffset(1.55);
     h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
     h3->GetYaxis()->SetLabelSize(15);
  
     h3->GetXaxis()->SetTitle("m_#gamma#gamma/GeV");
     h3->GetXaxis()->SetTitleSize(20);
     h3->GetXaxis()->SetTitleFont(43);
     h3->GetXaxis()->SetTitleOffset(4.);
     h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
     h3->GetXaxis()->SetLabelSize(15);
  
     //c->SaveAs(name+"_"+bin.first+".png");
     c->SaveAs(name+"_"+cat.first+".png");

     delete lg;
     delete axis;
     delete pad1;
     delete pad2;
     delete c;
  }
}
