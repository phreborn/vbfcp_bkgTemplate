#include "utils.h"

#ifdef __CLING__
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasLabels.C"
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasUtils.C"
#endif

void ratioPlot(TH1F *h1, TH1F *h2, TF1 *func, TString name){
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

   TH1F *h2_noErr = (TH1F*) h2->Clone(Form("%s_noErr", h2->GetName()));
   for(int i=1; i<=h2_noErr->GetNbinsX(); i++){
     h2_noErr->SetBinError(i, 0.);
   }

   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   h2->SetStats(0);          // No statistics on upper plot
   h2->SetMinimum(0.);
   h2->Draw("ep");               // Draw h1
   h1->Draw("same hist");         // Draw h2 on top of h1

   ATLASLabel(0.22,0.90,"Internal");
   myText(0.22, 0.85, 1, "#sqrt{s}= 13 TeV, 139 fb^{-1}");
   myText(0.22, 0.8, 1, "id:TightAntitight,iso:LooseAntiloose");

   h2->GetYaxis()->SetLabelSize(.04);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(15);
   axis->Draw();

   TLegend* lg = new TLegend(0.7, 0.75, 0.95, 0.9);
   if(h2!=nullptr) lg->AddEntry(h2, "Data", "l");
   if(h1!=nullptr) lg->AddEntry(h1, "MC #gamma#gamma", "f");
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
   h3->SetMinimum(0.);  // Define Y ..
   h3->SetMaximum(.1); // .. range
   //h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(h2_noErr);
   h3->SetMarkerStyle(20);
   h3->Draw("hist");       // Draw the ratio plot

   //h3->Fit(func, "r");
   //func->SetLineColor(6);
   //func->Draw("same");

   h2->SetLineColor(kBlack);
   h2->SetLineWidth(2);

   //h2->GetYaxis()->SetTitle("normalized to 1");
   h2->GetYaxis()->SetTitleSize(20);
   h2->GetYaxis()->SetTitleFont(43);
   h2->GetYaxis()->SetTitleOffset(1.55);

   //h1->SetLineColor(7);
   h1->SetLineWidth(0);
   //h1->SetFillColor(TColor::GetColor(255,215,0));
   h1->SetFillColor(kRed+2);

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

   c->SaveAs(name+".png");
   c->SaveAs(name+".pdf");

   delete axis;
   delete lg;
   delete pad1;
   delete pad2;
   delete c;
}

double polynomial2(double *x, double *par){
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

void invCut_hist(){
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  TString dirname = "subtractPlots";

  int nBins = 110;

  char *cf_cats = (char*)"cats.cfg";
  map<TString, string> catCuts;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;

  map<TString, pair<float, float>> bins;
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

  string config = "config";
  string blindCut = "";
  readConfigFile(config.data(), "blindSel", blindCut);
  string TTCut = "";
  readConfigFile(config.data(), "TTSel", TTCut);
  string invIDCut = "";
  readConfigFile(config.data(), "invID", invIDCut);
  string invID_invIsoCut = "";
  readConfigFile(config.data(), "invID_invIso", invID_invIsoCut);

  TChain ch_data("output", "output");
  TChain ch_yy("output", "output");

  ch_data.Add("h026_data_v20.root");
  ch_yy.Add("pre_slim/fullrun2/h026_364352.diphoton_v20_AF2_slim.root");

  ROOT::RDataFrame df_data(ch_data, {"m_yy"});
  ROOT::RDataFrame df_yy(ch_yy, {"m_yy"});

  //TFile *f_data = new TFile("/scratchfs/bes/chenhr/atlaswork/VBF_CP/ntuples/data17/data17_slim.root","read");
  //TFile *f_yy = new TFile("/scratchfs/bes/chenhr/atlaswork/VBF_CP/ntuples/mc16d/mc16d.364352.diphoton_AF2_slim.root","read");
  TFile *f_data = new TFile("h026_data_v20.root","read");
  TFile *f_yy = new TFile("pre_slim/fullrun2/h026_364352.diphoton_v20_AF2_slim.root","read");

  TTree *t_data = (TTree*) f_data->Get("output");
  TTree *t_yy = (TTree*) f_yy->Get("output");
  
  Int_t N_j_30,N_photon,cutflow,catCoup_XGBoost_ttH;
  Float_t m_yy,pT_y1,pT_y2,m_jj_30,DeltaEta_jj,Zepp,oo1,oo2,WeightDtilde1,WeightDtilde2,weight,xsec_kF_eff,total_weights, wt;
  bool isDalitz,passBasic,passPresel,passTrigMatch,passRelPt,passMass,passIso,passPID,isPassed,passPID_y1,passPID_y2,passIso_y1,passIso_y2;
  
  t_data->SetBranchAddress("passPID_y1", &passPID_y1);
  t_data->SetBranchAddress("passPID_y2", &passPID_y2);
  t_data->SetBranchAddress("passIso_y1", &passIso_y1);
  t_data->SetBranchAddress("passIso_y2", &passIso_y2);
  t_data->SetBranchAddress("passPID", &passPID);
  t_data->SetBranchAddress("passIso", &passIso);
  t_data->SetBranchAddress("isPassed", &isPassed);
  t_data->SetBranchAddress("catCoup_XGBoost_ttH", &catCoup_XGBoost_ttH);
  
  t_data->SetBranchAddress("m_yy", &m_yy);
  t_data->SetBranchAddress("oo1", &oo1);
  t_data->SetBranchAddress("oo2", &oo2);
  t_data->SetBranchAddress("wt", &weight);
 
  t_yy->SetBranchAddress("passPID_y1", &passPID_y1);
  t_yy->SetBranchAddress("passPID_y2", &passPID_y2);
  t_yy->SetBranchAddress("passIso_y1", &passIso_y1);
  t_yy->SetBranchAddress("passIso_y2", &passIso_y2);
  t_yy->SetBranchAddress("passPID", &passPID);
  t_yy->SetBranchAddress("passIso", &passIso);
  t_yy->SetBranchAddress("isPassed", &isPassed);
  t_yy->SetBranchAddress("catCoup_XGBoost_ttH", &catCoup_XGBoost_ttH);
 
  t_yy->SetBranchAddress("m_yy", &m_yy);
  t_yy->SetBranchAddress("oo1", &oo1);
  t_yy->SetBranchAddress("oo2", &oo2);
  t_yy->SetBranchAddress("wt", &weight);

  //for(auto bin : bins){
  for(auto cat : catCuts){

    //cout<<"========  "<<bin.first<<"  ========"<<endl;
    cout<<"========  "<<cat.first<<"  ========"<<endl;

    //float b_l = bin.second.first;
    //float b_r = bin.second.second;
  
    TH1F *h_data_nom = new TH1F("data_nom", "", nBins, 105,160);
    TH1F *h_data_invID = new TH1F("data_invID", "", nBins, 105,160);
    TH1F *h_data_invID_invIso = new TH1F("data_invID_invIso", "", nBins, 105,160);
  
    TH1F *h_yy_nom = new TH1F("yy_nom", "", nBins, 105,160);
    TH1F *h_yy_invID = new TH1F("yy_invID", "", nBins, 105,160);
    TH1F *h_yy_invID_invIso = new TH1F("yy_invID_invIso", "", nBins, 105,160);
  
    int n_data = t_data->GetEntries();

    string catCut = cat.second;

    string cat_TT = "";
    cat_TT = Form("%s && %s", catCut.data(), TTCut.data()); cout<<"cat_TT: "<<cat_TT<<endl;

    string cat_invID = "";
    cat_invID = Form("%s && %s", catCut.data(), invIDCut.data()); cout<<"cat_invID: "<<cat_invID<<endl;

    string cat_invID_invIso = "";
    cat_invID_invIso = Form("%s && %s", catCut.data(), invID_invIsoCut.data()); cout<<"cat_invID_invIso: "<<cat_invID_invIso<<endl;

    df_data.Filter(cat_TT).Filter(blindCut).Foreach([&h_data_nom](float m_yy){ h_data_nom->Fill(m_yy/1000); }, {"m_yy"});

    df_data.Filter(cat_invID).Foreach([&h_data_invID](float m_yy, float wRatio){ h_data_invID->Fill(m_yy/1000, wRatio); }, {"m_yy", "wRatio"});

    df_data.Filter(cat_invID_invIso).Foreach([&h_data_invID_invIso](float m_yy, float wRatio){ h_data_invID_invIso->Fill(m_yy/1000, wRatio); }, {"m_yy", "wRatio"});

    //for (int i = 0; i < n_data; i++){
    //  t_data->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if (! isPassed) continue;
    //  h_data_nom->Fill(m_yy/1000); 
    //}
  
    //for (int i = 0; i < n_data; i++){
    //  t_data->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if (! (((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso)) continue;
    //  h_data_invID->Fill(m_yy/1000); 
    //}
  
    //for (int i = 0; i < n_data; i++){
    //  t_data->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if (! ((!passPID_y1&&!passIso_y1&&passPID_y2&&passIso_y2)||(!passPID_y2&&!passIso_y2&&passPID_y1&&passIso_y1))) continue;
    //  h_data_invID_invIso->Fill(m_yy/1000);
    //}
  
    int n_yy = t_yy->GetEntries();

    df_yy.Filter(cat_TT).Foreach([&h_yy_nom](float m_yy, float wt, float wRatio){ h_yy_nom->Fill(m_yy/1000, wt*wRatio); }, {"m_yy", "wt", "wRatio"});

    df_yy.Filter(cat_invID).Foreach([&h_yy_invID](float m_yy, float wt, float wRatio){ h_yy_invID->Fill(m_yy/1000, wt*wRatio); }, {"m_yy", "wt", "wRatio"});

    df_yy.Filter(cat_invID_invIso).Foreach([&h_yy_invID_invIso](float m_yy, float wt, float wRatio){ h_yy_invID_invIso->Fill(m_yy/1000, wt*wRatio); }, {"m_yy", "wt", "wRatio"});

    //for(int i = 0; i < n_yy; i++){
    //  t_yy->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if(!isPassed) continue;
    //  h_yy_nom->Fill(m_yy/1000,weight);
    //}
    //
    //for(int i = 0; i < n_yy; i++){
    //  t_yy->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if(! (((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso)) continue;
    //  h_yy_invID->Fill(m_yy/1000,weight);
    //}
  
    //for(int i = 0; i < n_yy; i++){
    //  t_yy->GetEntry(i);
    //  if(oo1 < b_l || oo1 >= b_r) continue;
    //  if(! ((!passPID_y1&&!passIso_y1&&passPID_y2&&passIso_y2)||(!passPID_y2&&!passIso_y2&&passPID_y1&&passIso_y1))) continue;
    //  h_yy_invID_invIso->Fill(m_yy/1000,weight);
    //}
  
    //h_data_invID->Sumw2();
    //h_data_invID_invIso->Sumw2();
    //h_yy_invID->Sumw2();
    //h_yy_invID_invIso->Sumw2();
  
    TH1F *h_yj_invID = (TH1F*)h_data_invID->Clone("yj_invID");
    h_yj_invID->Add(h_yy_invID, -1);
    TH1F *h_yj_invID_invIso = (TH1F*)h_data_invID_invIso->Clone("yj_invID_invIso");
    h_yj_invID_invIso->Add(h_yy_invID_invIso, -1);
  
    cout<<"data nominal: "<<h_data_nom->Integral()<<endl;
    cout<<"data inverted ID: "<<h_data_invID->Integral()<<endl;
    cout<<"data inverted ID Iso: "<<h_data_invID_invIso->Integral()<<endl;
    cout<<"yy nominal: "<<h_yy_nom->Integral()<<endl;
    cout<<"yy inverted ID: "<<h_yy_invID->Integral()<<endl;
    cout<<"yy inverted ID Iso: "<<h_yy_invID_invIso->Integral()<<endl;
    cout<<"yj inverted ID: "<<h_yj_invID->Integral()<<endl;
    cout<<"yj inverted ID Iso: "<<h_yj_invID_invIso->Integral()<<endl;
  
    //h_yj_invID->Draw();
    //h_yj_invID_invIso->Draw();
  
    //TF1 *poly2 = new TF1("poly2", polynomial2, 105, 160, 3);
    TF1 *poly2 = new TF1("poly2", "pol1", 105, 160);
  
    //ratioPlot(h_yy_invID, h_data_invID, poly2, dirname+"/invID_"+bin.first);
    //ratioPlot(h_yy_invID_invIso, h_data_invID_invIso, poly2, dirname+"/invID_invIso_"+bin.first);
    ratioPlot(h_yy_invID, h_data_invID, poly2, dirname+"/invID_"+cat.first);
    ratioPlot(h_yy_invID_invIso, h_data_invID_invIso, poly2, dirname+"/invID_invIso_"+cat.first);
  
    delete h_data_nom;
    delete h_data_invID;
    delete h_data_invID_invIso;
    delete h_yy_nom;
    delete h_yy_invID;
    delete h_yy_invID_invIso;
    delete poly2;
  }
}
