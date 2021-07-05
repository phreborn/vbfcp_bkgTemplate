#ifdef __CLING__
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasLabels.C"
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasUtils.C"
#endif

#include "utils.h"

bool getFrac(TString file, TString bin, double &frac, double &uncer){
  ifstream inf;
  inf.open(file.Data());

  TString bin_tmp;
  double frac_tmp;
  double uncer_tmp;

  while(1){
    inf>>bin_tmp>>frac_tmp>>uncer_tmp;
    if(!inf.good()) break;
    if(bin_tmp!=bin) continue;
    frac = frac_tmp;
    uncer = uncer_tmp;
    return 1;
  }
  return 0;
}

void ratioPlot(TH1F *h1, TH1F *h2, TF1 *func_nom, TF1 *func_alt, TH1F *h4, TString name){
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
   h1->SetStats(0);          // No statistics on upper plot
   h1->Draw();               // Draw h1
   h4->Draw("same");
   h2->Draw("same hist");         // Draw h2 on top of h1

   ATLASLabel(0.22,0.85,"Internal");
   myText(0.22, 0.8, 1, "#gamma+jet background");

   h1->GetYaxis()->SetLabelSize(0.03);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(15);
   axis->Draw();

   TLegend* lg = new TLegend(0.65, 0.75, 0.95, 0.9);
   if(h1!=nullptr) lg->AddEntry(h1, "Data #gammaj(inv. id)", "lp");
   if(h4!=nullptr) lg->AddEntry(h4, "Data #gammaj(inv. id, inv. iso)", "lp");
   if(h2!=nullptr) lg->AddEntry(h2, "MC #gamma#gamma(nominal sel.)", "l");
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
   h3->SetMinimum(-1.5);  // Define Y ..
   h3->SetMaximum(5.5); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(h2_noErr);
   h3->SetMarkerStyle(20);
   h3->Draw("ep");       // Draw the ratio plot

   //h3->Fit(func_nom, "r");
   //func_nom->SetLineColor(6);
   //func_nom->Draw("same");

   //h3->Fit(func_alt, "r");
   //func_alt->SetLineColor(6-1);
   //func_alt->SetLineStyle(kDashed);
   //func_alt->Draw("same");

   TH1F *h5 = (TH1F*)h4->Clone("h5");
   h5->SetLineColor(kBlack);
   h5->Sumw2();
   h5->Divide(h2_noErr);
   h5->SetMarkerStyle(24);
   h5->Draw("same");

   TFile *ft = new TFile(name+".root", "recreate");
   ft->cd();
   h5->Write();

   h5->Fit(func_nom, "rWW");
   func_nom->SetLineColor(6);
   func_nom->Draw("same");

   h5->Fit(func_alt, "rWW");
   func_alt->SetLineColor(6-1);
   func_alt->SetLineStyle(kDashed);
   func_alt->Draw("same");

   h1->SetLineColor(kBlack);
   h1->SetLineWidth(2);

   h1->GetYaxis()->SetTitle("a.u.");
   h1->GetYaxis()->SetTitleSize(20);
   h1->GetYaxis()->SetTitleFont(43);
   h1->GetYaxis()->SetTitleOffset(1.55);

   h4->SetLineColor(kBlack);
   h4->SetLineWidth(2);
   h4->SetMarkerStyle(24);

   h2->SetLineColor(7);
   h2->SetLineWidth(2);

   h3->SetTitle(""); // Remove the ratio title

   h3->GetYaxis()->SetTitle("ratio to #gamma#gamma ");
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

   delete lg;
   delete axis;
   delete pad1;
   delete pad2;
   delete c;
}

void multiplyRatio(TH1F *hist, TH1F *h_ratio, int nBins){
  for(int i = 1; i <= nBins; i++){
    double bc = hist->GetBinContent(i);
    double be = hist->GetBinError(i);
    double ratio = h_ratio->GetBinContent(i);

    bc *= ratio;
    hist->SetBinContent(i, bc);
    hist->SetBinError(i, 0.);
  }
}

double polynomial2(double *x, double *par){
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

void bkgTemplate(){
  SetAtlasStyle();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  TString dirname = "smoothPlots/";
  TString fracfile = "fraction.csv";

  int nBins = 55;
  int binToMerge = 5;

  char *cf_cats = (char*)"cats.cfg";
  map<TString, string> catCuts;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;

  map<TString, pair<float, float>> bins;
  bins["OOincl"] = make_pair(-999999999, 999999999);
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

  TFile *f_out = new TFile("template.root", "recreate");

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

  ch_data.Add("h026_data.root");
  ch_yy.Add("h026_364352.diphoton_AF2_slim.root");

  ROOT::RDataFrame df_data(ch_data, {"m_yy"});
  ROOT::RDataFrame df_yy(ch_yy, {"m_yy"});

  //TFile *f_data = new TFile("/scratchfs/bes/chenhr/atlaswork/VBF_CP/ntuples/data17/data17_slim.root","read");
  //TFile *f_yy = new TFile("/scratchfs/bes/chenhr/atlaswork/VBF_CP/ntuples/mc16d/mc16d.364352.diphoton_AF2_slim.root","read");
  TFile *f_data = new TFile("h026_data.root","read");
  TFile *f_yy = new TFile("h026_364352.diphoton_AF2_slim.root","read");

  TTree *t_data = (TTree*) f_data->Get("output");
  TTree *t_yy = (TTree*) f_yy->Get("output");

  Int_t N_j_30,N_photon,cutflow,catCoup_XGBoost_ttH;
  Float_t m_yy,pT_y1,pT_y2,m_jj_30,DeltaEta_jj,Zepp,oo1,oo2,WeightDtilde1,WeightDtilde2,weight,xsec_kF_eff,total_weights;
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
    float b_l = bins[cat.first].first;
    float b_r = bins[cat.first].second;

    //TH1F *h_data_nom = new TH1F("data_nom_"+bin.first, "", nBins, 105,160);
    TH1F *h_data_nom = new TH1F("data_nom_"+cat.first, "", nBins, 105,160);
    TH1F *h_data_invID = new TH1F("data_invID", "", nBins, 105,160);
    TH1F *h_data_invID_invIso = new TH1F("data_invID_invIso", "", nBins, 105,160);
  
    TH1F *h_yy_nom = new TH1F("yy_nom", "", nBins, 105,160);
    TH1F *h_yy_invID = new TH1F("yy_invID", "", nBins, 105,160);
    TH1F *h_yy_invID_invIso = new TH1F("yy_invID_invIso", "", nBins, 105,160);

    TH1F *th_data_nom = new TH1F("tdata_nom_"+cat.first, "", nBins, 105,160);
    TH1F *th_data_invID = new TH1F("tdata_invID", "", nBins, 105,160);
    TH1F *th_data_invID_invIso = new TH1F("tdata_invID_invIso", "", nBins, 105,160);
    TH1F *th_yy_nom = new TH1F("tyy_nom", "", nBins, 105,160);
    TH1F *th_yy_invID = new TH1F("tyy_invID", "", nBins, 105,160);
    TH1F *th_yy_invID_invIso = new TH1F("tyy_invID_invIso", "", nBins, 105,160);
  
    // filling data
    int n_data = t_data->GetEntries();

    //string catCut = catCuts[bin.first];
    string catCut = cat.second;

    string cat_TT = "";
    cat_TT = Form("%s && %s", catCut.data(), TTCut.data()); cout<<"cat_TT: "<<cat_TT<<endl;

    string cat_invID = "";
    cat_invID = Form("%s && %s", catCut.data(), invIDCut.data()); cout<<"cat_invID: "<<cat_invID<<endl;

    string cat_invID_invIso = "";
    cat_invID_invIso = Form("%s && %s", catCut.data(), invID_invIsoCut.data()); cout<<"cat_invID_invIso: "<<cat_invID_invIso<<endl;

    df_data.Filter(cat_TT).Filter(blindCut).Foreach([&h_data_nom](float m_yy){ h_data_nom->Fill(m_yy/1000); }, {"m_yy"});
    TH1F * oh_data_nom = (TH1F*) h_data_nom->Clone("odata_nom");

    df_data.Filter(cat_invID).Foreach([&h_data_invID](float m_yy){ h_data_invID->Fill(m_yy/1000); }, {"m_yy"});
    TH1F * oh_data_invID = (TH1F*) h_data_invID->Clone("odata_invID");

    df_data.Filter(cat_invID_invIso).Foreach([&h_data_invID_invIso](float m_yy){ h_data_invID_invIso->Fill(m_yy/1000); }, {"m_yy"});
    TH1F * oh_data_invID_invIso = (TH1F*) h_data_invID_invIso->Clone("odata_invID_invIso");
  
    for (int i = 0; i < n_data; i++){
      t_data->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if (! isPassed) continue;
      if (m_yy/1000>120&&m_yy/1000<130) continue;
      //h_data_nom->Fill(m_yy/1000); 
      th_data_nom->Fill(m_yy/1000); 
    }

    for (int i = 0; i < n_data; i++){
      t_data->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if (! (((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso)) continue;
      //h_data_invID->Fill(m_yy/1000); 
      th_data_invID->Fill(m_yy/1000); 
    }
  
    for (int i = 0; i < n_data; i++){
      t_data->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if (! ((!passPID_y1&&!passIso_y1&&passPID_y2&&passIso_y2)||(!passPID_y2&&!passIso_y2&&passPID_y1&&passIso_y1))) continue;
      //h_data_invID_invIso->Fill(m_yy/1000);
      th_data_invID_invIso->Fill(m_yy/1000);
    }
  
    // filling MC yy
    int n_yy = t_yy->GetEntries();

    df_yy.Filter(cat_TT).Foreach([&h_yy_nom](float m_yy, float wt){ h_yy_nom->Fill(m_yy/1000, wt); }, {"m_yy", "wt"});
    TH1F * oh_yy_nom = (TH1F*) h_yy_nom->Clone("oyy_nom");

    df_yy.Filter(cat_invID).Foreach([&h_yy_invID](float m_yy, float wt){ h_yy_invID->Fill(m_yy/1000, wt); }, {"m_yy", "wt"});
    TH1F * oh_yy_invID = (TH1F*) h_yy_invID->Clone("oyy_invID");

    df_yy.Filter(cat_invID_invIso).Foreach([&h_yy_invID_invIso](float m_yy, float wt){ h_yy_invID_invIso->Fill(m_yy/1000, wt); }, {"m_yy", "wt"}); cout<<"df h_yy_invID_invIso: "<<h_yy_invID_invIso->Integral()<<endl;
    TH1F * oh_yy_invID_invIso = (TH1F*) h_yy_invID_invIso->Clone("oyy_invID_invIso");
  
    for(int i = 0; i < n_yy; i++){
      t_yy->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if(!isPassed) continue;
      //h_yy_nom->Fill(m_yy/1000,weight);
      th_yy_nom->Fill(m_yy/1000,weight);
    }
    
    for(int i = 0; i < n_yy; i++){
      t_yy->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if(! (((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso)) continue;
      //h_yy_invID->Fill(m_yy/1000,weight);
      th_yy_invID->Fill(m_yy/1000,weight);
    }
  
    for(int i = 0; i < n_yy; i++){
      t_yy->GetEntry(i);
      if(oo1 < b_l || oo1 >= b_r) continue;
      if(! ((!passPID_y1&&!passIso_y1&&passPID_y2&&passIso_y2)||(!passPID_y2&&!passIso_y2&&passPID_y1&&passIso_y1))) continue;
      //h_yy_invID_invIso->Fill(m_yy/1000,weight);
      th_yy_invID_invIso->Fill(m_yy/1000,weight);
    }
  
    //h_data_invID->Sumw2();
    //h_data_invID_invIso->Sumw2();
    //h_yy_invID->Sumw2();
    //h_yy_invID_invIso->Sumw2();
  
  
    // yj = data - yy
    TH1F *h_yj_invID = (TH1F*) h_data_invID->Clone("yj_invID");
    h_yj_invID->Add(h_yy_invID, -1);
    TH1F *h_yj_invID_invIso = (TH1F*) h_data_invID_invIso->Clone("yj_invID_invIso");
    h_yj_invID_invIso->Add(h_yy_invID_invIso, -1);
  
    cout<<"sideband data nominal: "<<h_data_nom->Integral()<<endl;
    cout<<"data inverted ID: "<<h_data_invID->Integral()<<endl;
    cout<<"data inverted ID Iso: "<<h_data_invID_invIso->Integral()<<endl;
    cout<<"yy nominal: "<<h_yy_nom->Integral()<<endl;
    cout<<"yy inverted ID: "<<h_yy_invID->Integral()<<endl;
    cout<<"yy inverted ID Iso: "<<h_yy_invID_invIso->Integral()<<endl;
    cout<<"yj inverted ID: "<<h_yj_invID->Integral()<<endl;
    cout<<"yj inverted ID Iso: "<<h_yj_invID_invIso->Integral()<<endl;

    TH1F *cl_yy_nom = (TH1F*)h_yy_nom->Clone(); // save for stat uncertainty before normalization

    h_yy_nom->Scale(1./h_yy_nom->Integral()); // yy shape, nom. sel.
    h_yj_invID->Scale(1./h_yj_invID->Integral()); // yj shape, inv. id.
    h_yj_invID_invIso->Scale(1./h_yj_invID_invIso->Integral()); // yj shape, inv. id. and iso.
  
    TF1 *poly1 = new TF1("poly1", "pol1", 105, 160);
    //TF1 *poly2 = new TF1("poly2", polynomial2, 105, 160, 3);
    TF1 *poly2 = new TF1("poly2", "pol2", 105, 160);

    TH1F *h_yj_invID_clone = (TH1F*) h_yj_invID->Clone("yj_invID_clone");
    TH1F *h_yy_nom_clone = (TH1F*) h_yy_nom->Clone("yy_nom_clone");
    TH1F *h_yj_invID_invIso_clone = (TH1F*) h_yj_invID_invIso->Clone("yj_invID_invIso_clone");

    h_yj_invID_clone->Rebin(binToMerge);
    h_yy_nom_clone->Rebin(binToMerge);
    h_yj_invID_invIso_clone->Rebin(binToMerge);

    //ratioPlot(h_yj_invID_clone, h_yy_nom_clone, poly2, poly1, h_yj_invID_invIso_clone, dirname+"invID_invIso_smooth_ratio_"+bin.first); // draw and fit
    ratioPlot(h_yj_invID_clone, h_yy_nom_clone, poly2, poly1, h_yj_invID_invIso_clone, dirname+"invID_invIso_smooth_ratio_"+cat.first); // draw and fit
    //ratioPlot(h_yj_invID_clone, h_yy_nom_clone, poly2, poly1, h_yj_invID_invIso_clone, dirname+"invID_smooth_ratio_"+bin.first); // draw and fit
//    ratioPlot(h_yj_invID_invIso_clone, h_yy_nom_clone, poly2, poly1, h_yj_invID_clone, dirname+"invID_invIso_smooth_ratio_"+bin.first); // draw and fit
  
    double p0, p1, p2;
    p0 = poly2->GetParameter(0); //cout<<"poly2 p0:"<<p0<<endl;
    p1 = poly2->GetParameter(1); //cout<<"poly2 p1:"<<p1<<endl;
    p2 = poly2->GetParameter(2); //cout<<"poly2 p2:"<<p2<<endl;
  
    TH1F *h_ratio = new TH1F("ratio", "", nBins, 105, 160); // nominal 2nd polynomial ratio
  
    double interval = (160.-105)/nBins; // !!! no all integer in division
    for(int i = 0; i < nBins; i++){
      double x = 105+i*interval+interval/2; //cout<<"x: "<<x<<endl;
      double val = p0 + p1*x +p2*x*x; //cout<<"val: "<<val<<endl;
      h_ratio->Fill(x, val);
    }

    p0 = poly1->GetParameter(0); //cout<<"poly1 p0: "<<p0<<endl;
    p1 = poly1->GetParameter(1); //cout<<"poly1 p1: "<<p1<<endl;
  
    TH1F *h_ratio_pol1 = new TH1F("ratio_pol1", "", nBins, 105, 160); // linear ratio to derive uncertainty
  
    for(int i = 0; i < nBins; i++){
      double x = 105+i*interval+interval/2; //cout<<"x: "<<x<<endl;
      double val = p0 + p1*x; //cout<<"val: "<<val<<endl;
      h_ratio_pol1->Fill(x, val);
    }

    TH1F *h_yj_reweight = (TH1F*) h_yy_nom->Clone("yj_reweight");
    h_yj_reweight->Multiply(h_ratio);
    h_yj_reweight->Scale(1./h_yj_reweight->Integral()); // reweighted yj shape, 2nd polynominal ratio
  
    TH1F *h_yj_reweight_pol1 = (TH1F*) h_yy_nom->Clone("yj_reweight_pol1");
    h_yj_reweight_pol1->Multiply(h_ratio_pol1);
    h_yj_reweight_pol1->Scale(1./h_yj_reweight_pol1->Integral()); // reweighted yj shape, linear ratio
  
    //double frac_yy = 0.78875;
    //double uncer = 0.05*frac_yy;
    double frac_yy = -1.;
    double uncer = -1.;

    //if(!getFrac(fracfile, bin.first, frac_yy, uncer)){
    if(!getFrac(fracfile, cat.first, frac_yy, uncer)){
      //cout<<"WARNING!!! can't get fraction of bin: "<<bin.first<<endl;
      cout<<"WARNING!!! can't get fraction of cat: "<<cat.first<<endl;
    }

    cout<<"yy fraction from 2x2D sideband: "<<frac_yy<<" +- "<<uncer<<endl;

    double frac_yy_u = frac_yy+uncer;
    double frac_yy_d = frac_yy-uncer;
  
    //double frac_yj_yy = 0.26783;
    double frac_yj_yy = (1-frac_yy)/frac_yy; // frac_yj = 1 - frac_yy
    double frac_yj_yy_u = (1-frac_yy_u)/frac_yy_u;
    double frac_yj_yy_d = (1-frac_yy_d)/frac_yy_d;
  
    double n_SB_data = h_data_nom->Integral(); // number of sideband data
  
    //TH1F *h_template = (TH1F*)h_yy_nom->Clone("bkg_template_"+bin.first);
    TH1F *h_template = (TH1F*)h_yy_nom->Clone("bkg_template_"+cat.first);
    h_template->Add(h_yj_reweight, frac_yj_yy); // 1*yy_shape + (frac_yj/frac_yy)*yj_shape
    h_template->Scale(1./h_template->Integral()); // template shape

    double frac_SB_yy = 1-h_template->Integral(h_template->GetXaxis()->FindBin(120),h_template->GetXaxis()->FindBin(130));cout<<"nom/2nd frac SB yy:"<<frac_SB_yy<<endl;
    h_template->Scale(n_SB_data/frac_SB_yy); // template normalized to sideband data
  
    //TH1F *h_template_pol1 = (TH1F*)h_yy_nom->Clone("bkg_template_pol1_"+bin.first);
    TH1F *h_template_pol1 = (TH1F*)h_yy_nom->Clone("bkg_template_pol1_"+cat.first);
    h_template_pol1->Add(h_yj_reweight_pol1, frac_yj_yy);
    h_template_pol1->Scale(1./h_template_pol1->Integral()); // smoothing func variated template
  
    double frac_SB_yy_pol1 = 1-h_template_pol1->Integral(h_template_pol1->GetXaxis()->FindBin(120),h_template_pol1->GetXaxis()->FindBin(130));cout<<"poly1 frac SB yy:"<<frac_SB_yy_pol1<<endl;
    h_template_pol1->Scale(n_SB_data/frac_SB_yy_pol1);
  
    //TH1F *h_template_u = (TH1F*)h_yy_nom->Clone("bkg_template_u_"+bin.first);
    TH1F *h_template_u = (TH1F*)h_yy_nom->Clone("bkg_template_u_"+cat.first);
    h_template_u->Add(h_yj_reweight, frac_yj_yy_u);
    h_template_u->Scale(1./h_template_u->Integral()); // up yy fraction variated template
  
    double frac_SB_yy_u = 1-h_template_u->Integral(h_template_u->GetXaxis()->FindBin(120),h_template_u->GetXaxis()->FindBin(130));cout<<"up frac SB yy:"<<frac_SB_yy_u<<endl;
    h_template_u->Scale(n_SB_data/frac_SB_yy_u);
  
    //TH1F *h_template_d = (TH1F*)h_yy_nom->Clone("bkg_template_d_"+bin.first);
    TH1F *h_template_d = (TH1F*)h_yy_nom->Clone("bkg_template_d_"+cat.first);
    h_template_d->Add(h_yj_reweight, frac_yj_yy_d);
    h_template_d->Scale(1./h_template_d->Integral()); // down yy fraction variated template
  
    double frac_SB_yy_d = 1-h_template_d->Integral(h_template_d->GetXaxis()->FindBin(120),h_template_d->GetXaxis()->FindBin(130));cout<<"down frac SB yy:"<<frac_SB_yy_d<<endl;
    h_template_d->Scale(n_SB_data/frac_SB_yy_d);

    // set stat uncertainty for template
    for (int i = 1; i <= nBins; i++){
      h_template->SetBinError(i, cl_yy_nom->GetBinError(i)*h_template->GetBinContent(i)/cl_yy_nom->GetBinContent(i));
    }
  
    // all uncertainty: stat, frac_yy, smoothing
    //TH1F *h_uncer = (TH1F*) h_template->Clone("template_uncer_"+bin.first);
    TH1F *h_uncer = (TH1F*) h_template->Clone("template_uncer_"+cat.first);
  
    for(int i = 1; i < nBins+1; i++){
      double uncer_stat = h_template->GetBinError(i); //cout<<"uncer stat: "<<uncer_stat<<endl;

      double up_frac = std::abs(h_template_u->GetBinContent(i)-h_template->GetBinContent(i));
      double down_frac = std::abs(h_template_d->GetBinContent(i)-h_template->GetBinContent(i));
      double uncer_frac = up_frac > down_frac? up_frac : down_frac; //cout<<"uncer fraction: "<<uncer_frac<<endl;
  
      double uncer_smooth = std::abs(h_template_pol1->GetBinContent(i)-h_template->GetBinContent(i)); //cout<<"uncer smoothing: "<<uncer_smooth<<endl;
  
      double uncer = std::sqrt(uncer_stat*uncer_stat+uncer_frac*uncer_frac+uncer_smooth*uncer_smooth);
  
      h_uncer->SetBinError(i, uncer);
      //h_uncer->SetBinError(i, h_template->GetBinContent(i)*0.3);
    }
  
    //TH1F *h_yj = (TH1F*) h_yj_reweight->Clone("yj_component_"+bin.first);
    TH1F *h_yj = (TH1F*) h_yj_reweight->Clone("yj_component_"+cat.first);
    h_yj->Scale((1-frac_yy)*n_SB_data/frac_SB_yy); // frac_yj*h_template->Integral()
  
    h_data_nom->Draw("ep");
    h_template_pol1->SetLineColor(kYellow);
    h_template_pol1->SetLineWidth(1);
    h_template_pol1->Draw("same hist");
    h_template->SetLineColor(kBlack);
    h_template->SetLineWidth(1);
    //h_template_u->SetLineColor(kBlue);
    //h_template_u->SetLineWidth(1);
    //h_template_d->SetLineColor(kRed);
    //h_template_d->SetLineWidth(1);
    h_template->Draw("same hist");
    //h_template_u->Draw("same hist");
    //h_template_d->Draw("same hist");
    h_uncer->SetMarkerSize(0);
    h_uncer->SetFillStyle(3254);
    h_uncer->SetFillColor(kBlack);
    h_uncer->Draw("same e2");
    h_yj->Draw("same hist");
  
    f_out->cd();
    h_data_nom->Sumw2();
    h_data_nom->Write();
    h_template->Write();
    h_template_pol1->Write();
    h_template_u->Write();
    h_template_d->Write();
    h_uncer->Write();
    h_yj->Write();

    TFile * tf = new TFile("testfile_"+cat.first+".root", "recreate");
    tf->cd();
    //h_yy_invID->Write();
    //h_yy_invID_invIso->Write();
    oh_data_nom->Write();
    oh_data_invID->Write();
    oh_data_invID_invIso->Write();
    oh_yy_nom->Write();
    oh_yy_invID->Write();
    oh_yy_invID_invIso->Write();

    th_data_nom->Write();
    th_data_invID->Write();
    th_data_invID_invIso->Write();
    th_yy_nom->Write();
    th_yy_invID->Write();
    th_yy_invID_invIso->Write();

    delete h_data_nom;
    delete h_data_invID;
    delete h_data_invID_invIso;
    delete h_yy_nom;
    delete h_yy_invID;
    delete h_yy_invID_invIso;
    delete poly1;
    delete poly2;
    delete h_ratio;
    delete h_ratio_pol1;
  }
}
