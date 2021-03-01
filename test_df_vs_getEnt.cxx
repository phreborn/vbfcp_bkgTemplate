void test_df_vs_getEnt(){
  //TFile *f_data = new TFile("h026_data.root","read");
  TFile *f_yy = new TFile("h026_364352.diphoton_AF2_slim.root","read");

  TChain ch_yy("output", "output");
  ch_yy.Add("h026_364352.diphoton_AF2_slim.root");
  ROOT::RDataFrame df_yy(ch_yy, {"m_yy"});

  TTree * t_yy = (TTree*) f_yy->Get("output");

  float m_yy, oo1, weight;
  bool isPassed, passPID_y1, passPID_y2, passIso_y1, passIso_y2, passIso;

  t_yy->SetBranchAddress("m_yy", &m_yy);
  t_yy->SetBranchAddress("oo1", &oo1);
  t_yy->SetBranchAddress("wt", &weight);
  t_yy->SetBranchAddress("passPID_y1", &passPID_y1);
  t_yy->SetBranchAddress("passPID_y2", &passPID_y2);
  t_yy->SetBranchAddress("passIso_y1", &passIso_y1);
  t_yy->SetBranchAddress("passIso_y2", &passIso_y2);
  t_yy->SetBranchAddress("passIso", &passIso);
  t_yy->SetBranchAddress("isPassed", &isPassed);

  TH1F *th_yy = new TH1F("th_yy", "", 55, 105, 160);
  TH1F *dfh_yy = new TH1F("dfh_yy", "", 55, 105, 160);

  int n_yy = t_yy->GetEntries();
  for(int i = 0; i<n_yy; i++){
    t_yy->GetEntry(i);
    //if(oo1 < -999999999 || oo1 >= 999999999) continue;
    if(! (oo1>=-999999999 && oo1<999999999)) continue;
    //if(! (((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso)) continue;
    if(! ((!passPID_y1&&!passIso_y1&&passPID_y2&&passIso_y2)||(!passPID_y2&&!passIso_y2&&passPID_y1&&passIso_y1))) continue;
    th_yy->Fill(m_yy/1000, weight);
  }

  //df_yy.Filter("((!passPID_y1&&passPID_y2)||(passPID_y1&&!passPID_y2))&&passIso").Foreach([&dfh_yy](float m_yy, float wt){ dfh_yy->Fill(m_yy/1000, wt); }, {"m_yy", "wt"});
  df_yy.Filter("(oo1>=-999999999 && oo1<999999999) && ((!passPID_y1&&!passIso_y1 && passPID_y2&&passIso_y2) || (!passPID_y2&&!passIso_y2 && passPID_y1&&passIso_y1))").Foreach([&dfh_yy](float m_yy, float wt){ dfh_yy->Fill(m_yy/1000, wt); }, {"m_yy", "wt"});


  cout<<"tree: "<<th_yy->Integral()<<endl;
  cout<<"df: "<<dfh_yy->Integral()<<endl;

  th_yy->SetMarkerColor(kRed);
  dfh_yy->SetLineColor(kBlue);
  th_yy->Draw("e");
  dfh_yy->Draw("same hist");
}
