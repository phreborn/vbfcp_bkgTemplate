#include <iterator>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <dirent.h>

#include <sstream>
#include <fstream>

using namespace std;

int fileNameFilter(const struct dirent *cur) {
    std::string str(cur->d_name);
    if (str.find("..") != std::string::npos) {
        return 0;
    }
    return 1;
}

std::vector<std::string> getDirBinsSortedPath(std::string dirPath) {
    struct dirent **namelist;
    std::vector<std::string> ret;
    int n = scandir(dirPath.c_str(), &namelist, fileNameFilter, alphasort);
    if (n < 0) {
        return ret;
    }
    for (int i = 0; i < n; ++i) {
        std::string filePath(namelist[i]->d_name);
        ret.push_back(filePath);
        free(namelist[i]);
    };
    free(namelist);
    return ret;
}

void getFilesList(std::vector<std::string> &files, TString dirpath, TString kw){
  std::string path_str = dirpath.Data();
  std::vector<std::string> sub_dirs = getDirBinsSortedPath(path_str);

  for(auto d : sub_dirs){
    if(d==".") continue;
    if(d.find(kw.Data()) == std::string::npos) continue;
    cout<<"d: "<<path_str+d<<endl;
    std::vector<std::string> fs = getDirBinsSortedPath(path_str+d+"/");
    for(auto f : fs){
      if(f==".") continue;
      cout<<"f: "<<path_str+"/"+d+"/"+f<<endl;
      files.push_back(path_str+d+"/"+f);

    }
  }
}

static bool readConfigFile(const char * cfgfilepath, const string & key, string & value)
{
    fstream cfgFile;
    cfgFile.open(cfgfilepath);
    if( ! cfgFile.is_open())
    {
        cout<<"can not open cfg file!"<<endl;
        return false;
    }
    char tmp[1000];
    while(!cfgFile.eof())
    {
        cfgFile.getline(tmp,1000);
        string line(tmp);
        size_t pos = line.find(':');
        if(pos==string::npos) continue;
        string tmpKey = line.substr(0,pos);
        if(key==tmpKey)
        {   
            value = line.substr(pos+1);
        }
    }
    return false;
}

TString getMCSampleName(int mcID){
  string name;
  readConfigFile("../MCSamples.config", Form("SampleName.%d", mcID), name);
  while(name.find(" ")!=std::string::npos) { name.replace(name.find(" "), 1, ""); }
  return name.data();
}

TH1F *getCutFlowHist(int mcID, TFile* file){
  TString suffix = "_noDalitz_weighted";
  TString cutFlowName = Form("CutFlow_%s%s", getMCSampleName(mcID).Data(), suffix.Data());
  TH1F *cutFlow = (TH1F*) file->Get(cutFlowName);
  return cutFlow;
}

double getSumOfWeights(int mcID, TFile* file){
  double NxAOD = getCutFlowHist(mcID, file)->GetBinContent(1);
  double NDxAOD = getCutFlowHist(mcID, file)->GetBinContent(2);
  double WDxAOD = getCutFlowHist(mcID, file)->GetBinContent(3);

  double weightSum = WDxAOD*NxAOD/NDxAOD;
  cout<<"xAOD, DxAOD, allEvt: "<<NxAOD<<", "<<NDxAOD<<", "<<WDxAOD<<endl;
  return weightSum;
}

void slim(int ini, int fin){
  map<TString,float> lumi;
  lumi["mc16a"] = 36207.66;
  lumi["mc16d"] = 44307.4;
  lumi["mc16e"] = 58450.1;

  TString camp = "mc16a";
  int id = 364352;

  TFile *fout = new TFile(Form(camp+".364352.diphoton_AF2_slim_%i_%i.root", ini, fin), "recreate");
  TTree *tout = new TTree("output","");

  Int_t N_j_30,N_photon,cutflow,catCoup_HybridICHEP,catCoup_XGBoost_ttH,catCoup_Moriond2017BDT;
  Float_t m_yy,pT_y1,pT_y2,m_jj_30,DeltaEta_jj,Zepp,oo1,oo2,WeightDtilde1,WeightDtilde2,weight,xsec_kF_eff,total_weights,wt,weightJvt_30,weightJvt_50,weightFJvt,weightFJvt_30,weight_catCoup_HybridICHEP,weight_catCoup_XGBoost_ttH,weight_catCoup_Moriond2017BDT;
  bool isDalitz,passBasic,passPresel,passTrigMatch,passRelPt,passMass,passIso,passPID,isPassed,passPID_y1,passIso_y1,passPID_y2,passIso_y2;

  Float_t wRatio;
  Float_t minDeltaR_y_j;
  tout->Branch("wRatio", &wRatio, "wRatio/F");
  tree->SetBranchAddress("minDeltaR_y_j", &minDeltaR_y_j);

  float BDTout_ggH, BDTout_yy; // to change to int
  int cat_BDTggH_BDTyy;
  tout->Branch("cat_BDTggH_BDTyy",&cat_BDTggH_BDTyy, "cat_BDTggH_BDTyy/I");
  tout->Branch("BDTout_ggH", &BDTout_ggH, "BDTout_ggH/F");
  tout->Branch("BDTout_yy", &BDTout_yy, "BDTout_yy/F");

  tout->Branch("m_yy", &m_yy, "m_yy/F");
  tout->Branch("catCoup_HybridICHEP",&catCoup_HybridICHEP, "catCoup_HybridICHEP/I");
  tout->Branch("catCoup_XGBoost_ttH",&catCoup_XGBoost_ttH, "catCoup_XGBoost_ttH/I");
  tout->Branch("catCoup_Moriond2017BDT",&catCoup_Moriond2017BDT, "catCoup_Moriond2017BDT/I");
  tout->Branch("passPID", &passPID, "passPID/O");
  tout->Branch("passIso", &passIso, "passIso/O");
  tout->Branch("isPassed", &isPassed, "isPassed/O");
  tout->Branch("passPID_y1", &passPID_y1, "passPID_y1/O");
  tout->Branch("passIso_y1", &passIso_y1, "passIso_y1/O");
  tout->Branch("passPID_y2", &passPID_y2, "passPID_y2/O");
  tout->Branch("passIso_y2", &passIso_y2, "passIso_y2/O");
  tout->Branch("oo1", &oo1, "oo1/F");
  tout->Branch("oo2", &oo2, "oo2/F");
  tout->Branch("WeightDtilde1", &WeightDtilde1, "WeightDtilde1/F");
  tout->Branch("WeightDtilde2", &WeightDtilde2, "WeightDtilde2/F");
  tout->Branch("weightJvt_30", &weightJvt_30, "weightJvt_30/F");
  tout->Branch("weightJvt_50", &weightJvt_50, "weightJvt_50/F");
  tout->Branch("weightFJvt", &weightFJvt, "weightFJvt/F");
  tout->Branch("weightFJvt_30", &weightFJvt_30, "weightFJvt_30/F");
  tout->Branch("weight_catCoup_HybridICHEP", &weight_catCoup_HybridICHEP, "weight_catCoup_HybridICHEP/F");
  tout->Branch("weight_catCoup_XGBoost_ttH", &weight_catCoup_XGBoost_ttH, "weight_catCoup_XGBoost_ttH/F");
  tout->Branch("weight_catCoup_Moriond2017BDT", &weight_catCoup_Moriond2017BDT, "weight_catCoup_Moriond2017BDT/F");
  tout->Branch("wt", &wt, "wt/F");

  TFile *f = new TFile("/eos/user/h/huirun/vbf_cp/h026/nominal/"+camp+"/VBFCP_v8_mc16a.diphoton_ext.root","read");

  double sumOfWeights = getSumOfWeights(id, f);

  TTree *tree = (TTree*) f->Get("output");

  tree->SetBranchAddress("Nominal.isPassedBasic", &passBasic);
  tree->SetBranchAddress("Nominal.isPassedTriggerMatch", &passTrigMatch);
  tree->SetBranchAddress("Nominal.isPassedRelPtCuts", &passRelPt);
  tree->SetBranchAddress("Nominal.isPassedMassCut", &passMass);
  tree->SetBranchAddress("Nominal.N_j_30", &N_j_30);
  tree->SetBranchAddress("Nominal.N_photon", &N_photon);
  tree->SetBranchAddress("Nominal.m_jj_30", &m_jj_30);
  tree->SetBranchAddress("Nominal.DeltaEta_jj", &DeltaEta_jj);
  tree->SetBranchAddress("Nominal.Zepp", &Zepp);
  tree->SetBranchAddress("isDalitz", &isDalitz);

  tree->SetBranchAddress("Nominal.weight", &weight);
  tree->SetBranchAddress("xsec_kF_eff", &xsec_kF_eff);

  tree->SetBranchAddress("Nominal.m_yy", &m_yy);
  tree->SetBranchAddress("Nominal.isPassedIsolation", &passIso);
  tree->SetBranchAddress("Nominal.isPassedPID", &passPID);
  tree->SetBranchAddress("Nominal.isPassed", &isPassed);
  tree->SetBranchAddress("Nominal.isPassedPID_y1", &passPID_y1);
  tree->SetBranchAddress("Nominal.isPassedIso_y1", &passIso_y1);
  tree->SetBranchAddress("Nominal.isPassedPID_y2", &passPID_y2);
  tree->SetBranchAddress("Nominal.isPassedIso_y2", &passIso_y2);
  tree->SetBranchAddress("Nominal.oo1", &oo1);
  tree->SetBranchAddress("Nominal.oo2", &oo2);
  tree->SetBranchAddress("WeightDtilde1", &WeightDtilde1);
  tree->SetBranchAddress("WeightDtilde2", &WeightDtilde2);
  tree->SetBranchAddress("Nominal.weightJvt_30", &weightJvt_30);
  tree->SetBranchAddress("Nominal.weightJvt_50", &weightJvt_50);
  tree->SetBranchAddress("Nominal.weightFJvt", &weightFJvt);
  tree->SetBranchAddress("Nominal.weightFJvt_30", &weightFJvt_30);
  tree->SetBranchAddress("Nominal.catCoup_HybridICHEP", &catCoup_HybridICHEP);
  tree->SetBranchAddress("Nominal.weight_catCoup_HybridICHEP", &weight_catCoup_HybridICHEP);
  tree->SetBranchAddress("Nominal.catCoup_XGBoost_ttH", &catCoup_XGBoost_ttH);
  tree->SetBranchAddress("Nominal.weight_catCoup_XGBoost_ttH", &weight_catCoup_XGBoost_ttH);
  tree->SetBranchAddress("Nominal.catCoup_Moriond2017BDT", &catCoup_Moriond2017BDT);
  tree->SetBranchAddress("Nominal.weight_catCoup_Moriond2017BDT", &weight_catCoup_Moriond2017BDT);
  tree->SetBranchAddress("Nominal.BDTout_ggH", &BDTout_ggH);
  tree->SetBranchAddress("Nominal.BDTout_yy", &BDTout_yy);
  tree->SetBranchAddress("Nominal.cat_BDTggH_BDTyy",&cat_BDTggH_BDTyy);

  int n_entries = tree->GetEntries();

  for (int i = 0; i < n_entries; i++){
    if(i < ini || i > fin) continue;
    if(i%1000000==0) std::cout<<i<<" events"<<std::endl;

    tree->GetEntry(i);
    if(isDalitz==1) continue;
    if(!passBasic||N_photon<2||!passTrigMatch||!passRelPt||!passMass) continue;
    if(N_j_30<2) continue;
    //if(m_jj_30/1000<400) continue;
    if(DeltaEta_jj>-2&&DeltaEta_jj<2) continue;
    if(Zepp>5||Zepp<-5) continue;
    //if(catCoup_XGBoost_ttH<11||catCoup_XGBoost_ttH>14) continue;
    //if(catCoup_Moriond2017BDT<11||catCoup_Moriond2017BDT>14) continue;

    wt = lumi[camp]*weight_catCoup_XGBoost_ttH*xsec_kF_eff/sumOfWeights;

    float ratio_Zepp = 0.807723 + 0.15129*Zepp - 0.0174526*Zepp*Zepp;
    float ratio_minDR = 1.21246 - 0.155163*minDeltaR_y_j;
    wRatio = 1./(ratio_Zepp*ratio_minDR);

    tout->Fill();
  }

  fout->cd();
  tout->Write();

}
