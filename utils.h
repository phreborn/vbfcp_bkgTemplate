#include <iterator>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <dirent.h>

#include <sstream>
#include <fstream>

#include <string.h>

using namespace std;

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

vector<string> readInLines(const char * cfgfilepath){
    vector<string> readin;

    fstream cfgFile;
    cfgFile.open(cfgfilepath);
    if( ! cfgFile.is_open())
    {   
        cout<<"can not open cfg file!"<<endl;
    }
    char tmp[1000];
    while(!cfgFile.eof())
    {   
        cfgFile.getline(tmp,1000);
        string line(tmp); 
        if(line.find("#") != string::npos) continue;
        if(line=="") continue;
        readin.push_back(line);
    }

    return readin;
}

map<TString, string> sepKeyValue(string cfg){
  map<TString, string> cats;

  vector<string> lines = readInLines(cfg.data());
  for(auto l : lines){
    int pos = l.find(":");
    string sKey = l.substr(0,pos);
    string sValues = l.substr(pos+1);
    TString tsKey = sKey.data();
    tsKey.ReplaceAll(" ", "");

    cats[tsKey] = sValues;
  }

  return cats;
}

void getCatCuts(string cfg, map<TString, string> &catCuts){
  catCuts = sepKeyValue(cfg);
}

