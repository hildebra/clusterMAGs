#pragma once
//#include "IO.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <cstring>
#include <map>
//#include <list>
#include <stdlib.h>
//#include <algorithm>
#include <math.h>
//#include <cmath>
#include <time.h>
//#include <random>
#include <assert.h>
#include <unordered_map>
//#include <numeric>
#include <future>
#include <mutex>
#include <chrono>
//#include <random>
#include "include/robin_hood.h"
#include <algorithm>    
#include <sys/stat.h>




#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#pragma warning(disable:4996)
#else
#define _gzipread
#endif
#define notRpackage


#ifdef _gzipread
#include "gzstream.h"
#endif
typedef double mat_fl;
typedef float smat_fl;

bool isGZfile(const std::string fi);

using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;



inline bool file_exists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}



struct options
{
public:
	options(int argc, char** argv);
//	options(std::string, std::string , int repeats, std::vector<double> depth, 
//		int NoOfMatrices, 
//		bool verbose, unsigned int threads);
	void print_details();
	int maxTier() { return  std::min((int)complTiers.size(), maxTierRds); }
	void show_quals();

	//~options();

	//vars
  std::string map = "";
  std::string outDir = "";
  std::string tmp  = "";
  std::string LCAdir = "";
  std::string canopyFile = "";
  std::string GTDBfile = "";
 // std::string taxMGSreport = "";
  string geneCatIdx = "";
  string path2Bins = "Binning/SB/";
  string path2Assmbl = "assemblies/metag/";
  string file2Assmbl;
  //std::string map = "";
  uint threads = 1;
  bool useCheckm1;
  bool useCheckm2;

  vector<float> complTiers; vector<float> contaTiers;
  vector<float> LCAcomplTier; 
  vector<int> allowNovelMGSTiers;
  vector<float> match2MGS_Tier; 
  vector<float> globalUniq4MGS_Tier; vector<float> maxMatch4NovelMGS_Tier;
  int maxTierRds;
  string MGSprefix;
  string outTag;
  string qualSuffix;
  bool processOnlyHQ;
  int minMGassigned2tax;//min num of MGs that have a tax assignment (not "?") to accept this tax for MAG
};




//testcalls:
//-CMsuffix .cm2 -FILEtag SBx -MGStag MM2 -geneCatIdx C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/compl.incompl.95.fna.clstr.idx -LCAdir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/MGs/ -outDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/out/ -map C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/map.0.txt,C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/map.1.txt -canopyDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/Cano/clusters.txt.filt -MGfile C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/GTDBmg.subset.cats
//-CMsuffix .cm2 -FILEtag DI -MGStag di. -geneCatIdx C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/compl.incompl.95.fna.clstr.idx -LCAdir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/MGs/ -outDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/out/ -map C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/map.0.txt -canopyDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/Cano/clusters.txt.filt -MGfile C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsDIME/GTDBmg.subset.cats
//-show_quals