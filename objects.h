#pragma once
#include "options.h"
#include <filesystem>
#include <regex>



//typedef robin_hood::unordered_map <uint, uint> assoMap;
typedef unsigned long genetype;
typedef robin_hood::unordered_map <genetype, string> gene2stringM;
typedef robin_hood::unordered_map <genetype, vector<string>> gene2strvecM;

typedef robin_hood::unordered_map < genetype, int> geneCnter;
typedef robin_hood::unordered_map < genetype, float> geneCnterFloat;

namespace fs = std::filesystem;



typedef std::unordered_map<string, string > string2string;
typedef std::unordered_map<string, int> SmplOccur;
typedef std::unordered_map<string, int >::iterator SmplOccurIT;
typedef std::unordered_map<string, vector<int> > SmplOccurMult;


vector<pair<genetype, size_t>> srtGeneCnt(geneCnter& in);

void printMatrixNice(vector<vector<int>> mat, vector<float> thrs, string Rlab = "");
void printMatrixNice(vector<vector<int>> mat, vector<float> thrs, vector<string> Rlab, size_t start=0);
void printVectorNice(vector<int> vec, vector<float> thrs, string descr="");
vector<int> createHist(vector<float>uniqq, vector<float> histThr);



class MFmap {
public:
	MFmap(string mapF);
	~MFmap() {}
	vector<string> getAllSmpls(bool onlyAssmbled, bool assmGrpID=false);
	string getAssmblGrp(uint i) { return assmblGr[i]; }
	string getMapGrp(uint i) { return mapGr[i]; }
	size_t AssmblGrpSize(string s) { return (CntAssGrps)[s].size(); }
	string SmplPath(uint i) { return smplLoc[i]; }
	string SmplIds(uint i) { return smplIDs[i]; }
	string SmplIdsAlt(uint i) { return smplAltId[smplIDs[i]]; }
	size_t totalSmpls() { return smplIDs.size(); }
	SmplOccurMult& getAssmGrpCnts() {return CntAssGrps;}
private:
	void read_map(const string mapF, string& baseP);

	string2string smplRid;//reverse track from XXM7 to XX ID
	string2string smplAltId;//track from XX to XXM7 ID
	SmplOccurMult CntMapGrps;
	SmplOccurMult smpls;
	vector<string> mapGr;
	vector<string> assmblGr;
	vector<string> smplLoc; // only needed in readmap
	vector<string> smplIDs;

	SmplOccurMult CntAssGrps;

	size_t smplN;
	int curr;
	bool oldFolderStructure;
};

class MG_LCA {
public:
	MG_LCA(string inD);
	void addMGs(string);
	~MG_LCA() {}
	void report();
	bool hasTax(genetype,int lvl=6);
	bool getTax(genetype,vector<string>& out, string& cls);
	const vector<string>& getCats(bool srt) { if (srt) { sort(geneCats.begin(), geneCats.end()); }return geneCats; }
	string getGeneCat(genetype x) {
		auto xx = gene2Cat.find(x); if (xx == gene2Cat.end()) { return""; }
		else { return xx->second; }
	}
	size_t maxTaxLvl() { return TaxLvls.size(); }
	vector<string> getTaxLvls() { return TaxLvls; }
private:
	//functions
	bool addLCA(string in, string cat);

	//objects
	gene2strvecM gene2Tax;
	gene2stringM gene2Cat;
	vector<string> geneCats; //stores the discovered gene categories
	vector<string> TaxLvls;//header of tax names (eg 7 entries)
};

class MGS;
//complex map that stores: [contig] [num_on_ctg, geneID_in_GC]
typedef robin_hood::unordered_map<string, vector<pair<int,genetype>>> ctg2geneLST;

class MAG {
public:
	MAG(string bin, string smpl, vector<string>& ctgIn ):
		TaxAssignments(0), LCAscores(0), majorityTax(0),
		contigs(ctgIn), markerGenes(0),
		Completeness(0.f), Contamination(100.f),  GCcontent(0.f), codingDens(0.f),//LCAcompl(-1.f),
		Ngenes(0), GenomeSize(0), ctg_N50(0), centreScore(0.f), compndScore(0.f),
		binName(bin), sampleID(smpl), foundGenes(0), foundMGs(0),
		selfTier(-1), mgsAssign(nullptr), mgsRank(1000), isCanopy(false),
		bestMatchFrac(-1.f), uniqFrac(-1.f), bestMGS("")
	{}
	~MAG() {}
	string getUniqName() { return sampleID+"__"+ binName; }
	size_t size() { return contigs.size(); }
	const string& getSmplID() {return sampleID;}
	void setGC(float x) { GCcontent = x; }
	void setN50(int x) { ctg_N50 = x; }
	void setCompl(float x) { Completeness = x; }
	void setConta(float x) { Contamination = x; }
	void setNgenes(int x) { Ngenes = x; }
	void setGenomeSize(int x) { GenomeSize = x; }
	void setcodingDens(float x) { codingDens = x; }
	vector<string>& getContigs() {return contigs;}
	void addGene(genetype, string, int);
	//data exchange
	void setGenes(vector<genetype> x);
	void setCanopy(bool b) { isCanopy = b; }
	void clearCtg() { contigs.clear(); }
	robin_hood::unordered_map <genetype, bool>& getGenes() { return genes; }
	MGS* getAssignedMGS() { return mgsAssign; }
	string getAssignedMGSname();
	void setAssignedMGS(MGS* t) { mgsAssign = t; }
	geneCnter& getMGgenes() {return markerGenes;}
	string formatMAG(vector<string>&, bool,bool);
	//quality of MAG
	float getConta() { return Contamination; }
	float getComple() { return Completeness; }
	int getN50() { return ctg_N50; }
	float getGC() { return GCcontent; }
	float getCodingDens() { return codingDens; }
	void setCentreScore(float f) { centreScore = f; }
	float getCentreScore() { return centreScore; }
	void assignTier(uint T) { selfTier = T; }
	void assignTier(options* op);
	int qualTier() {return selfTier;}
	int NContigs() { return contigs.size(); }
	int NGenes() { return genes.size(); }
	float compoundScore(int maxGenes, int maxN50);
	float getCompoundScore() { return compndScore; }
	void setMgsRank(uint x) { mgsRank = x; }
	void storeBestHit(float bestMtchFrac, MGS* matchMGS, float uniq);
	//marker genes /tax related
	void calcLCAscore(MG_LCA*, options* opts);
	//6: species, 5:genus
	float getLCAscore(int lvl = 6) { if (lvl > LCAscores.size()) { return -1.f; }return LCAscores[lvl]; }
	float getLCAscoreWei(int lvl = 6) { if (lvl > LCAscores.size()) { return -1.f; }return LCAscores[lvl] * ((float)mkgene2type.size()/120.f); }
	void getTopTaxa(int lvl, float& frac, string& nm);
	vector<string>& getMajorityTax() { return majorityTax; }

	int getFoundGenes(int mode = 0) {
		if (mode == 0) { return foundGenes; }//actual gene cat genes
		else if (mode == 1) { return Ngenes; }//predicted from checkM2 file
		else if (mode == 2) { return foundMGs; }//marker gene
		return 0;
	}
	float MGoverlap(MAG* other);
	int hitGenes(geneCnter& x);
	size_t getNcats() { return mkgene2type.size(); }
	//provides a comma separated list of all genes in MAG
private:
	//functions
	float LCAscore2(SmplOccur&,string&, int minAssigns=10);
	string printGenes(bool includeMGs = false);
	string printMGbyCat(vector<string>&);

	
	//objects
	robin_hood::unordered_map <genetype, bool> genes;//stored if gene is marker gene
	//SmplOccur speciesAssignments;
	//taxonomy related
	vector<SmplOccur> TaxAssignments;
	vector<float> LCAscores;
	vector<string> majorityTax;
	
	ctg2geneLST ctg2Gene;
	vector<string> contigs;
	//vector< genetype> markerGenes;
	geneCnter markerGenes;
	unordered_map<string, vector<genetype>> mkgene2type;
	//quality of MAG/MGS
	float Completeness, Contamination,  GCcontent, codingDens;//LCAcompl,
	int Ngenes,GenomeSize,ctg_N50;
	float centreScore, compndScore;
	int foundGenes, foundMGs;
	int selfTier;
	MGS* mgsAssign;
	uint mgsRank;
	//identity
	string binName;
	string sampleID;
	bool isCanopy;

	//log best hit etc. has to be kept separate as for logging purposes only..
	float bestMatchFrac, uniqFrac;
	string bestMGS;
};


typedef robin_hood::unordered_map<string, MAG*> MAGlst;

class MGS {
public:
	MGS(string x):MGSid(x), uniqnessMG(0.f), uniqnessMGweigh(0.f), uniqnessALL(0.f),
		centre(nullptr),majorityTax(0), taxConsistency(0), MAGtaxConsistency(0), MAGtaxAssigns(0){
	}
	~MGS() {}
	bool addMAG(MAG* m);
	string getID() { return MGSid; }
	float matchGenes(geneCnter& tar, vector<float>& hits);
	size_t size() { return members.size(); }
	size_t g_size() { return geneCntAll.size(); }
	float getUniqness(int lvl = 0) {
		if (lvl == 0) { return uniqnessMG; }
		else if (lvl == 2) { return uniqnessALL; }
		else if (lvl == 1) { return uniqnessMGweigh; }
		return 0.f;
		}
	void setUniqness(float x, int lvl ) { if (lvl==0) { uniqnessMG = x;} else if (lvl == 2) { uniqnessALL = x; } else if (lvl == 1) { uniqnessMGweigh = x; }
}
	geneCnter& getMarkerGenes() { return geneCnt; }
	geneCnter& getAllGenes() { return geneCntAll; }
	void makeMGScentreMAG();
	MAG* getCentreMAG() {return centre;}
	MAGlst& getMAGs() { return members; }
	void createAllGeneCnts(); //accum all genes used
	string writeAllGenesSorted(geneCnter& allGene);
	string writeMGSscore();
	string writeMGSobservations();

	int totalMG() { return geneCnt.size(); }
	int totalMGassigned(MG_LCA* lca,int lvl=6);
	
	//taxonomy related..
	void calcLCAconsistency(MG_LCA* locLCA);//return LCA consistency at different levels
	vector<float> getMGStaxConsis() { return taxConsistency; }
	vector<float> getMAGtaxConsis() { return MAGtaxConsistency; }
	string getTaxStr();
	string getTaxConsisStr();
	string getTaxMAGconsisStr();
	vector<string>& getTax() { return majorityTax; }
	string getTax(int wh) { assert(wh < majorityTax.size()); return majorityTax[wh]; }
	float taxaHit(SmplOccur&,int lvl=6);
	//void setTaxConsistency(vector<float> in) { taxConsistency = in; }
private:
	//functions
	int getMultiBin(genetype g, geneCnter& allGene);

	//object
	MAGlst members;
	string MGSid;
	geneCnter geneCnt;//counts only marker genes
	geneCnter geneCntAll; //counts all marker + non-marker genes
	float uniqnessMG, uniqnessMGweigh, uniqnessALL;
	MAG* centre;

	//tax related functions.. calculated after MGS is fixed
	vector<string> majorityTax;
	vector<float> taxConsistency;
	vector<float> MAGtaxConsistency;
	vector<SmplOccur> MAGtaxAssigns;

}; 


//typedef robin_hood::unordered_map<string, MGS*> MGSlst;

class MAGs {
public:
	MAGs(MFmap* maps, options* opts);
	~MAGs();
	void readCanopies();
	bool assignTiers();
	void transform2GC(MFmap* maps);
	void createMAGsByCtg(MFmap* maps);
	bool readBinnerOut(string bFile, string qFile, string sample);
	void addLCAscores(MG_LCA*);
	void report();
	void reportMGS();
	void removeCtgInfo();
	void cluster2MGS();
	void makeMGScentreMAG();
	void calcLCAconsistency();
	//output
	void writeClusterCentres();
	void writeAllMAGs( MG_LCA* lca);
	void writeClusters();
	void writeClusterScores();
	void writeClusterObservations();
	void evalRemainderMAGs();

private:
	//functions
	bool passTier(MAG* m) { if (m->qualTier() < 0 || m->qualTier() >= maxTier) { return false; }return true; }
	void match2MGS(MAG* curMAG, MGS*& matchMGS, float& bestMtchFrac, float& uniqFrac);
	void calcMGSuni(bool markerGeneOnly=true);
	void add2allGenes(geneCnter&);
	bool readMAGqual(const string& qFile, string& sample, MAGlst& locMAGlist);



	//objects
	size_t smplN;
	MAGlst MAGlist;
	unordered_map<string, MAGlst> MAGsbyCtg;
	vector<MGS*> MGSs;

	geneCnter allGenes; //marker + non marker genes
	//geneCnter MGGenes; //marker genes

	int maxTier;

	options* opts;
	MG_LCA* lcaO;
};


