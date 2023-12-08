#include "objects.h"


MG_LCA::MG_LCA( string inD):
	gene2Tax(0), gene2Cat(0), geneCats(0, "") , TaxLvls(0,"")
{
    if (inD == "") { cerr << "path to MG_LCA is not correct!"; exit(33); }
    //read in dir with .LCA files..
    std::string path(inD);
    std::string ext(".LCA");

	struct stat sb;
	if (stat(path.c_str(), &sb) != 0) {
		cerr << "Path " << path << " does not exist!\n";
		exit(88);
	}

    for (auto& p : fs::recursive_directory_iterator(path))
    {
        if (p.path().extension() != ext) {
            continue;
        }
        std::cout << p.path().stem().string() << ' ';
        //std::cout << p.path() << '\n';
        geneCats.push_back(p.path().stem().string());

        addLCA(p.path().string(), p.path().stem().string());

    }
    report();
}

void MG_LCA::addMGs(string in) {
	std::ifstream infile(in);
	if (!infile) { cerr << "Could not find -MGfile " << in << endl; return; }
	std::string line;
	int cnt = -1;
	int fndGenes(0), newGenes(0);

	vector<string> emptyTax(TaxLvls.size(), "?");

	SmplOccur catCheck;
	for (auto ct : geneCats) { catCheck[ct] = 1; }

	while (std::getline(infile, line)) {
		if (line.substr(0, 1) == "#" || line.length() == 0) {
			continue;
		}
		//TIGR03723       165     1273130,1273131,1273132,1273133,
		stringstream ss(line); string category, geneCnt, genesCSV;
		getline(ss, category, '\t');
		getline(ss, geneCnt, '\t');
		getline(ss, genesCSV, '\t');
		stringstream ss2(genesCSV);

		vector<genetype> curLine(0);		string segments;
		while (getline(ss2, segments, ',')) {
			curLine.push_back((genetype)atoi(segments.c_str()));
		}
		if (atoi(geneCnt.c_str()) != curLine.size()) {
			cerr << category<<": Expected "<< atoi(geneCnt.c_str()) << " genes, found " << curLine.size() << endl;
		}

		//now transfer new info

		auto catFnd = catCheck.find(category);
		if (catFnd == catCheck.end()) {
			cerr << "Found new MG category (unexpected)! :" << category << " .. adding\n";
			geneCats.push_back(category);
		}
		
		for (auto gn : curLine) {
			auto found = gene2Tax.find(gn);
			if (found == gene2Tax.end()){//new gene!! add this
				newGenes++;
				gene2Tax[gn] = emptyTax;
				gene2Cat[gn] = category;
			}else {
				fndGenes++;
			}

		}

	}
	infile.close();
	cout << "Found " << fndGenes << " known MGs, added " << newGenes <<" marker genes"<< endl;
}

bool MG_LCA::hasTax(genetype gene,int lvl) {
	assert(lvl < 7);
	auto g2t = gene2Tax.find(gene);
	if (g2t != gene2Tax.end() && g2t->second[lvl] != "?") {
		return true;
	}
	return false;

}

bool MG_LCA::getTax(genetype gene, vector<string>& out, string& cls) {
	auto g2t = gene2Tax.find(gene);
	if (g2t != gene2Tax.end()) {
		out = g2t->second;
		cls = gene2Cat.find(gene)->second;
		return true;
	}
	return false;
}
void MG_LCA::report() {
    cout << "\nMarker Gene LCAs read.. Found " << gene2Tax.size() << " genes representing " << geneCats.size();
    cout << " MG categories with " << TaxLvls.size() << " taxonomic levels.\n";
}

bool MG_LCA::addLCA(string in, string cat) {
    std::ifstream infile(in);
	if (!infile) { return false; }
    std::string line;
    int cnt = -1;
    while (std::getline(infile, line))
    {
        cnt++;
        std::istringstream iss(line);
		//genetype a;
        string catcher; vector<string> tmp;
        if (cnt == 0) {
            if (TaxLvls.size() == 0) {
				getline(iss, catcher, '\t');
                while (getline(iss ,catcher,'\t')) {
                    tmp.push_back(catcher);
                }
                TaxLvls = tmp;
            }
            continue;
        }
        if (!getline(iss, catcher, '\t')) {
            cerr << "Couldn't read line: " << line<<endl;
            continue; 
        } // error
		genetype gene = (genetype)stoi(catcher);
        while (getline(iss, catcher, '\t')) {
            tmp.push_back(catcher);
        }
        gene2Tax[gene] = tmp;
        gene2Cat[gene] = cat;

        // process pair (a,b)
    }
    infile.close();
    return true;
}



MFmap::MFmap(string mapF):smplRid(0), smpls(0), mapGr(0),
smplLoc(0),smplN(0),curr(0),oldFolderStructure(false)
{
    stringstream ss2(mapF.c_str());  string segments;
    while (getline(ss2, segments, ',')) {
        string baseP = ""; //SmplOccurMult CntMapGrps ;
        read_map(segments, baseP);
        //read_abundances(CntMapGrps, baseP, opts->calcCoverage, opts->calcCovMedian, opts->oldMapStyle);

    }
}

void MFmap::read_map(const string mapF, string& baseP) {
	ifstream in;
	int map2folderIdx = 0; if (oldFolderStructure) { map2folderIdx = 1; }
	curr++;//keep track of different maps and inPaths

	in.open(mapF.c_str());
	if (!in) {
		cerr << "Couldn't open mapping file " << mapF << endl;
		exit(56);
	}
	cout << "Reading map " << mapF << endl;// " on path " << baseP[curr] << endl;
	
	string line(""); int cnt(-1); int assGrpN(-1);
	//mapping group params
	int mapGrpN(-1); //bool fillMapGrp(false); 
	int artiCntAssGrps(0); int skSmplCol(-1);

	string baseP1 = ""; string baseP2 = "";//#OutPath #RunID
	//string baseP = "";


	//1st part: just read headers
	while (getline(in, line)) {
		cnt++; int sbcnt(-1);
		stringstream ss(line); string segments;
		if (line.substr(0, 1) == "#") { //read column position of smpl id, path & assGrps
			//if (cnt > 0) { continue; }
			if (cnt == 0) {
				while (getline(ss, segments, '\t')) {
					sbcnt++;
					if (sbcnt == 0 && segments != "#SmplID") {
						cerr << "Map has to start with tag \"#SmplID\"\n"; exit(83);
					}
					if (sbcnt == 1 && !(segments == "Path" || segments == "SmplPrefix")) {
						cerr << "Map has to have tag \"Path\" || \"SmplPrefix\" as second entry\n";
						exit(83);
					}
					if (segments == "AssmblGrps") {
						assGrpN = sbcnt;
						cout << "Found Assembly groups in map\n";
					}
					if (segments == "MapGrps") {
						mapGrpN = sbcnt;
						//fillMapGrp = true;
						cout << "Found Mapping groups in map\n";
					}
					if (segments == "ExcludeAssembly") {
						skSmplCol = sbcnt;
						cout << "Samples can be excluded from assembly\n";
					}
				}
			}
			while (getline(ss, segments, '\t')) {
				if (segments == "#OutPath") {
					getline(ss, segments, '\t'); baseP1 = segments;
				}
				if (segments == "#RunID") {
					getline(ss, segments, '\t'); baseP2 = segments;
				}
			}

			continue;
		}
	}
	if (baseP2 == "") { cerr << "No #RunID specified in map!?"; exit(123); }
	if (baseP1 == "") { cerr << "No #OutPath specified in map!?"; exit(123); }
	baseP = baseP1 + "/" + baseP2 + "/";

	//2nd part: read samples
	
	in.clear();                 // clear fail and eof bits
	in.seekg(0, std::ios::beg); // back to the start!
	//cnt = 0;
	while (getline(in, line)) {
		//cnt++;
		if (line.substr(0, 1) == "#" || line.length()==0) {
			continue;
		}
		stringstream ss(line); string segments;
		vector<string> curLine(0);
		while (getline(ss, segments, '\t')) {
			curLine.push_back(segments);
		}
		if (curLine.size() == 0) { continue; }
		if (skSmplCol > -1 && curLine[skSmplCol] == "1") { continue; }
		string smpID = curLine[0];
		if (smpID == "") { continue; }//empty sample.. skip

		if (assGrpN >= 0 && curLine.size() <= assGrpN) {
			curLine.resize(assGrpN + 1, "");
		}
		//getline(ss, segments, '\t');
		//idx 1 for old folder structure, 0 for new folder structure
		string subDir = curLine[map2folderIdx];


		//assembly groups
		string assGrp("");
		if (assGrpN != -1) {
			//handles assembly groups from here
			assGrp = curLine[assGrpN];
		}
		else {//simulate CntAssGrps
			assGrp = to_string(artiCntAssGrps);
			artiCntAssGrps++;
		}
		if (CntAssGrps.find(assGrp) != CntAssGrps.end()) {
			CntAssGrps[assGrp].push_back((int)smplLoc.size());
		}
		else {
			CntAssGrps[assGrp] = vector<int>(1, (int)smplLoc.size());
		}

		if (assGrp != "" && CntAssGrps[assGrp].size() > 1) {
			string nsmpID = smpID + "M" + std::to_string(CntAssGrps[assGrp].size());
			if (smpls.find(nsmpID) != smpls.end()) {
				cerr << "Double sample ID: " << nsmpID << endl;
				exit(12);
			}
			smpls[nsmpID] = CntAssGrps[assGrp];//(int)smplLoc.size();
			smplRid[nsmpID] = smpID;
			smplAltId[smpID] = nsmpID;
		} else {
			if (smpls.find(smpID) != smpls.end()) {
				cerr << "Double sample ID: " << smpID << endl;
				exit(12);
			}
			smpls[smpID] = vector<int>(1, (int)smplLoc.size());
			smplRid[smpID] = smpID;
			smplAltId[smpID] = smpID;
		}
		assmblGr.push_back(assGrp);

		//mapping groups
		string mapGrp("");
		if (mapGrpN != -1) {
			mapGrp = curLine[mapGrpN];
		}
		if (mapGrp != "" && CntMapGrps.find(mapGrp) != CntMapGrps.end()) {
			(CntMapGrps)[mapGrp].push_back((int)smplLoc.size());
		}
		else if (mapGrp != "") {
			(CntMapGrps)[mapGrp] = vector<int>(1, (int)smplLoc.size());
		}


		mapGr.push_back(mapGrp);
		//useSmpl.push_back(true);
		smplLoc.push_back(baseP + "/" + subDir + "/");
		smplIDs.push_back(smpID);

	}
	in.close();
	smplN = smplLoc.size();

	cout << "read Map with " << smplN << " samples\n";


	//getAllSmpls(true, true);	int x = 0;//DEBUG

}

float MAG::LCAscore2(SmplOccur& taxs, string& bester, int minAssigns) {
	float hiSc (0), totSc(0), qCnt(0);
	bester = "";
	float ret (0.f);
	for (auto tx : taxs) {
		//completely ignore "?" cats..
		if (tx.first == "?") { continue; }// qCnt += (float)tx.second;}
		totSc +=(float) tx.second;
		if (tx.second > (float)hiSc) {
			hiSc = (float)tx.second;
			bester = tx.first;
		}
	}
	//needs to have enough evidence..
	if (totSc >= (float)minAssigns) {
		ret = hiSc / (totSc- qCnt);
	}
	if (qCnt > 5.f && totSc < float(minAssigns)) {
		ret = 1.f;//well at least consistent for "unknown"..
		bester = "?";
	}
	return ret;
}

void MAG::getTopTaxa(int lvl, float& frac, string& nm) {
	frac = LCAscores[lvl];
	nm = majorityTax[lvl];
	
	/*SmplOccur tmp;
	tmp = TaxAssignments[lvl];
	
	float hiCnt(0.f); float totalC(0.f); string hiTax("");
	float qCnt(0.f);
	for (auto tax : tmp) {
		totalC += (float)tax.second;
		if (tax.first == "?") {
			qCnt += (float)tax.second;
		} else if (tax.second > hiCnt) {
			hiCnt = (float) tax.second; hiTax = tax.first;
		}
	}
	nm = hiTax;
	frac = hiCnt / (totalC - qCnt);
	*/
}

void MAG::calcLCAscore(MG_LCA* lca,options* opts) {
	vector<string> dummy; string cls;
	if (TaxAssignments.size() == 0) {
		TaxAssignments.resize(lca->maxTaxLvl());
		majorityTax.resize(lca->maxTaxLvl(),"");
		LCAscores.resize(lca->maxTaxLvl(),-1.f);
	}
	for (auto gn : genes) {
		genetype gene = gn.first;
		if (lca->getTax(gene, dummy, cls)) {
			//genus level assignments of gene
			for (size_t k = 0; k < lca->maxTaxLvl(); k++) {
				auto genf = TaxAssignments[k].find(dummy[k]);
				if (genf != TaxAssignments[k].end()) {
					genf->second++;
				}else { TaxAssignments[k][dummy[k]] = 1; }

			}
			gn.second = true;
			foundMGs++;
			auto mgfnd = markerGenes.find(gene);
			if (mgfnd == markerGenes.end()) { markerGenes[gene] = 1; }
			else { mgfnd->second++; }
			auto g2cl = mkgene2type.find(cls);
			if (g2cl == mkgene2type.end()) {
				mkgene2type[cls] = vector<genetype>(1, gene);
			}
			else {
				g2cl->second.push_back(gene);
			}
		}
	}
	int minAssigns = opts->minMGassigned2tax;
	//calc actual score by best hit
	for (size_t k = 0; k < lca->maxTaxLvl(); k++) {
		string tmp("");
		LCAscores[k] = LCAscore2(TaxAssignments[k], tmp, minAssigns);
		majorityTax[k] = tmp;
	}
	//genLCAscore = LCAscore2(TaxAssignments[5]);
}



vector<string> MFmap::getAllSmpls(bool onlyAssmbled, bool assmGrpID ) {
	if (!onlyAssmbled && !assmGrpID) {
		return smplIDs;
	}
	SmplOccur currCntAsGr;// = map->getAssmGrpCnts();
	vector<string> ret;
	for (size_t i = 0; i < smplIDs.size(); i++) {
		string assmblGrp = getAssmblGrp(i);
		if (onlyAssmbled && assmblGrp != "") {
			SmplOccurIT cMGcnts = currCntAsGr.find(assmblGrp);
			if (cMGcnts == currCntAsGr.end()) {
				currCntAsGr[assmblGrp] = 1;
			}
			else { (*cMGcnts).second++; }
			if (AssmblGrpSize(assmblGrp) != (uint)currCntAsGr[assmblGrp]) {
				continue;
			}
		}

		if (assmGrpID) {
			ret.push_back(smplAltId[smplIDs[i]]);
		}
	}

	//DEBUG
	//for (auto smpl : ret) { cout << "\""<< smpl << "\"" << endl; }//exit(0);

	return ret;
}







void MAG::storeBestHit(float bestMtchFrac, MGS* matchMGS, float uniq) {
	bestMatchFrac = bestMtchFrac;
	uniqFrac = uniq;
	if (matchMGS==nullptr){
		bestMGS = "";
	} else {
		bestMGS = matchMGS->getID();
	}
}




int MAG::hitGenes(geneCnter& x) {
	int ret = 0;
	for (auto zz : x) {
		if (genes.find(zz.first) != genes.end()) {
			ret++;
		}
	}
	return ret;
}

void MAG::assignTier(options* op) {
	int maxTierRds = op->maxTier();
	selfTier = 100;
	for (uint T = 0; T < maxTierRds; T++) {
		if (getComple() > op->complTiers[T] && getConta() < op->contaTiers[T]
			&& (LCAscores.size()==0 || LCAscores[6] < 0 || LCAscores[6] > op->LCAcomplTier[T])
			) {
			selfTier = T;
			break;
		}
	}
}

float MAG::MGoverlap(MAG* other) {
	if (markerGenes.size() == 0) { return 0.f; }
	int shrdGenes = other->hitGenes(markerGenes);
	return shrdGenes / (float)markerGenes.size();
}


//case for Canopies: we don't have info on ctgs..
void MAG::setGenes(vector<genetype> x) {
	assert(genes.size() == 0);
	for (auto g : x) {
		genes[g] = false;
	}
	vector < pair<int, genetype> > yy (x.size());
	for (int i = 0; i < (int)x.size(); i++) {
		yy[i] = pair<int, genetype>(i, x[i]);
	}
	ctg2Gene["Ctg???"] = yy;
	foundGenes = genes.size();
	
}
void MAG::addGene(genetype Gn, string Ctg,int num){
	auto found = ctg2Gene.find(Ctg);
	if (found == ctg2Gene.end()) {
		ctg2Gene[Ctg] = vector< pair<int, genetype >>(1, pair<int,genetype>(num,Gn));
	} else {
		found->second.push_back(pair<int, genetype>(num, Gn));
	}
	auto f2 = genes.find(Gn);
	if (f2 == genes.end()){
		genes[Gn] = false;
	}
	foundGenes++;
}


float MAG::compoundScore(int maxGenes, int maxN50) {
	compndScore = 0.f; float cnter(0.f);
	compndScore += (getComple()/100.f); cnter+=1.f;//
	if (getLCAscore() >= 0.f) {
		compndScore += ((getLCAscore() + getLCAscore(5)) / 2.f); cnter += 1.f;
	}
	if (maxGenes > 0 && maxN50 > 0) {
		compndScore += float(foundGenes / maxGenes) + float(ctg_N50 / maxN50) / 2.f; cnter += 1.f;
	}
	compndScore -= (2.f * cnter * getConta() / 100.f);
	compndScore /= cnter;
	return compndScore;
}


string MAG::getAssignedMGSname(bool retAlt) {
	if (retAlt) {
		if (altMGShits == nullptr) { return ""; }
		return altMGShits->getID();
	} else {
		if (mgsAssign == nullptr) { return ""; }
		return mgsAssign->getID();
	}
}
string MAG::printMGbyCat(vector<string>& cats) {
	string str;
	for (string cat : cats) {
		auto found = mkgene2type.find(cat);
		if (found == mkgene2type.end()) {
			str += "\t";
		} else {
			str += to_string(found->second[0]);
			for (size_t ii = 1; ii < found->second.size(); ii++) {
				str += "," + to_string(found->second[ii]);
			}
			str += "\t";
		}
	}
	return str;
}

string MAG::printGenes(bool includeMGs ) {
	string outStr = "";
	for (auto ctg : ctg2Gene) {
		vector< pair<int, genetype>> gens = ctg.second;
		sort(gens.begin(), gens.end(), [](const pair<int, genetype>& a, const pair<int, genetype>& b) { return(a.first < b.first); });

		int lastG(0);
		for (auto gg:gens){
			int diff = gg.first- lastG;
			for (int x = 0; x < (diff-1); x++) {
				outStr += "?,";
			}
			if (!includeMGs && genes.find(gg.second)->second) {
				outStr += "?,";//do nothing..
			}else {
				outStr += to_string(gg.second) + ",";
			}
			lastG = gg.first;
			//break;//DEBUG
		}
		outStr += ",";
		//break;
	}
	outStr = outStr.substr(0, (outStr.length() - 2));
	return outStr ;
}

string MAG::formatMAG(vector<string>& cats, MG_LCA* lcaL, bool isCtr, bool MGs_with_others) {
	stringstream of;
	//MAG\tMGS\tRepresentative4MGS\tMatch2MGS\tUniqueness\tAssociatedMGS\tCompleteness
	//bool MGs_with_others(true);if (cats.size() > 0) { MGs_with_others = false; }
	of << getUniqName() << "\t" << getAssignedMGSname(false) << "\t";
	if (isCtr) { of << "*\t"; }else { of << "\t"; }//centre of MAG?
	of << MGSmatch << "\t" << uniqScore << "\t" << getAssignedMGSname(true) << "\t";
	of << getComple() << "\t" << getConta() << "\t"  << getLCAscore() << "\t" << getN50();
	of << "\t" << getFoundGenes(0) << "\t" << getGC() << "\t" << getCodingDens();
	of << "\t" << getCentreScore() << "\t" << getCompoundScore() << "\t";

	//add tax here
	size_t sizTax = lcaL->getTaxLvls().size();
	if (majorityTax.size() == sizTax) {
		for (size_t x = 0; x < majorityTax.size(); x++) {
			of << majorityTax[x] << "\t";
		}
	}
	else {
		for (size_t x = 0; x < sizTax; x++) {//just inserrt empty
			of << "\t";
		}
	}

	if (cats.size()>0) {
		of<<printMGbyCat(cats);
	}
	of << printGenes(MGs_with_others) ;
	return of.str();
}

MAGs::MAGs(MFmap* map,options* opt):
	smplN(0), MGSs(0), maxTier(opt->maxTier()), opts(opt), lcaO(nullptr)
{
	string path2Bins = opts->path2Bins;
	string path2Assmbl = opts->path2Assmbl;
	string file2Assmbl = opts->file2Assmbl;
	smplN = map->totalSmpls();
	string qualSuffix = opts->qualSuffix;
	uint geneN = 0;
	uint preMapSize(smplN);
	cout << "Reading Bins from " << smplN << " samples\n";
	//read the gene abundances sample-wise in
	SmplOccur currCntAsGr;// = map->getAssmGrpCnts();
	uint smplsIncl = 0;
	for (uint i = 0; i < smplN; i++) {
		smplsIncl++;
		//only include last sample of mapping group..
		string assmblGrp = map->getAssmblGrp(i);
		if (assmblGrp != "") {
			SmplOccurIT cMGcnts = currCntAsGr.find(assmblGrp);
			if (cMGcnts == currCntAsGr.end()) {currCntAsGr[assmblGrp] = 1;
			}else {(*cMGcnts).second++;}

			if (map->AssmblGrpSize(assmblGrp) == (uint)currCntAsGr[assmblGrp]) {
				string gFile(map->SmplPath(i) + path2Assmbl + path2Bins + map->SmplIds(i));
				if (!file_exists(gFile)) {
					gFile = map->SmplPath(i) + path2Assmbl + file2Assmbl;
					if (file_exists(gFile)) {
						ifstream inf(gFile); getline(inf, gFile);
						gFile += path2Bins + map->SmplIds(i);
					}
				}
				//cout << "Assmgrp:"<<gFile << endl;
				cout << "Assmgrp:";
				readBinnerOut(gFile, gFile + qualSuffix, map->SmplIdsAlt(i));

			}
		}
		else {
			string gFile(map->SmplPath(i) + path2Assmbl + path2Bins + map->SmplIds(i));
			//cout << gFile << endl;
			readBinnerOut(gFile, gFile + qualSuffix, map->SmplIds(i));

		}

	}
	cout << "Including " << smplsIncl << " samples\n";// , using "<< geneN<<" genes\n";
}
MAGs::~MAGs() {
	for (auto mag : MAGlist) {
		delete mag.second;
	}
}


bool MAGs::readBinnerOut(string bFile, string qFile, string sample) {
	std::ifstream infile(bFile);
	if (!infile) { 
		cerr << "Could not find bin file " << bFile << endl;
		return false; 
	}
	std::string line;
	int cnt = 0; int binCnt(0);
	string lastBin = ""; vector<string> contigs;

	MAGlst locMAGlist;
	//vector<string> bins;//DEBUG

	while (std::getline(infile, line)) {
		cnt++;
		std::istringstream iss(line);
		string ctg, bin;
		if (!(iss >> ctg >> bin)) { continue; }
		if (bin != lastBin) {
			if (lastBin != "") {//create this bin..
				string uniID = sample + "__" + lastBin;
				auto res = locMAGlist.find(uniID);
				if (res == locMAGlist.end()) {
					locMAGlist[uniID] = new MAG(lastBin, sample, contigs);
					binCnt++;
				}
				else {
					cerr << "MAG uniq ID: " << uniID << " was found twice!!";
				}
				//bins.push_back(lastBin);//DEBUG
			}
			lastBin = bin; contigs.clear();
		}
		size_t pos = ctg.find("__");
		if (pos != std::string::npos) {
			contigs.push_back(ctg.substr(pos + 2));
		}
		else {
			cerr << "Wrongly formatted contig:" << ctg << endl;
		}
		//contigs.push_back(ctg);
	}
	//remaining Bin..
	if (lastBin != "") {
		locMAGlist[(sample + "__" + lastBin)] = new MAG(lastBin, sample, contigs);
		binCnt++;
	}

	infile.close();

	readMAGqual(qFile, sample, locMAGlist);


	//now transfer into big list..
	int passed(0); int passedContig(0);
	for (auto MG : locMAGlist) {
		if (opts->processOnlyHQ) { MG.second->assignTier(opts); }
		else { MG.second->assignTier((uint)0); }
		if ( passTier(MG.second)) { passed++; passedContig += MG.second->NContigs(); }
		MAGlist[MG.first] = MG.second;
	}

	cout << sample << ":: passed " << passed<< "/" << binCnt << " MAGs with " << passedContig <<"/" << cnt << " contigs total.\n";

	return true;

}

bool MAGs::readMAGqual(const string& qFile, string& sample, MAGlst& locMAGlist){
	ifstream infile(qFile);
	if (!infile) { cerr << "Could not open " << qFile << " !"; }
	vector<string> inTmp(20, "");
	bool isCM2(true);
	float Conta, Compl, GC, codingDens;
	int Ngenes, GenomeSize, N50;
	string bin; string line; int cnt(0);
	while (std::getline(infile, line)) {
		cnt++;
		std::istringstream iss(line);
		string tmp; size_t cc = 0;
		while (getline(iss, tmp,'\t')) {
			inTmp[cc] = tmp;
			cc++;
			if (cc > 19) { break; }
		}
		if (inTmp[0] == "Name") { continue; }//header in .cm2 file..
		if (isCM2) {
			if (cc < 2) { cerr << "On .cm2 line:\n" << line << "\nless than 2 elements were found!\nAborting\n"; exit(61); }
			bin = inTmp[0];
			string uniID = sample + "__" + bin;
			auto res = locMAGlist.find(uniID);
			if (res == locMAGlist.end()) {
				cerr << "Could not find MAG " << uniID << " in ref Bin! Cannot add qual info.\n";
				continue;
			}
			Compl = atof(inTmp[1].c_str());
			Conta = atof(inTmp[2].c_str());
			res->second->setCompl(Compl); res->second->setConta(Conta);
			if (cc > 5) {
				GC = atof(inTmp[9].c_str());
				Ngenes = atoi(inTmp[10].c_str());
				GenomeSize = atoi(inTmp[8].c_str());
				codingDens = atof(inTmp[5].c_str());
				N50 = atoi(inTmp[6].c_str());
				res->second->setGC(GC); res->second->setN50(N50);
				res->second->setNgenes(Ngenes); res->second->setGenomeSize(GenomeSize);
				res->second->setcodingDens(codingDens);
			}
		}
	}

	infile.close();

	return true;
}

void  MAGs::readCanopies() {
	string inF1 = opts->canopyFile;
	if (inF1 == "") { return; }
	string qFile = inF1 + ".cm2";
	ifstream in(inF1.c_str());
	if (!in) { cerr << "Could not open canpy file : " << inF1; exit(523); }

	std::string line;
	int cnt = -1; int binCnt(0);
	string lastBin = ""; vector<genetype> genes;

	MAGlst locMAGlist;
	//vector<string> bins;//DEBUG
	string sample("Cano");
	vector<string>emptyVec(0, "");

	while (std::getline(in, line)) {
		cnt++;
		std::istringstream iss(line);
		string  bin;
		genetype gene;
		if (!(iss >> bin >> gene)) { continue; }
		if (bin != lastBin) {
			if (lastBin != "") {//create this bin..
				string uniID = sample + "__" + lastBin;
				auto res = locMAGlist.find(uniID);
				if (res == locMAGlist.end()) {
					MAG* cano = new MAG(lastBin, sample, emptyVec);
					cano->setCanopy(true); cano->setGenes(genes);
					cano->setGC(-1);
					cano->setN50(1000);//some average gene size
					locMAGlist[uniID] = cano;
					binCnt++;
				}
				else {
					cerr << "MAG uniq ID: " << uniID << " was found twice!!";
				}
				//bins.push_back(lastBin);//DEBUG
			}
			lastBin = bin; genes.clear();
		}
		genes.push_back(gene);
		//contigs.push_back(ctg);
	}
	//remaining Bin..
	MAG* cano = new MAG(lastBin, sample, emptyVec);
	cano->setCanopy(true); cano->setGenes(genes);
	locMAGlist[(sample + "__" + lastBin)] = cano;
	//bins.push_back(lastBin);//DEBUG
	binCnt++;
	in.close();

	readMAGqual(qFile, sample, locMAGlist);


	//now transfer into big list..
	int passed(0), passedGenes(0);
	for (auto MG : locMAGlist) {
		if (opts->processOnlyHQ) { MG.second->assignTier(opts); }
		else { MG.second->assignTier((uint)0); }

		if (passTier(MG.second)) { passed++; passedGenes += MG.second->NGenes(); }

		MAGlist[MG.first] = MG.second;
	}

	cout << sample << ":: passed " << passed<< "/" << binCnt << " MGS with " << passedGenes<<"/" <<cnt << " genes total.\n";

}



void MAGs::reportMGS() {
	float avgMAGsInMGS(0.f);
	float avgUni(0.f), avgUniAll(0.f);
	for (auto mgs : MGSs) {
		avgMAGsInMGS+=(float)mgs->size();
		avgUni += mgs->getUniqness(0);
		avgUniAll += mgs->getUniqness(2);
	}
	int totMAGs = (int)avgMAGsInMGS;
	avgMAGsInMGS /= (float)MGSs.size();
	avgUni /= (float)MGSs.size();
	avgUniAll /= (float)MGSs.size();
	//calc how many MAGs are without match
	int mgsMAG(0),  badMAG(0);
	vector<int> goodMAG(maxTier, 0);
	for (auto mag : MAGlist) {
		if (mag.second == nullptr || mag.second->getAssignedMGS() != nullptr) { 
			mgsMAG++;
			continue; 
		}
		if (passTier(mag.second)) {
			goodMAG[mag.second->qualTier()]++;
		} else {
			badMAG++;
		}
	}
	cout << "\n\n-----------------------------------------------------------------\n";
	cout << "Found " << MGSs.size() << " MGS with " << totMAGs << " MAGs (avg " << round(avgMAGsInMGS * 100) / 100 << " MAGs/MGS).\n";
	cout << round(avgUni * 100) / 100 << "/" << round(avgUniAll * 100) / 100 << " avg MG/ALL Uniqueness.\n";
	cout << mgsMAG << " MAGs in MGS, unassigned ";
	string goodMAGstr1(to_string(goodMAG[0])); string goodMAGstr2("T0"); int goodSum(0);
	for (int i = 1; i < maxTier; i++) {
		goodMAGstr1 += "," + to_string(goodMAG[i]);
		goodMAGstr2 += ",T" + to_string(i);
		goodSum += goodMAG[i];
	}
	cout << goodSum << " ("<<goodMAGstr1 << " in " << goodMAGstr2 << ") high qual MAGs";
	cout<< " and " << badMAG << " low qual MAGs.\n";
	cout << "-----------------------------------------------------------------\n\n";
}

void MAGs::report() {
	float avgGenesFnd(0.f), avgGenesExpec(0.f), avgCont(0.f),avgComp(0.f);
	float avgMGgenes(0.f), avgLCA(0.f), avgLCAg(0.f);
	float cnt(0);
	int LCAmiss(0);
	for (auto mag : MAGlist) {
		avgGenesFnd += (float)mag.second->getFoundGenes();
		avgGenesExpec += (float)mag.second->getFoundGenes(false);
		avgCont += (float)mag.second->getConta();
		avgComp += (float)mag.second->getComple();
		float MGgenes = (float)mag.second->getFoundGenes(2);
		float lcaloc = (float)mag.second->getLCAscore();
		if (MGgenes == 0.f) { LCAmiss++; 
		} else {
			avgLCA += lcaloc;
			avgLCAg += (float)mag.second->getLCAscore(5);
			avgMGgenes += MGgenes;
		}
		cnt += 1.f;
	}
	avgGenesFnd /= cnt; avgCont /= cnt; avgGenesExpec /= cnt;
	avgComp /= cnt; 
	avgLCA /= (cnt- LCAmiss); avgMGgenes /= (cnt - LCAmiss); avgLCAg /= (cnt - LCAmiss);

	cout << "Found "<< MAGlist.size() << " total MAGs across " << smplN << " Samples.\n";
	cout << "These had on average:\n";
	cout << "  - " << avgGenesFnd << " genes asigned (" << avgGenesExpec << " expected, "<< avgMGgenes<< " MG genes)\n";
	cout << "  - " << avgCont << " Contamination\n";
	cout << "  - " << avgComp << " Completeness\n";
	cout << "  - " << avgLCA << " species, " << avgLCAg << " genus LCA ("<< LCAmiss<<" missing)\n";
}

void MAGs::createMAGsByCtg(MFmap* maps) {

	assert(MAGsbyCtg.size() == 0);

	cout << "Creating MAG by contig list..";
	maxTier = opts->maxTier();
	uint skippedMAGs(0), MAGkept(0);

	//first create dummy samples
	for (string smpl : maps->getAllSmpls(true,true)) {
		MAGlst local;
		MAGsbyCtg[smpl] = local;
	}

	//string ctg = "";
	for (auto MG : MAGlist) {
		if (!passTier(MG.second)) {
			//cerr << MG.second->qualTier() << " ";//DEBUG
			skippedMAGs++; continue; 
		}
		MAGkept++;
		string smplID = MG.second->getSmplID();
		auto smplFnd = MAGsbyCtg.find(smplID);
		MAGlst local;
		if (smplFnd != MAGsbyCtg.end()){
			local = smplFnd->second;
		}
		for (string ctg : MG.second->getContigs()) {
			local[ctg] = MG.second;
		}
		//create anew
		if (smplFnd == MAGsbyCtg.end()) {
			MAGsbyCtg[smplID] = local;
		} else {
			smplFnd->second = local;
		}
	}

	cout << "done with " <<MAGkept<< " MAGs.";
	if (skippedMAGs) {
		cout << " Skipped " << skippedMAGs << " MAGs due to low qual.";
	}
	cout << endl;
}

void MAGs::addLCAscores(MG_LCA* lcas ) {
	lcaO = lcas;
	cout << "Calculating LCA completeness scores..\n";
	size_t lcacnt(0); int LCAassigned = 0;
	for (auto mags : MAGlist) {
		if (!passTier(mags.second)) { continue; }
		mags.second->calcLCAscore(lcaO, opts);
		if (mags.second->getLCAscore() > 0.f) { LCAassigned++; }
		lcacnt++;
	}
	cout << "Done calculating. "<< LCAassigned<<"/"<< lcacnt << " MAGs with LCA scores. \n";
}

void MAGs::removeCtgInfo(){
	MAGsbyCtg.clear();
	for (auto mag : MAGlist) {
		mag.second->clearCtg();
	}
}

void MAGs::cluster2MGS() {
	vector<vector<pair<float,MAG*>> > MAGinTiers;
	//vector<vector<float>> ContaInTiers;//records contamination of MAGs in MAGinTiers
	vector<pair<float, MAG*>> nullTV(0, pair<float, MAG*>(-2.f,nullptr));
	//create 3 Tiers at which MAGs could be..
//	my @ComplTiers = (98, 95, 80, 60); my @ContaTiers = (1, 3, 5, 10);
//	my @LCAcompl = (0.97, 0.95, 0.92, 0.9); my @AllowNovelMGS = (1, 1, 1, 0);
//	my $maxRounds = 3; my @OverlapTiers = (0.8, 0.8, 0.9, 0.9); my @unknwnTiers = (0.75, 0.75, 0.8, 1.1);

	vector<float> complTiers = opts->complTiers;
	vector<float> contaTiers = opts->contaTiers;
	
	vector<float> LCAcomplTier = opts->LCAcomplTier; vector<int> allowNovelMGSTiers = opts->allowNovelMGSTiers;
	vector<float> match2MGS_Tier = opts->match2MGS_Tier; 
	vector<float> globalUniq4MGS_Tier = opts->globalUniq4MGS_Tier;
	vector<float> maxMatch4NovelMGS_Tier = opts->maxMatch4NovelMGS_Tier;
	int maxTierRds = opts->maxTier();

	MAGinTiers.resize(maxTierRds, nullTV);
	//ContaInTiers.resize(maxTierRds, vector<float>(0));


	vector<uint> TierCnts(complTiers.size(), 0); size_t addedCnt(0);

	for (auto mag : MAGlist) {
		//reassign since now LCA scores were added..
		mag.second->assignTier(opts);
		int curT = mag.second->qualTier();
		if (curT>=0 && curT < maxTierRds) {//do something with this..
			TierCnts[curT]++;addedCnt++;
			MAGinTiers[curT].push_back(pair<float, MAG*>(mag.second->getLCAscoreWei(), mag.second));
			//ContaInTiers[curT].push_back(mag.second->getConta());
		}
	}
	cout << "MAG in Tiers: "; int cntt(0);
	for (auto x : TierCnts) { cout << "T" << cntt << "->" << x << " "; cntt++; }
	cout << "out of "<< MAGlist.size()<< " MAGs." << endl;

	vector<MAG*> srtList(addedCnt,nullptr);//figure out later how to get the ordered vector..
	size_t cnter(0);
	for (uint T = 0; T < maxTierRds; T++) {
		sort(MAGinTiers[T].begin(), MAGinTiers[T].end(), [](const pair<float, MAG* >& a, const pair<float, MAG*>& b) { return(a.first > b.first); });
		for (auto curMAG : MAGinTiers[T]) {
			srtList[cnter] = curMAG.second;
			cnter++;
		}
	}
	string MGSprefix = opts->MGSprefix;
	for (uint T = 0; T < maxTierRds; T++) {
		uint addMAGcnt(0), matchMAGcnt(0), ambMAGcnt(0);
		cout <<"-------------------------- At Create Tier "<<T<<" --------------------------\n";
		for (MAG* curMAG : srtList) {
			if (curMAG==nullptr || curMAG->getAssignedMGS() != nullptr ) { continue; }
			int curT = curMAG->qualTier();
			if(curT<0 || curT > maxTierRds) { continue; }
			bool makeMGS(false);
			//try matching to existing MGS (core genes)
			MGS* matchMGS(nullptr); float bestMtchFrac(0.f);float uniqFrac(0.f);
			
			match2MGS(curMAG, matchMGS, bestMtchFrac, uniqFrac);
			
			//testing..
			//uniqFrac = 1 - bestMtchFrac;
			
			curMAG->storeBestHit(bestMtchFrac, matchMGS, uniqFrac);
			if (uniqFrac > globalUniq4MGS_Tier[T]
					&& bestMtchFrac < maxMatch4NovelMGS_Tier[T]
					&& curMAG->qualTier() <= allowNovelMGSTiers[T]
				) {
				makeMGS = true;
			
			} else if(matchMGS != nullptr && bestMtchFrac >= match2MGS_Tier[T]) {
				matchMGS->addMAG(curMAG); curMAG->setAssignedMGS(matchMGS, bestMtchFrac, uniqFrac);
				matchMAGcnt++;
				//cerr << "Found match " << matchMGS->getID() << "<->";
				//cerr<< curMAG->getUniqName() << " (" << bestMtchFrac << ", " << uniqFrac << ")\n";
			
			} else {//worst case, MAG can't be assigned..
				//cout << "A ";
				ambMAGcnt++;
				//log this..
			}
			
			if (makeMGS) {
				string mgsName = MGSprefix + to_string(MGSs.size());
				matchMGS = new MGS(mgsName);
				MGSs.push_back(matchMGS);
				matchMGS->addMAG(curMAG);
				curMAG->setAssignedMGS(matchMGS, bestMtchFrac, uniqFrac);
				addMAGcnt++;
				//cerr << "Creating " << mgsName << " (" << bestMtchFrac << ", " << uniqFrac << "; Compl:"<< curMAG->getComple() << " Conta:"<< curMAG->getConta() <<" LCAspec:" << curMAG->getLCAscore() << ")\n";
			}
		}
		cout << "Created " << addMAGcnt << " MGS, matched " << matchMAGcnt << " MAGs to these, " << ambMAGcnt << " ambigous\n";
	}

	
	//postprocessing / reporting on created MGS
	calcMGSuni();


	for (auto mgs : MGSs) {
		mgs->createAllGeneCnts();
	}
	calcMGSuni(false);
	reportMGS();
}


void MAGs::evalRemainderMAGs() {
	//function that check how many MAGs could be potentially assigned to MGS

	vector<float> mtThr({.9,.8,.7,.6,.5,.4,.3,.2,.1,.0});
	vector<int> emptVec(mtThr.size());
	vector<vector<int>> matchMatr(maxTier, emptVec);
	vector<vector<int>> uniqMatr(maxTier, emptVec);
	for (auto mag : MAGlist) {
		MAG* curMAG = mag.second;
		if (curMAG == nullptr || curMAG->getAssignedMGS() != nullptr || !passTier(curMAG)) { continue; }
		MGS* matchMGS(nullptr); float bestMtchFrac(0.f); float uniqFrac(0.f);
		match2MGS(curMAG, matchMGS, bestMtchFrac, uniqFrac);
		//catalog best matches
		for (size_t x = 0; x < mtThr.size(); x++) {
			if (bestMtchFrac > mtThr[x]) {
				matchMatr[curMAG->qualTier()][x]++;
				break;
			}
		}
		//catalog how unique these MAGs really are..
		for (size_t x = 0; x < mtThr.size(); x++) {
			if (uniqFrac > mtThr[x]) {
				uniqMatr[curMAG->qualTier()][x]++;
				break;
			}
		}

		if (matchMGS != nullptr) {
			curMAG->setNextBestMGS(matchMGS, bestMtchFrac, uniqFrac);
		}

		
	}

	cout << "Unassigned MAGs could match at the following threshholds:\n";
	printMatrixNice(matchMatr, mtThr,"Tier");
	cout << "Unassigned MAGs are unique at the following threshholds:\n";
	printMatrixNice(uniqMatr, mtThr,"Tier");
	cout << endl;


}

void MAGs::match2MGS(MAG* curMAG, MGS*& matchMGS, float& bestMtchFrac, float& uniqFrac) {
	geneCnter tarSet = curMAG->getMGgenes();
	if (tarSet.size() == 0) { return; }
	float bestHits = 0.f;
	vector<float> matchedTars(tarSet.size(), 0.f);//records which of the marker genes were so far matched
	for (auto mgs : MGSs) {
		float hits = mgs->matchGenes(tarSet, matchedTars);
		if (hits > bestHits) {
			bestHits = hits;
			matchMGS = mgs;
		}
	}
	bestMtchFrac = float(bestHits) / float(tarSet.size());
	int uniHit(0);
	for (float mm : matchedTars) {
//		if (mm) { uniHit++; }
		if (mm>0.1) { uniHit++; }
	}
	//uniqFrac = 1.f - (float(uniHit) / float(tarSet.size()));
	float div = float(tarSet.size());// *(MGSs.size() - 1));
	if (div == 0) { uniqFrac = 1.f; }
	else { uniqFrac = 1.f - (uniHit / div); }

	return;
}


void MAGs::add2allGenes(geneCnter& lst) {
	size_t lstS = lst.size();
	size_t allGS = allGenes.size();

	for (auto gn : lst) {
		auto found = allGenes.find(gn.first);
		if (found == allGenes.end()) {
			allGenes[gn.first] = 1;
		}
		else {
			found->second++;
		}
	}	
}


void MAGs::calcMGSuni(bool markerGeneOnly) {
	int genesTot(0), genesMulti(0);
	geneCnterFloat glbGenesW;
	geneCnter MGGenes;//create anew each time..
	float Nmgs = (float)MGSs.size();
	for (auto mgs : MGSs) {
		float lsize = float(mgs->size());
		geneCnter locals;
		if (markerGeneOnly) { locals = mgs->getMarkerGenes(); }
		else { locals = mgs->getAllGenes(); }
		for (auto gn : locals) {
			auto found = MGGenes.find(gn.first);genesTot++;
			if (found == MGGenes.end()) {MGGenes[gn.first] = 1;
			}else {found->second++;genesMulti++;}
			auto fnd2 = glbGenesW.find(gn.first);
			if (fnd2 == glbGenesW.end()) {glbGenesW[gn.first] = gn.second/ lsize;
			}else { fnd2->second+= gn.second/ lsize;  }
		}
	}
	//just get max value for a global occurring gene (weighted)
	float maxGlbGenesW(0.f);
	for (auto X : glbGenesW) { if (X.second > maxGlbGenesW) { maxGlbGenesW = X.second; } }

	vector<float> uniqq; vector<float> uniqqWei;
	for (auto mgs : MGSs) {
		geneCnter locals;
		float lsize = float(mgs->size());
		assert(lsize > 0);
		if (markerGeneOnly) { locals = mgs->getMarkerGenes(); }
		else { locals = mgs->getAllGenes(); }
		int uni(0), multi(0);  float uniWei(0.f), multiWei(0.f);
		for (auto gn : locals) {
			auto found = MGGenes.find(gn.first);
			if (found->second == 1) { uni++; }
			else { multi++; }
			auto fnd2 = MGGenes.find(gn.first);
			float w_MGS = (gn.second / lsize);
			float w_GLB = (fnd2->second- w_MGS) / maxGlbGenesW;//substract itself from glb gene weight..
			uniWei += w_MGS * w_GLB;
		}
		float locUni = (float)uni / ((float)uni + (float)multi);
		float locUniWei = 1.f - (uniWei / float(locals.size()));
		uniqq.push_back(locUni);
		uniqqWei.push_back(locUniWei);
		if (markerGeneOnly) {
			mgs->setUniqness(locUni, 0); mgs->setUniqness(locUniWei, 1);
		} else{ mgs->setUniqness(locUni, 2); }

	}
	if (markerGeneOnly) {
		cout << "Marker genes in MGS: " << genesTot << " unique, " << genesMulti << " shared among MGS.\n";
	}
	else {
		cout << "All genes in MGS:    " << genesTot << " unique, " << genesMulti << " shared among MGS.\n";
	}

	sort(uniqq.begin(), uniqq.end());
	if (markerGeneOnly) {
		vector<float> histThr({ 1.,0.98,.95,.9,.85,.8,.7,.6,.5,.25,.1,.0 });
		vector<int> histo = createHist(uniqq, histThr);
		vector<int> histoWei = createHist(uniqqWei, histThr);
		vector<vector<int>> hist2; hist2.push_back(histo); hist2.push_back(histoWei);

		cout << "\nUniqueness of MGS (based on marker genes): \n";
		//printVectorNice(histo, histThr, "Uniqueness");
		//cout << "Weighted uniqueness of MGS (based on marker genes): \n";
		//printVectorNice(histoWei, histThr, "Uniqueness");
		//for (auto locUni : uniqq) {cout << int(locUni * 100) << " ";}
		printMatrixNice(hist2, histThr, vector<string>({ "Uniqueness","Weighted uniq" }));
		cout << endl;
	}
}

void MAGs::calcLCAconsistency() {

	vector<float> LCAcons(lcaO->maxTaxLvl(), 0.f);
	vector<float> LCAconsMAG(lcaO->maxTaxLvl(), 0.f);
	vector<float> empt(MGSs.size(), 0.f);
	vector<vector<float>> LCAconsMatr(lcaO->maxTaxLvl(), empt);
	vector<vector<float>> LCAconsMatrMAG(lcaO->maxTaxLvl(), empt);
	int cnt(0);
	for (auto mgs : MGSs) {
		mgs->calcLCAconsistency(lcaO);
		vector<float> tmp = mgs->getMGStaxConsis();
		for (size_t i = 0; i < tmp.size(); i++) {
			LCAcons[i] += tmp[i]; LCAconsMatr[i][cnt] = tmp[i];
		}
		tmp = mgs->getMAGtaxConsis();
		for (size_t i = 0; i < tmp.size(); i++) {
			LCAconsMAG[i] += tmp[i]; LCAconsMatrMAG[i][cnt] = tmp[i];
		}
		cnt++;
	}

	for (size_t i = 0; i < LCAcons.size(); i++) {
		LCAcons[i] /= (float) cnt;
	}

	vector<float> histThr({ 1.,0.98,.95,.9,.8,.7,.5,.25,.0 });
	vector<vector<int>> histoM(LCAconsMatr.size());
	for (size_t i = 0; i < LCAconsMatr.size(); i++) {
		histoM[i] = createHist(LCAconsMatr[i], histThr);
	}
	vector<vector<int>> histoMAG(LCAconsMatrMAG.size());
	for (size_t i = 0; i < LCAconsMatrMAG.size(); i++) {
		histoMAG[i] = createHist(LCAconsMatrMAG[i], histThr);
	}


	//also check if taxa names are shared among MGS
	SmplOccur allSpecs; float nonUni(0.),uni(0.);
	for (auto mgs : MGSs) {
		string sp = mgs->getTax(6);
		if (sp == "?") { continue; }
		auto found = allSpecs.find(sp);
		if (found == allSpecs.end()) { allSpecs[sp] = 1; uni++; }
		else { found->second++; nonUni++; }
	}
	float MAGsHitOther(0.);
	for (auto mgs : MGSs) {
		MAGsHitOther += (float)mgs->taxaHit(allSpecs);
	}

	//show taxonomic consisitency (at MGS level)
	vector<string> taxLs = lcaO->getTaxLvls();
	cout << round(100 * (nonUni / (nonUni + uni))) / 100 << "% of species shared amongst MGS," <<
		float(MAGsHitOther) / float(MGSs.size()) << " average MAGs <-> other MGS (1.0 is good).\n";
	cout << "MGS Taxonomic consistency (MAGs are same tax) avg (genus, species): " << round(LCAcons[5] * 100) / 100 << ", " << round(LCAcons[6] * 100) / 100 << "\n";
	int printMStart = 0;
	for (printMStart = 0; printMStart < (histoM.size()-1); printMStart++) {if (histoM[printMStart][0] < MGSs.size()) { break; }}
	printMatrixNice(histoM, histThr, taxLs, printMStart);
	cout << "MAGs (in MGS) tax consistency (MG in MAGs same tax, averaged/MGS):" << endl;
	printMStart = 0; for (printMStart = 0; printMStart < (histoM.size() - 1); printMStart++) { if (histoMAG[printMStart][0] < MGSs.size()) { break; } }
	printMatrixNice(histoMAG, histThr, taxLs, printMStart);
	
	cout << endl;
	
	//write a log file about this
	string outFile = opts->outDir + opts->outTag + ".clusters.tax";
	cout << "Writing tax info to " << outFile << endl;
	ofstream of(outFile);
	if (!of) { cerr << "Couldn't open cluster file " << outFile << endl; exit(71); }
	of << "MGS";
	for (auto tl : taxLs) { of << "\t" << tl; }
	for (auto tl : taxLs) { of << "\tConsis_" << tl; }
	for (auto tl : taxLs) { of << "\tMAGcons_" << tl; }
	of << "\tUniquenessMG\tweightedUniquenessMG\tUniquenessAll\tN_MG\t";
	for (auto tl : taxLs) { of << "\t" << "N_MG_Assi2" <<  tl; }
	
	of << "\n";
	for (auto mgs : MGSs) {
		of << mgs->getID()<< "\t"<< mgs->getTaxStr() << "\t"<< mgs->getTaxConsisStr() << "\t" << mgs->getTaxMAGconsisStr();
		of << "\t" << mgs->getUniqness(0) << "\t" << mgs->getUniqness(1) << "\t" << mgs->getUniqness(2);
		of << "\t" << mgs->totalMG();
		for (size_t i = 0; i < taxLs.size(); i++) {
			of<< "\t" << mgs->totalMGassigned(lcaO,i) ;
		}
		of << "\n";
	}
	of.close();


}


void MAGs::makeMGScentreMAG() {
	int cnt(0);
	for (auto mgs : MGSs) {
		mgs->makeMGScentreMAG();
		cnt++;
	}
	cout << "Found "<< cnt<<" MGS centres\n";
}


bool MAGs::assignTiers() {
	int cnt(0);
	for (auto mag : MAGlist) {
		mag.second->assignTier(opts);
		cnt++;
	}
	cout << "Assigned pre-Tiers to " << cnt<< " MAGs\n";

	return true;
}


void MAGs::transform2GC(MFmap* maps) {

	cout << "Converting MAG contings -> gene cat ids..\n";

	ifstream inf(opts->geneCatIdx);
	if (!inf) {cerr << "Could not open gene cat index file: " << opts->geneCatIdx << "\nAborting\n"; exit(183);}

	if (MAGsbyCtg.size() == 0) { this->createMAGsByCtg(maps); }

	//for (auto mm : MAGsbyCtg) { cerr << mm.second.size() << " "; }cerr << endl;

	std::string line; int cnt (0); int gcnt (0); int gfnd(0);
	std::regex Mremove("M\\d+$");
	SmplOccur FailSamples;
	int FailCnt(0);
	while (std::getline(inf, line)) {

		cnt++;
		if (line.substr(0, 1) == "#") { continue; }
		std::istringstream iss(line);
		string tmp; size_t cc = 0;
		size_t gene(0);
		if (getline(iss, tmp, '\t')) {
			gene = (genetype) stoul(tmp);
		} else { cerr << "Could not read line " << cnt << ": " << line; continue; }
		while (getline(iss, tmp, ',')) {
			tmp = tmp.substr(1);
			size_t found = tmp.find("__");
			if (found == string::npos) {
				cerr << "Gene " << tmp << " is in wrong format! "; continue;
			}
			string smpl = tmp.substr(0, found);
			//size_t xx = tmp.find_last_of("_");
			string ctg = tmp.substr(found + 2, tmp.find_last_of("_") - found - 2);

			//gene & contig & sample found.. now translate this into MAG
			auto SmplFound = MAGsbyCtg.find(smpl);
			if (SmplFound == MAGsbyCtg.end()) {
				//smpl = regex_replace(smpl, Mremove, "");SmplFound = MAGsbyCtg.find(smpl);
				//if (SmplFound == MAGsbyCtg.end()) {
					if (FailCnt < 50) { 
						cerr << "Could not find sample " << smpl << " from gene " << tmp; 
						if (FailCnt == 49) { cerr << "\n..\n"; }
					}
					auto fafi = FailSamples.find(smpl);
					if (fafi == FailSamples.end()) { FailSamples[smpl] = 1; }
					else { fafi->second++; }
					FailCnt++;
					continue;
				//}
			}

			//int dbcnt(0); for (auto c : SmplFound->second) {cerr << c.first << " "; dbcnt++; if (dbcnt > 10) { break; }}
			auto CtgFound = SmplFound->second.find(ctg);
			if (SmplFound->second.end() != CtgFound){
				//all good to log the gene, now get gene number
				int GeneNumOnCtg = stoi(tmp.substr(tmp.find_last_of("_") + 1));
				CtgFound->second->addGene(gene, ctg, GeneNumOnCtg);
				gfnd++;
			}
			gcnt++;
		}
		

	}

	cout << "Read " << cnt << " lines from index file. Selected "<< gfnd<<"/"<< gcnt << " genes\n";
	if (FailCnt) {
		cerr << "\nFailed to find " << FailCnt << " genes->samples links. Samples not found:\n";
		for (auto fafi : FailSamples) {
			cerr << " - " << fafi.first << ": " << fafi.second << " genes\n";
		}
		cerr << "\n\n";
	}
}


float MGS::matchGenes(geneCnter& tar, vector<float>& hits) {

	assert(tar.size() == hits.size());

	float ret(0.f);
	float totalMAGs = (float) members.size();
	int cnt(0);
	for (auto gen:tar) {
		auto found = geneCnt.find(gen.first);
		if (found != geneCnt.end()){//ok no hit..
			//weighted version..
			float weiVal = float(found->second) / totalMAGs;
			ret += weiVal;
			if (weiVal > hits[cnt]) {
				hits[cnt] = weiVal;
			}
		}
		cnt++;
	}
	return ret;
}

void MGS::makeMGScentreMAG() {
	if (centre != nullptr) { return; }
	float totalMAGs = (float)members.size();
	if (totalMAGs == 0) { cerr << "makeMGScentreMAG:::MGS was empty!"; exit(35); }
	if (totalMAGs == 1) { //early abort
		for (auto mag : members) { centre = mag.second; break; }
		centre->setMgsRank(0);
		return; 
	}
	//get some max values for later eval..

	int maxGenes(0); int maxN50(0); int maxCats(0);
	for (auto mag : members) {
		assert(mag.second != nullptr);
		if (mag.second->getFoundGenes() > maxGenes) {
			maxGenes = mag.second->getFoundGenes() ;
		}
		if ((int)mag.second->getN50() > maxN50) {
			maxN50 = (int)mag.second->getN50();
		}
		if ((int)mag.second->getNcats() > maxCats) {
			maxCats = (int) mag.second->getNcats();
		}
	}
	if (maxCats == 0) {maxCats = 1;}
	float centreQ(0.f); float bestSco = 0.f;

	//sort(mag_vec.begin(), mag_vec.end(), [](const MAG* a, const MAG* b) { return(a->compoundScore() < b->compoundScore()); }));
	for ( auto mag : members) {
		geneCnter locGenes = mag.second->getMGgenes();
		assert(locGenes.size() > 0);
		float currSco(0.f); float qSco(0.f);
		for (auto gene : locGenes) {
			auto found = geneCnt.find(gene.first);
			currSco += found->second;
		}
		currSco /= (totalMAGs* maxCats);
		qSco = mag.second->compoundScore(maxGenes, maxN50);
		mag.second->setCentreScore(currSco);
		if ( (currSco+ qSco) > (centreQ+bestSco)) {//new MAG found!
			centre = mag.second;
			bestSco = currSco;
			centreQ = centre->compoundScore(maxGenes, maxN50);
		}
	}
	centre->setMgsRank(0);
}

void MGS::createAllGeneCnts() {

	if (geneCntAll.size() > 0) {return;}

	for (auto mag : members) {
		auto lst = mag.second->getGenes();
		for (auto gn : lst) {
			auto found = geneCntAll.find(gn.first);
			if (found == geneCntAll.end() ) {
				geneCntAll[gn.first] = 1;
			} else {
				found->second++;
			}
		}
	}
}

bool MGS::addMAG(MAG* m) { 
	string id = m->getUniqName();
	auto fnd = members.find(id); 
	if (fnd == members.end()) {
		members[id] = m;
	} else {
		cerr << "MAG added twice to MGS " << MGSid << endl;
		return false;
	}
	//create the jelly core..
	for (auto g : m->getMGgenes()) {
		auto ff = geneCnt.find(g.first);
		if (ff == geneCnt.end()) {
			geneCnt[g.first] = 1;
		} else {
			ff->second++;
		}
	}

	//establish basic tax
	if (majorityTax.size() == 0) {
		majorityTax = m->getMajorityTax();
	} else {
		vector<string> tt = m->getMajorityTax();
		int misses(0);
		for (size_t T = 0; T < tt.size();T++) {
			if (tt[T] != majorityTax[T] && tt[T] != "?") {
				misses++;
			}
		}
		if (misses > 0) {
			cerr << "Adding ill fitting MAG "<< m->getUniqName()<< " to MGS " << this->getID()<<endl;
			for (auto T : tt) { cerr << T << " "; } cerr << " != ";
			for (auto T : majorityTax) { cerr << T << " "; } cerr << endl;
		}
	}

	return true;
}


vector<pair<genetype, size_t>> srtGeneCnt(geneCnter& in) {
	vector<pair<genetype, size_t>> ret; ret.reserve(in.size());
	for (auto& [key, value] : in) {
		ret.emplace_back(std::pair<genetype, size_t>{key, value});
	}
	sort(ret.begin(), ret.end(), [](const pair<genetype, size_t >& a, const pair<genetype, size_t>& b) { return(a.second > b.second); });

	return ret;
}

void printMatrixNice(vector<vector<int>> mat, vector<float> thrs, string Rlab) {
	vector<string> rowL(mat.size(), "");
	for (size_t T = 0; T < mat.size(); T++) {rowL[T] = Rlab + to_string(T);}
	printMatrixNice(mat, thrs, rowL, 0);
}
void printMatrixNice(vector<vector<int>> mat, vector<float> thrs, vector<string> Rlab, size_t start ){
	assert(Rlab.size() == mat.size());
	string tabs; tabs.insert(0, 1 + int(Rlab[0].size() / 8), '\t');
	assert(start < mat.size());
	string spacer; spacer.insert(0, size_t(8 * thrs.size()), '-');
	cout << spacer << "\n"<< tabs;
	for (auto thr : thrs) {
		cout << ">=" << thr << "\t";
	}
	cout << "\n";
	for (size_t T = start; T < mat.size(); T++) {
		cout <<Rlab[T] << "\t";
		for (size_t x = 0; x < mat[T].size(); x++) {
			cout << mat[T][x] << "\t";
		}
		cout << endl;
	}
}

void printVectorNice(vector<int> vec, vector<float> thrs, string descr) {
	bool transformFreq(false);
	string tabs; tabs.insert(0, 1 + int(descr.size() / 8), '\t');
	string spacer; spacer.insert(0, size_t(8 * thrs.size()), '-');
	cout << spacer << "\n"<< tabs;
	for (auto thr : thrs) {
		if (transformFreq) {
			cout << ">=" << round(thr*100) << "%\t";
		}else { cout << ">=" << thr << "\t"; }
	}
	cout << "\n";
	cout << descr << "\t";
	for (size_t T = 0; T < vec.size(); T++) {
		cout << vec[T] << "\t";
	}
	cout << endl;
}

vector<int> createHist(vector<float>uniqq, vector<float> histThr){
	vector<int> histo(histThr.size(), 0);
	for (size_t x = 0; x < uniqq.size(); x++) {
		for (size_t U = 0; U < histThr.size(); U++) {
			if (uniqq[x] >= histThr[U]) {
				histo[U]++; break;
			}
		}
	}
	return histo;
}




int MGS::getMultiBin(genetype g, geneCnter& allGene) {
	int ret(1);
	auto found = allGene.find(g);
	if (found == allGene.end()) {
		return 1;
	} else {
		return found->second;
	}
	return ret;
}

string MGS::writeAllGenesSorted(geneCnter& allGene) {
	//TODO
	stringstream str;
	//first print marker genes
	//	of << "MGS\tGene\tOcc\tMultiCopy\tMultiBin\tisMarkerGene\n";
	vector<pair<genetype, size_t>> MGgenes = srtGeneCnt(geneCnt);
	for (auto gn : MGgenes) {
		int mutlBin = getMultiBin(gn.first, allGene);
		str << MGSid << "\t" << gn.first << "\t" << gn.second << "\t0\t" << mutlBin << "\t1\n";
	}

	//overwrite old object, sort rest of genes..
	MGgenes = srtGeneCnt(geneCntAll);
	for (auto gn : MGgenes) {
		int mutlBin = getMultiBin(gn.first, allGene);
		str << MGSid << "\t" << gn.first << "\t" << gn.second << "\t0\t" << mutlBin << "\t0\n";
	}

	return str.str();
}
float MGS::taxaHit(SmplOccur& OT,int lvl) {
	float cnter(0.f);
	if (MAGtaxAssigns.size() < lvl) { return cnter; }
	for (auto species : MAGtaxAssigns[lvl]) {
		auto found = OT.find(species.first);
		if (OT.end() != found) { cnter++; }
	}
	return cnter;
}
string MGS::getTaxConsisStr() {
	stringstream ret("");
	assert(taxConsistency.size() > 0);
	ret << taxConsistency[0];
	for (size_t i = 1; i < taxConsistency.size(); i++) {
		ret << "\t" << taxConsistency[i];
	}
	return ret.str();
}

string MGS::getTaxMAGconsisStr() {
	stringstream ret("");
	assert(MAGtaxConsistency.size() > 0);
	ret << MAGtaxConsistency[0];
	for (size_t i = 1; i < MAGtaxConsistency.size(); i++) {
		ret << "\t" << MAGtaxConsistency[i];
	}
	return ret.str();
}


string 	MGS::getTaxStr() {
	string ret("");
	assert(majorityTax.size()>0);
	ret = majorityTax[0];
	for (size_t i = 1; i < majorityTax.size();i++) {
		ret += "\t" + majorityTax[i];
	}
	return ret;
}
int MGS::totalMGassigned(MG_LCA* lca,int lvl) { 
	int cnt = 0;
	for (auto gn : geneCnt) {
		if (lca->hasTax(gn.first, lvl)) {
			cnt++;
		}
	}
	return cnt;
}


void MGS::calcLCAconsistency(MG_LCA* locLCA) {
	//vector<float> ret(locLCA->maxTaxLvl(),0.f);
	if (taxConsistency.size() == 0) { taxConsistency.resize(locLCA->maxTaxLvl(),0.f); }
	if (MAGtaxConsistency.size() == 0) { MAGtaxConsistency.resize(locLCA->maxTaxLvl(),0.f); }
	if (majorityTax.size() == 0) { majorityTax.resize(locLCA->maxTaxLvl(), "?"); }
	if (MAGtaxAssigns.size() == 0) { MAGtaxAssigns.resize(locLCA->maxTaxLvl()); }

	
	float TfracG(0.f), TfracS(0.f); string Tgenus(""), Tspec("");
	for (auto mag : members) {
		//at genus lvl
		for (size_t k = 0; k < taxConsistency.size(); k++) {
			string tmp(""); float tmp2(0.f);
			mag.second->getTopTaxa(k, tmp2, tmp);
			MAGtaxConsistency[k] += tmp2;
			auto fnd = MAGtaxAssigns[k].find(tmp);
			if (fnd == MAGtaxAssigns[k].end()) { MAGtaxAssigns[k][tmp] = 1; }
			else { fnd->second++; }
		}
		//at species lvl
		//mag.second->getTopTaxa(6, TfracS, Tspec);
	}
	/*ret[5] = TfracG; ret[6] = TfracS;
	majorityTax[5] = Tgenus; majorityTax[6] = Tspec;
	*/
	for (size_t k = 0; k < taxConsistency.size(); k++) {
		MAGtaxConsistency[k] /= (float)members.size();
	}
	for (size_t k = 0; k < MAGtaxAssigns.size(); k++) {
		string bestT; int bestV(0),total(0);
		for (auto tax : MAGtaxAssigns[k]) {
			if (tax.first == "?") { continue; }
			if (tax.second > bestV) {
				bestV = tax.second; bestT = tax.first;
			}
			total+= tax.second;
		}
		majorityTax[k] = bestT;
		taxConsistency[k] = bestV / total;
	}

	//return taxConsistency;
}

string MGS::writeMGSobservations() {
	string str(""); 
	assert(centre != nullptr);
	str += this->MGSid + "\t" + to_string(members.size()) +
		"\t" + to_string(centre->qualTier()) + "\t"  ;
	for (auto mag : members) {
		str += mag.second->getUniqName() + ",";
	}
	str = str.substr(0, (str.length() - 1));
	return str;

}

string MGS::writeMGSscore() {
	string str(""); string cmModel(".cm2");
	str += this->MGSid + "\t" + to_string(centre->getComple()) +
		"\t" + to_string(centre->getConta()) + "\t" + cmModel +
		"\t" + centre->getUniqName();
	return str;
}


void MAGs::writeClusters() {
	string outFile = opts->outDir + opts->outTag+".clusters";
	cout << "Writing clusters to " << outFile << endl;
	for (auto mgs : MGSs) {
		add2allGenes(mgs->getAllGenes());
	}
	cout << "Created all gene list\n";
	ofstream of(outFile);
	if (!of) { cerr << "Couldn't open cluster file " << outFile << endl; exit(71); }
	of << "MGS\tGene\tOcc\tMultiCopy\tMultiBin\tisMarkerGene\n";
	for (auto mgs : MGSs) {
		of << mgs->writeAllGenesSorted(allGenes);
	}
	of.close();
}

void MAGs::writeClusterObservations(){
	string outFile = opts->outDir + opts->outTag + ".clusters.obs";
	cout << "Writing clusters observations to " << outFile << endl;
	ofstream of(outFile);
	if (!of) { cerr << "Couldn't open cluster file " << outFile << endl; exit(71); }
	of << "Bin\tObservations\tQualTier\tMembers\n";
	int cnt(0);
	for (auto mgs : MGSs) {
		if (mgs == nullptr) { cerr << "\n\nEmpty MGS found!!\n\n"; continue; }
		of << mgs->writeMGSobservations() << endl;
		cnt++;
	}
	of.close();
}

void MAGs::writeClusterScores() {
	string outFile = opts->outDir + opts->outTag + ".clusters.cm2";
	cout << "Writing clusters scores to " << outFile << endl;
	ofstream of(outFile);
	if (!of) { cerr << "Couldn't open cluster file " << outFile << endl; exit(71); }
	of << "Name\tCompleteness\tContamination\tCheckMmodel\tOrigin\n";
	for (auto mgs : MGSs) {
		of << mgs->writeMGSscore() << endl;
	}
	of.close();
}


void MAGs::writeAllMAGs(MG_LCA* lca) {
	string outFile = opts->outDir + "MAGvsGC.txt";
	cout << "Writing MAGs to " << outFile << endl;

	vector<string> taxLs = lcaO->getTaxLvls();


	ofstream of(outFile);
	if (!of) { cerr << "Couldn't open MAG_MGS_GC file " << outFile << endl; exit(71); }
	of << "MAG\tMGS\tRepresentative4MGS\tMatch2MGS\tUniqueness\tAssociatedMGS\tCompleteness\tContamination\tLCAcompleteness\tN50\t";
	of << "N_Genes\tGC\tCodingDensity\tCentreScore\tCompoundScore";
	for (auto tl : taxLs) { of << "\t" << tl; }

	vector<string> cats = lca->getCats(true);
	for (auto cat : cats) {
		of << "\t" << cat;
		//cats.push_back(cat);
	}
	of <<"\tother_Genes"<< endl;

	for (auto mgs : MGSs) {
		MAG* cntr = mgs->getCentreMAG();
		of << cntr->formatMAG(cats, lcaO,true,true) << endl;
		MAGlst mags = mgs->getMAGs();

		//get sorted list
		std::vector<std::pair<MAG*, float>> mag_vec;mag_vec.reserve(mags.size());
		for (auto& [key, value] : mags) {
			mag_vec.emplace_back(std::pair<MAG*, float>{value, (value->getCompoundScore() + value->getCentreScore())});
		}
		sort(mag_vec.begin(), mag_vec.end(), [](const pair<MAG*, float>& a, const pair<MAG*, float>& b) { return(a.second > b.second); });


		for (auto mag : mag_vec) {
			if (mag.first==nullptr || mag.first == cntr){continue;}
			of << mag.first->formatMAG(cats, lcaO,false,true) << endl;
		}

	}

	for (auto mag : MAGlist) {
		MAG* mg = mag.second;
		if (mg == nullptr || mg->getAssignedMGS() != nullptr) { continue; }
		of << mg->formatMAG(cats, lcaO,false,true)<<endl;
	}

	of.close();
}

// output of NGS centres... 
void MAGs::writeClusterCentres() {

}
