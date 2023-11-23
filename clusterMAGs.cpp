#include "clusterMAGs.h"
//5.11.23: started clusterMAGs v0.1 @Falk Hildebrand
//11.11.23: first release v0.2
//23.11.23: v.21: add association to "best" MGS file

const char* version="0.21";


void stateVersion(){
    printf("clusterMAGs v%s\n", version);
}

void versionMsg() {
	stateVersion();
	exit(0);
}

void announce() {
	cout << "===========================\n";
	stateVersion();
	cout << "===========================\n" << endl;
}

void helpMsg(){
	stateVersion();
	printf("\n");
	printf("USAGE\n");
	printf("    clusterMAGs -geneCatIdx <file> -map <file> -BinDir <dir> -MGdir <dir> -tmp <dir> [options] \n");
	printf("\n");
	printf("OPTIONS\n");
	printf("Files, naming, operations conventions ---------------\n");
	printf("    -map                 Path to gene catalog\n");
	printf("    -outDir              Path to a output directory\n");
	printf("    -tmp                 Depth or multiple comma seperated depths to rarefy to. Default is 0.95 times the minimal column sum.\n");
	printf("    -LCAdir              Path to Marker Gene (MG) LCAs.\n");
	printf("    -MGfile              File listing MG_category\\tnum_genes_in_cat\\tcsv_MGs.\n");
	printf("    -canopyDir           Path to reference Canopy dir.\n");
	printf("    -geneCatIdx          Path to gene catalog index file (transcribing assembly genes to gene cat genes).\n");
	
	printf("    -CMsuffix            Suffix of checkM2 files Default: \".cm2\".\n");
	printf("    -FILEtag             Prefix for output file Default: \"SB.\".\n");
	printf("    -MGStag              Tag for clustered metageonmic species (MGS). Default: \"MGS.\"\n");
	printf("    -t                   Number of threads to use. Default: 1\n");

	printf("Quality criteria ---------------\n");
	printf("    -minMGassigned2tax   Minimal number of genes assigned to same tax to accept MAG taxonomy. Default: 10\n");
	printf("    -processLowQual      Process MAGs that are below minimal inclusion criteria? Default: false\n");

	printf("    -useCheckm1          Use CheckM1 MAG quals? Default: false\n");
	printf("    -useCheckm2          Use CheckM1 MAG quals? Default: true\n");
	printf("Misc ---------------\n");
	printf("    -h                   Show help\n");
	printf("    -v                   Show version\n");
	printf("    -help_output         Help file regarding output files from clusterMAGs\n");
	printf("    -show_quals          Print used quality criteria.\n");
	printf("\n");
	printf("Author: Falk Hildebrand, QIB/EI, (c) 2023\n");

	printf("\n");

	std::exit(0);
}


bool isGZfile(const std::string fi) {
	std::string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}





vector<double> parseDepths(string a){
    std::vector<double> vect;
    std::stringstream ss(a);

    float i;

    while (ss >> i)
    {
        vect.push_back(i);
        cout << i << " ";
        if (ss.peek() == ',')
            ss.ignore();
    }


    return vect;
}



options::options(int argc, char** argv) :map(""), outDir(""), tmp(""),
LCAdir(""), canopyFile(""), geneCatIdx(""), 
path2Bins("Binning/SB/"), path2Assmbl("assemblies/metag/"), file2Assmbl("assembly.txt"),
threads(10), useCheckm1(false), useCheckm2(true),
complTiers(0), contaTiers(0), LCAcomplTier(0), allowNovelMGSTiers(0), match2MGS_Tier(0), globalUniq4MGS_Tier(0),
maxTierRds(4), MGSprefix("MGS."), outTag("SB"), qualSuffix(".cm2"), processOnlyHQ(true), minMGassigned2tax(10){


    bool hasErr = false;
	bool showQs(false);


	if (argc == 0) { 
		cerr << "Not enough options given to clusterMAGs"; exit(23); 
	}//

	if (!strcmp(argv[1] ,"-h") || !strcmp(argv[1] ,"--help")) {
		helpMsg();
	}else if (!strcmp(argv[1] , "version") || !strcmp(argv[1] , "-version") || !strcmp(argv[1], "-v") || !strcmp(argv[1] ,"--version")) {
		versionMsg();

	} else if (!strcmp(argv[1], "-show_quals")) {
		showQs = true;
	} else {
		announce();

	}

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-map"))
			map = argv[++i];
		else if (!strcmp(argv[i], "-outDir"))
			outDir = argv[++i];
		///else if (!strcmp(argv[i], "-m"))
		else if (!strcmp(argv[i], "-tmp")) {
			tmp = argv[++i];
		}else if (!strcmp(argv[i], "-LCAdir")) {
			LCAdir = argv[++i];
		}else if (!strcmp(argv[i], "-MGfile")) {
			GTDBfile = argv[++i];
		}else if (!strcmp(argv[i], "-FILEtag")) {
			outTag = argv[++i];
		}else if (!strcmp(argv[i], "-CMsuffix")) {
			qualSuffix = argv[++i];
		}else if (!strcmp(argv[i], "-MGStag")) {
			MGSprefix = argv[++i];
		}else if (!strcmp(argv[i], "-geneCatIdx")) {
			geneCatIdx = argv[++i];
		} else if (!strcmp(argv[i], "-canopyDir")) {
			canopyFile = argv[++i];
		}else if (!strcmp(argv[i], "-t")) {
			threads = atoi(argv[++i]);
		}else if (!strcmp(argv[i], "-minMGassigned2tax")) {
			minMGassigned2tax = atoi(argv[++i]);
		}else if (!strcmp(argv[i], "-processLowQual")) {  // no swap
			processOnlyHQ = false;
		}else if (!strcmp(argv[i], "-useCheckm1")) {  // no swap
			useCheckm1 = true; useCheckm2 = false;
		}else if (!strcmp(argv[i], "-useCheckm2")) {  // no swap
			useCheckm2 = true; useCheckm1 = false;
		}

    }
	

    // sanity checks
    // we need input
	if (!showQs) {

		if (map == "") {//just set some defaults
			cerr << "-map must be specified\n";
			hasErr = true;
		}
		if (outDir == "") {//just set some defaults
			cerr << "-outDir must be specified\n";
			hasErr = true;
		}
		if (LCAdir == "") {//just set some defaults
			cerr << "-LCAdir must be specified\n";
			hasErr = true;
		}


		if (hasErr) {
			cerr << "Error in option parsing.\nUse \"clusterMags -h\" to get full help.\n";
			exit(98);
		}
	}

/*	complTiers.resize(maxTierRds, 0.f);
	contaTiers.resize(maxTierRds, 0.f);
	LCAcomplTier.resize(maxTierRds, 0.f);
	allowNovelMGSTiers.resize(maxTierRds, false);
	overlapTiers.resize(maxTierRds, 0.f);
	unknwnTiers.resize(maxTierRds, 0.f);*/

	complTiers.insert(complTiers.end(), { 98,95,80,60 });
	contaTiers.insert(contaTiers.end(), { 1,3,5,10 });
	LCAcomplTier.insert(LCAcomplTier.end(), { .9,.8,.6,.5 });


	//in the cluster rounds, which Tier can make a new MGS?
	allowNovelMGSTiers.insert(allowNovelMGSTiers.end(), { 0,1,2,2});
	//how much of marker genes should match in order to merge MAGs into MGS?
	match2MGS_Tier.insert(match2MGS_Tier.end(), { 0.9,0.85,0.8,0.75 });
	//how much of marker genes should not be in other MAGs?
	globalUniq4MGS_Tier.insert(globalUniq4MGS_Tier.end(), { 0.9,0.85,0.8, 0.8 });

	maxMatch4NovelMGS_Tier.insert(maxMatch4NovelMGS_Tier.end(), { 0.3,0.3,0.3,0.5 });

	if (showQs) {
		show_quals();
	}


}

void options::print_details(){

    stateVersion();
    // print run mode:
    cout << "------------------------------------ "  << std::endl;
    cout << "Run information:" << std::endl;
    cout << "output dir:    " << outDir << std::endl;
    cout << "map:           " << map << std::endl;
	cout << "MG dir:     " << LCAdir << std::endl;
	cout << "canopy dir:     " << canopyFile << std::endl;
	cout << "threads:     " << threads << std::endl;
	cout << "------------------------------------ " << std::endl;

    //cout << "mode:           " << mode  << std::endl;
    cout << std::endl;

}


void options::show_quals() {
	stateVersion();
	cout << "------------------------------------ " << std::endl;
	cout << "Quality criteria for processing MAGs:\n";
	cout << " - complTiers: "; for (auto x : complTiers) { cout << x << " "; }cout << endl;
	cout << " - contaTiers: "; for (auto x : contaTiers) { cout << x << " "; }cout << endl;
	cout << " - LCAcomplTier: "; for (auto x : LCAcomplTier) { cout << x << " "; }cout << endl;
	cout << " - allowNovelMGSTiers: "; for (auto x : allowNovelMGSTiers) { cout << x << " "; }cout << endl;
	cout << " - match2MGS_Tier: "; for (auto x : match2MGS_Tier) { cout << x << " "; }cout << endl;
	cout << " - globalUniq4MGS_Tier: "; for (auto x : globalUniq4MGS_Tier) { cout << x << " "; }cout << endl;
	cout << " - maxMatch4NovelMGS_Tier: "; for (auto x : maxMatch4NovelMGS_Tier) { cout << x << " "; }cout << endl;
	cout << "------------------------------------ " << std::endl;

	cout << std::endl;

	exit(0);

}





int main(int argc, char* argv[])
{

	if (argc < 2) { cerr << "Not enough arguments. Use \"clusterMAGs -h\" for getting started.\n"; exit(3); }

	//clock_t tStart = clock();
	Benchmark* benchmark = new Benchmark("Time clusterMAGs: ");
	benchmark->start();

	options* opts = new options(argc, argv);

	MG_LCA* mglca = nullptr;
	mglca = new MG_LCA(opts->LCAdir);
	mglca->addMGs(opts->GTDBfile);
	benchmark->now_total_time();

	MFmap* inmaps = new MFmap(opts->map);
	benchmark->now_total_time();

	MAGs* mags = new MAGs(inmaps, opts);
	//benchmark->now_total_time();
	//read additional MAGs from canopies
	mags->readCanopies();

	//mags->assignTiers();
	//benchmark->now_total_time();

	mags->createMAGsByCtg(inmaps);
	benchmark->now_total_time();


	mags->transform2GC(inmaps);
	benchmark->now_total_time();

	mags->addLCAscores(mglca);
	benchmark->now_total_time();

	mags->report();
	//mags->removeCtgInfo();

	//---------------------------------------------
	mags->cluster2MGS();
	benchmark->now_total_time();
	//---------------------------------------------

	mags->makeMGScentreMAG();
	benchmark->now_total_time();

	mags->evalRemainderMAGs();
	mags->calcLCAconsistency();
	benchmark->now_total_time();

	//MAGvsGC.txt
	mags->writeAllMAGs( mglca);
	benchmark->now_total_time();

	//SB.clusters
	mags->writeClusters();
	// SB.clusters.cm2
	mags->writeClusterScores();
	//SB.clusters.obs
	mags->writeClusterObservations();

	delete mags, inmaps, mglca, opts;

	// all done.. time for the benchmark
	
	cout << "Finished clusterMAGs algorithm\n";
	benchmark->stop();
	benchmark->printResults();
	delete benchmark;

	return 0;
}
