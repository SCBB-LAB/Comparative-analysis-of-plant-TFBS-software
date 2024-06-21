//g++ cmf.cpp -o cmf 
//./cmf -w 7 -m 2 -d 1 -l 5 -u 20 -x 3 -t 50 -i1 ../MotifScan/Seqs/Oct4ESSeqs4Scanning.txt -i2 ../MotifScan/Seqs/ControlSeqs500.txt -o Oct4test.txt -f Oct4ES/
// basic file operations
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "cmf.h" 
#include <sys/stat.h>
using namespace std;

// some global variables
double X2cutoff      = 16.27; // chiSquared values for 0.001 alpha level (degrees of freedom = 3)
int    lowerBound    = 5;     // motif length variables
int    upperBound    = 20;    //   '    '
double fdrCut        = 0.667; // fdr cutoff
 
int main (int argc, char **argv) 
{
  int w   = 7;
  int m   = 2;
  int MAX = 10;
  int d   = 1;
  int checkBaseLR = 0; 
  FILE * pFile;
  const char *fastaFile1;
  const char *fastaFile2;
  string pwmFolder = "defaultCMFfolder";
  string seedFile  = "defaultCMFseedfile.txt";

  // get cmd line parameters
  if(argc == 1)
  {
    cout << "-w,  the length of the motif seed" << endl;
    cout << "-m,  numer of mismatches in seed" << endl;
    cout << "-F,  fdr cutoff" << endl;
    cout << "-t,  the number of top seeds to test" << endl;
    cout << "-i1, first set of sequences" << endl;
    cout << "-i2, second set of sequences" << endl;
    cout << "-d,  1: enrichment only in i1, 2: enrichment in both datasets" << endl;
    cout << "-l,  lower bound on length of motifs" << endl;
    cout << "-u,  upper bound on length of motifs" << endl;
    cout << "-o,  output of seed statistics" << endl;
    cout << "-f,  folder for all other output (must exist before running)" << endl;
    cout << "-c,  if > 0 filter out seeds based on cg content (ex -c 4, filter out seeds with more than 4 C or G's)" << endl;

    return(0);
  }

  for(int i=0;i<argc;i++) 
    { 
    if(argv[i][0]=='-') 
      { 
	if(argv[i][1]=='w')      { w             = atoi(argv[i+1]);         if(w < 5 || w > 20){cout << "w must be 5-20" << endl;  return(0);}
                                                                       else if(w > 8){cout << "WARNING: w > 8 may be computationally infeasible" << endl;}} 
	else if(argv[i][1]=='m') { m             = atoi(argv[i+1]);         if(m < 1 || m > 2){cout << "m must be 1 or 2" << endl; return(0);}} 
	else if(argv[i][1]=='F') { fdrCut        = atof(argv[i+1]);         if(fdrCut < .001 || fdrCut > .999){cout << "unreasonable FDR cutoff" << endl; return(0);}} 
	else if(argv[i][1]=='t') { MAX           = atoi(argv[i+1]);         if(MAX > 50){cout << "WARNING: large t may be infeasible" << endl; }} 
	else if(argv[i][1]=='d') { d             = atoi(argv[i+1]);         if(d < 1 || d > 2){cout << "d must be 1 or 2" << endl; return(0);}} 
        else if(argv[i][1]=='f') { pwmFolder     = argv[i+1];               mkdir(pwmFolder.c_str(), 77777);} 
        else if(argv[i][1]=='l') { lowerBound    = max(atoi(argv[i+1]),5);  } 
        else if(argv[i][1]=='u') { upperBound    = min(atoi(argv[i+1]),20); } 
        else if(argv[i][1]=='o') { seedFile      = argv[i+1];               } 
        else if(argv[i][1]=='c') { checkBaseLR   = atoi(argv[i+1]);         }
        else if(argv[i][1]=='i') 
        { 
	if(argv[i][2]=='1')     { fastaFile1 = argv[i+1]; }
	else if(argv[i][2]=='2'){ fastaFile2 = argv[i+1]; }
	} 
      } 
    } 
  // check last few inputs
  if(lowerBound > upperBound || lowerBound > w){cout << "conflicting value for motif length lower bound (-l)" << endl; return(0);}
  if(lowerBound > upperBound || upperBound < w){cout << "conflicting value for motif length upper bound (-u)" << endl; return(0);}
  if(checkBaseLR < 0 | checkBaseLR > w){cout << "-c must be 0-"<< w << endl; return(0);}
  pFile = fopen(seedFile.c_str(),"w");  

  // get sequences    
  mapTypeString fastaSeqs1 = getFastaSeqs(fastaFile1); if(fastaSeqs1.empty()){return(0);}
  mapTypeString fastaSeqs2 = getFastaSeqs(fastaFile2); if(fastaSeqs2.empty()){return(0);}
  
  // get background counts (just like wmer but w=1 for marginals and w=2 for one degree Markov Chain )
  mapTypeInt marginalCounts1 = countWmers(fastaSeqs1,1,0,(int)fastaSeqs1.size(),0);
  mapTypeInt marginalCounts2 = countWmers(fastaSeqs2,1,0,(int)fastaSeqs2.size(),0);

  double CGdif = singleDifOfProps(marginalCounts1["C"] + marginalCounts1["G"], marginalCounts2["C"] + marginalCounts2["G"],  marginalCounts1["total"],  marginalCounts2["total"]); 

  cout << "t-score for difference of CG content is " << CGdif << ", with p1 = " << (marginalCounts1["C"] + marginalCounts1["G"]) / (double)marginalCounts1["total"] << " and p2 = "  << (marginalCounts2["C"] + marginalCounts2["G"]) / (double)marginalCounts2["total"] << endl;

  mapTypeInt oneDegMCCounts1 = countWmers(fastaSeqs1,2,0,(int)fastaSeqs1.size(),0);
  mapTypeInt oneDegMCCounts2 = countWmers(fastaSeqs2,2,0,(int)fastaSeqs2.size(),0);
  mapTypeInt marginalCountsT; 
  mapTypeInt oneDegMCCountsT; 
  for(mapTypeInt::const_iterator it = marginalCounts1.begin(); it != marginalCounts1.end(); ++it){marginalCountsT[it->first] = marginalCounts1[it->first] + marginalCounts2[it->first];}
  for(mapTypeInt::const_iterator it = oneDegMCCounts2.begin(); it != oneDegMCCounts2.end(); ++it){oneDegMCCountsT[it->first] = oneDegMCCounts2[it->first] + oneDegMCCounts1[it->first];}  

  // make marginal and 1deg MC freqs from counts
  struct4mapTypeDouble freqs, freqs1, freqs2, temp;
  for(mapTypeInt::const_iterator it = marginalCountsT.begin(); it != marginalCountsT.end(); ++it){ if(it->first != "total"){freqs.structMap[it->first] = log(marginalCountsT[it->first]/(double)marginalCountsT["total"]);} }  
  for(mapTypeInt::const_iterator it = oneDegMCCountsT.begin(); it != oneDegMCCountsT.end(); ++it){ if(it->first != "total"){freqs.structMap[it->first] = log(oneDegMCCountsT[it->first]/(double)marginalCountsT[it->first.substr(0,1)]); }}  

  for(mapTypeInt::const_iterator it = marginalCounts1.begin(); it != marginalCounts1.end(); ++it){ if(it->first != "total"){freqs1.structMap[it->first] = log(marginalCounts1[it->first]/(double)marginalCounts1["total"]);} }  
  for(mapTypeInt::const_iterator it = oneDegMCCounts1.begin(); it != oneDegMCCounts1.end(); ++it){ if(it->first != "total"){freqs1.structMap[it->first] = log(oneDegMCCounts1[it->first]/(double)marginalCounts1[it->first.substr(0,1)]); }}  

  for(mapTypeInt::const_iterator it = marginalCounts2.begin(); it != marginalCounts2.end(); ++it){ if(it->first != "total"){freqs2.structMap[it->first] = log(marginalCounts2[it->first]/(double)marginalCounts2["total"]);} }  
  for(mapTypeInt::const_iterator it = oneDegMCCounts2.begin(); it != oneDegMCCounts2.end(); ++it){ if(it->first != "total"){freqs2.structMap[it->first] = log(oneDegMCCounts2[it->first]/(double)marginalCounts2[it->first.substr(0,1)]); }}  

  // count wmers in each set of sequences 
  mapTypeInt wmerCounts1     = countWmers(fastaSeqs1,w,m,(int)fastaSeqs1.size(),1); cout << "total w-mer count in grp1 is " << wmerCounts1["total"] << endl;
  mapTypeInt wmerCounts2     = countWmers(fastaSeqs2,w,m,(int)fastaSeqs2.size(),1); cout << "total w-mer count in grp2 is " << wmerCounts2["total"] << endl;
  mapTypeInt wmerCountsT; 
  for(mapTypeInt::const_iterator it = wmerCounts1.begin(); it != wmerCounts1.end(); ++it){wmerCountsT[it->first] = wmerCounts1[it->first] + wmerCounts2[it->first];}
  for(mapTypeInt::const_iterator it = wmerCounts2.begin(); it != wmerCounts2.end(); ++it){wmerCountsT[it->first] = wmerCounts2[it->first] + wmerCounts1[it->first];}  
  mapTypeDouble wmerDiffProp = diffOfProps(wmerCounts1, wmerCounts2, wmerCountsT, wmerCounts1["total"], wmerCounts2["total"]);

  mapTypeVectorInt allPosSites = findWmerSites(fastaSeqs1, w);
  cout << "done finding all " << w << "-mers in group 1" << endl;
  mapTypeVectorInt allNegSites = findWmerSites(fastaSeqs2, w);
  cout << "done finding all " << w << "-mers in group 2" << endl;

  int tL1 = wmerCounts1["total"]/2;
  int tL2 = wmerCounts2["total"]/2;

  ///////////////////////////////////////////////////////
  // step 1 create seeds using w choose m neighbors    //
  ///////////////////////////////////////////////////////
  cout << "begin seeding" << endl;
  struct4mapTypeWmerNeighborhood neighborhoods;
  struct4mapTypeInt              wmerCnts1;     wmerCnts1.structMap         = wmerCounts1; 
  wmerNeighborhood               holder;        wmerNeighborhood *ptrholder = &holder;

  int CGcount;
  for(mapTypeDouble::const_iterator it = wmerDiffProp.begin(); it != wmerDiffProp.end(); ++it)
  {
    CGcount =0;
    for(int i=0;i<w;i++){if(it->first.substr(i,1) == "C" || it->first.substr(i,1) == "G"){CGcount++;}}
    if(checkBaseLR > 0 && CGcount > checkBaseLR){wmerCounts1[it->first] = 0; wmerCounts2[it->first] = 0; wmerDiffProp[it->first] = 0;}
  }
  if(checkBaseLR == 0){freqs1 = freqs; freqs2 = freqs;}

  struct4mapTypeInt              wmerCnts2;     wmerCnts2.structMap         = wmerCounts2; 
  struct4mapTypeDouble           wmerDifProps;  wmerDifProps.structMap      = wmerDiffProp; 

  int counter = 0; double lastPerc = 0;
  fprintf(pFile,"wmer\tcnts1\tcnts2\twmerDiffProp\twmerEnrichment\tnhoodCnts1\tnhoodCnts2\tnhoodDiffProp\tadjCnts1\tadjCnts2\tadjLen1\tadjLen2\tnhoodAdjustedDiffProp\tnhoodEnrichment\tnhoodMaxFlag\tneighborCnts\n");

  for(mapTypeInt::const_iterator it = wmerCountsT.begin(); it != wmerCountsT.end(); it++)
  {
    if(it->first != "total")
    {
      string key = it->first;

      if(abs(wmerDiffProp[key]) > 0)
      {
        int *positions = new int[m]; 
        ptrholder->name                  = key;
        ptrholder->width                 = w;
        ptrholder->difPropBestMisMatches = wmerDiffProp[key]/10;
        ptrholder->difPropBestMMadjusted = wmerDiffProp[key]/10;
        ptrholder->difPropExactMatch     = wmerDiffProp[key];
        ptrholder->difPropUpdated        = 0; ptrholder->threshold = log((double)2000); ptrholder->pb = (double)fastaSeqs1.size() / ((double)tL1 * 2); 
        ptrholder->counter               = 0; 
        ptrholder->maxFlag               = 1;
        ptrholder->exactEnrichment       = log2(((double)wmerCounts1[key]/wmerCounts1["total"]) / ((double)wmerCounts2[key]/wmerCounts2["total"]));; 
        ptrholder->bestEnrichment        = 0; ptrholder->tempEnrichment = 0;
	ptrholder->totLen1               = tL1;
	ptrholder->totLen2               = tL2;

        ptrholder->tempCnts1.clear();     ptrholder->tempCnts1["tot"] = 0;
        ptrholder->tempCnts2.clear();     ptrholder->tempCnts2["tot"] = 0; 
        ptrholder->neighborCnts1.clear(); ptrholder->neighborCnts1["tot"] = 0;
        ptrholder->neighborCnts2.clear(); ptrholder->neighborCnts2["tot"] = 0; 
  
        cycleThruWchooseM(positions, key, 0, m, m, &wmerCnts1, &wmerCnts2, ptrholder, &wmerDifProps, &allPosSites, &allNegSites);  
        ptrholder->convergedFlag         = 0;
        ptrholder->indexOffset           = 0;

        for(int j = 0; j < 0; j++){ for(int k = 0; k<4; k++){ptrholder->pwm[j][k] = 0;} }

        neighborhoods.structMap[key]     = *ptrholder;
        delete [] positions;
      }

      double perc = (double)counter/wmerCountsT.size(); if(perc - lastPerc > .1){lastPerc = perc; cout << perc << endl;}  
      counter++;
    }
  }

  double *exactDiffP = new double[(int)   neighborhoods.structMap.size()];
  double *nhoodDiffP = new double[(int)   neighborhoods.structMap.size()];
  string *keys1      = new string[(int)  neighborhoods.structMap.size()];
  string *keys2      = new string[(int)  neighborhoods.structMap.size()];
  int    index       = 0;

  for(mapTypeWmerNeighborhood::const_iterator it = neighborhoods.structMap.begin(); it != neighborhoods.structMap.end(); ++it)
  {
    string key = it->first;
    ptrholder = &neighborhoods.structMap[key]; 

    for(mapTypeInt::const_iterator it2 = ptrholder->neighborCnts1.begin(); it2 != ptrholder->neighborCnts1.end(); it2++)
    { 
      if(neighborhoods.structMap[it2->first].difPropBestMMadjusted == ptrholder->difPropBestMMadjusted && abs(wmerDiffProp[it2->first]) > abs(wmerDiffProp[key])){ptrholder->maxFlag = 0;}
      if(neighborhoods.structMap[reverseStrand(it2->first)].difPropBestMMadjusted == ptrholder->difPropBestMMadjusted && abs(wmerDiffProp[reverseStrand(it2->first)]) > abs(wmerDiffProp[key])){ptrholder->maxFlag = 0;}
    }
    
    if(it->second.maxFlag == 1)
    { 
      keys1[index] = it->first; keys2[index] = it->first; 
      exactDiffP[index] = wmerDiffProp[key]; 
      nhoodDiffP[index] = it->second.difPropBestMMadjusted; 
      index++; 

      if(it->first != reverseStrand(it->first)){neighborhoods.structMap[reverseStrand(it->first)].maxFlag =0;}
    }

    fprintf(pFile,"%s\t%i\t%i\t%5.5f\t%5.5f\t%i\t%i\t%5.5f\t%5.5f\t%5.5f\t%i\t%i\t%5.5f\t%5.5f\t%i\t",key.c_str(),wmerCounts1[key],wmerCounts2[key],wmerDiffProp[key],ptrholder->exactEnrichment,ptrholder->neighborCnts1["tot"],ptrholder->neighborCnts2["tot"],ptrholder->difPropBestMisMatches,ptrholder->adjCnts1,ptrholder->adjCnts2,ptrholder->nhoodPosSites.size(),ptrholder->nhoodNegSites.size(),ptrholder->difPropBestMMadjusted,ptrholder->bestEnrichment,ptrholder->maxFlag); 
    for(mapTypeInt::const_iterator it2 = ptrholder->neighborCnts1.begin(); it2 != ptrholder->neighborCnts1.end(); it2++){ fprintf(pFile,"%s(%i, %i), ",it2->first.c_str(),it2->second,ptrholder->neighborCnts2[it2->first]); } fprintf(pFile,"\n"); 
  }
  fclose(pFile);

  index--;  cout << "done with neighbor seeding (" << index << " seeds)" << endl;

  vector<string> seedsV;
  int flag, j;   

  cout << "end seeding" << endl << "update motifs" << endl;;

  //////////////////////////////////////////////////
  // step 2, update best MAX seeds into motifs    //
  //////////////////////////////////////////////////
  bubble_sort(nhoodDiffP, keys1, index);
  bubble_sort(exactDiffP, keys2, index);

  for(j = 0; j < MAX; j++) // positive group ranked by nhood diffProp
  { 
    flag = updateSeed(&neighborhoods, fastaSeqs1, fastaSeqs2, &neighborhoods.structMap[keys1[j]], &freqs1, &freqs2, tL1, tL2, seedsV, 0); if(flag == 1){seedsV.push_back(keys1[j]);} cout << j+1 << endl; 
  }
  if(d == 2) // negative group ranked by nhood diffProp
  {
    for(j=index;j>(index-MAX);j--)
    { 
      flag = updateSeed(&neighborhoods, fastaSeqs2, fastaSeqs1, &neighborhoods.structMap[keys1[j]], &freqs2, &freqs1, tL2, tL1, seedsV, 1); if(flag == 1){seedsV.push_back(keys1[j]);} cout << j+1 << endl; 
    }
  } 
  for(j=0;j<MAX;j++) // positive group ranked by exact diffProp
  {
    if(find(seedsV.begin(), seedsV.end(), keys2[j]) == seedsV.end())
    { 
      flag = updateSeed(&neighborhoods, fastaSeqs1, fastaSeqs2, &neighborhoods.structMap[keys2[j]], &freqs1, &freqs2, tL1, tL2, seedsV, 0); if(flag == 1){seedsV.push_back(keys2[j]);} cout << j+1 << endl; 
    }
  }
  if(d == 2) // negative group ranked by exact diffProp
  {
    for(j = index; j > (index-MAX); j--)
    {
      if(find(seedsV.begin(), seedsV.end(), keys2[j]) == seedsV.end())
      { 
        flag = updateSeed(&neighborhoods, fastaSeqs2, fastaSeqs1, &neighborhoods.structMap[keys2[j]], &freqs2, &freqs1, tL2, tL1, seedsV, 1); if(flag == 1){seedsV.push_back(keys2[j]);} cout << j+1 << endl; 
      }
    }
  } 
  string *seeds = new string[seedsV.size()]; for(int j=0; j<seedsV.size(); j++){seeds[j]=seedsV.at(j);}

  // write out results
  if(seedsV.size() > d*MAX*2 || seedsV.size() < 1){cout << "Error, No seeds converged. Check cg content " <<  pwmFolder << " "<< seedsV.size() << endl;}
  else
  {
    writeOutMotifs(seeds, &neighborhoods,  pwmFolder, seedsV.size(), fastaSeqs1, fastaSeqs2, &freqs1, &freqs2, 0.05, m);
  }

  delete [] exactDiffP;
  delete [] nhoodDiffP;
  delete [] keys1;
  delete [] keys2;
  delete [] seeds;
}

/////////////////////
//    functions    //
/////////////////////
mapTypeString getFastaSeqs(const char* fastaFileName)
{
  string line;
  ifstream fastaFile;
  mapTypeString seqs;
  
  fastaFile.open(fastaFileName);
  if (fastaFile.is_open())
    {
      string tempName;
      string tempSeq;
      int    cntr = 0;
    
      while (! fastaFile.eof() )
	{
	  getline(fastaFile,line);
          if(line[line.size() - 1] == '\r'){line = line.substr(0,line.size()-1);}

          if(line[0] == '>')
          {
 	    size_t tabSpot  = line.find("\t");
	    if(tabSpot != -1){line.replace(tabSpot,1,",");}
            line = line.substr(1,line.size());

	    if(cntr > 0){seqs[tempName] = tempSeq;}  
            tempName = line; 
	    tempSeq  = "";
            cntr++;
	  }
	  else
	  {
            for(int j = 0; j<line.size();j++)
	    {
	      if(     line[j] == 'a'){line.replace(j,1,"A");}
	      else if(line[j] == 'c'){line.replace(j,1,"C");}
	      else if(line[j] == 'g'){line.replace(j,1,"G");}
	      else if(line[j] == 't'){line.replace(j,1,"T");}
	      else if(line[j] == 'n'){line.replace(j,1,"N");}
	    }

            tempSeq = tempSeq+line; 
          }
	}
      seqs[tempName] = tempSeq; 
      fastaFile.close();
      return(seqs);
    }
  else{ cout << "error opening " << fastaFileName << endl;}    

}

mapTypeInt countWmers(mapTypeString seqs, int w, int m, int t, int revFlag)
{
  int cntr = 0;
  mapTypeInt wmerCounts;

  for(mapTypeString::const_iterator it = seqs.begin(); it != seqs.end() && cntr < t; ++it)
    {
      string tempSeq = it->second; 

      for(int j = 0; j<(tempSeq.size()-w+1);j++)
	{
          string wmer = tempSeq.substr(j,w);
	  if(wmer.find("N",0) > w)
          {
            if(wmer == ""){cout << "blank read in " << it->first << " at " << j << endl;}
            else
	    {
  	      ++wmerCounts[wmer];
	      ++wmerCounts["total"];

              if(revFlag==1)
	      {
	        string revWmer = reverseStrand(wmer);
                ++wmerCounts[revWmer];
                ++wmerCounts["total"];
	      }
	    }
	  }
	}
      cntr++;
    }
  return(wmerCounts);
}


mapTypeVectorInt findWmerSites(mapTypeString seqs, int w)
{
  mapTypeVectorInt sites;
  signed int seqInd = 0;
 
  for(mapTypeString::const_iterator it = seqs.begin(); it != seqs.end(); ++it)
  {
    string tempSeq = it->second; 

    for(signed int j = 0; j<(tempSeq.size()-w+1);j++)
    {
      string wmer = tempSeq.substr(j,w);
      if(wmer.find("N",0) > w)
      {
        sites[wmer].push_back(j+seqInd);

        string rWmer = reverseStrand(wmer);
        if(rWmer != wmer){ sites[rWmer].push_back(j+seqInd); }
      }
    }
    seqInd += 100000;
  }
  return(sites);
}

string reverseStrand(string seq)
{
  mapTypeString mirror;
  mirror.insert(pair<string,string>("A", "T"));
  mirror.insert(pair<string,string>("C", "G"));
  mirror.insert(pair<string,string>("G", "C"));
  mirror.insert(pair<string,string>("T", "A"));

  string revSeq; 
  for(int i = seq.size()-1; i >= 0; i--){ revSeq.append(mirror[seq.substr(i,1)]); } 
  return(revSeq);
}

int getNeighborhoodCounts(mapTypeInt counts, string wmer, int M, int start, int matchFlag)
{
  string bases    = "ACGT";
  int wmerCounts  = 0;

  if(start == 0){ wmerCounts += counts[wmer]; }

  for(int i = start; i < wmer.size() - M + 1; i++)
    {
      string tempWmer = wmer;

      for(int j = 0; j < bases.size(); j++)
	{
	  tempWmer.replace(i,1,bases.substr(j,1));
          int matchFlag2 = 0;
          if(matchFlag ==1 && wmer.substr(i,1) == bases.substr(j,1)){matchFlag2 = 1;}

	  if(M>1){     wmerCounts += getNeighborhoodCounts(counts, tempWmer, M-1,i+1,matchFlag2); }
	  else if(M==1)
	    {
	      if(matchFlag2 ==0 ){wmerCounts += counts[tempWmer];}
	    }
	  else {cout << "problem with m in getNeighborhoodCounts, m = " << M << endl;}
	}
    }
  //cout << wmerCounts << endl;
  return(wmerCounts);
}

// recurssion
void cycleThruWchooseM(int *pos,string wmer,int start, int M, int origM, struct4mapTypeInt *cnts1, struct4mapTypeInt *cnts2, wmerNeighborhood *ptrholder, struct4mapTypeDouble *test,mapTypeVectorInt *posSites, mapTypeVectorInt *negSites)
{
  for(int i = start; i < wmer.size()-M+1; i++)
  { 
    pos[origM-M] = i;
    if(M>1){ cycleThruWchooseM(pos, wmer, i+1, M-1, origM, cnts1, cnts2, ptrholder, test, posSites, negSites); }
    if(M==1)
    {
      ptrholder->tempCnts1.clear(); ptrholder->tempCnts1["tot"] = 0;
      ptrholder->tempCnts2.clear(); ptrholder->tempCnts2["tot"] = 0; 
      ptrholder->tempPosSites.clear();  ptrholder->tempNegSites.clear();
      cycleThruNucleotides(wmer, pos, origM, origM, 0, cnts1, cnts2, ptrholder, test, posSites, negSites);

      double adjustedDifProp = computeAdjustedDifProp(ptrholder,cnts1->structMap["total"], cnts2->structMap["total"]);  
      double difProp4Wmer    = singleDifOfProps(ptrholder->tempCnts1["tot"], ptrholder->tempCnts2["tot"], cnts1->structMap["total"], cnts2->structMap["total"]);
      double tempEnrich      = log2((double(ptrholder->tempCnts1["tot"])/ cnts1->structMap["total"]) / ((double)ptrholder->tempCnts2["tot"]/ cnts2->structMap["total"]));

      if(abs(adjustedDifProp) >= abs(ptrholder->difPropBestMMadjusted))
      { 
        ptrholder->difPropBestMisMatches = difProp4Wmer; 
        ptrholder->difPropBestMMadjusted = adjustedDifProp; 
        ptrholder->neighborCnts1         = ptrholder->tempCnts1; 
        ptrholder->neighborCnts2         = ptrholder->tempCnts2; 
        ptrholder->nhoodPosSites         = ptrholder->tempPosSites;
        ptrholder->nhoodNegSites         = ptrholder->tempNegSites;
        ptrholder->adjCnts1              = ptrholder->tempAdjCnts1;
        ptrholder->adjCnts2              = ptrholder->tempAdjCnts2 ;
        ptrholder->bestEnrichment        = ptrholder->tempEnrichment;
        ptrholder->maxFlag               = 1;

        for(int j = 0; j < origM; j++){ptrholder->flexiblePositions[j] = pos[j];}

        makePWM(ptrholder,test); 
        updateLogPWM(ptrholder);
      }
    }
  }
}

void cycleThruNucleotides(string wmer, int *pos, int M, int origM, int start, struct4mapTypeInt *cnts1, struct4mapTypeInt *cnts2, wmerNeighborhood *ptrholder, struct4mapTypeDouble *test,mapTypeVectorInt *posSites, mapTypeVectorInt *negSites)
{
  string bases = "ACGT";
  for(int i = start; i < origM; i++)
  {
    string tempWmer = wmer;
    for(int j = 0; j < bases.size(); j++)
    {
      tempWmer.replace(pos[i],1,bases.substr(j,1));
      cycleThruNucleotides(tempWmer, pos, M-1, origM, i+1, cnts1, cnts2, ptrholder, test, posSites, negSites);
      if(M==1)
      {
        ptrholder->counter++;

	if((ptrholder->difPropBestMisMatches > 0 && test->structMap[tempWmer] > 0) || (ptrholder->difPropBestMisMatches < 0 && test->structMap[tempWmer] < 0))
        {
          ptrholder->tempCnts1[tempWmer] = cnts1->structMap[tempWmer]; ptrholder->tempCnts1["tot"] += cnts1->structMap[tempWmer];
          ptrholder->tempCnts2[tempWmer] = cnts2->structMap[tempWmer]; ptrholder->tempCnts2["tot"] += cnts2->structMap[tempWmer];

          ptrholder->tempPosSites.insert(ptrholder->tempPosSites.end(),(*posSites)[tempWmer].begin(),(*posSites)[tempWmer].end()); 
          ptrholder->tempNegSites.insert(ptrholder->tempNegSites.end(),(*negSites)[tempWmer].begin(),(*negSites)[tempWmer].end()); 
        }
      }
    }
  }
}

double computeAdjustedDifProp(wmerNeighborhood *ptrholder, int n1, int n2)  
{
  sort(ptrholder->tempPosSites.begin(), ptrholder->tempPosSites.end()); //unique(ptrholder->tempPosSites.begin(), ptrholder->tempPosSites.end());
  sort(ptrholder->tempNegSites.begin(), ptrholder->tempNegSites.end()); //unique(ptrholder->tempNegSites.begin(), ptrholder->tempNegSites.end());

  double posTotBP = 0; 
  for(vector<signed int>::iterator it = ptrholder->tempPosSites.begin(); it != ptrholder->tempPosSites.end(); ++it)
  {
    if(it+1 == ptrholder->tempPosSites.end()){ posTotBP += ptrholder->width; }else{ posTotBP += min(*(it+1) - *it, ptrholder->width); }
  }

  double negTotBP = 0; 
  for(vector<signed int>::iterator it = ptrholder->tempNegSites.begin(); it != ptrholder->tempNegSites.end(); ++it)
  {
    if(it+1 == ptrholder->tempNegSites.end()){ negTotBP += ptrholder->width; }else{ negTotBP += min(*(it+1) - *it, ptrholder->width); }
  }

  double cnts1 = posTotBP/ptrholder->width;
  double cnts2 = negTotBP/ptrholder->width;

  double p1 = cnts1 / n1;
  double p2 = cnts2 / n2;
  double p  = (cnts1 + cnts2) / (n1 + n2);
  double dp = (p1 - p2) / sqrt(p*(1-p)*(1/((double)n1) + 1/((double)n2)));

  ptrholder->tempAdjCnts1 = cnts1;
  ptrholder->tempAdjCnts2 = cnts2;
  ptrholder->tempEnrichment = log2(p1/p2);
  return(dp);
}

double singleDifOfProps(long double cnts1, long double cnts2, int n1, int n2)
{
  double p1 = (double) cnts1 / n1;
  double p2 = (double) cnts2 / n2;
  double p  = (double) (cnts1 + cnts2) / (n1 + n2);
  double dp = (p1 -p2) / sqrt(p*(1-p)*(1/((double)n1) + 1/((double)n2)));

  return(dp);
}

mapTypeDouble diffOfProps(mapTypeInt gCnts1, mapTypeInt gCnts2, mapTypeInt gCntsTot, int n1, int n2)
{
  mapTypeDouble dOfProps;
  for(mapTypeInt::const_iterator it = gCntsTot.begin(); it != gCntsTot.end(); ++it)
    {
      if(it->first != "total")
	{
	  string key = it->first;
          dOfProps[key] = singleDifOfProps(gCnts1[key], gCnts2[key], n1, n2);
	}
    }
  return(dOfProps);
}

void bubble_sort(double *array, string *names, int len)
{
  int i, j, flag = 1;    // set flag to 1 to begin initial pass
  double  temp;          // holding variable
  string  temp2;         // holding variable

  for(i = 1; (i <= len) && flag; i++)
  {
    flag = 0;
    for (j=0; j < (len -1); j++)
    {
      if(array[j+1] > array[j])     
      { 
	temp  = array[j]; array[j] = array[j+1]; array[j+1] = temp;   // swap values
	temp2 = names[j]; names[j] = names[j+1]; names[j+1] = temp2;  // swap keys
	flag  = 1;        // indicates that a swap occurred.
      }
    }
  }
  return;
}


void add2pwm(double mat[][4], string wmer, double count)
{
  mapTypeCharInt nucleotide;
  nucleotide['A'] = 0;
  nucleotide['C'] = 1;
  nucleotide['G'] = 2;
  nucleotide['T'] = 3;

  for(int i=0; i<wmer.size();i++){mat[i][nucleotide[wmer.at(i)]] += count;}
  return;
}

void makePWM(wmerNeighborhood *ptrholder, struct4mapTypeDouble *test)
{
  for(int j = 0; j < ptrholder->width; j++){ for(int k = 0; k<4; k++){ptrholder->pwm[j][k] = 0;} }

  for(mapTypeInt::const_iterator it = ptrholder->neighborCnts1.begin(); it != ptrholder->neighborCnts1.end(); it++)
  { 
    if(it->first != "tot")
    { 
      if(ptrholder->difPropExactMatch > 0){ add2pwm(ptrholder->pwm, it->first, it->second); }
      else{ add2pwm(ptrholder->pwm, it->first, -1*it->second); }
    } 
  }

  for(mapTypeInt::const_iterator it = ptrholder->neighborCnts2.begin(); it != ptrholder->neighborCnts2.end(); it++)
  { 
    if(it->first != "tot")
    { 
      if(ptrholder->difPropExactMatch > 0){ add2pwm(ptrholder->pwm, it->first, -1*(double)it->second * ptrholder->totLen1 / (double)ptrholder->totLen2 ); }
      else{ add2pwm(ptrholder->pwm, it->first, (double)it->second * ptrholder->totLen1 / (double)ptrholder->totLen2 ); }
    } 
  }

  double total = 0;
  for(int j = 0; j < ptrholder->width; j++){ for(int k = 0; k<4; k++){ptrholder->pwm[j][k] = max(ptrholder->pwm[j][k], 0.0);}}
  for(int k = 0; k<4; k++){total +=  ptrholder->pwm[0][k];}
  for(int j = 0; j < ptrholder->width; j++){total=0; for(int k=0; k<4; k++){total += ptrholder->pwm[j][k];} for(int k=0; k<4; k++){ptrholder->pwm[j][k] += 0.05 * total;} }
}

void updateLogPWM2(wmerNeighborhood *ptrholder)
{ 
  double pwmTotal;
  double prevPWM[ptrholder->width][4]; for(int j = 0; j < ptrholder->width; j++){for(int i=0; i<4; i++){prevPWM[j][i] = exp(ptrholder->prevlogPWM[j][i]);}}
  double currentPWM[ptrholder->width][4]; for(int j = 0; j < ptrholder->width; j++){for(int i=0; i<4; i++){currentPWM[j][i] = exp(ptrholder->logPWM[j][i]);}}
  double prevdist = 0;
  double dist     = 0;

  for(int j = 0; j < ptrholder->width; j++)
  { 
    pwmTotal  = 0; for(int i=0; i<4; i++){pwmTotal+=ptrholder->pwm[j][i];}
    for(int k = 0; k<4; k++)
    { 
      ptrholder->logPWM[j][k]     = log(ptrholder->pwm[j][k]/(double)pwmTotal); 
      dist                        = max(dist,    abs(currentPWM[j][k] - exp(ptrholder->logPWM[j][k]))); 
      prevdist                    = max(prevdist,abs(prevPWM[j][k] - exp(ptrholder->logPWM[j][k]))); 
      ptrholder->prevlogPWM[j][k] = log(currentPWM[j][k]);
    }
  }

  if(dist <= 0.01 || prevdist <= 0.01){ptrholder->convergedFlag = 1;}
}

void updateLogPWM(wmerNeighborhood *ptrholder)
{
  double pwmTotal;

  for(int j = 0; j < ptrholder->width; j++)
  { 
    pwmTotal  = 0; for(int i=0; i<4; i++){pwmTotal+=ptrholder->pwm[j][i];}
    for(int k = 0; k<4; k++){ ptrholder->logPWM[j][k] = log(ptrholder->pwm[j][k]/(double)pwmTotal); }
  }
}

double computeLogLikelihoodRatio(string wmer, double pwm[][4], struct4mapTypeDouble *freqs, int revFlag)
{
  mapTypeCharInt nucleotide;
  nucleotide['A'] = 0;
  nucleotide['C'] = 1;
  nucleotide['G'] = 2;
  nucleotide['T'] = 3;
 
  double numerator   = pwm[0][nucleotide[wmer.at(0)]]; 
  for(int i=1; i<wmer.size();i++){numerator = numerator + pwm[i][nucleotide[wmer.at(i)]];} 

  if(revFlag == 1){wmer =  reverseStrand(wmer);}
  double denominator = freqs->structMap[wmer.substr(0,1)];
  for(int i=0; i<wmer.size()-1;i++){denominator = denominator + freqs->structMap[wmer.substr(i,2)];}

  double llr = numerator-denominator;
  return(llr);
}

void updateMotif(wmerNeighborhood *seed, mapTypeString seqs1, mapTypeString seqs2, int w,  struct4mapTypeDouble *freqs1, struct4mapTypeDouble *freqs2, int totLen1, int totLen2, mapTypeVectorInt inds1, mapTypeVectorInt inds2, double extra,int updateSitesFlag)
{
  mapTypeVectorInt posInds;
  mapTypeVectorInt negInds;

  double   llr, rllr, weighted, rweighted, LR, rLR;
  string   tempSeq, wmer, rwmer;
  long double thresh = seed->threshold; 
  int      l,j;

  mapTypeCharInt nucleotide;
  nucleotide['A'] = 0; nucleotide['C'] = 1; nucleotide['G'] = 2; nucleotide['T'] = 3;
 
  double newPWM[w][4]; for(int i = 0; i < w; i++){for(j = 0; j < 4; j++){newPWM[i][j] =0 ;}}
  double N0[4][4];     for(int i = 0; i < 4; i++){for(j = 0; j < 4; j++){N0[i][j]     =0 ;}}
  double N1[4][4];     for(int i = 0; i < 4; i++){for(j = 0; j < 4; j++){N1[i][j]     =0 ;}}
  double Nw[4][4];     for(int i = 0; i < 4; i++){for(j = 0; j < 4; j++){Nw[i][j]     =0 ;}}
  double Nwp1[4][4];   for(int i = 0; i < 4; i++){for(j = 0; j < 4; j++){Nwp1[i][j]   =0 ;}}

  for(mapTypeString::const_iterator it = seqs1.begin(); it != seqs1.end(); ++it)
  {
    vector<signed int> posV; 

    for(l = 0; l < inds1[it->first].size(); l++)
    {
      j     = inds1[it->first].at(l) + seed->indexOffset;
 
      if(j>=0 && j<(it->second.size() - w +1))
      {  
        wmer  = it->second.substr(j,w);
        if(wmer.find("N",0) >= w)
        {
          rwmer = reverseStrand(wmer);
          LR    = computeLogLikelihoodRatio(wmer, seed->logPWM, freqs1, 0); 
          rLR   = computeLogLikelihoodRatio(rwmer,seed->logPWM, freqs1, 1);  
  
          if(LR > thresh && LR >= rLR)
          {       
            for(int i=0; i<w;i++){newPWM[i][nucleotide[wmer.at(i)]]++;} 
            posV.push_back(j);

            if(j>1 && j<(it->second.size() - w -1))
            {
              wmer  = it->second.substr(j-2,w+4); //extend 2 bp in front and behind
              if(wmer.find("N",0) >= w+4)
              {
                N0[nucleotide[wmer.at(0)]][nucleotide[wmer.at(1)]]++;
                N1[nucleotide[wmer.at(1)]][nucleotide[wmer.at(2)]]++;
                Nw[nucleotide[wmer.at(w)]][nucleotide[wmer.at(w+1)]]++;
                Nwp1[nucleotide[wmer.at(w+1)]][nucleotide[wmer.at(w+2)]]++;
              }
            }
          }	
          else if(rLR > thresh && rLR > LR)
          { 

            for(int i=0; i<w;i++){newPWM[i][nucleotide[rwmer.at(i)]]++;} 
            posV.push_back(j);

            if(j>1 && j<(it->second.size() - w -1))
            {
              wmer  = it->second.substr(j-2,w+4); //extend 2 bp in front and behind
              rwmer = reverseStrand(wmer);

              if(wmer.find("N",0) >= w+4)
              {

                N0[nucleotide[rwmer.at(0)]][nucleotide[rwmer.at(1)]]++;
                N1[nucleotide[rwmer.at(1)]][nucleotide[rwmer.at(2)]]++;
                Nw[nucleotide[rwmer.at(w)]][nucleotide[rwmer.at(w+1)]]++;
                Nwp1[nucleotide[wmer.at(w+1)]][nucleotide[rwmer.at(w+2)]]++;
	      }
	    }
	  }
        posInds[it->first] = posV;
	}
      }
    }
  }

  for(mapTypeString::const_iterator it = seqs2.begin(); it != seqs2.end(); ++it)
  {
    vector<signed int> negV;
 
    for(l = 0; l < inds2[it->first].size(); l++)
    {
      j     = inds2[it->first].at(l) + seed->indexOffset;

      if(j>=0 && j<(it->second.size() - w +1))
      {  
        wmer  = it->second.substr(j,w); 
        if(wmer.find("N",0) >= w)
        {
          rwmer = reverseStrand(wmer);
          LR    = computeLogLikelihoodRatio(wmer.substr(2,w), seed->logPWM, freqs2, 0); 
          rLR   = computeLogLikelihoodRatio(rwmer.substr(2,w),seed->logPWM, freqs2, 1);  
  
          if(LR > thresh && LR >= rLR)
          {
            for(int i=0; i<w;i++){ newPWM[i][nucleotide[wmer.at(i)]]  -= (double)totLen1/totLen2;} 
            negV.push_back(j);

            if(j>1 && j<(it->second.size() - w -1))
            {
              wmer  = it->second.substr(j-2,w+4);  //extend 2bp extra in each direction 
              if(wmer.find("N",0) >= w+4)
              {
                N0[nucleotide[wmer.at(0)]][nucleotide[wmer.at(1)]]              -= (double)totLen1/totLen2;
                N1[nucleotide[wmer.at(1)]][nucleotide[wmer.at(2)]]              -= (double)totLen1/totLen2;
                Nw[nucleotide[wmer.at(w)]][nucleotide[wmer.at(w+1)]]            -= (double)totLen1/totLen2;
                Nwp1[nucleotide[wmer.at(w+1)]][nucleotide[wmer.at(w+2)]]        -= (double)totLen1/totLen2;
	      }
	    }
          }
          else if(rLR > thresh && rLR > LR)
          {   
            for(int i=0; i<w;i++){ newPWM[i][nucleotide[rwmer.at(i)]]   -= (double)totLen1/totLen2;} 
            negV.push_back(j);

            if(j>1 && j<(it->second.size() - w -1))
            {
              wmer  = it->second.substr(j-2,w+4);  //extend 2bp extra in each direction 
              rwmer = reverseStrand(wmer);
              if(wmer.find("N",0) >= w+4)
              {
                N0[nucleotide[rwmer.at(0)]][nucleotide[rwmer.at(1)]]              -= (double)totLen1/totLen2;
                N1[nucleotide[rwmer.at(1)]][nucleotide[rwmer.at(2)]]              -= (double)totLen1/totLen2;
                Nw[nucleotide[rwmer.at(w)]][nucleotide[rwmer.at(w+1)]]            -= (double)totLen1/totLen2;
                Nwp1[nucleotide[rwmer.at(w+1)]][nucleotide[rwmer.at(w+2)]]        -= (double)totLen1/totLen2;
	      }
	    }  
  	  }
          negInds[it->first] = negV;
	}
      }
    }
  }

  for(int i = 0; i < w; i++){for(j = 0; j < 4; j++){newPWM[i][j] = max(newPWM[i][j], (double)0); }}
  for(int i = 0; i < w; i++){int j; double total = 0; for(j = 0; j < 4; j++){total+=newPWM[i][j];} for(j = 0; j < 4; j++){seed->pwm[i][j] = newPWM[i][j] + extra * total; }}

  if(seed->width < 20 && updateSitesFlag == 0){ determineMotifLength(seed, freqs1 , N0, N1, Nw, Nwp1,lowerBound, upperBound); }
  if(updateSitesFlag == 1){seed->posSites=posInds; seed->negSites=negInds;}

  updateLogPWM2(seed);
}

void writeOutMotifs(string *seeds, struct4mapTypeWmerNeighborhood *nhoods, string folder, int maxNumSeeds, mapTypeString seqs1, mapTypeString seqs2, struct4mapTypeDouble *freqs1, struct4mapTypeDouble *freqs2, double extra, int origM)
{
  FILE * pwmFile;
  string file = folder + "/output.txt";
  pwmFile     = fopen(file.c_str(),"w");
  string consensus, wmer, rwmer;;
  mapTypeString seqsI;
  mapTypeString seqsII;
  vector<signed int> v;
  double LR1, rLR1;
  wmerNeighborhood *ptrholder;

  double *bestDiffP = new double[maxNumSeeds];
  for(int i = 0; i<maxNumSeeds; i++){bestDiffP[i] = nhoods->structMap[seeds[i]].bestEnrichment;}
  bubble_sort(bestDiffP, seeds, maxNumSeeds-1);

  for(int i = 0; i<maxNumSeeds; i++)  
  {
    ptrholder = &nhoods->structMap[seeds[i]];
    consensus =  makeConsensus(seeds[i], nhoods);
    fprintf(pwmFile,"MOTIF:\t%s\nInitial Seed:\t%s\n", consensus.c_str(), seeds[i].c_str());
    fprintf(pwmFile, "Flexible Seed Positions: %i", nhoods->structMap[seeds[i]].flexiblePositions[0] + 1);
    for(int j = 1; j < origM; j++){fprintf(pwmFile, ", %i", nhoods->structMap[seeds[i]].flexiblePositions[j] + 1);}
    fprintf(pwmFile,"\nLikelihood Threshold:\t%5.5f\nt-score:\t%5.5f\nEnrichment:\t%5.5f(log2)\n\nPWM:\n", exp((double)nhoods->structMap[seeds[i]].threshold), nhoods->structMap[seeds[i]].difPropUpdated, nhoods->structMap[seeds[i]].bestEnrichment);

    fprintf(pwmFile,"A\tC\tG\tT\n");
    for(int i = 0; i < ptrholder->width; i++)
    {
      int j; double total = 0; 
      for(j = 0; j < 4; j++){total+=ptrholder->pwm[i][j];} 
      for(j = 0; j < 4; j++){fprintf(pwmFile,"%5.2f\t", ptrholder->pwm[i][j] / total); ptrholder->pwm[i][j] = ptrholder->pwm[i][j] + extra * total; }
      fprintf(pwmFile,"\n");
    }

    updateLogPWM2(ptrholder);

    if(ptrholder->difPropUpdated < 0){seqsI = seqs2; seqsII = seqs1;}else{seqsI = seqs1; seqsII = seqs2;}

    fprintf(pwmFile,"\n\nPositive Sites:\n");
    for(mapTypeString::const_iterator it = seqsI.begin(); it != seqsI.end(); ++it) 
    { 
      v = ptrholder->posSites[it->first];  

      for(int l = 0; l < v.size(); l++) 
      { 
        wmer  = it->second.substr(v.at(l), ptrholder->width); rwmer = reverseStrand(wmer); 
        LR1   = computeLogLikelihoodRatio(wmer, ptrholder->logPWM, freqs1, 0);     rLR1 = computeLogLikelihoodRatio(rwmer,ptrholder->logPWM, freqs1, 1);   
        fprintf(pwmFile,">%s\t%i\t%7.2f\n",it->first.c_str(),v.at(l),exp(max(LR1,rLR1))); if(rLR1 > LR1){fprintf(pwmFile,"%s\n",rwmer.c_str());}else{fprintf(pwmFile,"%s\n",wmer.c_str());}
      }
    }

    fprintf(pwmFile,"\n\nNegative Sites:\n");
    for(mapTypeString::const_iterator it = seqsII.begin(); it != seqsII.end(); ++it) 
    { 
      v = ptrholder->negSites[it->first];  

      for(int l = 0; l < v.size(); l++) 
      { 
        wmer  = it->second.substr(v.at(l), ptrholder->width); rwmer = reverseStrand(wmer); 
        LR1   = computeLogLikelihoodRatio(wmer, ptrholder->logPWM, freqs1, 0);     rLR1 = computeLogLikelihoodRatio(rwmer,ptrholder->logPWM, freqs1, 1);   
        fprintf(pwmFile,">%s\t%i\t%7.2f\n",it->first.c_str(),v.at(l),exp(max(LR1,rLR1))); if(rLR1 > LR1){fprintf(pwmFile,"%s\n",rwmer.c_str());}else{fprintf(pwmFile,"%s\n",wmer.c_str());}
      }
    }
    fprintf(pwmFile,"\n####################\n\n");
  }
  fclose(pwmFile);
}

string makeConsensus(string seed, struct4mapTypeWmerNeighborhood *nhoods)
{
  string nucleotides = "ACGT";
  double best;
  string bestChar; 
  string returnStr;

  for(int i = 0; i < nhoods->structMap[seed].width; i++)
  {
    best = 0;
    for(int j = 0; j < 4; j++)
    {
      if(nhoods->structMap[seed].pwm[i][j] > best){bestChar = nucleotides.substr(j,1); best = nhoods->structMap[seed].pwm[i][j];}
    }
    returnStr.append(bestChar);
  }
  return(returnStr);
}

mapTypeVectorInt scanNfindInds(mapTypeString seqs, wmerNeighborhood* seed, struct4mapTypeDouble* freqs, double cutoff)
{
  mapTypeVectorInt inds;
  int w = seed->width;
  string tempSeq, wmer, rwmer;
  double LR, rLR;
  int j;

  for(mapTypeString::const_iterator it = seqs.begin(); it != seqs.end(); ++it)
  {
    vector<signed int> v; vector<double> val; vector<signed int> flag;
    tempSeq   = it->second; 
    int cntr  = 0; int cntr2 = 0;

    for(j = 0; j<(tempSeq.size()-w+1);j++)
    {
      wmer    = tempSeq.substr(j,w); 
      rwmer   = reverseStrand(wmer);
      if(wmer.find("N",0) >= w)
      {
        LR    = exp(computeLogLikelihoodRatio(wmer, seed->logPWM, freqs, 0)); rLR   = exp(computeLogLikelihoodRatio(rwmer,seed->logPWM, freqs, 1));  
	if((LR > cutoff) || (rLR > cutoff))
        {
          v.push_back(cntr); val.push_back(max(LR,rLR)); flag.push_back(0); 
	  if(cntr2 > 0){ if(v.at(cntr2) - v.at(cntr2-1) < 2*w){ flag.at(cntr2)   = 1; flag.at(cntr2-1) = 1;} } 
          cntr2++; 
        }
      }
      cntr++;
    }

    for(j = 0; j < v.size(); j++)
    {
      if(flag.at(j) == 1)
      {
        int    maxloc = j;
        double maxval = val.at(j);
 
        while(j+1 < v.size() && flag.at(j+1) == 1)
        {
          if(val.at(j+1) > maxval){ maxval = val.at(j+1); maxloc = j+1; }
          j++;
        }
        flag.at(maxloc) = 0;
      }
    }

    vector<signed int> v2; for(j=0; j < v.size(); j++){if(flag.at(j) == 0){v2.push_back(v.at(j));}}
    inds[it->first] = v2;
  }
  return(inds);
}

mapTypeVectorInt scanNfindIndsFancy(mapTypeString seqs, wmerNeighborhood* seed, struct4mapTypeDouble* freqs,mapTypeVectorInt lowerCutoffInds, double cutoff)
{
  mapTypeVectorInt inds;
  int w = seed->width;
  string tempSeq, wmer, rwmer;
  double LR, rLR;
  int j, l;

  for(mapTypeString::const_iterator it = seqs.begin(); it != seqs.end(); ++it)
  {
    vector<signed int> v; 
    tempSeq   = it->second; 
    int cntr  = 0; int cntr2 = 0;

    for(l = 0; l < lowerCutoffInds[it->first].size(); l++)
    {
      j     = lowerCutoffInds[it->first].at(l) + seed->indexOffset;

      if(j>=0 && j<(it->second.size() - w +1))
      {
        wmer  = it->second.substr(j,w);
      
        if(wmer.find("N",0) >= w)
        {
          rwmer = reverseStrand(wmer);
          LR    = exp(computeLogLikelihoodRatio(wmer, seed->logPWM, freqs, 0)); 
          rLR   = exp(computeLogLikelihoodRatio(rwmer,seed->logPWM, freqs, 1));  

          if(LR > cutoff || rLR > cutoff)
	  { 
            v.push_back(j - seed->indexOffset);  
          }
        }
      }
    }
    inds[it->first] = v;
  }
  return(inds);
}

int assessOverlaps(struct4mapTypeWmerNeighborhood *nhoods, wmerNeighborhood *ptrholder1, vector<string>& seedsV)
{
  int total1  = 0;
  for(mapTypeVectorInt::const_iterator it = ptrholder1->posSites.begin(); it != ptrholder1->posSites.end(); ++it){ total1 += ptrholder1->posSites[it->first].size(); }
  for(mapTypeVectorInt::const_iterator it = ptrholder1->negSites.begin(); it != ptrholder1->negSites.end(); ++it){ total1 += ptrholder1->negSites[it->first].size(); }
  
  int total2, j, k, l, overLap, dist;
  wmerNeighborhood holder2;

  for(j=0; j<seedsV.size(); j++)
  {
    holder2 = nhoods->structMap[seedsV.at(j)]; 
    total2  = 0;
    overLap = 0;

    for(mapTypeVectorInt::const_iterator it = holder2.posSites.begin(); it != holder2.posSites.end(); ++it)
    { 
      total2 += holder2.posSites[it->first].size(); 
      for(k = 0; k < ptrholder1->posSites[it->first].size(); k++)
      { 
        for(l = 0; l < holder2.posSites[it->first].size(); l++)
        { 
          dist = abs(ptrholder1->posSites[it->first].at(k) - holder2.posSites[it->first].at(l)); 
          if(dist <= ptrholder1->width || dist <= holder2.width){overLap++;} }
        }
    }

    for(mapTypeVectorInt::const_iterator it = holder2.negSites.begin(); it != holder2.negSites.end(); ++it)
    { 
      total2 += holder2.negSites[it->first].size(); 
      for(k = 0; k < ptrholder1->negSites[it->first].size(); k++)
      { 
        for(l = 0; l < holder2.negSites[it->first].size(); l++)
        { 
          dist = abs(ptrholder1->negSites[it->first].at(k) - holder2.negSites[it->first].at(l)); 
          if(dist <= ptrholder1->width || dist <= holder2.width){overLap++;} }
        }
    }
    if(min(overLap/(double)total1,overLap/(double)total2) > .66){cout << "too much overlap with previous motif" << endl; return(1);}
  }
  return(0);
}

int updateSeed(struct4mapTypeWmerNeighborhood *nhoods, mapTypeString seqs1, mapTypeString seqs2, wmerNeighborhood *ptrholder,struct4mapTypeDouble *freqs1,struct4mapTypeDouble *freqs2, int totLen1, int totLen2, vector<string>& seedsV, int directionFlag)
{
  int flag               = 1; 
  int overlapFlag        = 0;
  ptrholder->indexOffset = 0;
  mapTypeVectorInt basicPossibleSites1 = scanNfindInds(seqs1, ptrholder, freqs1,2);        
  mapTypeVectorInt basicPossibleSites2 = scanNfindInds(seqs2, ptrholder, freqs2,2);        

  mapTypeVectorInt possibleSites1;
  mapTypeVectorInt possibleSites2;

  for(int l=0; l < 30; l++)
  {
    if(((double)l/5 - floor((double)l/5)) == 0 )
    { 
      possibleSites1      = scanNfindIndsFancy(seqs1, ptrholder, freqs1, basicPossibleSites1, 5); 
      possibleSites2      = scanNfindIndsFancy(seqs2, ptrholder, freqs2, basicPossibleSites2, 5); 
      ptrholder->posSites = scanNfindIndsFancy(seqs1, ptrholder, freqs1, basicPossibleSites1, 1000); 
      ptrholder->negSites = scanNfindIndsFancy(seqs2, ptrholder, freqs2, basicPossibleSites2, 1000); 
    }
    flag = updateSingleMotif(possibleSites1, possibleSites2, seqs1, seqs2, ptrholder, freqs1, freqs2, totLen1, totLen2);

    if(ptrholder->convergedFlag == 1){break;}
    if(flag == 0){l = 10000; break;}
  }
  overlapFlag = assessOverlaps(nhoods, ptrholder, seedsV);

  updateMotif(ptrholder,seqs1,seqs2,ptrholder->width,freqs1,freqs2,totLen1,totLen2,possibleSites1,possibleSites2,0.0,1); 
  double total=0; for(int i=0; i < 4; i++){total += ptrholder->pwm[0][i];}

  if(directionFlag == 1){ ptrholder->difPropUpdated = -ptrholder->difPropUpdated; ptrholder->bestEnrichment = -ptrholder->bestEnrichment; }

  if(flag == 1 && overlapFlag == 0 && total > 0 ){  return(1); }
  else if(flag == 1 && overlapFlag == 1 && total > 0){return(2);}
  else{return(0);}
}

int updateSingleMotif(mapTypeVectorInt inds1, mapTypeVectorInt inds2, mapTypeString seqs1, mapTypeString seqs2, wmerNeighborhood *ptrholder,struct4mapTypeDouble *freqs1,struct4mapTypeDouble *freqs2, int totLen1, int totLen2)
{ 
  int           l, j, k, o;
  int           w = ptrholder->width;
  double        pwm[w][4];
  long double   llr, rllr, weight, LR, rLR, C;
  string        tempSeq; string   wmer; string   Rwmer;
  int           llrCutN = 20;
  double        *thresholds = new double[llrCutN]; for(j = 0; j < llrCutN; j++){thresholds[j] = log((double)100*(j+1));} 
  int           llrCuts = 0;
  double        maxDifProp = 0;
  signed int    motifCnts1[seqs1.size()][llrCutN];
  signed int    motifCnts2[seqs2.size()][llrCutN];

  ptrholder->FDR = 1;
  for(j = 0;j<seqs1.size();j++){ for(k = 0; k<llrCutN; k++){ motifCnts1[j][k]= 0; } }
  for(j = 0;j<seqs2.size();j++){ for(k = 0; k<llrCutN; k++){ motifCnts2[j][k]= 0; } }
  for(j = 0;j< w;          j++){for( k = 0 ;k < 4;     k++){pwm[j][k] = ptrholder->logPWM[j][k];}}

  // first set of seqs 
  int cnt1    = 0; int *tot1 = new int[llrCutN]; for(int n = 0; n < llrCutN; n++){tot1[n] =0;} 
  for(mapTypeString::const_iterator it = seqs1.begin(); it != seqs1.end(); ++it)
  {
    for(l = 0; l < inds1[it->first].size(); l++)
    {
      j    = inds1[it->first].at(l) + ptrholder->indexOffset; 
      if(j > 0)
      {
        wmer = it->second.substr(j,w); 
        if(wmer.find("N",0) >= w)
        {
          LR = computeLogLikelihoodRatio(wmer, pwm, freqs1, 0); rLR = computeLogLikelihoodRatio(reverseStrand(wmer), pwm, freqs1, 1);
          for(k = 0; k<llrCutN; k++)
          { 
            if(LR > thresholds[k] || rLR > thresholds[k]){motifCnts1[cnt1][k]++;}
          }
        } 
      }
    }
    cnt1++;
  }
  for(l = 0; l < cnt1; l++){for(int m = 0; m < llrCutN; m++){tot1[m] += motifCnts1[l][m];}}
  
  // second set of seqs 
  int cnt2    = 0; int *tot2 = new int[llrCutN]; for(int n = 0; n < llrCutN; n++){tot2[n] =0;}
  for(mapTypeString::const_iterator it = seqs2.begin(); it != seqs2.end(); ++it)
  {
    for(l = 0; l < inds2[it->first].size(); l++)
    {
      j    = inds2[it->first].at(l) + ptrholder->indexOffset; 
      if(j > 0)
      {
        wmer = it->second.substr(j,w); 
        if(wmer.find("N",0) >= w)
        {
          LR   = computeLogLikelihoodRatio(wmer, pwm, freqs2, 0); rLR = computeLogLikelihoodRatio(reverseStrand(wmer), pwm, freqs2, 1);
          for(k = 0; k<llrCutN; k++)
          { 
            if(LR > thresholds[k] || rLR > thresholds[k]){motifCnts2[cnt2][k]++;}
          }
        } 
      }
    }
    cnt2++;
  }
  for(l = 0; l < cnt2; l++){for(int m = 0; m < llrCutN; m++){tot2[m] += motifCnts2[l][m];}}
  
  // now do diffofproportions
  double p1                 = 0; double p2     = 0;
  ptrholder->difPropUpdated = 0;
  ptrholder->bestEnrichment = 0;
  ptrholder->logLratio      = 0;
  int totI = 0;   int totII = 0;  
  
  double difProp4Wmer, binP1, enrichment, logLR, fdr; 
  int dirFlag;
  
  for(o = 0; o < llrCutN; o++)
  {
    difProp4Wmer = singleDifOfProps(tot1[o], tot2[o], totLen1, totLen2); p1 = tot1[o]/(double)totLen1; p2 = tot2[o]/(double)totLen2; 
    dirFlag      = 1;
    fdr          = 1;
  
    if(difProp4Wmer > 0){ fdr = tot2[o] * (double)totLen1/(double)totLen2/tot1[o]; }
    else{dirFlag = 0; if( o == (llrCutN - 1)){ptrholder->convergedFlag = 1;}}
  
    if(fdr < fdrCut && dirFlag == 1) 
    {
      int O        = max(o-1, 0);
      enrichment   = log2(p1 / p2);
  
      ptrholder->FDR            = fdr;
      ptrholder->bestEnrichment = enrichment;
      ptrholder->difPropUpdated = difProp4Wmer; 
      ptrholder->threshold      = thresholds[o]; 
      ptrholder->p1             = p1; 
      ptrholder->p2             = p2; 
      ptrholder->seqP1          = tot1[o]/(double)seqs1.size(); 
      ptrholder->seqP2          = tot2[o]/(double)seqs2.size(); 
      totI = tot1[o];     totII = tot2[o];
      break;
    }
  }

  updateMotif(ptrholder,seqs1,seqs2,w,freqs1,freqs2,totLen1,totLen2,inds1,inds2,0.05,0);

  double rowSum;
  for(int i = 0; i < ptrholder->width; i++){rowSum = 0; int j; for(j = 0; j < 4; j++){rowSum += ptrholder->pwm[i][j];} if(rowSum == 0){return(0);}}
  delete [] thresholds;
  delete [] tot1;
  delete [] tot2;

  if(o == 20){ return(0); }else{ return(1); }
}

void determineMotifLength(wmerNeighborhood *seed, struct4mapTypeDouble *freqs ,double N0[][4],double N1[][4],double Nw[][4],double Nwp1[][4], int lBound, int uBound) 
{ 
  string maxPos   = "asis";
  double maxLBF   = 1;
  double marginals[4];    
  double total    = 0;
 
  double lBF0   = computeBayesFactor(freqs, N0);   if(lBF0   > maxLBF){maxLBF = lBF0;   maxPos = "0";}
  double lBF1   = -computeBayesFactor(freqs, N1);  if(lBF1   > maxLBF){maxLBF = lBF1;   maxPos = "1";}
  double lBFw   = -computeBayesFactor(freqs, Nw);  if(lBFw   > maxLBF){maxLBF = lBFw;   maxPos = "w";}
  double lBFwp1 = computeBayesFactor(freqs, Nwp1); if(lBFwp1 > maxLBF){maxLBF = lBFwp1; maxPos = "wp1";}
  double pwm[20][4]; for(int i=0; i < 20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = 0; }}
  double IC     = 0; // information content

  if(seed->width < uBound)
  {
    if(maxPos == "0")
    {
      for(int i = 0; i<4; i++){ marginals[i] = 0; for(int j=0; j<4; j++){ marginals[i] += max(N0[j][i],0.001); } total += marginals[i]; }
      for(int i = 0; i<4; i++){ IC += (double)marginals[i]/total * log2((double)marginals[i]/total);} IC += 2;
      double X2 = chiSquared(freqs, N0);

      if(X2 > X2cutoff && IC > .1)
      {
        for(int i=0; i<20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = 0;}}
  
        pwm[0][0]   = marginals[0]; pwm[0][1] = marginals[1]; pwm[0][2] = marginals[2]; pwm[0][3] = marginals[3]; 
        for(int i   = 0; i < seed->width; i++){ for(int j=0; j<4; j++){ pwm[i+1][j] = seed->pwm[i][j]; }}
        seed->width = seed->width + 1; seed->indexOffset -= 1;
      }
      else{ for(int i=0; i < seed->width; i++){ for(int j=0; j<4; j++){ pwm[i][j] = seed->pwm[i][j];}} }
    }
    else if(maxPos == "1" && seed->width > lBound)
    {
      for(int i=0; i<20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = 0;}}
      for(int i   = 1; i < seed->width; i++){ for(int j=0; j<4; j++){ pwm[i-1][j] = seed->pwm[i][j]; }}
      seed->width = seed->width - 1; seed->indexOffset += 1;
    }
    else if(maxPos == "w" && seed->width > lBound)
    {
      for(int i=0; i<20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = seed->pwm[i][j];}}
      pwm[seed->width-1][0] = 0; pwm[seed->width-1][1] = 0; pwm[seed->width-1][2] = 0; pwm[seed->width-1][3] = 0; 
      seed->width = seed->width - 1; 
    }
    else if(maxPos == "wp1")
    {
      for(int i = 0; i<4;  i++){ marginals[i] = 0; for(int j=0; j<4; j++){ marginals[i] += max(Nwp1[j][i],0.0001); } total += marginals[i]; } 
      for(int i = 0; i<4; i++){ IC += (double)marginals[i]/total * log2((double)marginals[i]/total);} IC += 2;
      double X2 = chiSquared(freqs, N0);

      if(X2 > X2cutoff && IC > .1)
      {
        for(int i=0; i<20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = seed->pwm[i][j];}}
        pwm[seed->width][0] = marginals[0]; pwm[seed->width][1] = marginals[1]; pwm[seed->width][2] = marginals[2]; pwm[seed->width][3] = marginals[3]; 
        seed->width = seed->width + 1; 
      }
      else{ for(int i=0; i < seed->width; i++){ for(int j=0; j<4; j++){ pwm[i][j] = seed->pwm[i][j];}} } 
    }
    else
    {
      for(int i=0; i<20; i++){ for(int j=0; j<4; j++){ pwm[i][j] = seed->pwm[i][j];}}
    }

    for(int i=0; i < seed->width; i++){ for(int j=0; j<4; j++){ seed->pwm[i][j] = pwm[i][j]; }} 
  }
}

double computeBayesFactor(struct4mapTypeDouble* freqs, double n[][4]) 
{
  string nucs         = "ACGT"; mapTypeCharInt nucleotide; nucleotide['A'] = 0; nucleotide['C'] = 1; nucleotide['G'] = 2; nucleotide['T'] = 3;
  double lnumerator   = 0;
  double ldenominator = 0;
  string dinuc;
  double marginals[4];    
  double total = 0;
  double N[4][4];     for(int i = 0; i < 4; i++){for(int j = 0; j < 4; j++){N[i][j] = max(0.0,n[i][j])+1;}}
  
  // set N to integers
  //for(int i=0; i<4; i++){ marginals[i] = 0; for(int j=0; j<4; j++){ N[j][i] = max(0.0,floor(N[j][i])); marginals[i] += N[j][i]; } total += marginals[i]; } 
  for(int i=0; i<4; i++){ marginals[i] = 0; for(int j=0; j<4; j++){ marginals[i] += N[j][i]; total += N[j][i]; }} 

  // lnumerator: assume alpha (a_A , a_C, a_G, and a_T) = 1 from the dirichlet  
  double a   = 1;
  lnumerator = log((double)3*2*1) - log((double)1); // just log(6) 
  for(int i=0; i<4; i++){for(double j=(marginals[i]+a-1); j>0; j--){ lnumerator += log(j); }}
  for(double i=(total+4-1); i>0; i--){lnumerator -= log(i);}

  // ldenominator
  for(int i = 0; i < 4; i++){for(int j = 0; j < 4; j++){ string dinuc = nucs.substr(i,1); dinuc.append(nucs.substr(j,1)); ldenominator += freqs->structMap[dinuc] * N[i][j]; }}

  double lBF = lnumerator - ldenominator;

  return(lBF);
} 

double chiSquared(struct4mapTypeDouble* freqs, double N[][4])
{ 
  string nucs         = "ACGT"; 
  double observed[4];
  double expected[4];
  double Nprime[4];
  double X2           = 0;

  int i, j;

  for(i=0; i<4; i++){observed[i]=0; expected[i]=0; Nprime[i]=0;}
  for(i=0; i<4; i++){for(j=0; j<4; j++){observed[i] += N[j][i]; Nprime[i] += N[i][j];}}
  for(i=0; i<4; i++){for(j=0; j<4; j++){ string dinuc = nucs.substr(j,1); dinuc.append(nucs.substr(i,1)); expected[i] += exp(freqs->structMap[dinuc]) * Nprime[i]; }}

  for(i=0; i<4; i++){ X2 += (pow(observed[i] - expected[i],2))/expected[i];}

  return(X2);
}

