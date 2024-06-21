#include <string>
#include <map>
using namespace std;

// Data types 
typedef map<string, string>              mapTypeString;
typedef map<string, int>                 mapTypeInt;
typedef map<string, double>              mapTypeDouble;
typedef map<char,   char>                mapTypeChar;
typedef map<char,   int>                 mapTypeCharInt;
typedef map<string, vector<signed int> > mapTypeVectorInt;
typedef map<string, mapTypeVectorInt>    nestedMapTypeVectorInt;

struct  wmerNeighborhood 
{
  string     name;

  int        width;
  int        counter;
  int        maxFlag, convergedFlag;
  int        totLen1, totLen2;
  int        indexOffset;
  int        flexiblePositions[20];

  double     difPropExactMatch;
  double     difPropBestMisMatches, difPropBestMMadjusted;
  double     difPropUpdated;
  double     logPWM[20][4];
  double     prevlogPWM[20][4];
  double     p1, p2, exactEnrichment, bestEnrichment, tempEnrichment, seqP1, seqP2; 
  double     exactLogLR, logLratio, FDR;
  double     threshold, pb;
  double     pwm[20][4];
  double     basalLLR;
  
  mapTypeInt neighborCnts1;
  mapTypeInt neighborCnts2;
  mapTypeInt tempCnts1;
  mapTypeInt tempCnts2;
  double     adjCnts1, adjCnts2, tempAdjCnts1, tempAdjCnts2;

  vector<signed int> nhoodPosSites;
  vector<signed int> nhoodNegSites;
  vector<signed int> tempPosSites;
  vector<signed int> tempNegSites;

  mapTypeVectorInt posSites; // these get used for exact matches at first then for sites
  mapTypeVectorInt negSites; // with high LR later
};

typedef map<string, wmerNeighborhood> mapTypeWmerNeighborhood;

struct  struct4mapTypeString{           mapTypeString           structMap; };
struct  struct4mapTypeInt{              mapTypeInt              structMap; };
struct  struct4mapTypeDouble{           mapTypeDouble           structMap; };
struct  struct4mapTypeWmerNeighborhood{ mapTypeWmerNeighborhood structMap; };


// functions
mapTypeString           getFastaSeqs(const char*);
mapTypeInt              countWmers(mapTypeString, int, int, int, int);
mapTypeVectorInt        findWmerSites(mapTypeString, int);
string                  reverseStrand(string);
int                     getNeighborhoodCounts(mapTypeInt, string, int, int, int);
mapTypeDouble           diffOfProps(mapTypeInt, mapTypeInt, mapTypeInt, int, int);
double                  singleDifOfProps(long double, long double, int, int);
void                    bubble_sort(double*, string*, int);
void                    add2pwm(double mat[][4], string, double);
double                  computeLogLikelihoodRatio(string, double pwm[][4], struct4mapTypeDouble*, int);
void                    cycleThruWchooseM(int*, string, int, int, int, struct4mapTypeInt*, struct4mapTypeInt*, wmerNeighborhood*, struct4mapTypeDouble*, mapTypeVectorInt*, mapTypeVectorInt*);  
void                    cycleThruNucleotides(string, int*, int, int, int, struct4mapTypeInt*, struct4mapTypeInt*, wmerNeighborhood*, struct4mapTypeDouble*, mapTypeVectorInt*, mapTypeVectorInt*);
void                    makePWM(wmerNeighborhood*, struct4mapTypeDouble*);
void                    updateMotif(wmerNeighborhood*, mapTypeString, mapTypeString, int, struct4mapTypeDouble*,struct4mapTypeDouble*,int,int,mapTypeVectorInt,mapTypeVectorInt, double,int);
void                    updateLogPWM(wmerNeighborhood*);
void                    updateLogPWM2(wmerNeighborhood*);
void                    writeOutMotifs(string*, struct4mapTypeWmerNeighborhood*, string, int, mapTypeString, mapTypeString, struct4mapTypeDouble*, struct4mapTypeDouble*,double,int);
string                  makeConsensus(string, struct4mapTypeWmerNeighborhood*);
mapTypeVectorInt        scanNfindInds(mapTypeString, wmerNeighborhood*,struct4mapTypeDouble*, double);
mapTypeVectorInt        scanNfindIndsFancy(mapTypeString, wmerNeighborhood*,struct4mapTypeDouble*,mapTypeVectorInt, double);
int                     updateSingleMotif(mapTypeVectorInt, mapTypeVectorInt,mapTypeString,mapTypeString,wmerNeighborhood*,struct4mapTypeDouble*,struct4mapTypeDouble*,int,int);
void                    determineMotifLength(wmerNeighborhood*,struct4mapTypeDouble*,double N0[][4],double N1[][4],double Nw[][4],double Nwp1[][4], int, int);
double                  computeBayesFactor(struct4mapTypeDouble*, double N[][4]);
double                  chiSquared(struct4mapTypeDouble*, double N[][4]);
int                     updateSeed(struct4mapTypeWmerNeighborhood*, mapTypeString, mapTypeString, wmerNeighborhood*,struct4mapTypeDouble*,struct4mapTypeDouble*, int, int, vector<string>&, int);
int                     assessOverlaps(struct4mapTypeWmerNeighborhood*, wmerNeighborhood*, vector<string>&);
double                  computeAdjustedDifProp(wmerNeighborhood*,int,int);  
