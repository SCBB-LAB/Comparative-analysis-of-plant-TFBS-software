/********************************************************************
 * FILE: iForm.c
 * AUTHOR: Chen Hebing
 * CREATE DATE: 4/27/3013
 * PROJECT: Combine fimo consensus RSAT Storm Homer
 * COPYRIGHT: 2013-2016
 ********************************************************************/
#define DEFINE_GLOBALS
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include "lib/matrix.h"
#include "lib/alphabet.h"
#include "lib/heap.h"
#include "lib/dir.h"
#include "lib/fasta-io.h"
#include "lib/hash_alph.h"
#include "lib/io.h"
#include "lib/meme-io.h"
#include "lib/projrel.h"
#include "lib/pssm.h"
#include "lib/prior-reader-from-psp.h"
#include "lib/reservoir.h"
#include "lib/seq-reader-from-fasta.h"
#include "lib/simple-getopt.h"
#include "lib/string-list.h"
#include "lib/utils.h"

#define ALPHABET_SIZE 4
#define BASE 2.71828182846   /* The base for taking logarithms. */
#define INFINITY_LARGE 1.0e100 
#define ROUND_TO_INT(n) (((n) > 0.0) ? ((int)((n) + 0.5)) : ((int)((n) - 0.5)))
#define SUM_LN(n1, n2) (((n1) >= (n2)) ? ((n1) + log(exp((n2) - (n1)) + 1.0)) : ((n2) + log(exp((n1) - (n2)) + 1.0)))
#define YES 1
#define NO 0	
#define con_Fudge 0.25
#define LN_PROB_SUM 0.0

#define MOTIF_PROBABILITY_MINIMUM 0.001

#define FASTA_READER_BUFSIZE 4096
#define ALPHABET "ACGT"
#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 
#endif
#define ASSERT(assertion, msg)                          \
{                                                       \
    if (!(assertion)) {                                 \
        fprintf(stderr, "Assertion error (%s)\n", msg); \
        exit(1);                                        \
    }                                                   \
} 
#define ERROR(...)                   \
{                                    \
    fprintf(stderr, "ERROR ");       \
    fprintf(stderr, __VA_ARGS__);    \
    fprintf(stderr, "\n");           \
    exit(1);                         \
}

double Multiple; //Consensus main value;
int scale_storm=1000;
double translation_storm;	
int total_discretized_scores=0;
  	
char* program_name = "iForm";
VERBOSE_T verbosity = NORMAL_VERBOSE;
const int SEQUENCE_INCREMENT_my = 50;

double *my_base;
	
// Structure for tracking fimo command line parameters.
typedef struct options {

  BOOLEAN_T allow_clobber;      // Allow overwritting of files in output directory.
  BOOLEAN_T compute_qvalues;    // Compute q-values
  BOOLEAN_T output_pthresh_set; // p-value threshold has been set from command line.
  BOOLEAN_T output_qthresh_set; // q-value threshold has been set from command line.
  BOOLEAN_T text_only;          // Generate only plain text output
  BOOLEAN_T scan_both_strands;  // Scan forward and reverse strands

  char* bg_filename;            // Name of file file containg background freq.
  char* command_line;           // Full command line
  char* meme_filename;          // Name of file containg motifs.
  char* motif_name;             // Use this motif name in the output.
  char* output_dirname;         // Name of the output directory
  char* seq_filename;           // Name of file containg sequences.
  char* seq_name;               // Use this sequence name in the output.

  int max_seq_length;    // Maximum allowed sequence length.
  int max_stored_scores; // Maximum number of matches to store per pattern.

  double alpha;       // Non-motif specific scale factor.
  double beta;        // Expected number of sites
  double pseudocount; // Pseudocount added to Motif PSFM.
  double pthresh;     // Maximum p-value to report.
  double qthresh;     // Maximum q-value to report.

  ALPH_T alphabet;    // Alphabet specified by MEME file.
  STRING_LIST_T* selected_motifs; // Indices of requested motifs.

  char *prior_distribution_filename; // Path to file containing prior distribution
  char *pval_lookup_filename;   // Print p-value lookup table.
  char *html_stylesheet_path; // Path to master copy of HTML stlesheet for CisML
  char *html_stylesheet_local_path; // Path to working copy of HTML stylesheet for CisML
  char *css_stylesheet_path; // Path to master copy of CSS style sheet for CisML
  char *css_stylesheet_local_path; // Path to working copy of CSS stylesheet for CisML
  char *text_stylesheet_path; // Path to plain-text stylesheet for CisML
  char *gff_stylesheet_path; // Path to GFF stylesheeet for CisML
  char *wiggle_stylesheet_path; // Path to wiggle stylesheeet for CisML
  char *html_path; // Path to FIMO HTML output file
  char *text_path; // Path to FIMO plain-text output file
  char *gff_path; // Path to FIMO GFF output file
  char *wiggle_path; // Path to FIMO wiggle output file
  char *fimo_path; // Pathe to FIMO XML output file
  char *cisml_path; // Path to CisML XML output file

  const char* HTML_STYLESHEET;  // Name of HTML XSLT stylesheet.
  const char* CSS_STYLESHEET;   // Name of CSS stylesheet.
  const char* GFF_STYLESHEET;   // Name of GFF XSLT stylesheet.
  const char* WIGGLE_STYLESHEET;   // Name of wiggle XSLT stylesheet.
  const char* TEXT_STYLESHEET;  // Name of plain-text XSLT stylesheet.
  const char* HTML_FILENAME;    // Name of HTML output file.
  const char* TEXT_FILENAME;    // Name of plain-text output file.
  const char* GFF_FILENAME;     // Name of GFF output file.
  const char* WIGGLE_FILENAME;  // Name of wiggle output file.
  const char* FIMO_FILENAME;    // Name of FIMO XML output file.
  const char* CISML_FILENAME;   // Name of CisML XML output file.

  const char* usage; // Usage statment

} FIMO_OPTIONS_T;
typedef struct pattern PATTERN_T;
typedef struct scanned_sequence SCANNED_SEQUENCE_T;
typedef struct matched_element MATCHED_ELEMENT_T;

struct pattern {

  char *accession; // Required.
  char *name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  char *db; // May be NULL.
  char *lsid; // May be NULL.

  int num_allocated_sequences; // Size of child scanned-sequence array.
  int num_allocated_elements; // Size of matched-element array.
  int num_sequences; // Count of child scanned-sequences.
  long num_scanned_positions; // Count of all positions scanned for
                            // matched_element.
  int num_stored_matches; // Count of the matched-elements stored in the elements array.
  int max_stored_matches;   // maximum number of matches to store for this pattern.
  double max_pvalue_retained; // Largest pvalue of retained records.

  BOOLEAN_T has_all_pvalues;  // Has retained all matched elements, not just those
                              // with best p-values.
  BOOLEAN_T qvalues_computed; // Have q-values been calcuatled for matched-elements
  BOOLEAN_T is_complete;     /// All matched elements have been added to pattern
  SCANNED_SEQUENCE_T **sequences; // Array of child scanned-sequence pointers.
  HEAP *element_heap; // Heap of matched elements ordered by ascending score.
  MATCHED_ELEMENT_T **elements; // Array of matched element pointers, ordered by p-value
};

struct scanned_sequence {

  char* accession; // Required.
  char* name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  int *length; // May be NULL.
  char* db; // May be NULL.
  char* lsid; // May be NULL.

  long num_scanned_positions; // Count of all positions scanned for
                            // matched_elements.
  int num_matched_elements; // Count of all matched_elements
                            // <= num_scanned_elements because of filtering.
  int num_allocated_elements; // Number of elements in elements array.

  MATCHED_ELEMENT_T **elements; // Array of all matched elements for this sequence
  PATTERN_T *parent_pattern; // Pointer to containing pattern.

};

struct matched_element {

  int start; // Required.
  int stop; // Required.

  double score;
  BOOLEAN_T has_score;
  double pvalue;
  BOOLEAN_T has_pvalue;
  double qvalue;
  BOOLEAN_T has_qvalue;
  char* clusterid; // May be NULL.
  char* sequence; // May be NULL.
  char strand;

  SCANNED_SEQUENCE_T *parent_sequence; // Pointer to containing scanned-sequence.

};
typedef struct pvalues_t Pvalues_t;
struct pvalues_t
{
    double *data;
    int size;
    double w_min;
    double w_max;
} ;
typedef struct Array ARRAY;
struct  Array
{
    double **data;
    int I;
    int J;
    double p;
    double pseudo;
    
};
typedef struct Markov MARKOV;
struct Markov 
{
    int order;
    double pseudo;     // pseudo-count (used for smoothing)
    double **T;        // transition matrix
    double *S;         // stationary vector
    int msize;         // transition matrix size (number of rows)
    double *priori;    // priori vector (pA, pC, pG, pT)
    double *logpriori; // log priori vector (pA, pC, pG, pT)

    int alphabet_size;
    
    double p;
    int prefix;
    int suffix;
    int idx;
};	
typedef struct values_t Values_t;
struct values_t
{
    int *data;
    int size;
    double e;
    double min;
    double max;
} ;
typedef struct fasta_reader_t Fasta_reader_t;
struct fasta_reader_t
{
    char *buffer;
    int bufsize;
    int pos;
    FILE *fp;
};
typedef struct seq_t Seq_t;
struct seq_t
{
    char *data;
    int size;
    int msize;
    char *name; // limited to 1024 chars
};


Values_t *new_values(double min, double max, double e);

void free_values(Values_t *values);
void free_Pvalues(Pvalues_t *pvalues);
void values_add(Values_t *values, double value);
void values_print(Pvalues_t *table,Values_t *values);
ARRAY Array_alloc(int dim1, int dim2);
void Array_dealloc(ARRAY AAA);
MARKOV bernoulli(double *priori);
ARRAY read_matrix_RSAT(char *filename, double pseudo);
MARKOV Markov_alloc(int o);
void Markov_dealloc(MARKOV m);  
void transform2logfreq(ARRAY AAA,MARKOV m);
Fasta_reader_t *new_fasta_reader(FILE *fp);
Seq_t *new_seq(int size);
void free_seq_RSAT(Seq_t *seq);
Seq_t *fasta_reader_next(Fasta_reader_t *reader,double *base_comp);
Pvalues_t *new_pvalues();
int scan_seq_pvalue(Seq_t *seq, ARRAY matrix, MARKOV bg, Values_t *values, double **homer_mat,
           Values_t *values_homer,int rc);
Seq_t *new_seq_rc(Seq_t *seq);
void seq_append_c(Seq_t *seq, char c);
double *scan_seq_RSAT(char *faline, ARRAY matrix, MARKOV bg, Pvalues_t *pvalues,int rc);
double *scan_seq_homer(char *faline, double **homer_mat,int L, Pvalues_t *pvalues,int rc);	
double *scan_seq_storm(char *faline, double **homer_mat,int L,double *index_pvalue,double *base_comp,int *alphabase,int rc);
		
SCANNED_SEQUENCE_T *allocate_scanned_sequence_my(
  char *accession,
  char *name,
  PATTERN_T *parent
);
MATCHED_ELEMENT_T *allocate_matched_element_with_score_my(
  int start,
  int stop,
  double score,
  double pvalue,
  SCANNED_SEQUENCE_T *parent
);
PATTERN_T *allocate_pattern_my(char *accession, char *name);
MATCHED_ELEMENT_T *allocate_matched_element_my(
  int start,
  int stop,
  SCANNED_SEQUENCE_T *parent
);

int get_pattern_num_stored_matches_my(PATTERN_T *pattern); 
MATCHED_ELEMENT_T **get_pattern_elements(PATTERN_T *pattern);
int get_elements_start(MATCHED_ELEMENT_T *elements);
int get_elements_stop(MATCHED_ELEMENT_T *elements);
char get_elements_strand(MATCHED_ELEMENT_T *elements);
double get_elements_score(MATCHED_ELEMENT_T *elements);
double get_elements_pvalue(MATCHED_ELEMENT_T *elements);
double get_elements_qvalue(MATCHED_ELEMENT_T *elements);
char *get_elements_sequence(MATCHED_ELEMENT_T *elements);
SCANNED_SEQUENCE_T *get_elements_parent_sequence(MATCHED_ELEMENT_T *elements);
char *get_scanned_sequence_name_my(SCANNED_SEQUENCE_T *scanned_sequence);
char *get_pattern_name_my(PATTERN_T* pattern);
double Array_logP(ARRAY AAA,char *word);
int Markov_word2index(MARKOV m,char *word, int len);
double Markov_logP(MARKOV m,char *word, int len);	   
double Markov_logPBernoulli(MARKOV m,char *word, int l);
char fasta_reader_getc(Fasta_reader_t *reader);  
int word_is_valid(Seq_t *seq, int start, int l);
double score2pvalue_RSAT(Pvalues_t *pvalues, double score);
double fisher_combine_pvalue(double *inpvalue,double chisq_k);
double Stouffer_combine_pvalue(double *inpvalue,double chisq_k,double *score);
//#######Consensus #######
//double inf_mat(double **align_mat,int A_size,int L);
//double adjust_inf(double **align_mat, double inf, int A_size ,int L);
double **weight_matrix(double **align_mat, int A_size ,int L);
double *det_marginal_prob(double **weight_mat, int T_SIZE,int width);
double *scan_seq_Consensus(char *faline,double **Align_mat,double *Marginal_prob,int rc,int L,int *alphabase); 

double *cal_storm_index_pvalue(double **storm_mat,int Align_matI,int Align_matJ,double *base_comp);
/**********************************************************************
  copy_matched_element_my

  Copy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
void *copy_matched_element_my(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  MATCHED_ELEMENT_T *new_e = allocate_matched_element_with_score_my(
      e->start,
      e->stop,
      e->score,
      e->pvalue,
      e->parent_sequence
    );
  new_e->qvalue = e->qvalue;
  new_e->clusterid = e->clusterid;
  new_e->sequence = e->sequence;
  new_e->strand = e->strand;
  return new_e;

}
/**********************************************************************
  compare_matched_elements_my
  Compare two objects, treating them as matched elements
**********************************************************************/
int compare_matched_elements_my(void *p1, void *p2) {
  MATCHED_ELEMENT_T *e1 = (MATCHED_ELEMENT_T *) p1;
  MATCHED_ELEMENT_T *e2 = (MATCHED_ELEMENT_T *) p2;

  if (e1->pvalue < e2->pvalue) {
    return 1;
  }
  else if (e1->pvalue > e2->pvalue) {
    return -1;
  }
  else {
    // If p-values are equal, compare staring postions
    // to break the time.
    return e1->start < e2->start ? 1: -1;
  }

}
/**********************************************************************
  destroy_matched_element_my

  Destroy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
void destroy_matched_element_my(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  free_matched_element(e);

}

/**********************************************************************
  add_pattern_scanned_sequence

  Adds a pointer to a scanned_sequence to the array of pointers to
  scanned_sequences in a cisml pattern object.
**********************************************************************/
static void add_pattern_scanned_sequence(
  PATTERN_T *pattern,
  SCANNED_SEQUENCE_T *sequence
) {

  assert(pattern != NULL);
  assert(sequence != NULL);
  assert(pattern->num_sequences <= pattern->num_allocated_sequences);

  sequence->parent_pattern = pattern;

  if (pattern->num_sequences == pattern->num_allocated_sequences) {
    pattern->num_allocated_sequences += SEQUENCE_INCREMENT_my;
    pattern->sequences = mm_realloc(
      pattern->sequences,
      pattern->num_allocated_sequences * sizeof(SCANNED_SEQUENCE_T *)
    );
  }
  pattern->sequences[pattern->num_sequences] = sequence;
  pattern->num_sequences++;

}


/***********************************************************************
  Print plain text record for motif site to standard output.
 ***********************************************************************/
static void print_site_as_text(
  const char *motif_id,
  const char *seq_name,
  const char *raw_seq,
  int start,
  int stop,
  char strand,
  double score,
  double pvalue,
  double qvalue
);

/***********************************************************************
  Free memory allocated in options processing
 ***********************************************************************/
/***********************************************************************
  Process command line options
 ***********************************************************************/
static void process_command_line(
  int argc,
  char* argv[],
  FIMO_OPTIONS_T *options
) {

  // Define command line options.
  const int num_options = 17;
  cmdoption const fimo_options[] = {
    {"alpha", REQUIRED_VALUE},
    {"beta", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"max-seq-length", REQUIRED_VALUE},
    {"max-stored-scores", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"motif-pseudo", REQUIRED_VALUE},
    {"norc", NO_VALUE},
    {"o", REQUIRED_VALUE},
    {"method", REQUIRED_VALUE},
    {"output-pthresh", REQUIRED_VALUE},
    {"output-qthresh", REQUIRED_VALUE},
    {"no-qvalue", NO_VALUE},
    {"prior-dist", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"pval-lookup", REQUIRED_VALUE} // This option is hidden from users.
  };

  // Define the usage message.
  options->usage =
    "USAGE: fimo [options] <motif file> <sequence file>\n"
    "\n"
    "   Options:\n"
    "     --alpha <double> (default 1.0)\n"
    "     --beta <double> (default 1.0)\n"
    "     --bgfile <background> (default from NR sequence database)\n"
    "     --max-seq-length <int> (default=2.5e8)\n"
    "     --max-stored-scores <int> (default=100000)\n"
    "     --motif <id> (default=all)\n"
    "     --motif-pseudo <float> (default=0.1)\n"
    "     --norc\n"
    "     --o <output dir> (default=fimo_out)\n"
    "     --method (default=Fisher)\n"
    "     --output-pthresh <float> (default 1e-4)\n"
    "     --output-qthresh <float> (default 1.0)\n"
    "     --no-qvalue\n"
    "     --prior-dist <PSP distribution filename> (default none)\n"
    "     --text\n"
    "     --verbosity [1|2|3|4] (default 2)\n"
    "\n"
    "   Use \'-\' for <sequence file> to read the database from standard input.\n"
    "   Use \'--bgfile motif-file\' to read the background from the motif file.\n"
    "\n";

  int option_index = 0;

  /* Make sure various options are set to NULL or defaults. */
  options->allow_clobber = TRUE;
  options->compute_qvalues = TRUE;
  options->output_pthresh_set = FALSE;
  options->output_qthresh_set = FALSE;
  options->text_only = FALSE;
  options->scan_both_strands = TRUE;

  options->bg_filename = NULL;
  options->command_line = NULL;
  options->meme_filename = NULL;
  options->motif_name = "motif";
  options->output_dirname = "Fisher";  //not dir Change to use for combine pvale method
  options->prior_distribution_filename = NULL;
  options->seq_filename = NULL;

  options->max_seq_length = MAX_SEQ;
  options->max_stored_scores = 100000;


  options->alpha = 1.0;
  options->beta = 1.0;
  options->pseudocount = 0.049;
  options->pthresh = 1e-4;
  options->qthresh = 1.0;

  options->selected_motifs = new_string_list();
  options->pval_lookup_filename = NULL;
  verbosity = 2;

  simple_setopt(argc, argv, num_options, fimo_options);

  // Parse the command line.
  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    }
    else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }
    if (strcmp(option_name, "bgfile") == 0){
      options->bg_filename = option_value;
    }
    else if (strcmp(option_name, "prior-dist") == 0){
      options->prior_distribution_filename = option_value;
    }
    else if (strcmp(option_name, "max-seq-length") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options->max_seq_length = (int)atof(option_value);
    }
    else if (strcmp(option_name, "alpha") == 0) {
      options->alpha = atof(option_value);
    }
    else if (strcmp(option_name, "beta") == 0) {
      options->beta = atof(option_value);
    }
    else if (strcmp(option_name, "max-seq-length") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options->max_seq_length = (int)atof(option_value);
    }
    else if (strcmp(option_name, "max-stored-scores") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options->max_stored_scores = (int)atof(option_value);
    }
    else if (strcmp(option_name, "motif") == 0){
      if (options->selected_motifs == NULL) {
        options->selected_motifs = new_string_list();
      }
      add_string(option_value, options->selected_motifs);
    }
    else if (strcmp(option_name, "motif-pseudo") == 0){
      options->pseudocount = atof(option_value);
    }
    else if (strcmp(option_name, "norc") == 0){
      options->scan_both_strands = FALSE;
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      //options->output_dirname = option_value;
      options->allow_clobber = FALSE;
    }
    else if (strcmp(option_name, "method") == 0){
      // Set output directory with clobber
      options->output_dirname = option_value;
      //options->allow_clobber = TRUE;
    }
    else if (strcmp(option_name, "output-pthresh") == 0){
      options->pthresh = atof(option_value);
      options->qthresh = 1.0;
      options->output_pthresh_set = TRUE;
    }
    else if (strcmp(option_name, "output-qthresh") == 0){
      options->qthresh = atof(option_value);
      options->pthresh = 1.0;
      options->output_qthresh_set = TRUE;
    }
    else if (strcmp(option_name, "no-qvalue") == 0){
      options->compute_qvalues = FALSE;
    }
    else if (strcmp(option_name, "text") == 0){
      options->text_only = TRUE;
    }
    else if (strcmp(option_name, "verbosity") == 0){
      verbosity = atoi(option_value);
    }
    else if (strcmp(option_name, "pval-lookup") == 0) {
      options->pval_lookup_filename = option_value;
    }
  }


  // Check that qvalue options are consistent
  if (options->compute_qvalues == FALSE && options->output_qthresh_set == TRUE) {
    die("The --no-qvalue option cannot be used with the --output-qthresh options");
  }

  // Turn off q-values if text only.
  if (options->text_only == TRUE) {
    if (options->compute_qvalues) {
      fprintf(stderr, "Warning: text mode turns off computation of q-values\n");
    }
    options->compute_qvalues = FALSE;
  }

  // Must have sequence and motif file names
  if (argc != option_index + 2) {
    fprintf(stderr, "%s", options->usage);
    exit(EXIT_FAILURE);
  }
  // Record the command line
  options->command_line = get_command_line(argc, argv);

  // Record the input file names
  options->meme_filename = argv[option_index];
  option_index++;
  options->seq_filename = argv[option_index];
  option_index++;

}

/**********************************************************psp***************
 * Calculate the log odds score for a single motif-sized window.
 *************************************************************************/
static inline BOOLEAN_T score_motif_site(
  char *seq,
  PSSM_T *pssm,
  double *pvalue, // OUT
  double *score // OUT
) {

  ARRAY_T* pv_lookup = pssm->pv;
  MATRIX_T* pssm_matrix = pssm->matrix;
  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);
  BOOLEAN_T scorable_site = TRUE;
  double scaled_log_odds = 0.0;

  // For each position in the site
  int motif_position;
  for (motif_position = 0; motif_position < pssm->w; motif_position++) {

    char c = seq[motif_position];
    int alph_index = alphabet_index(c, alphabet);

    // Check for gaps and ambiguity codes at this site
    if(c == '-' || c == '.' || alph_index >= alph_size) {
        scorable_site = FALSE;
        break;
    }

    scaled_log_odds += get_matrix_cell(motif_position, alph_index, pssm_matrix);
  }

  if (scorable_site == TRUE) {

    int w = pssm->w;
    *score = get_unscaled_pssm_score(scaled_log_odds, pssm);


    // Handle scores that are out of range
    if ((int) scaled_log_odds >= get_array_length(pv_lookup)) {
      scaled_log_odds = (float)(get_array_length(pv_lookup) - 1);
      *score = scaled_to_raw(scaled_log_odds, w, pssm->scale, pssm->offset);
    }
    *pvalue = get_array_item((int) scaled_log_odds, pv_lookup);

  }

  return scorable_site;

}


static void print_site_as_text(
  const char *motif_id,
  const char *seq_name,
  const char *raw_seq,
  int start,
  int stop,
  char strand,
  double score,
  double pvalue,
  double qvalue
) {
  fprintf(
    stdout, 
    "%s\t%s\t%d\t%d\t%c\t%g\t%.3g\t%.3g\t%s\n", 
    motif_id,
    seq_name,
    start,
    stop,
    strand,
    score,
    pvalue,
    qvalue,
    raw_seq
  );
}


/*************************************************************************
 * Entry point for fimo
 *************************************************************************/
int main(int argc, char *argv[]) {

  FIMO_OPTIONS_T options;
	
  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/
  process_command_line(argc, argv, &options);
  
//  double caoc=gsl_cdf_gaussian_P(4,1);   //pnorm in R
//  double kaoc=gsl_cdf_gaussian_Pinv(0.2,1);   //qnorm in R
//  printf("caoc\t%g\t%g\n",caoc,kaoc);
    	int i,j;
	//############RSAT input#####################
	char *seqfile = options.seq_filename;
	char *matfile = options.meme_filename;
    int rc = options.scan_both_strands;
    double precision = 0.495;
    double pseudo = 0.97;
//////////////my_base///////////////
	my_base=(double *)malloc(5*sizeof(double));
	for (i = 0; i < 5; i++)
		my_base[i]=0.0;
    FILE *fpbase;
    fpbase = fopen(seqfile, "r");
    if (fpbase == NULL)
        ERROR("unable to open '%s'", seqfile);	
	//while (!feof(fpbase))
	while (1)
	{
		char *infabase=(char *)malloc(1000 * sizeof(char));		
		if(!(fgets(infabase,1000,fpbase))){break;}
		//printf("%s",infabase);
		if(strstr(infabase,">")){continue;}
		int inlength=strlen(infabase);
		//printf("ca%d\n",inlength);
		for(i=0;i<inlength;i++)
		{
		     switch (infabase[i])
		    {
		        case 'a':
		        case 'A':
		            my_base[0]+=1.0;
		            break;
		        case 'c':
		        case 'C':
		            my_base[1]+=1.0;
		            break;
		        case 'g':
		        case 'G':
		            my_base[2]+=1.0;
		            break;
		        case 't':
		        case 'T':
		            my_base[3]+=1.0;
		            break;
		    }    			
		}
		free(infabase);
    }
	fclose(fpbase);
   //base_compsum=0.0;
   for (i = 0; i < 4; i++)
    	my_base[4]+=my_base[i];  //sum of base_comp number
   if(!(my_base[4]>0)){printf("fa file error!\n");exit(1);} //if fa file is empty
   for (i = 0; i < 4; i++)
   {
    	my_base[i]=my_base[i]/my_base[4];
		if(my_base[i]==0){my_base[i]=0.00001;}
    	//printf("%g\n",my_base[i]);
   }
          
          
////////////////////////my_base
    // set bg model
    MARKOV markov_my;
	double priori[4] = {0.25, 0.25, 0.25, 0.25};
	priori[0]=my_base[0];priori[1]=my_base[1];priori[2]=my_base[2];priori[3]=my_base[3];
    markov_my=bernoulli(priori);
    //printf("start bernoulli\n");
    // matrix
    ARRAY matrix_my;
    matrix_my=read_matrix_RSAT(matfile, pseudo);
    	 //
    //printf("%d\t%d\n", matrix_my.I,matrix_my.J);  I= 4  J= 21
    
   	//##################matrix for Consensus
   	double **Align_mat;  //Consensus input matrix
  	Align_mat = (double **)malloc(matrix_my.I * sizeof(double *));
	//int indexATCG[4]={0,3,1,2}; //ACGT  to ATCG

    for ( i= 0; i < matrix_my.I; i++)
    {
    	Align_mat[i] = (double *)malloc(matrix_my.J * sizeof(double));
        for (j = 0; j < matrix_my.J; j++)
        {
        	//Align_mat[i][j]=matrix_my.data[indexATCG[i]][j];
        	Align_mat[i][j]=matrix_my.data[i][j];
        	//printf("%g\t",Align_mat[i][j]);
        }
        //printf("\n");       
    } 
    
   	//##################matrix for Homer#######################
  	double *homer_total;
	homer_total = (double *)malloc(matrix_my.J * sizeof(double));  
    for (j = 0; j < matrix_my.J; j++)
    {  	homer_total[j]=0.0;
        for ( i= 0; i < matrix_my.I; i++)
			homer_total[j]+=matrix_my.data[i][j];
    }  //sum of each line   	
   	double **homer_mat;  //Consensus input matrix
  	homer_mat = (double **)malloc(matrix_my.I * sizeof(double *));
    for ( i= 0; i < matrix_my.I; i++)
    {
    	homer_mat[i] = (double *)malloc(matrix_my.J * sizeof(double));
        for (j = 0; j < matrix_my.J; j++)
	    {
        	homer_mat[i][j]=matrix_my.data[i][j]/homer_total[j];
				if (homer_mat[i][j] < MOTIF_PROBABILITY_MINIMUM)
					homer_mat[i][j] = MOTIF_PROBABILITY_MINIMUM;
				homer_mat[i][j] = log(homer_mat[i][j]/con_Fudge);
				//printf("%g\t",homer_mat[i][j]);
        }
		//printf("\n");    	
    }  
  	free(homer_total);    
      	
    transform2logfreq(matrix_my,markov_my);
    // values (distrib)
    /*int i,j;
        for ( j= 0; j < matrix_my.J; j++)
        {
            for (i = 0; i < matrix_my.I; i++)
                printf("%g\t",matrix_my.data[i][j]);
            printf("\n");
        } */
    Values_t *values_my = NULL;
    Values_t *values_homer = NULL;
    //if (distrib)
    values_my = new_values(-1000, 10000.0, precision);
    values_homer = new_values(-1000, 10000.0, precision);    
    // sequences
    FILE *fp;
    fp = fopen(seqfile, "r");
    if (fp == NULL)
        ERROR("unable to open '%s'", seqfile);
    Fasta_reader_t *reader = new_fasta_reader(fp);
	double *base_comp;	//for base_comp for storm
	//double base_compsum;
	base_comp=(double *)malloc(5*sizeof(double));
	
	for (i = 0; i < 5; i++)
		base_comp[i]=0.0;
	
    while (1)
    {
        Seq_t *seq = fasta_reader_next(reader,base_comp);
        if (seq == NULL)
            break;
        scan_seq_pvalue(seq, matrix_my, markov_my, values_my, homer_mat, values_homer,rc);
        free_seq_RSAT(seq);
    } 
   //base_compsum=0.0;
   for (i = 0; i < 4; i++)
    	base_comp[4]+=base_comp[i];  //sum of base_comp number
   if(!(base_comp[4]>0)){printf("fa file error!\n");exit(1);} //if fa file is empty
   for (i = 0; i < 4; i++)
   {
    	base_comp[i]=base_comp[i]/base_comp[4];
    	if(base_comp[i]==0){base_comp[i]=0.00001;}
    	//printf("%g\n",base_comp[i]);
   }
          
    Pvalues_t *pvalues_my = NULL; 
    pvalues_my = new_pvalues();
    values_print(pvalues_my,values_my);
   
    Pvalues_t *pvalues_homer = NULL;
    pvalues_homer = new_pvalues();          
    values_print(pvalues_homer,values_homer);
    
    //for (i = 0; i < pvalues_my->size; i++)
       //printf("%g\n",pvalues_my->data[i]);
    
    //for (i = 0; i < pvalues_homer->size; i++)
       //printf("%g\n",pvalues_homer->data[i]);    
    free_values(values_homer);    
    free_values(values_my);
    fclose(fp);
   

	int Align_matI=matrix_my.I;
	int Align_matJ=matrix_my.J;  
	 	
   //############Storm input##########################
   	double **storm_mat;  //Consensus input matrix
  	storm_mat = (double **)malloc(Align_matI * sizeof(double *));
  	double *storm_mat_sum;
  	storm_mat_sum = (double *)malloc(Align_matJ * sizeof(double));
  	for (j = 0; j < Align_matJ; j++)
  	{
  		storm_mat_sum[j]=0.0;
  		 for ( i= 0; i < Align_matI; i++)
  		 	 storm_mat_sum[j]+=Align_mat[i][j];
  	}
    for ( i= 0; i < Align_matI; i++)
    {
    	storm_mat[i] = (double *)malloc(Align_matJ * sizeof(double));
        for (j = 0; j < Align_matJ; j++)
        {
        	storm_mat[i][j]=log((Align_mat[i][j]+base_comp[i])/base_comp[i]/(storm_mat_sum[j]+1))/log(2);
        	//printf("%g\t",storm_mat[i][j]);
        }
        //printf("\n");       
    } 
    free(storm_mat_sum);
    double *storm_index_pvalue;
    storm_index_pvalue=cal_storm_index_pvalue(storm_mat,Align_matI,Align_matJ,base_comp);
   
   
   //############Consensus input##########################
for (i = 0; i < Align_matI; ++i){for (j = 0; j < Align_matJ; ++j){Align_mat[i][j] += my_base[i];}}
	//for (i = 0; i < Align_matI; ++i){for (j = 0; j < Align_matJ; ++j){printf("%.3f\t",Align_mat[i][j]);}printf("\n");}printf("\n");
	//double Inf_fudge = inf_mat(Align_mat,Align_matI,Align_matJ);
	Align_mat = weight_matrix(Align_mat,Align_matI,Align_matJ);	 
	
	
	//for (i = 0; i < Align_matI; ++i){for (j = 0; j < Align_matJ; ++j){printf("%.3f\t",Align_mat[i][j]);}printf("\n");}printf("\n");
	
	double *Marginal_prob;
	Marginal_prob=det_marginal_prob(Align_mat,Align_matI ,Align_matJ);
		

	//Consensus  matrix read ok  Pvalue OK####################
	  
	 int *alphabase;
	 alphabase=(int *)malloc(128*sizeof(int));
	 for(i=0;i<127;i++){alphabase[i]=-1;}
	 int alphaidx;
	 alphaidx=(int)'A';alphabase[alphaidx]=0;
	 alphaidx=(int)'T';alphabase[alphaidx]=3;
	 alphaidx=(int)'C';alphabase[alphaidx]=1;
	 alphaidx=(int)'G';alphabase[alphaidx]=2;
	 alphaidx=(int)'a';alphabase[alphaidx]=0;
	 alphaidx=(int)'t';alphabase[alphaidx]=3;
	 alphaidx=(int)'c';alphabase[alphaidx]=1;
	 alphaidx=(int)'g';alphabase[alphaidx]=2;
/*	  
	  alphabase[(int)'A']=1;
	  alphabase[(int)'T']=1;  	  
	  alphabase[(int)'C']=2;  	  
	  alphabase[(int)'G']=3; 	  
	  alphabase[(int)'a']=0;
	  alphabase[(int)'t']=1;  	   	  
	  alphabase[(int)'c']=2; 
	  int temp=(int)'c';
	  	 printf("%d\n",temp);  	   	  
	  alphabase[(int)'g']=3;     
*/		   	    	
     //printf("%d\t%d\t%d\t%d\n",alphabase[65],alphabase[84],alphabase[67],alphabase[71]);

  /**********************************************
   * Read the motifs.
   **********************************************/
  BOOLEAN_T has_reverse_strand;
  ARRAY_T* bg_freqs = NULL;
  FILE *cisml_file = NULL;

/////////////////////change meme
 if((options.bg_filename)&&(strcmp(options.bg_filename,"motif-file") == 0))
 {
	  FILE *fpmeme = fopen(options.meme_filename, "r");
	  if (fpmeme == NULL)
	    ERROR("unable to open '%s'", options.meme_filename);
	  char **memedata;
	  memedata = (char **)malloc(50 * sizeof(char *));
	  int memei=0;
	  while (1)
	  {
	  	memedata[memei]=(char *)malloc(100 * sizeof(char));
	    if(!(fgets(memedata[memei],100,fpmeme))){break;}
	    memei++;
	  }
	  fclose(fpmeme);
	  
	  FILE *fpmemenew = fopen(options.meme_filename, "w");
	  if (fpmemenew == NULL)
	    ERROR("unable to open '%s'", options.meme_filename);
	  for(i=0;i<memei-1;i++)
	  {
	  	  if(strstr(memedata[i],"A "))
	  	  {
	  	  	 fprintf(fpmemenew,"A %f C %f G %f T %f \n",base_comp[0],base_comp[1],base_comp[2],base_comp[3]); 
	  	  }
	  	  else
	  	  {
	  	  	  fprintf(fpmemenew,"%s",memedata[i]);
	  	  }
	  	  
	  }
	  fprintf(fpmemenew,"\n");
	  fclose(fpmemenew);
	 
	  for(i=0;i<memei;i++)
		free(memedata[i]);
	  free(memedata);   	 
 }

/////////////////////change meme  

  ARRAYLST_T *motifs = arraylst_create();
  read_meme_file2(
    options.meme_filename,
    options.bg_filename,
    options.pseudocount,
    REQUIRE_PSPM,  //need PSPMs
    motifs,
    NULL, //motif occurrences, not used
    &has_reverse_strand,
    &bg_freqs
  );

  // Doesn't include rev comp motifs
  int num_motif_names = arraylst_size(motifs);

  DATA_BLOCK_READER_T *fasta_reader 
    = new_seq_reader_from_fasta(options.seq_filename);


  PRIOR_DIST_T *prior_dist = NULL;
  if (options.prior_distribution_filename) {
    prior_dist = new_prior_dist(options.prior_distribution_filename);
  }

  // If motifs use protein alphabet we will not scan both strands
  ALPH_T alphabet_type = which_alphabet();
  if (alphabet_type == PROTEIN_ALPH) {
    options.alphabet = PROTEIN_ALPH;
    options.scan_both_strands = FALSE;
  }
  else {
    options.alphabet = DNA_ALPH;
  }

  if (options.scan_both_strands == TRUE) {
    // Set up hash tables for computing reverse complement
    setup_hash_alph(DNAB);
    setalph(0);
    // Correct background by averaging on freq. for both strands.
    average_freq_with_complement(bg_freqs);
    int alph_size = get_alph_size(ALPH_SIZE);
    normalize_subarray(0, alph_size, 0.0, bg_freqs);
    fill_in_ambiguous_chars(FALSE, bg_freqs);
    // Make reverse complement motifs.
    add_reverse_complements2(motifs);
  }

  // Print the text-only header line.
	 
    fprintf(stdout, "Motif");
    fprintf(stdout, "\tSeq");
    fprintf(stdout, "\tStart");
    fprintf(stdout, "\tStop");
    fprintf(stdout, "\tStrand");
    fprintf(stdout, "\tLog-odds");
    fprintf(stdout, "\tp-value");
    fprintf(stdout, "\tq-value");
    fprintf(stdout, "\tSite\n");

 
  // Create p-value sampling reservoir
  RESERVOIR_SAMPLER_T *reservoir = NULL;
  if (!options.text_only) {
    reservoir = new_reservoir_sampler(10000);
  }

  /**************************************************************
   * Score each of the sites for each of the selected motifs.
   **************************************************************/
  int motif_index = 0;
  int num_scanned_sequences = 0;
  long num_scanned_positions = 0;

  // Count of all motifs including rev. comp. motifs

	motif_index = 0;
    MOTIF_T* motif = (MOTIF_T *) arraylst_get(motif_index, motifs);
 
    char* motif_id = get_motif_id(motif);
    char* bare_motif_id = get_bare_motif_id(motif);

    if ((get_num_strings(options.selected_motifs) == 0)
      || (have_string(bare_motif_id, options.selected_motifs) == TRUE)) {

      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using motif %s of width %d.\n",
          motif_id,
          motif->length
        );
      }

      // Build PSSM for motif and tables for p-value calculation.
      // FIXME: the non-averaged freqs should be used for p-values
      PSSM_T* pos_pssm 
        = build_motif_pssm(
            motif, 
            bg_freqs, 
            bg_freqs, 
            prior_dist, 
            options.alpha,
            PSSM_RANGE, 
            0,    // no GC bins
            FALSE // make log-likelihood pssm
          );



      // If required, do the same for the reverse complement motif.
      MOTIF_T* rev_motif = NULL;
      PSSM_T* rev_pssm = NULL;
      if (options.scan_both_strands) {
        ++motif_index;
        rev_motif = (MOTIF_T *) arraylst_get(motif_index, motifs);
        motif_id = get_motif_id(rev_motif);
        // FIXME: the non-averaged freqs should be used for p-values
        rev_pssm 
          = build_motif_pssm(
              rev_motif, 
              bg_freqs, 
              bg_freqs, 
              prior_dist, 
              options.alpha,
              PSSM_RANGE, 
              0, // GC bins
              FALSE
            );
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(
            stderr,
            "Using motif %s of width %d.\n",
            motif_id,
            motif->length
          );
        }
      }
			 		       
      char strand = (which_alphabet() == PROTEIN_ALPH ? '.' : '+');

      // Create cisml pattern for this motif.
      PATTERN_T* pattern = NULL;
      if (!options.text_only) {
        pattern = allocate_pattern_my(bare_motif_id, bare_motif_id);
        set_pattern_max_stored_matches(pattern, options.max_stored_scores);
        set_pattern_max_pvalue_retained(pattern, options.pthresh);
      }

      // Read the FASTA file one sequence at a time.
      num_scanned_positions = 0L;
      num_scanned_sequences = 0L;
      while (
        go_to_next_sequence_in_data_block_reader(fasta_reader) != FALSE
      ) {

        char *fasta_seq_name = NULL;

        BOOLEAN_T fasta_result 
          = get_seq_name_from_data_block_reader(fasta_reader, &fasta_seq_name);

		  long num_positions = 0L;
		  // Create a scanned_sequence record and record it in pattern.
		  SCANNED_SEQUENCE_T* scanned_seq = NULL;
		  if (!options.text_only) {
		    scanned_seq =
		      allocate_scanned_sequence_my(fasta_seq_name, fasta_seq_name, pattern);
		  }
		  
			DATA_BLOCK_T *total_block = new_sequence_block(FASTA_READER_BUFSIZE);
			get_next_block_from_data_block_reader(fasta_reader, total_block);
			char *total_seq = get_sequence_from_data_block(total_block);
			
			//printf("%s\n",total_seq);
			
		  // Score and record each possible motif site in the sequence
		  //DATA_BLOCK_T *seq_block = new_sequence_block(pssm->w);
			int total_len;
			total_len=strlen(total_seq);
			//printf("%d\n",total_len);
		  //while (get_next_block_from_data_block_reader(fasta_reader, seq_block) != FALSE) 
		  	 
		 
		  
		  //####################RSAT cal###################
		  //add pvalue here
		  int rc;
		  if (rev_pssm != NULL){rc =1;} //both
		 
		  double *pvalue_RSAT;
		  pvalue_RSAT=scan_seq_RSAT(total_seq, matrix_my, markov_my,pvalues_my,rc); 
		  //printf("%g\n",pvalue_RSAT[0]);
		  //printf("%g\n",pvalue_RSAT[10]);
		  
		  //####################homer cal###################
		  double *pvalue_homer;
		  //printf("%d\t%d\n",Align_matJ,pvalues_homer->size);
		  pvalue_homer=scan_seq_homer(total_seq, homer_mat,Align_matJ,pvalues_homer,rc); 
			 

		 
		  //####################Consensus cal###################
		  //add pvalue here		  
		  double *pvalue_Con;
		  //printf("cao%d\n",Align_matJ);
		  pvalue_Con=scan_seq_Consensus(total_seq, Align_mat, Marginal_prob,rc,Align_matJ,alphabase); 
		  //printf("%g\n",pvalue_RSAT[0]);
		  //printf("%g\n",pvalue_RSAT[10]);		  

		  //####################storm cal###################
		  double *pvalue_storm;
		  //printf("%d\t%d\n",Align_matJ,pvalues_homer->size);
		  pvalue_storm=scan_seq_storm(total_seq, storm_mat,Align_matJ,storm_index_pvalue,base_comp,alphabase,rc); 
			 

		  char *raw_seq;
		  raw_seq=(char *)malloc((pos_pssm->w+1) * sizeof(char));
		  raw_seq[pos_pssm->w]='\0';
		  //char raw_seq[pos_pssm->w];
		  //printf("%s\n",raw_seq);
		  //####################fimo cal###################
		  	double *repvalue;
			repvalue=(double *)malloc(5*sizeof(double));
			double *inpvalue;
			inpvalue=(double *)malloc(5*sizeof(double));
		  int tti;
		  for(tti =0; tti<total_len - pos_pssm->w +1;tti++)
		  {
		    double pvalue = 0.0;
		    double score = 0.0;


		    // Track number of positions in sequence
		    ++num_positions;

		    // Get corresponding prior.
		    //char *raw_seq = get_sequence_from_data_block(seq_block);
		    //int start = get_start_pos_for_data_block(seq_block);
		    //char raw_seq[pos_pssm->w];
		    //raw_seq=mm_malloc((pos_pssm->w) * sizeof(char));
		    	//printf("%s\n",fimo_total_seq+tti);
		 	//char *raw_seq;
		    //raw_seq=(char *)malloc((pos_pssm->w) * sizeof(char));
		    int ti;	
		    for(ti=0;ti<pos_pssm->w;ti++)
		    {
		    	raw_seq[ti]=total_seq[tti+ti];
		    }
		    //memcpy(raw_seq, total_seq+tti, pos_pssm->w);
		    //fimo_total_seq++;
		    	//printf("%s\n",raw_seq);
		    int start=tti+1;
		    int stop = start + motif->length - 1;

		    // Score and record forward strand
		    BOOLEAN_T scoreable_site = score_motif_site( raw_seq,  pos_pssm, &pvalue, &score);

		    if (scoreable_site == TRUE) {
		      if (options.text_only != TRUE) {
		        add_scanned_sequence_scanned_element(scanned_seq);
		      }
		      //here combine pvalue
		      //score=pvalue_RSAT[tti];
		      //score=pvalue_Con[tti];
		      //score=pvalue_homer[tti];
		      //score=pvalue_storm[tti];
		      //score= pvalue;
		      //##############ok add RSAT here!######
			 inpvalue[0]=pvalue;
			 inpvalue[1]=pvalue_RSAT[tti];
			 inpvalue[2]=pvalue_Con[tti];
			 inpvalue[3]=pvalue_homer[tti];
			 inpvalue[4]=pvalue_storm[tti];		
			 
			 if(strcmp(options.output_dirname,"Stouffer")==0)
			 {
			 	 int caoi=0;
//for(caoi=0;caoi<5;caoi++){printf("%g\t",inpvalue[caoi]);}printf("\n");
			 	 pvalue=Stouffer_combine_pvalue(inpvalue,5.0,&score);
			 }
			 else
			 {
			 	 pvalue=fisher_combine_pvalue(inpvalue,5.0);
			 }
	 		 //pvalue=pvalue_homer[tti];		
	 		 //printf("%g\t%g\t%g\t%g\t%g\t%g\n",inpvalue[0],inpvalue[1],inpvalue[2],inpvalue[3],inpvalue[4],pvalue);
		      
		      if (reservoir != NULL) {
		        reservoir_sample(reservoir, pvalue);
		      }
		      if (pvalue > options.pthresh) {
		        if (options.text_only != TRUE) {
		          set_pattern_has_all_pvalues(pattern, FALSE);
		        }
		      }
		      else {
		      	  if (options.text_only == TRUE) {
		            print_site_as_text(
		              get_bare_motif_id(motif),
		              fasta_seq_name, 
		              raw_seq,
		              start,
		              stop,
		              '+',
		              score,
		              pvalue,
		                0
		            );
		       	 } else
		       	 {
			          MATCHED_ELEMENT_T *element = allocate_matched_element_with_score_my(
			            start,
			            stop,
			            score,
			            pvalue,
			            scanned_seq
			          );
			          BOOLEAN_T added = add_pattern_matched_element(pattern, element);
			          // If already have dropped elements of greater p-value
			          // we won't be keeping this one.
			          if (added == TRUE) {
			            set_matched_element_sequence(element, raw_seq);
			            set_matched_element_strand(element, strand);
			          }
			          else {
			            free_matched_element(element);
			          }		        	
		       	 }
		      	}
		    }

		    // Score and record reverse strand if appropriate.
		    if (rev_pssm != NULL) {

		      scoreable_site = score_motif_site(raw_seq, rev_pssm, &pvalue, &score);

		      if (scoreable_site == TRUE) {
		        if (options.text_only != TRUE) {
		          add_scanned_sequence_scanned_element(scanned_seq);
		        }
		        //here combine pvalue
		      //score=pvalue_RSAT[total_len - pos_pssm->w +1+tti];
		      //score=pvalue_Con[total_len - pos_pssm->w +1+tti];
		      //score=pvalue_homer[total_len - pos_pssm->w +1+tti];
		      //score=pvalue_storm[total_len - pos_pssm->w +1+tti];
		      //score= pvalue;
		      
		      //##############ok add RSAT here!######
               repvalue[0]=pvalue;
               repvalue[1]=pvalue_RSAT[total_len - pos_pssm->w +1+tti];
               repvalue[2]=pvalue_Con[total_len - pos_pssm->w +1+tti];
               repvalue[3]=pvalue_homer[total_len - pos_pssm->w +1+tti];
               repvalue[4]=pvalue_storm[total_len - pos_pssm->w +1+tti];
               
			 if(strcmp(options.output_dirname,"Stouffer")==0)
			 {
			 	 int caoi=0;
//for(caoi=0;caoi<5;caoi++){printf("%g\t",repvalue[caoi]);}printf("\n");
			 	 pvalue=Stouffer_combine_pvalue(repvalue,5.0,&score);
			 }
			 else
			 {
			 	 pvalue=fisher_combine_pvalue(repvalue,5.0);
			 }
			  //pvalue =pvalue_homer[total_len - pos_pssm->w +1+tti];             
	        
		       //printf("%g\t%g\t%g\t%g\t%g\t%g\n",repvalue[0],repvalue[1],repvalue[2],repvalue[3],repvalue[4],pvalue); 
		        if (reservoir != NULL) {
		          reservoir_sample(reservoir, pvalue);
		        }
		        if (pvalue > options.pthresh) {
		          if (options.text_only != TRUE) {
		            set_pattern_has_all_pvalues(pattern, FALSE);
		          }
		        }
		        else {
		          // Since we're using the reverse complemment motif
		          // convert sequence to reverse complment for output.
		          char *invcomp_seq = strdup(raw_seq);
		          if (invcomp_seq == NULL) {
		            die("Unable to allocate memory for RC sequence string while printing\n");
		          }
		          invcomp_dna(invcomp_seq, motif->length);	
		          	        	
		        if (options.text_only == TRUE) {
		            print_site_as_text(
		              get_bare_motif_id(motif),
		              fasta_seq_name, 
		              invcomp_seq,
		              start,
		              stop,
		              '-',
		              score,
		              pvalue,
		                0
		            );
		       	 } else
		       	 {
		            MATCHED_ELEMENT_T *element = allocate_matched_element_with_score_my(
		              stop,
		              start,
		              score,
		              pvalue,
		              scanned_seq
		            );
		            BOOLEAN_T added = add_pattern_matched_element(pattern, element);
		            // If already have dropped elements of greater p-value
		            // we won't be keeping this one.
		            if (added == TRUE) {
		               set_matched_element_sequence(element, invcomp_seq);
		               set_matched_element_strand(element, '-');
		            }
		            else {
		              free_matched_element(element);
		            }		       	 	 
		       	 }

 

		          myfree(invcomp_seq);
		        }
		      }
		    }
		   
		  }
		  free(pvalue_homer);
		  free(pvalue_Con);
		  free(pvalue_RSAT);
		  free(pvalue_storm);
		  free(inpvalue);
		  free(repvalue);
		  free(raw_seq);
		  free_data_block(total_block);

        myfree(fasta_seq_name);

      }  // All sequences parsed

      // The pattern is complete.
      if (!options.text_only) {
        set_pattern_is_complete(pattern);
      }

      // Compute q-values, if requested.
      if (options.compute_qvalues) {
        int num_samples = get_reservoir_num_samples_retained(reservoir);
        ARRAY_T *sampled_pvalues = allocate_array(num_samples);
        fill_array(get_reservoir_samples(reservoir), sampled_pvalues);
        pattern_calculate_qvalues(pattern, sampled_pvalues);
        free_array(sampled_pvalues);
      }
    
      if (!options.text_only)
      {
	     int num_stored_matches = get_pattern_num_stored_matches_my(pattern);
	      int num_ii,start,stop;
	      char ostrand;
	      double score,pvalue,qvalue;
	      char *motif_id=NULL;
	      char *seq_name=NULL;
	      char *raw_seq=NULL;
	      SCANNED_SEQUENCE_T *my_scanned_sequence;
	      motif_id=get_pattern_name_my(pattern);  
	      MATCHED_ELEMENT_T **o_element=get_pattern_elements(pattern); 
	      for(num_ii=0;num_ii<num_stored_matches;num_ii++)
	      {
	      	  //motif_id=pattern->elements[num_ii]->start;
	      	  //seq_name=pattern->elements[num_ii]->stop; 
	      	  //raw_seq=pattern->elements[num_ii]->stop;  
	      	   	       
	      	  start=get_elements_start(o_element[num_ii]);
	      	  stop=get_elements_stop(o_element[num_ii]); 
	      	  ostrand=get_elements_strand(o_element[num_ii]);
	      	  score=get_elements_score(o_element[num_ii]); 	  
	      	  pvalue=get_elements_pvalue(o_element[num_ii]);
	      	  qvalue=get_elements_qvalue(o_element[num_ii]);
	      	  raw_seq=get_elements_sequence(o_element[num_ii]);  	  
  
	      	  my_scanned_sequence=get_elements_parent_sequence(o_element[num_ii]);
	      	  seq_name=get_scanned_sequence_name_my(my_scanned_sequence);
	      	  
	             print_site_as_text(
	              motif_id,
	              seq_name, 
	              raw_seq,
	              start,
	              stop,
	              ostrand,
	              score,
	              pvalue,
	              qvalue  
	            );     	  
	      }	  
      }
 

      // Free memory associated with this motif.
      free_pssm(pos_pssm);
      free_pssm(rev_pssm);
    }
    else {
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Skipping motif %s.\n", motif_id);
      }
    }

    if (!options.text_only) {
      clear_reservoir(reservoir);
    }




  // Clean up.

  // Close input readers
  close_data_block_reader(fasta_reader);

  free_data_block_reader(fasta_reader);

  if (reservoir != NULL) {
    free_reservoir(reservoir);
  }

  if (prior_dist != NULL) {
    free_prior_dist(prior_dist);
  }

  free(base_comp);
  free_motifs(motifs);
  free_array(bg_freqs);
  free_string_list(options.selected_motifs);
  free_Pvalues(pvalues_my);
  free_Pvalues(pvalues_homer);
  free(storm_index_pvalue);
  for(i=0;i<Align_matI;i++)
	free(Align_mat[i]);
  free(Align_mat);
  for(i=0;i<Align_matI;i++)
	free(homer_mat[i]);
  free(homer_mat);
  for(i=0;i<Align_matI;i++)
	free(storm_mat[i]);
  free(storm_mat);  
  Array_dealloc(matrix_my);
  Markov_dealloc(markov_my);
  

  return 0;

}

/**********************************************************************
  allocate_scanned_sequence_my

  Constructor for the cisml scanned_sequence data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
SCANNED_SEQUENCE_T *allocate_scanned_sequence_my(
  char *accession,
  char *name,
  PATTERN_T *parent_pattern
) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  SCANNED_SEQUENCE_T *scanned_sequence = mm_malloc(sizeof(SCANNED_SEQUENCE_T));
  scanned_sequence->accession = NULL;
  scanned_sequence->name = NULL;
  scanned_sequence->pvalue = NULL;
  scanned_sequence->score = NULL;
  scanned_sequence->length = NULL;
  scanned_sequence->db = NULL;
  scanned_sequence->lsid = NULL;
  scanned_sequence->num_scanned_positions = 0L;
  scanned_sequence->num_matched_elements = 0;
  scanned_sequence->num_allocated_elements = 0;

  // Set required fields
  int length = strlen(accession) + 1;
  scanned_sequence->accession = mm_malloc(length * sizeof(char));
  strncpy(scanned_sequence->accession, accession, length);
  length = strlen(name) + 1;
  scanned_sequence->name = mm_malloc(length * sizeof(char));
  strncpy(scanned_sequence->name, name, length);
  add_pattern_scanned_sequence(parent_pattern, scanned_sequence);
  scanned_sequence->elements = NULL;

  return scanned_sequence;
}
/**********************************************************************
  allocate_matched_element_with_score_my

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Sets score and pvalue.
  Other fields set to NULL.

**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_with_score_my(
  int start,
  int stop,
  double score,
  double pvalue,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = allocate_matched_element_my(start, stop, parent_sequence);

  set_matched_element_score(element, score);
  set_matched_element_pvalue(element, pvalue);

  return element;
}
/**********************************************************************
  allocate_pattern_my

  Constructor for the cisml pattern data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
PATTERN_T *allocate_pattern_my(char *accession, char *name) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  PATTERN_T *pattern = mm_malloc(sizeof(PATTERN_T));
  pattern->accession = NULL;
  pattern->name = NULL;
  pattern->pvalue = NULL;
  pattern->score = NULL;
  pattern->db = NULL;
  pattern->lsid = NULL;
  pattern->sequences = NULL;
  pattern->elements = NULL;

  pattern->num_allocated_sequences = 0;
  pattern->num_allocated_elements = 0;

  pattern->num_sequences = 0;
  pattern->num_scanned_positions = 0L;
  pattern->max_stored_matches = 100000;
  pattern->num_stored_matches = 0;
  pattern->max_pvalue_retained = 1.0;
  pattern->qvalues_computed = FALSE;
  pattern->has_all_pvalues = TRUE;
  pattern->is_complete = FALSE;

  // Set required fields
  int length = strlen(accession) + 1;
  pattern->accession = mm_malloc(length * sizeof(char));
  strncpy(pattern->accession, accession, length);
  length = strlen(name) + 1;
  pattern->name = mm_malloc(length * sizeof(char));
  strncpy(pattern->name, name, length);
  pattern->element_heap = create_heap(
    pattern->max_stored_matches,
    compare_matched_elements_my,
    copy_matched_element_my,
    destroy_matched_element_my,
    NULL, // Key function
    NULL  // Print function
  );
  pattern->elements = NULL;

  return pattern;
}

/**********************************************************************
  allocate_matched_element_my

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.
**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_my(
  int start,
  int stop,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = mm_malloc(sizeof(MATCHED_ELEMENT_T));
  element->start = start;
  element->stop = stop;
  element->parent_sequence = parent_sequence;

  // Initialze optional fields
  element->score = 0.0;
  element->has_score = FALSE;
  element->pvalue = 0.0;
  element->has_pvalue = FALSE;
  element->qvalue = 0.0;
  element->has_qvalue = FALSE;
  element->clusterid = NULL;
  element->sequence = NULL;
  element->strand = '\0';

  return element;
}
/**********************************************************************
  get_pattern_num_stored_matches_my

  Gets the total number of matched element objects stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_num_stored_matches_my(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->num_stored_matches;
}
MATCHED_ELEMENT_T **get_pattern_elements(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->elements;
};

int get_elements_start(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->start;	
};
int get_elements_stop(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->stop;	
};	
char get_elements_strand(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->strand;	
};
double get_elements_score(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->score;	
};
double get_elements_pvalue(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->pvalue;	
};
double get_elements_qvalue(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->qvalue;	
};
char *get_elements_sequence(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->sequence;	
};
SCANNED_SEQUENCE_T *get_elements_parent_sequence(MATCHED_ELEMENT_T *elements)
{
  assert(elements != NULL);
  return elements->parent_sequence;	
};
/**********************************************************************
  get_scanned_sequence_name_my

  Gets the program_name member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_name_my(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->name;

}
/**********************************************************************
  get_pattern_name_my
  Gets the program_name member from a cisml pattern object.
**********************************************************************/
char *get_pattern_name_my(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->name;
}


Values_t *new_values(double min, double max, double e)
{
    Values_t *values = (Values_t *) malloc(sizeof(Values_t));
    values->min = min;
    values->max = max;
    values->e = e;
    values->size = (int) (values->e + (values->max - values->min) / values->e);
    values->data = (int *) malloc(sizeof(int) * values->size);
    int i;
    for (i = 0; i < values->size; i++)
        values->data[i] = 0;

    return values;
}
Pvalues_t *new_pvalues()
{
    Pvalues_t *pvalues = (Pvalues_t *) malloc(sizeof(Pvalues_t));
    pvalues->w_min = 1E300;
    pvalues->w_max = -1E300;
    pvalues->size = 0;
    pvalues->data = (double *) malloc(sizeof(double) * 1);
    return pvalues;
}

void free_values(Values_t *values)
{
    free(values->data);
    free(values);
}
void free_Pvalues(Pvalues_t *pvalues)
{
    free(pvalues->data);
    free(pvalues);
}
void values_add(Values_t *values, double value)
{
    ASSERT(value > values->min && value < values->max, "score out of range");
    int i = (int) ((value - values->min) / values->e);
    ASSERT(i >= 0 && i < values->size, "score out of table range");
    values->data[i] += 1;
}

void values_print(Pvalues_t *table,Values_t *values)
{
    int i;
    // find [a,b]
    i = 0;
    while (i < values->size - 1 && values->data[i] == 0)
        i++;
    int a = i;
    i = values->size-1;
    while (i > 0 && values->data[i] == 0)
        i--;
    int b = i;
    
    //printf("%d\t%d\t%d\n",a,b,values->size);
    assert(a >= 0 && b >= 0 && a < values->size && b < values->size);

    // total
    int total = 0;
    for (i = a; i <= b; i++)
        total += values->data[i];
    //printf("%d\t%d\t%d\n",a,b,total);
    // print
    //fprintf(fout, "#score\tocc\tco\tcco\tdCDF\n");
    float weight, Pval;
    //Pvalues_t *table = new_pvalues();   
    
    int cum = 0;
    for (i = a; i <= b; i++)
    {
        cum += values->data[i];
 
        weight= values->min + i * values->e;
        
        Pval= (total - cum) / (double) total;
        if(Pval == 0)
        {
        	Pval= 0.1 / (double) total;
        }
        table->w_min = MIN(table->w_min, weight);
        table->w_max = MAX(table->w_max, weight);
        table->data = (double *) realloc(table->data, sizeof(double) * (table->size + 1));
        table->data[table->size] = Pval;
        table->size += 1;        
    }
   // return table;
}
    ARRAY Array_alloc(int dim1, int dim2)
    {
    	ARRAY AAA;
    	double val=0.0;
        AAA.I = dim1;
        AAA.J = dim2;
        if (AAA.I <= 0 || AAA.J <= 0)
            return;

        // alloc
        AAA.data = (double **)malloc((AAA.I) * sizeof(double *));
        int i,j;
        for (i = 0; i < AAA.I; i++)
            AAA.data[i] = (double *)malloc((AAA.J) * sizeof(double));
        
        // init
        for (i = 0; i < AAA.I; i++){
            for (j = 0; j < AAA.J; j++)
                AAA.data[i][j] = val;
        }
        return AAA;
    }

    MARKOV Markov_alloc(int o)
    {
    	MARKOV m;
    	double p=1.0;
        m.pseudo = 1.0;
        m.order = o;
        m.msize = (int) pow((double) ALPHABET_SIZE, m.order);
        m.priori = (double *)malloc((ALPHABET_SIZE) * sizeof(double));
        m.logpriori = (double *)malloc((ALPHABET_SIZE) * sizeof(double));

        m.S = (double *)malloc((m.msize) * sizeof(double));
        m.T = (double **)malloc((m.msize) * sizeof(double *));
        int i,j;
        for (i = 0; i < m.msize; i++)
            m.T[i] = (double *)malloc((ALPHABET_SIZE) * sizeof(double));
        // init values
        for (i = 0; i < m.msize; i++)
        {
            m.S[i] = 1.0;
            for (j = 0; j < ALPHABET_SIZE; j++)
                m.T[i][j] = 0.0;
        }
        return m;
    }	
	void Array_dealloc(ARRAY AAA)
    {
        if (AAA.data != NULL)
        {
        	int i;
            for (i = 0; i < AAA.I; i++)
                free(AAA.data[i]);
            free(AAA.data);
        }
    }
    
	void Markov_dealloc(MARKOV m)
    {
        if (m.order != -1)
        {
            free(m.priori);
            free(m.logpriori);
            int i;
             for (i = 0; i < m.msize; i++)
                free (m.T[i]);
            free(m.T);
            free(m.S);     
        }
    }
 
MARKOV bernoulli(double *priori)
{
    //Markov_dealloc(m);  
    MARKOV m; 
    m=Markov_alloc(0);
    m.S[0] = 1.0;
    int i;
    for (i = 0; i < ALPHABET_SIZE; i++)
    {
        m.priori[i] = priori[i];
        m.logpriori[i] = log(priori[i]);
        m.T[0][i] = priori[i];
    }
    return m;
}

ARRAY read_matrix_RSAT(char *filename, double pseudo)
{
    // open file
    ARRAY matrix_my;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
        ERROR("unable to open '%s'", filename);

    // read data
    float p[4][256];
    //char base[4];
    int l = 0;
    int c = -1;
    
    char mtemp[30];  //finish title
    while (!feof(fp))
    {
    	fscanf(fp, "%s", mtemp);
		if(strcmp(mtemp,"letter-probability") ==0){break;}
    }
    char mtemp2[100];
    fgets(mtemp2,100,fp);  //finish "letter-probability" line  go to next;
    //printf("%s\n",mtemp2);
    char *pnsites=strstr(mtemp2,"nsites=");
    char *pE=strstr(mtemp2,"E=");
    char nvalue[10];
    int ii;
    for(ii=0;ii<strlen(pnsites)-strlen(pE)-7;ii++)
    {
    	nvalue[ii]=pnsites[ii+7];
    }
    nvalue[ii]='\0';
    float nsites;
    nsites=atof(nvalue);
    //printf("%f\n",nsites);

    //printf("pre %s\n",mtemp2);
    //start to read matrix_my
    while (!feof(fp))
    {
        c = -1;
        while (++c < 4)
        {
            float v;
            if (fscanf(fp, "%f", &v) != 1)
                break;
            p[c][l] = v;
        }
        if(c==4)
        {
        	l++;	
        }
    }
    c=4;
    
    matrix_my=Array_alloc(c,l);
    matrix_my.pseudo = pseudo;

    int i,j;
    for (i = 0; i < l; i++)
    {
        for (j = 0; j < c; j++)
     	{
        	matrix_my.data[j][i] = p[j][i]*nsites;
        	//printf("%f\t",p[j][i]);
        }
         //printf("\n");   
    }

    // security check

    // close stream
    fclose(fp);
    
    return matrix_my;
}

    void transform2logfreq(ARRAY AAA,MARKOV m)
    {
		int i,j;
        for (j = 0; j < AAA.J; j++)
        {
            double N = 0;
            for (i = 0; i < AAA.I; i++)
                N += AAA.data[i][j];
            for (i = 0; i < AAA.I; i++)
            {
            	AAA.data[i][j] = log((AAA.data[i][j] + m.priori[i] * AAA.pseudo) / (N + AAA.pseudo));
            }
        }
    }


Fasta_reader_t *new_fasta_reader(FILE *fp)
{
    assert(fp != NULL);
    Fasta_reader_t *fasta_reader = (Fasta_reader_t *) malloc(sizeof(Fasta_reader_t));
    fasta_reader->bufsize = FASTA_READER_BUFSIZE;
    fasta_reader->pos = 0;
    fasta_reader->buffer = (char *) malloc(sizeof(char) * fasta_reader->bufsize);
    fasta_reader->fp = fp;
    return fasta_reader;
}

Seq_t *new_seq(int size)
{
    Seq_t *seq = (Seq_t *) malloc(sizeof(Seq_t));
    seq->size = 0;
    seq->msize = size;
    seq->data = (char *) malloc(sizeof(char) * seq->msize);
    seq->name = (char *) malloc(sizeof(char) * 1024);
    return seq;
}
void free_seq_RSAT(Seq_t *seq)
{
    free(seq->data);
    free(seq->name);
    free(seq);
}

Seq_t *fasta_reader_next(Fasta_reader_t *reader,double *base_comp)
{
    Seq_t *seq = new_seq(1000);

    // read header
    char c;
    int i = 0;
    int short_name_found = 0;
    do 
    {
        c = fasta_reader_getc(reader);
        if ((c == '\t' || c == ' ' || c == '\r') && i > 0)
            short_name_found = 1;
        if (!short_name_found && i < 1024 && \
            c!= EOF && c != '\n' && c != '\t' && c != ' ' && c != '\r' && c != '>')
            seq->name[i++] = c;
    } while (c != EOF && c != '\n');
    seq->name[i] = '\0';
    if (c == EOF)
        return NULL;

    // read sequence
    do 
    {
        c = fasta_reader_getc(reader);
        if (c == EOF)
          break;
        if (c != '\n' && c != '>')
       {
       	   seq_append_c(seq, c);
		     switch (c)
		    {
		        case 'a':
		        case 'A':
		            base_comp[0]+=1.0;
		            break;
		        case 'c':
		        case 'C':
		            base_comp[1]+=1.0;
		            break;
		        case 'g':
		        case 'G':
		            base_comp[2]+=1.0;
		            break;
		        case 't':
		        case 'T':
		            base_comp[3]+=1.0;
		            break;
		    }      	   
       }  
            

    } while (c != '>');

    return seq;
}


int scan_seq_pvalue(Seq_t *seq, ARRAY matrix, MARKOV bg, Values_t *values, double **homer_mat,
           Values_t *values_homer,int rc)
{
    int l = matrix.J;
    ASSERT(l < 256, "invalid matrix size");

    Seq_t *seqrc = NULL;
    if (rc)   //2str
        seqrc = new_seq_rc(seq);
    int maxpos = seq->size - l;  //maxpos=seq length - motif length
    //printf("%d\n", maxpos);
    int i,j;
    for (i = 0; i <= maxpos; i++)
    {
        if (!word_is_valid(seq, i, l))
            continue;
        double W;
        if (bg.order == 0)
            W = Array_logP(matrix ,&seq->data[i]) - Markov_logPBernoulli(bg,&seq->data[i], l);
        else
            W = Array_logP(matrix,&seq->data[i]) - Markov_logP(bg,&seq->data[i], l);
		
		//########for homer##########
		double h_score=0.0;
		char *h_word=&seq->data[i];
		for (j=0;j<l;j++) {
			h_score += homer_mat[(int) h_word[j]][j];
		}
		//printf("%g\t%g\n",h_score,W);
            if (values_homer != NULL)
            {
                values_add(values_homer, h_score);
            }
        //###################    		
		
		    if (values != NULL)
            {
                values_add(values, W);
            }


        if (!rc)
            continue;

        double Wrc;
        if (bg.order == 0)
            Wrc =Array_logP(matrix,&seqrc->data[maxpos - i]) - Markov_logPBernoulli(bg,&seqrc->data[maxpos - i], l);
        else
            Wrc = Array_logP(matrix,&seqrc->data[maxpos - i]) - Markov_logP(bg,&seqrc->data[maxpos - i], l);

		//########for homer##########
		double hr_score=0.0;
		char *hr_word=&seqrc->data[maxpos - i];
		for (j=0;j<l;j++) {
			hr_score += homer_mat[(int) hr_word[j]][j];
		}
		//printf("%g\t%g\n",hr_score,Wrc);
            if (values_homer != NULL)
            {
                values_add(values_homer, hr_score);
            }
        //################### 
        

            if (values != NULL)
            {
                values_add(values, Wrc);
            }
    }  
    if (rc)
        free_seq_RSAT(seqrc);
    
    return 1;
}

Seq_t *new_seq_rc(Seq_t *seq)
{
    Seq_t *rc = new_seq(seq->size);
    int i;
    for (i = 0; i < seq->size; i++)
    {
        rc->data[seq->size - i - 1] = 3 - seq->data[i];
    }
    return rc;
}
    double Array_logP(ARRAY AAA,char *word)
    {
        AAA.p = 0.0;
        int i;
        for (i = 0; i < AAA.J; i++)
            AAA.p += AAA.data[(int) word[i]][i];
        return AAA.p;
    }
    double Markov_logPBernoulli(MARKOV m,char *word, int l)
    {
        double s = 0.0;
        int i;
        for (i = 0; i < l; i++)
            s += m.logpriori[(int) word[i]];
        return s;
    }
    double Markov_logP(MARKOV m,char *word, int len)
    {
        if (m.order == 0)
        {
            double p = 0.0;
            int i;
            for (i = 0; i < len; i++)
                //p += log(priori[(int) word[i]]);
                p += m.logpriori[(int) word[i]];
            return p;
        }
        else
        {
            m.prefix = Markov_word2index(m,word, m.order);
            m.p = log(m.S[m.prefix]);
            int i;
            for (i = m.order; i < len; i++)
            {
                m.suffix = (int) word[i]; 
                m.prefix = Markov_word2index(m,&word[i-m.order], m.order);
                m.p += log(m.T[m.prefix][m.suffix]);
            }
            return m.p;
        }
    }    
    int Markov_word2index(MARKOV m,char *word, int len)
    {
        m.idx = 0;
        int P = 1;
        int i;
        for (i = 0; i < len; i++)
        {
            m.idx += (int) word[len-i-1] * P;
            P = P * ALPHABET_SIZE;
        }
        return m.idx;
    }    
char fasta_reader_getc(Fasta_reader_t *reader)
{
    return fgetc(reader->fp);
}
int word_is_valid(Seq_t *seq, int start, int l)
{
    int i;
    for (i = start; i < start + l; i++)
    {
        if (seq->data[i] < 0)
            return 0;
    }
    return 1;
}  

void seq_append_c(Seq_t *seq, char c)
{
    if (seq->size >= seq->msize)
    {
        seq->msize += 8;
        seq->data = (char *) realloc(seq->data, sizeof(char) * seq->msize);
    }
    switch (c)
    {
        case 'a':
        case 'A':
            seq->data[seq->size++] = 0;
            break;

        case 'c':
        case 'C':
            seq->data[seq->size++] = 1;
            break;

        case 'g':
        case 'G':
            seq->data[seq->size++] = 2;
            break;

        case 't':
        case 'T':
            seq->data[seq->size++] = 3;
            break;

        case 'n':
        case 'N':
            seq->data[seq->size++] = -1;
            break;
    }
} 
 	   
double *scan_seq_RSAT(char *faline, ARRAY matrix, MARKOV bg, Pvalues_t *pvalues,int rc)
{
	int faline_len=strlen(faline);
	//printf("%d\n",faline_len);
    int l = matrix.J;
    ASSERT(l < 256, "invalid matrix size");
    
    //make seq	
	Seq_t *seq = new_seq(faline_len); 
	int iii;
	char ccc;
	for(iii=0;iii<faline_len;iii++)
    {
        ccc = faline[iii];
        
        seq_append_c(seq, ccc);
    } 
    //make seqrc
    Seq_t *seqrc = NULL;

    int maxpos = faline_len - l;  //maxpos=seq length - motif length
    
    double *outpvalue;
    if (rc)   //2str
    {
    	seqrc = new_seq_rc(seq);    
    	outpvalue=(double *)malloc((2*maxpos+2)*sizeof(double));
    }
    else
    {
    	outpvalue=(double *)malloc((maxpos+1)*sizeof(double));
    }
    int i;
    for (i = 0; i <= maxpos; i++)
    {
        if (!word_is_valid(seq, i, l))
            continue;
            	
        double W;
        if (bg.order == 0)
            W = Array_logP(matrix ,&seq->data[i]) - Markov_logPBernoulli(bg,&seq->data[i], l);
        else
            W = Array_logP(matrix,&seq->data[i]) - Markov_logP(bg,&seq->data[i], l);

        double Pval = score2pvalue_RSAT(pvalues, W);

		outpvalue[i]=Pval;
		//printf("%g\t%g\n",W,Pval);
        if (!rc)
            continue;
        double Wrc;
        if (bg.order == 0)
            Wrc =Array_logP(matrix,&seqrc->data[maxpos - i]) - Markov_logPBernoulli(bg,&seqrc->data[maxpos - i], l);
        else
            Wrc = Array_logP(matrix,&seqrc->data[maxpos - i]) - Markov_logP(bg,&seqrc->data[maxpos - i], l);

        double Pval_rc = score2pvalue_RSAT(pvalues, Wrc);
        
        //outpvalue[2*maxpos+1-i]=Pval_rc;
        outpvalue[maxpos+1+i]=Pval_rc;
        //printf("%g\t%g\n",Wrc,Pval_rc);
    }
    free_seq_RSAT(seq);
    if (rc)
        free_seq_RSAT(seqrc);
    
    return outpvalue;
}

double *scan_seq_homer(char *faline, double **homer_mat,int L, Pvalues_t *pvalues,int rc)
{
	int faline_len=strlen(faline);
	//printf("%d\n",faline_len);
    int l = L;
    ASSERT(l < 256, "invalid matrix size");
    //make seq
    int i,j;	
	Seq_t *seq = new_seq(faline_len); 
	char ccc;
	for(i=0;i<faline_len;i++)
    {
        ccc = faline[i];
        seq_append_c(seq, ccc);
    } 
    //make seqrc
    Seq_t *seqrc = NULL;

    int maxpos = faline_len - l;  //maxpos=seq length - motif length
    
    double *outpvalue;
    if (rc)   //2str
    {
    	seqrc = new_seq_rc(seq);    
    	outpvalue=(double *)malloc((2*maxpos+2)*sizeof(double));
    }
    else
    {
    	outpvalue=(double *)malloc((maxpos+1)*sizeof(double));
    }
    
    for (i = 0; i <= maxpos; i++)
    {
        if (!word_is_valid(seq, i, l))
            continue;
            	
		double h_score=0.0;
		char *h_word=&seq->data[i];
		for (j=0;j<l;j++) {
			h_score += homer_mat[(int) h_word[j]][j];
		}
		
        double Pval = score2pvalue_RSAT(pvalues, h_score);
		//printf("%g\t%g\n",h_score,Pval);
		outpvalue[i]=Pval;
		
        if (!rc)
            continue;
		double hr_score=0.0;
		char *hr_word=&seqrc->data[maxpos - i];
		for (j=0;j<l;j++) {
			hr_score += homer_mat[(int) hr_word[j]][j];
		}
        double Pval_rc = score2pvalue_RSAT(pvalues, hr_score);
        //printf("%g\t%g\n",hr_score,Pval_rc);
        outpvalue[maxpos+1+i]=Pval_rc;
    }
    free_seq_RSAT(seq);
    if (rc)
        free_seq_RSAT(seqrc);
    
    return outpvalue;
}
	   
double score2pvalue_RSAT(Pvalues_t *pvalues, double score)
{
    if (pvalues == NULL)
        return 0.0;
    double index  = ((score - (double) pvalues->w_min) / (double) (pvalues->w_max - pvalues->w_min)) * ((double) pvalues->size);
    index = MIN(MAX(round(index) - 1, 0), pvalues->size - 1);
    return pvalues->data[(int) index];
}

double fisher_combine_pvalue(double *inpvalue,double chisq_k)
{
	double chisq_value=0.0;
	int i;
	int k;
	k=(int)chisq_k;
	for(i=0;i<k;i++)
	{
		if(inpvalue[i]>0)
		{
			chisq_value=chisq_value-2*log(inpvalue[i]);
		}
		else  //RSAT Homer pvale may be 0;
		{
			chisq_k=chisq_k-1;
		}
		
	}
	if(chisq_k<=0){printf("pvalue error!\n");exit(1);}
	
 	double gslra=gsl_cdf_chisq_Q(chisq_value, chisq_k*2);
 	//printf("%g\t%g\t%g\n",chisq_value,chisq_k,gslra);	
 	return gslra;
	
}
double Stouffer_combine_pvalue(double *inpvalue,double chisq_k,double *score)
{
   //double caoc=gsl_cdf_gaussian_P(0.1,1);   //pnorm in R
  //double kaoc=gsl_cdf_gaussian_Pinv(0.2,1);   //qnorm in R
	double sum=0.0;
	int i;
	int k;
	k=(int)chisq_k;
	for(i=0;i<k;i++)
	{
		if(inpvalue[i]>0 && inpvalue[i]<=1)
		{
			sum=sum+gsl_cdf_gaussian_Pinv((1-inpvalue[i]),1);
			//printf("%g\n",sum);
		}
		else  //RSAT Homer pvale may be 0;
		{
			chisq_k=chisq_k-1;
		}
		
	}
	if(chisq_k<=0){printf("pvalue error!\n");exit(1);}
	double zscore=sum/sqrt(chisq_k);
 	double gslra=1-gsl_cdf_gaussian_P(zscore,1);
	if(zscore<-10){zscore=-10;}
 	*score=zscore;
 	//printf("%g\t%g\t%g\n",chisq_value,chisq_k,gslra);	
 	return gslra;
	
}

double **weight_matrix(double **align_mat, int A_size ,int L)
{
  int i, j;
  double col_sum;              /* The sum of the current column. */
  double **matrix;             /* The weight matrix. */


  /* Allocate space for "matrix[][]". */
  matrix = (double **)malloc(A_size * sizeof(double *));
  for (i = 0; i < A_size; ++i)
	matrix[i] = (double *)malloc(L * sizeof(double));

  /* Manipulate each position of each row. */
  for (j = 0; j < L; ++j)
    {
      for (i = 0, col_sum = 0.0; i < A_size; ++i) col_sum += align_mat[i][j];

      /* Determine the weights. */
      for (i = 0; i < A_size; ++i)
	{
	  if (align_mat[i][j] == 0.0) matrix[i][j] = -INFINITY_LARGE;
	  else matrix[i][j] = log(align_mat[i][j] / col_sum / con_Fudge);
	  //printf("%3.3f\t%3.3f\t%3.3f\t%3.3f\n",align_mat[i][j],col_sum,con_Fudge,matrix[i][j]);
	}
   }
  return(matrix);
}

double *det_marginal_prob(double **weight_mat, int T_SIZE,int width)
{
  int Min_int_score = 0;
  double *Marginal_prob;
  /* Reset Marginal_prob[0] to the beginning of the allocated memory. */
  Marginal_prob += Min_int_score;  //all zeros
	//printf("%.3f\t%.3f\n",Marginal_prob,Min_int_score);
  /* Allocate space for STmax_score[] and STmin_score[]; determine Minus_inf,
   * Max_exact_score, Min_exact_score, Min_pseudo_score, Multiple,
   * STmax_score[], STmin_score[], Max_int_score, Min_int_score,
   * and STmax_int_range. */
  int i, j,k;
  char max_min_equality;  /* YES: maximum and minimum (non -INFINITY_LARGE)  *      scores are equal. */
  char Minus_inf;
  int min_int_score;      /* Minimum integral score greater than -INFINITY_LARGE. */
  double min_prob_score;  /* Minimum score for estimating probabilities. */
  double adj_minus_infinity(); /* Assign values to the -INFINITY_LARGE elements.
				* Determine Min_pseudo_score. */
  double Min_score= -INFINITY_LARGE;
  double Max_int_range = 1.0e4;
  
  int *STmax_score = (int *)malloc(width*sizeof(int));
  int *STmin_score = (int *)malloc(width*sizeof(int));
  double *max_score = (double *)malloc(width*sizeof(double));/* Value of column element with largest score. */
  double *min_score = (double *)malloc(width*sizeof(double)); /* Value of column element with lowest score. */
 

  /* Determine "Max_exact_score" and "Min_exact_score". */
  double Max_exact_score,Min_exact_score;
  int Max_int_score;

  for (j = 0, Max_exact_score = Min_exact_score = 0.0; j < width; ++j)
    {
      for (i = 0, max_score[j] = -INFINITY_LARGE, min_score[j] = INFINITY_LARGE;
	   i < T_SIZE; ++i)
	{
	  if (weight_mat[i][j] > -INFINITY_LARGE)
	    {
	      if (weight_mat[i][j] < min_score[j])
		min_score[j] = weight_mat[i][j];
	      if (weight_mat[i][j] > max_score[j])
		max_score[j] = weight_mat[i][j];
	    }
	  else Minus_inf = YES;
	}
      if (min_score[j] < INFINITY_LARGE)
	{
	  Max_exact_score += max_score[j];
	  Min_exact_score += min_score[j];
	}
    }
  //printf("Multiple\t%g\t%g\t%g\t%g\n",Multiple,Max_int_range,Max_exact_score,min_prob_score);

  /* Determine the minimum score for calculating p-values. */
  if (Min_score >= Min_exact_score) min_prob_score = Min_score;
  else min_prob_score = Min_exact_score;
  /* Make sure the minimum score for calculating p-values is
   * not greater than the maximum possible score. */
  if (min_prob_score > Max_exact_score)
    {
      fprintf(stderr, "The minimum score for calculating p-values is greater");
      fprintf(stderr, "\nthan the maximum possible score.\n");
      exit(1);
    }
  /* Adjustment for the special case when the minimum score for
   * calculating p-values is equal to the maximum possible score. */
  else if (min_prob_score == Max_exact_score) max_min_equality = YES;
  else max_min_equality = NO;
  
  /* Determine "Multiple". */
  if (max_min_equality == YES) Multiple = Max_int_range / Max_exact_score;
  else Multiple = Max_int_range / (Max_exact_score - min_prob_score);  //10000/Max_exact_score
  /* Determine "STmax_score[]", "STmin_score[]", and Max_int_score. */
  for (j = 0, min_int_score = Max_int_score = 0; j < width; ++j)
    {
      if (min_score[j] < INFINITY_LARGE)
	{
	  STmax_score[j] = ROUND_TO_INT(Multiple * max_score[j]);
	  STmin_score[j] = ROUND_TO_INT(Multiple * min_score[j]);

	  Max_int_score += STmax_score[j];
	  min_int_score += STmin_score[j];
	}
      else STmax_score[j] = STmin_score[j] = 0;
    }
	
  /* Determine "Min_int_score". */
  Min_int_score = ROUND_TO_INT(Multiple * min_prob_score);
  if ((Min_exact_score >= Min_score) && (min_int_score < Min_int_score))
    Min_int_score = min_int_score;
  --Min_int_score;
  if ((max_min_equality == YES) || (Min_int_score >= Max_int_score))
    Min_int_score = Max_int_score - 1;
 
  /* Determine "STmax_int_range". */
  int STmax_int_range = Max_int_score - Min_int_score;
  if (STmax_int_range < 1){printf("error!\n");exit(1);}
  free(max_score);
  free(min_score);


  /* Allocate space for "STarrayA[]", "STarrayB[]", "STarrayC[]",
   * and "STarray_index[]". */
  int max_score_1 = STmax_int_range + 1;  /* Maximum size of arrays. */
  double *STarrayA = (double *)malloc(max_score_1*sizeof(double));
  double *STarrayB = (double *)malloc(max_score_1*sizeof(double));
  double *STarrayC = (double *)malloc(max_score_1*sizeof(double));;
  int *STarray_index = (int *)malloc(T_SIZE*sizeof(int));
  for (i = 0; i < max_score_1; ++i)
  {
      STarrayA[i] = -INFINITY_LARGE;
      STarrayC[i] = -INFINITY_LARGE;
      STarrayB[i] = 0.0;
   }



  /* Determine the logarithm of the distribution of the integral scores of
   * the "weight_mat[][]" and place in ArrayA[]. */
  double ln_prob_AB;       /* ln(column_probability in STarrayB[] multiplied
			    *    by previous_probability in STarrayA[]) */
  double ln_prob_C;        /* The ln(probability currently in STarrayC[]). */
  double *temp_array;      /* Temporary pointer when exchanging pointers to
			    * "STarrayA" and "STarrayC". */
  int max_score_old;           /* Maximum score prior to current multiplication. */
  int new_max_score;       /* Maximum score after current multiplication. */
  int min_score_old;           /* Minimum score prior to current multiplication. */
  int new_min_score;       /* Minimum score after current multiplication. */
  int diff_sub;      /* Number of different substitution scores (<= A_size). */

  /* Pointers to multiplication arrays so maximum index = maximum score. */
  double *arrayA;       /* Pointer to STarrayA: product of previous columns. */
  double *arrayB;       /* Pointer to STarrayB: current column. */
  double *arrayC;       /* Pointer to STarrayC: product of current column
			 *                      with previous columns. */
  /* Minimum index in each of the multiplication arrays. */
  int min_idx_A;        /* arrayA */
  int min_idx_B;        /* arrayB */
  int min_idx_C;        /* arrayC */

  /* The score index for each multiplication array. */
  int idx_A;            /* arrayA */
  int idx_B;            /* arrayB */
  int idx_C;            /* arrayC */

  /* Initialize "arrayA", "max_score_old", and "min_score_old" for
   * the first matrix column.*/
  max_score_old = STmax_score[0];
  min_score_old = STmin_score[0];
  min_idx_A = -(STmax_int_range - max_score_old);
  if (min_score_old > min_idx_A) min_idx_A = min_score_old - 1;
  arrayA = -min_idx_A + STarrayA;
  for (i = 0; i < T_SIZE; ++i)
    {
      if (weight_mat[i][0] <= -INFINITY_LARGE) idx_A = min_idx_A;
      else
	{
	  idx_A = ROUND_TO_INT(Multiple * weight_mat[i][0]);
	  if (idx_A < min_idx_A) idx_A = min_idx_A;
	}

      if (arrayA[idx_A] > -INFINITY_LARGE)
	arrayA[idx_A] = SUM_LN(log(con_Fudge), arrayA[idx_A]);
      else arrayA[idx_A] = log(con_Fudge);
    }

  for (j = 1; j < width; ++j)
    {
      /* Insert the probabilities for the current column into "arrayB[]".
       * Set the "STarray_index[]" to indicate once which scores occur
       * in the current column. */
      min_idx_B = STmax_score[j] - STmax_int_range;
      arrayB = -min_idx_B + STarrayB;
      for (i = 0, diff_sub = 0; i < T_SIZE; ++i)
	{
	  if ((weight_mat[i][j] <= -INFINITY_LARGE) ||
	     ((idx_B = ROUND_TO_INT(Multiple * weight_mat[i][j])) < min_idx_B))
	    idx_B = min_idx_B;

	  if (arrayB[idx_B] == 0.0)
	    {
	      arrayB[idx_B] = con_Fudge;
	      STarray_index[diff_sub++] = idx_B;
	    }
	  else arrayB[idx_B] += con_Fudge;
	}
      /* Convert the probabilities in "arrayB[]" to logarithms. */
      for (i = 0; i < diff_sub; ++i)
	arrayB[STarray_index[i]] = log(arrayB[STarray_index[i]]);


      /* Multiply "arrayA[]" by "arrayB[]". */
      new_max_score = max_score_old + STmax_score[j];
      new_min_score = min_score_old + STmin_score[j];
      min_idx_C = -(STmax_int_range - new_max_score);
      if (new_min_score > min_idx_C) min_idx_C = new_min_score - 1;
      arrayC = -min_idx_C + STarrayC;
      if (arrayA[min_idx_A] > -INFINITY_LARGE)
	{
	  arrayC[min_idx_C] = arrayA[min_idx_A] + LN_PROB_SUM;
	  arrayA[min_idx_A] = -INFINITY_LARGE;
	}
      for (idx_A = max_score_old; idx_A > min_idx_A; --idx_A)
	{
	  if (arrayA[idx_A] > -INFINITY_LARGE)
	    {
	      for (i = 0; i < diff_sub; ++i)
		{
		  idx_B = STarray_index[i];

		  if ((idx_C = idx_A + idx_B) < min_idx_C) idx_C = min_idx_C;

		  ln_prob_AB = arrayA[idx_A] + arrayB[idx_B];
		  ln_prob_C = arrayC[idx_C];

		  if (ln_prob_C <= -INFINITY_LARGE) arrayC[idx_C] = ln_prob_AB;
		  else arrayC[idx_C] = SUM_LN(ln_prob_AB, ln_prob_C);
		}

	      /* Zero the current element of "arrayA[]". */
	      arrayA[idx_A] = -INFINITY_LARGE;
	    }
	}


      /* Zero "arrayB[]". */
      for (i = 0; i < diff_sub; ++i) arrayB[STarray_index[i]] = 0.0;

      /* Exchange pointers to "STarrayA[]" and "STarrayC[]". */
      temp_array = STarrayA;
      STarrayA = STarrayC;
      STarrayC = temp_array;
      arrayA = arrayC;
      min_idx_A = min_idx_C;

      /* Update "max_score_old". */
      max_score_old = new_max_score;
      min_score_old = new_min_score;
    }

  /* Transfer information on arrayA into external variables.
   * Max_int_score: the maximum index determined in determine_variables().
   * Min_int_score: the minimum index might change from value set
   *                in determine_variables() due to rounding error.
   *      ArrayA[]: holds final ln(probabilities) of integer scores. */
  Min_int_score = min_idx_A;    
  //ArrayA = arrayA;
  
  
  

  /* Determine the logarithm of the p-values and place in Marginal_prob[]. */
  //determine_marginal_prob();

  int max_score_last;    /* The greater of the maximum integral score determined
		     * from the maximum exact score and the integral matrix. */
  double delta_ln_p;/* The change in the ln(p-value) between
		     * two observed integral scores. */


  /* Allocate space for "Marginal_prob[]". */
  max_score_last = ROUND_TO_INT(Multiple * Max_exact_score);
  if (Max_int_score > max_score_last) max_score_last = Max_int_score;
  Marginal_prob = (double *)malloc((max_score_last - Min_int_score + 1)*sizeof(double));
  Marginal_prob -= Min_int_score;


  /* Determine the ln(p-value) of the Max_int_score. */
  Marginal_prob[Max_int_score] = arrayA[Max_int_score];

//printf("Max_int_score%d\t%d\n",Min_int_score,max_score_last);
  /* Copy ln(p-values) for scores greater than Max_int_score. */
  for (i = Max_int_score + 1; i <= max_score_last; ++i) //9998~10000
    Marginal_prob[i] = Marginal_prob[Max_int_score];

	
  /* Determine the main part of Marginal_prob[]. */
  for (i = Max_int_score - 1, j = Max_int_score;
       i > Min_int_score; --i)           //Min_int_score=-1
    {
      if (arrayA[i] > -INFINITY_LARGE)
	{
	  Marginal_prob[i] = SUM_LN(arrayA[i], Marginal_prob[j]);

	  /* Interpolate ln(p-values) for unobserved scores. */
	  delta_ln_p = (Marginal_prob[i] - Marginal_prob[j]) / (double)(j - i);
	  for (k = j--; j > i; k = j--)
	    Marginal_prob[j] = Marginal_prob[k] + delta_ln_p;
	}
    }

  /* Copy ln(p-values) for scores less than the lowest observed score j. */
  for (k = j - 1; k > Min_int_score; --k) Marginal_prob[k] = Marginal_prob[j];

  /* Determine the ln(p-value) for the Min_int_score. */
  Marginal_prob[i] = SUM_LN(arrayA[i], Marginal_prob[j]);

 // int i;
//for (i = 0; i < 10000; ++i){printf("%.3f\n",Marginal_prob[i]);}
//for (i = 0; i < 10000; ++i){printf("%.3f\n",arrayA[i]);}
//for (i = 0; i < 10000; ++i){printf("%.3f\n",STarrayA[i]);}
  /* Free arrays except for STarrayA. */
  free(STmax_score);
  free(STmin_score);
  free(STarrayA);
  free(STarrayB);
  free(STarrayC);
  free(STarray_index);	
  return Marginal_prob;
}

double *scan_seq_Consensus(char *faline,double **Align_mat,double *Marginal_prob,int rc,int L,int *alphabase)
{
  int i,j,idx,postt;  
  char *l_mer;           /* The L-mer currently being scored. */
  int letter;              /* The current letter of the current L-mer. */
  double score;            /* Score of current L-mer or its complement. */
  int int_score;
  double ln_prob;
  int faline_len=strlen(faline);          
  int maxpos = faline_len - L; 

  double *outpvalue;
  if (rc ==1 )   //2str
  {outpvalue=(double *)malloc((2*maxpos+2)*sizeof(double));}
  else{outpvalue=(double *)malloc((maxpos+1)*sizeof(double));}
    
    
  /* Process each L-mer of the Sequence[]. */
  int badvalue;
  for (idx= 0 ;idx <= maxpos; idx++)
  {
  	badvalue=0;  
	for (i = idx; i < idx + L; i++)
	{
		postt=(int)faline[i];
		if(alphabase[postt] < 0) badvalue=1;
	}	
    if(badvalue)continue;
	  /* Determine the score of the current L-mer. */
	for (i = 0, l_mer = faline + idx, score = 0; i < L; ++i)
	{
		postt=(int)(l_mer[i]);
	    letter = alphabase[postt];
	    //printf("%d\t%d\t%c\t%d\n",postt,alphabase[postt],l_mer[i],L);
	    score += Align_mat[letter][i];
	}
	int_score = ROUND_TO_INT(Multiple * score);  /* Integral score. */
	ln_prob = Marginal_prob[int_score];
	 //printf("pos\t%d\t%g\t%g\n",int_score,score,ln_prob);
	//printf("%g\n",ln_prob);
	//printf("pos\t%d\t%g\t%.7f\n",int_score,score,ln_prob);
	outpvalue[idx]=pow(BASE,ln_prob);
		
	if (rc == 1)   //Yes both stand
	{    
		for (i = L-1,score = 0; i >=0; i--)
		{
			postt=(int)(l_mer[i]);
		    letter = 3- alphabase[postt];
		    //printf("%d\t%d\t%c\t%g\n",postt,letter,l_mer[i],score);
		    score += Align_mat[letter][L-1-i];
		}
/*
		for (i = 0, j = L_1, score = 0; i < L; ++i, --j)
		{
		  letter = A_comp[(int)(l_mer[j])];
		  score += Matrix[letter][i];
		}*/
	 	int_score = ROUND_TO_INT(Multiple * score); 
		ln_prob = Marginal_prob[int_score];
		outpvalue[maxpos+1+idx]=pow(BASE,ln_prob);
		//printf("rev\t%d\t%g\t%g\n",int_score,score,ln_prob);
		//printf("%g\n",ln_prob);
	//printf("rev\t%d\t%g\t%.7f\n",int_score,score,ln_prob);
		
	 }
	 
  }
  return outpvalue; 
}

double *cal_storm_index_pvalue(double **storm_mat,int Align_matI,int Align_matJ,double *base_comp)
{
	int i,j,k;

	int **smi;
    
  	double minvalue=0.0;
    for ( i= 0; i < Align_matI; i++)
        for (j = 0; j < Align_matJ; j++)
    		if(minvalue>storm_mat[i][j]){minvalue=storm_mat[i][j];}
  	//printf("min %g\n",minvalue);
    translation_storm= -minvalue *  scale_storm;
    	
  	smi = (int **)malloc(Align_matI * sizeof(int *));
    for ( i= 0; i < Align_matI; i++)
    {
    	smi[i] = (int *)malloc(Align_matJ * sizeof(int));
        for (j = 0; j < Align_matJ; j++)
        {
        	smi[i][j]=(int)((storm_mat[i][j]-minvalue)*scale_storm + 0.5);
        	//printf("%d\t",smi[i][j]);
        }
        //printf("\n");       
    } 
    
/*	double total_temp=0.0;
    for (j = 0; j < Align_matJ; j++)
    {	
    	double maxline=0.0;
        for ( i= 0; i < Align_matI; i++)
    		if(maxline<storm_mat[i][j]){maxline=storm_mat[i][j];} 
    	total_temp= total_temp+(maxline-minvalue);
    	//printf("%g\t%g\t%g\n",total_temp,maxline,minvalue);
    }    
    total_discretized_scores=(int)(total_temp*scale_storm)+1;
    printf("total_discretized_scores %d\n",total_discretized_scores);
*/
    for (j = 0; j < Align_matJ; j++)
    {	
    	int maxline=0;
        for ( i= 0; i < Align_matI; i++)
    	if(maxline<smi[i][j]){maxline=smi[i][j];} 
    	total_discretized_scores+=maxline;
    	//printf("%g\t%g\t%g\n",total_temp,maxline,minvalue);
    } 
    total_discretized_scores++;   
    //printf("total_discretized_scores %d\n",total_discretized_scores);
    
        
    int *used;
    int *prev_used;
    double *counts;
    double *prev_counts;
    double *outpvalue;
    used= (int *)malloc(total_discretized_scores * sizeof(int));
    prev_used= (int *)malloc(total_discretized_scores * sizeof(int));
    counts= (double *)malloc(total_discretized_scores * sizeof(double));
    prev_counts= (double *)malloc(total_discretized_scores * sizeof(double));
    outpvalue= (double *)malloc(total_discretized_scores * sizeof(double));
    
 
	for(i=0;i<total_discretized_scores;i++)
		used[i]=0;
	for(i=0;i<total_discretized_scores;i++)
		counts[i]=0.0;	

	for (j = 0; j < ALPHABET_SIZE; ++j) {
		counts[smi[j][0]] += base_comp[j];
		used[smi[j][0]] = 1;
  	}
  	
  	for(i=0;i<total_discretized_scores;i++)
  		prev_used[i]=used[i];
  	for(i=0;i<total_discretized_scores;i++)
  		prev_counts[i]=counts[i];  

	double *tempcounts;
	int *tempused;
	int offset;
	for (i = 1; i < Align_matJ; ++i) { //ooooo! need to read again and again   PVALUE is here!
    //vector<size_t>::const_iterator smi_iter = smi[i].begin();
		for (j = 0; j < Align_matI; ++j) {
			for (k = 0; k < total_discretized_scores; ++k)
			if (prev_used[k] == i) {
		  		offset = k + smi[j][i];
		  	if (used[offset] < i + 1)
		    	counts[offset] = 0;
		  	counts[offset] += prev_counts[k]*base_comp[j];
		  	used[offset] = i + 1;
		  	//printf("%d\t%d\t%d\t%g\t%g\t%d\n",i,j,k,counts[offset],prev_counts[k],used[offset]);
		  	
			}
    	}
    	tempused=used;//used.swap(prev_used);
    	used=prev_used;
    	prev_used=tempused;
   		tempcounts=counts;//counts.swap(prev_counts);
   		counts=prev_counts;
    	prev_counts=tempcounts;
  	}
  	//printf("translation_storm %g\n",translation_storm);

	double m_threth = 0.0;
	double m_pvalue = 0.0;
  	for (i = 0; i < total_discretized_scores; ++i)
  	{
		if (prev_used[total_discretized_scores - i - 1] == Align_matJ) {
      		m_pvalue += prev_counts[total_discretized_scores - i - 1];
  	 	} 
  	 	outpvalue[i]= m_pvalue;
	 	m_threth=(total_discretized_scores - i - translation_storm*Align_matJ)/scale_storm;
	 	//printf("%g\t%g\n",m_threth,outpvalue[i]);

	}  
   
	for(i=0;i<Align_matI;i++)
		free(smi[i]);
	free(smi); 	
	free(used);
	free(prev_used);
	free(counts);
	free(prev_counts);
	
	return  outpvalue;
}

double *scan_seq_storm(char *faline, double **homer_mat,int L,double *index_pvalue,double *base_comp,int *alphabase,int rc)
{
  int i,j,idx,postt;  
  char *l_mer;           /* The L-mer currently being scored. */
  int letter;              /* The current letter of the current L-mer. */
  double score;            /* Score of current L-mer or its complement. */
  int int_score;
  double ln_prob;
  int faline_len=strlen(faline);          
  int maxpos = faline_len - L; 

  double *outpvalue;
  if (rc ==1 )   //2str
  {outpvalue=(double *)malloc((2*maxpos+2)*sizeof(double));}
  else{outpvalue=(double *)malloc((maxpos+1)*sizeof(double));}
    
    
  /* Process each L-mer of the Sequence[]. */
  int badvalue;
  for (idx= 0 ;idx <= maxpos; idx++)
  {
  	badvalue=0;  
	for (i = idx; i < idx + L; i++)
	{
		postt=(int)faline[i];
		if(alphabase[postt] < 0) badvalue=1;
	}	
    if(badvalue)continue;
	  /* Determine the score of the current L-mer. */
	for (i = 0, l_mer = faline + idx, score = 0; i < L; ++i)
	{
		postt=(int)(l_mer[i]);
	    letter = alphabase[postt];
	    score += homer_mat[letter][i];
	}
	int_score = (int)(total_discretized_scores - translation_storm*L - score*scale_storm-0.5);  
	if(int_score>total_discretized_scores-1){int_score=total_discretized_scores-1;}
	outpvalue[idx] = index_pvalue[int_score];
	
	printf("storm\t%d\t%g\t%g\n", int_score, score, outpvalue[idx]);//printf("pos\t%d\t%g\t%g\n",int_score,score,outpvalue[idx]);
		
	if (rc == 1)   //Yes both stand
	{    
		for (i = L-1,score = 0; i >=0; i--)
		{
			postt=(int)(l_mer[i]);
		    letter = 3- alphabase[postt];
		    score += homer_mat[letter][L-1-i];
		}
	 	int_score = (int)(total_discretized_scores - translation_storm*L - score*scale_storm-0.5);
		if(int_score>total_discretized_scores-1){int_score=total_discretized_scores-1;}
		outpvalue[maxpos+1+idx]=index_pvalue[int_score];
		printf("storm\t%d\t%g\t%g\n", int_score, score, outpvalue[maxpos + 1 + idx]);//printf("rev\t%d\t%g\t%g\n",int_score,score,outpvalue[maxpos+1+idx]);	
	 }
	 
  }
  return outpvalue; 
}