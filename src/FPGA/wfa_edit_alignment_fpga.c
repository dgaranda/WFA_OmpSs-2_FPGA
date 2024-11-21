/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2024 by Diego García Aranda <diego.garcia1@bsc.es>
 *
 *  This file is part of Wavefront Alignments Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Diego García Aranda <diego.garcia1@bsc.es>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>

double wall_time () {
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC,&ts);
   return (double) (ts.tv_sec) + (double) ts.tv_nsec * 1.0e-9;
}

/*
 * Translate k and offset to coordinates h,v
 */
#define EWAVEFRONT_V(k,offset) ((offset)-(k))
#define EWAVEFRONT_H(k,offset) (offset)

#define EWAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define EWAVEFRONT_OFFSET(h,v)   (h)

#define MAX(a,b) (((a)>=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))

#define PRINTF(format, ...) do { printf(format, ##__VA_ARGS__); } while(0)
#define PRINTF_COND(condition,format, ...) do { if (condition) { printf(format, ##__VA_ARGS__); } } while(0)
#define PRINTF_ERROR(format, ...) do { fprintf(stderr, format, ##__VA_ARGS__); } while(0);

#define LO_IDX(distance) (-distance)
#define HI_IDX(distance) (distance)
#define OFFSET_IDX(distance, k) (distance*(distance+1) + k)

#ifndef FPGA_EMU
#define FPGA(p) _Pragma(p)
#else
#define FPGA(...)
#endif            

typedef int16_t ewf_offset_t;  // Edit Wavefront Offset

/*
 * Wavefront for FPGA device
 */
typedef struct {
  ewf_offset_t* offsets;

  // CIGAR
  char* edit_cigar;
  int edit_cigar_length;

} edit_wavefronts_fpga_t;

int edit_wavefronts_init(
    edit_wavefronts_fpga_t* const wavefronts,
    const int pattern_length,
    const int text_length,
    const bool aligned,
    const size_t page_size) {

  const int max_distance = pattern_length + text_length;
  // Allocate wavefronts diagonals index
  if(aligned){

    // Allocate wavefronts offsets aligned to page size
    wavefronts->offsets= (ewf_offset_t*) aligned_alloc(page_size,max_distance*max_distance*sizeof(ewf_offset_t));
    if(wavefronts->offsets == NULL){
      PRINTF_ERROR("Aligned allocation of wavefronts offsets failed");
      return EXIT_FAILURE;
    }

    // Allocate CIGAR aligned to page size
    wavefronts->edit_cigar = (char*) aligned_alloc(page_size,max_distance);
    if(wavefronts->edit_cigar == NULL){
      PRINTF_ERROR("Aligned allocation of CIGAR failed");
      return EXIT_FAILURE;
    }
  
  }
  else{

    // Allocate wavefronts offsets
    wavefronts->offsets= calloc(max_distance*max_distance,sizeof(ewf_offset_t));
    if(wavefronts->offsets == NULL){
      PRINTF_ERROR("Allocation of wavefronts offsets failed");
      return EXIT_FAILURE;
    }

    // Allocate CIGAR
    wavefronts->edit_cigar = malloc(max_distance);
    if(wavefronts->offsets == NULL){
      PRINTF_ERROR("Allocation of CIGAR failed");
      return EXIT_FAILURE;
    }

  }

  return EXIT_SUCCESS;

}

void edit_wavefronts_clean(
    edit_wavefronts_fpga_t* const wavefronts) {
  free(wavefronts->offsets);
  free(wavefronts->edit_cigar);
}


/*
 * Edit Wavefront Backtrace
 */
int edit_wavefronts_backtrace(
    const ewf_offset_t* const offsets_wavefronts,
    char* const edit_cigar,
    const int target_k,
    const int target_distance) {
  // Parameters
  int edit_cigar_idx = 0;
  int k = target_k, distance = target_distance;
  ewf_offset_t offset = offsets_wavefronts[OFFSET_IDX(distance,k)];
  while (distance > 0) {
    // Fetch
    const ewf_offset_t* const offsets = offsets_wavefronts + OFFSET_IDX((distance-1),0);
    const int lo = LO_IDX(distance);
    const int hi = HI_IDX(distance);
    // Traceback operation
    if (lo <= k+1 && k+1 <= hi && offset == offsets[k+1]) {
      edit_cigar[edit_cigar_idx++] = 'D';
      ++k;
      --distance;
    } else if (lo <= k-1 && k-1 <= hi && offset == offsets[k-1] + 1) {
      edit_cigar[edit_cigar_idx++] = 'I';
      --k;
      --offset;
      --distance;
    } else if (lo <= k && k <= hi && offset == offsets[k] + 1) {
      edit_cigar[edit_cigar_idx++] = 'X';
      --distance;
      --offset;
    } else {
      edit_cigar[edit_cigar_idx++] = 'M';
      --offset;
    }
  }
  // Account for last offset of matches
  while (offset > 0) {
    edit_cigar[edit_cigar_idx++] = 'M';
    --offset;
  }
  // Return CIGAR length
  return edit_cigar_idx;
}


/*
 * Extend Wavefront
 */
void edit_wavefronts_extend_wavefront(
    ewf_offset_t* const offsets_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int distance) {
  // Parameters
  ewf_offset_t* const offsets = offsets_wavefronts + OFFSET_IDX(distance,0);
  const int k_min = LO_IDX(distance);
  const int k_max = HI_IDX(distance);
  // Extend diagonally each wavefront point
  int k;
  for (k=k_min;k<=k_max;++k) {
    int v = EWAVEFRONT_V(k,offsets[k]);
    int h = EWAVEFRONT_H(k,offsets[k]);
    while (v<pattern_length && h<text_length && pattern[v++]==text[h++]) {
      ++(offsets[k]);
    }
  }
}

/*
 * Edit Wavefront Compute
 */
void edit_wavefronts_compute_wavefront(
    ewf_offset_t* const offsets_wavefronts,
    const int distance) {
  const int distance_minus_one = distance-1;
  // Fetch lo and hi
  const int lo = LO_IDX(distance_minus_one);
  const int hi = HI_IDX(distance_minus_one);

  // Fetch offsets
  ewf_offset_t* const offsets = offsets_wavefronts + OFFSET_IDX(distance_minus_one,0);
  ewf_offset_t* const next_offsets = offsets_wavefronts + OFFSET_IDX(distance,0);

  // Loop peeling (k=lo-1)
  next_offsets[lo-1] = offsets[lo];

  // Loop peeling (k=lo)
  const ewf_offset_t bottom_upper_del = ((lo+1) <= hi) ? offsets[lo+1] : -1;
  next_offsets[lo] = MAX(offsets[lo]+1,bottom_upper_del);

  // Compute next wavefront starting point
  int k;
  for (k=lo+1;k<=hi-1;++k) {
    const ewf_offset_t max_ins_sub = MAX(offsets[k],offsets[k-1]) + 1;
    next_offsets[k] = MAX(max_ins_sub,offsets[k+1]);
  }

  // Loop peeling (k=hi)
  const ewf_offset_t top_lower_ins = (lo <= (hi-1)) ? offsets[hi-1] : -1;
  next_offsets[hi] = MAX(offsets[hi],top_lower_ins) + 1;

  // Loop peeling (k=hi+1)
  next_offsets[hi+1] = offsets[hi] + 1;
}


/*
 * Edit distance alignment using wavefronts
 */
FPGA("oss task device(fpga) in([pattern_length]pattern, [text_length]text) out([1]score, [1]edit_cigar_length, [max_distance]edit_cigar) inout([max_distance*max_distance]offsets_wavefronts)")
void edit_wavefronts_align(
    ewf_offset_t* offsets_wavefronts,
    char* edit_cigar,
    int* edit_cigar_length,
    const char* pattern,
    const int pattern_length,
    const char* text,
    const int text_length,
    const int max_distance,
    int* score) {
FPGA("HLS inline")
  // Parameters
  const int target_k = EWAVEFRONT_DIAGONAL(text_length,pattern_length);
  const int target_k_abs = ABS(target_k);
  const ewf_offset_t target_offset = EWAVEFRONT_OFFSET(text_length,pattern_length);

  // Init wavefronts
  int distance;
  offsets_wavefronts[0] = 0;

  // Compute wavefronts for increasing distance
  for (distance=0;distance<max_distance;++distance) {

    // Extend diagonally each wavefront point
    edit_wavefronts_extend_wavefront(offsets_wavefronts,
        pattern,pattern_length,
        text,text_length,distance);
    // Exit condition
    if (target_k_abs <= distance &&
        offsets_wavefronts[OFFSET_IDX(distance,target_k)] == target_offset) break;
    
    // Compute next wavefront starting point
    edit_wavefronts_compute_wavefront(
        offsets_wavefronts,distance+1);
  }
  (*score) = distance;

  // Backtrace
  (*edit_cigar_length) = edit_wavefronts_backtrace(offsets_wavefronts,edit_cigar,target_k,distance);

}

bool edit_wavefronts_check(
  const char* const cigar,
  const int cigar_length,
  const int score, 
  const char* const filename) {
  
  FILE* ref_file = fopen(filename, "r");
  if (ref_file == NULL) {
    PRINTF_ERROR("Error while opening check file %s\n", filename);
    return false;
  }

  // Read reference score
  int score_ref;
  if (fscanf(ref_file, "%d", &score_ref) != 1) {
    PRINTF_ERROR("Error while reading reference score in check file %s\n", filename);
    return false;
  }  

  // Check score
  if (score != score_ref) {
    PRINTF_ERROR("Check has failed: reference score != result score\n");
    PRINTF_ERROR("Reference score: %d\n", score_ref);
    PRINTF_ERROR("Result score: %d\n", score);
    return false;
  }

  fgetc(ref_file); // Skip newline

  char ref_ch;
  int i = 0;
  while ((ref_ch = fgetc(ref_file)) != EOF && ref_ch != '\n' && i < cigar_length) {
    if(ref_ch != cigar[i]){
      PRINTF_ERROR("Check has failed: reference CIGAR != result CIGAR at position %d\n", i);
      PRINTF_ERROR("Reference CIGAR: %c\n", ref_ch);
      PRINTF_ERROR("Result CIGAR: %c\n", cigar[i]);
      return false;
    }
  }

  // Check CIGAR length
  if (ref_ch != EOF || i != cigar_length) {
    PRINTF_ERROR("Check has failed: reference CIGAR length != result CIGAR length\n");
    return false;
  }

  // Return
  return true;

}

int edit_wavefronts_write_result(
  const char* const cigar,
  const int cigar_length,
  const int score,
  const char* const filename) {
  
  
  FILE* result_file = fopen(filename, "w");
  if (result_file == NULL) {
    PRINTF_ERROR("Error while opening result file %s\n", filename);
    return EXIT_FAILURE;
  }

  fprintf(result_file, "%d\n", score);
  fwrite(cigar, sizeof(char), cigar_length, result_file);
  fclose(result_file);

  return EXIT_SUCCESS;

}

// Display usage information
int usage(char* name){
  
  PRINTF_ERROR("Usage: %s\n", name);
  PRINTF_ERROR("Environment variables: \n");
  PRINTF_ERROR("\tUSAGE: print usage information\n");
  PRINTF_ERROR("\tREPS: number of reps to do the algorithm, value must be between 0 and %d, default (0) \n", INT32_MAX);
  PRINTF_ERROR("\tALIGNED: explicitly aligned data to page boundary, 0 -> inactive, 1 -> active, default (1) \n");
  PRINTF_ERROR("\tDEBUG: print debug information, 0 -> inactive, 1 -> active, default (0) \n");
  PRINTF_ERROR("\tTIMES: print timing information, 0 -> inactive, 1 -> active, default (0) \n");
  PRINTF_ERROR("\tCHECK: file to check the results\n");
  PRINTF_ERROR("\tWRITE_RESULT: file to write the results\n");

  return EXIT_FAILURE;

}


int main(int argc,char* argv[]) {

  PRINTF("\n");

  // Get program name
  char* name = argv[0];

  if(argc != 1){
    return usage(name);
  }


  // --------------------------------------------------------------------------------------------------------
  // Environment variables


  // USAGE variable
  const char* susage = getenv("USAGE");
  if (susage != NULL){
    return usage(name);
  }


  // Int REPS variable
  const char* sreps = getenv("REPS");
  int aux_reps = 0;
  if (sreps != NULL) {
    int aux = atoi(sreps);
    if (aux < 0){
      return usage(name);
    }
    aux_reps = aux;
  }
  const int reps = aux_reps;


  // Bool ALIGNED variable
  const char* saligned = getenv("ALIGNED");
  bool aux_aligned = true;
  if (saligned != NULL) {
    if (strcmp(saligned,"0")){
      aux_aligned = false;
    }
    else if(strcmp(saligned,"1")){
      aux_aligned = true;
    }
    else return usage(name);
  }
  const bool aligned = aux_aligned;


  // Bool DEBUG variable
  const char* sdebug = getenv("DEBUG");
  bool aux_debug = false;
  if (sdebug != NULL) {
    if (!strcmp(sdebug,"0")){
      aux_debug = false;
    }
    else if(!strcmp(sdebug,"1")){
      aux_debug = true;
    }
    else{
      PRINTF_ERROR("Invalid value for DEBUG\n");
      return usage(name);
    }
  }
  const bool debug = aux_debug;


  // Bool TIMES variable
  const char* stimes = getenv("TIMES");
  bool aux_times = false;
  if (stimes != NULL) {
    if (!strcmp(stimes,"0")){
      aux_times = false;
    }
    else if(!strcmp(stimes,"1")){
      aux_times = true;
    }
    else{
      PRINTF_ERROR("Invalid value for TIMES\n");
      return usage(name);
    }
  }
  const bool times = aux_times;


  // String CHECK variable
  const char* scheck = getenv("CHECK");
  if (scheck != NULL){
    if (access(scheck, F_OK) == 0){
      if (access(scheck, R_OK) != 0){
        PRINTF_ERROR("Check file %s is not readable\n", scheck);
        return EXIT_FAILURE;
      }
    }
    else{
      PRINTF_ERROR("Check file %s does not exist\n", scheck);
      return EXIT_FAILURE;
    }
  }
  const char* cfilename = scheck;
  const bool check = (cfilename != NULL);


  // String WRITE_RESULT variable
  const char* swrite_result = getenv("WRITE_RESULT");
  if (swrite_result != NULL){
    if (access(swrite_result, F_OK) == 0){
      if (access(swrite_result, W_OK) != 0){
        PRINTF_ERROR("File %s is not writable\n", swrite_result);
        return EXIT_FAILURE;
      }
    }
  }
  const char* rfilename = swrite_result;
  const bool write_result = (rfilename != NULL);

  

  // --------------------------------------------------------------------------------------------------------


  // Buffers
  char* pattern_mem_noalign =
      "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";
  char* text_mem_noalign    =
      "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

  // Buffers length
  const uint32_t pattern_length = strlen(pattern_mem_noalign);
  const uint32_t text_length = strlen(text_mem_noalign);
  const uint32_t max_distance = pattern_length + text_length;
  // Pattern & Text
  char* pattern = NULL;
  char* text = NULL;
  
  size_t page_size = 0;
  if (aligned){

    // Determine the page size
    page_size = sysconf(_SC_PAGESIZE);

    // Allocate memory aligned to page size for pattern string
    pattern = (char*) aligned_alloc(page_size,max_distance*sizeof(char));
    if (pattern == NULL) {
      PRINTF_ERROR("Aligned allocation of pattern failed");
      return EXIT_FAILURE;
    }

    // Allocate memory aligned to page size for text string 
    text = (char*) aligned_alloc(page_size,max_distance*sizeof(char));
    if (text == NULL) {
      PRINTF_ERROR("Aligned allocation of text failed");
      return EXIT_FAILURE;
    }

    memcpy(pattern, pattern_mem_noalign, pattern_length);
    memcpy(text,text_mem_noalign, text_length);

  }
  else{
    pattern = pattern_mem_noalign;
    text = text_mem_noalign;
  }

  //DEBUG(debug,"Edit wavefront data type size: %ld, position: %p\n",sizeof(edit_wavefront_t), pattern);
  //DEBUG(debug,"Edit wavefronts data type size: %ld, position: %p\n",sizeof(edit_wavefronts_t), text);

  PRINTF("#######################################################################################\n");
  PRINTF("Configuration summary:\n");

  PRINTF("\n\n");
  PRINTF("Environment variables\n");
  PRINTF("\tRepetitions: %d\n",reps);
  PRINTF("\tAligned: %d\n",aligned);
  PRINTF("\tDebug: %d\n",debug);
  PRINTF("\tTimes: %d\n",times);
  PRINTF_COND(check,"\tCheck results from filename: %s\n",cfilename);
  PRINTF_COND(write_result,"\tWrite result to filename: %s\n",rfilename);

  PRINTF("\n");
  PRINTF("Pattern length: %d\n",pattern_length);
  PRINTF("Text length: %d\n",text_length);
  PRINTF("\n");

  PRINTF("#######################################################################################\n");
  PRINTF("\n");

  edit_wavefronts_fpga_t wavefronts;

  // Initialize wavefronts
  PRINTF("\nInitializing wavefronts\n");
  const double tStartInit = wall_time();
  edit_wavefronts_init(&wavefronts,pattern_length,text_length,aligned,page_size);
  const double tEndInit = wall_time();
  PRINTF("Wavefronts initialized\n");
  PRINTF_COND(times,"Init time: %f\n", tEndInit-tStartInit);

  // Warm up execution
  /*PRINTF("\nWarming up...\n");
  const double tStartWarm = wall_time();
  edit_wavefronts_align(wavefronts.offsets,wavefronts.edit_cigar,&wavefronts.edit_cigar_length,pattern,pattern_length,text,text_length,max_distance,&wavefronts.edit_cigar_length);
  #pragma oss taskwait
  const double tEndWarm = wall_time();
  PRINTF("Warm up finished\n");
  PRINTF_COND(times,"Warm up time: %f\n", tEndWarm-tStartWarm);*/

  int i;
  int score = 0;
  for (i=0;i<reps;++i) {

    PRINTF("\n---------------------------------------------------------------------------------------\n");

    PRINTF("\nRepetition: %d\n",i);

    // Align Wavefronts
    PRINTF("\nAligning...\n");
    const double tStartAlign = wall_time();
    edit_wavefronts_align(wavefronts.offsets,wavefronts.edit_cigar,&wavefronts.edit_cigar_length,pattern,pattern_length,text,text_length,max_distance,&score);
    FPGA("oss taskwait")
    const double tEndAlign = wall_time();
    PRINTF("Alignment finished\n");
    PRINTF_COND(times,"WFA execution time: %f\n",tEndAlign-tStartAlign);

    // Check results
    PRINTF_COND(check,"\nChecking results...\n");
    const double tStartCheck = wall_time();
    if (check){
      if(edit_wavefronts_check(wavefronts.edit_cigar,wavefronts.edit_cigar_length,score,cfilename)) {
        return EXIT_FAILURE;
      }
    }
    PRINTF_COND(check,"Check finished\n");
    const double tEndCheck = wall_time();
    PRINTF_COND(check && times,"Check results time: %f\n", tEndCheck-tStartCheck);

    // Write results
    PRINTF_COND(write_result && !i,"\nWriting results...\n");
    const double tStartWrite = wall_time();
    if (write_result && !i){
      if(edit_wavefronts_write_result(wavefronts.edit_cigar,wavefronts.edit_cigar_length,score,rfilename)){
        return EXIT_FAILURE;
      }
    }
    PRINTF_COND(write_result && !i,"Results written\n");
    const double tEndWrite = wall_time();
    PRINTF_COND(write_result && !i && times,"Write results time: %f\n", tEndWrite-tStartWrite);
  }

  // Clean Wavefronts Offsets
  PRINTF("\nCleaning wavefronts offsets...\n");
  const double tStartClean = wall_time();
  edit_wavefronts_clean(&wavefronts);
  const double tEndClean = wall_time();
  PRINTF("Cleaning finished\n");
  PRINTF_COND(times,"Clean time: %f\n", tEndClean-tStartClean);

  if(aligned){
    PRINTF("\nFreeing memory space used for pattern and text...\n");
    const double tStartClean = wall_time();
    free(pattern);
    free(text);
    const double tEndClean = wall_time();
    PRINTF("Free finished\n");
    PRINTF_COND(times,"Free time: %f\n", tEndClean-tStartClean);
  }

}









