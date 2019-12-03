/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <iostream>
#include <string>
#include <cstring>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <atomic>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <map>
#include <thread>
#include <mutex>
#include <algorithm> 
#include "fasta.h"
#include "ntcoding.h"
#include "seed_pos_table.h"
#include "dsoft_hash_table.hpp"
#include "gact.h"
#include "ConfigFile.h"


// GACT scoring
int gact_sub_mat[25];
int gap_open;
int gap_extend;

// D-SOFT parameters
std::string seed_shape;
std::string seed_shape_str;
uint32_t bin_size;
int dsoft_threshold;
int num_seeds;
int seed_occurence_multiple;
int max_candidates;
int num_nz_bins;

#ifdef FPGA
// Function prototypes
bool init();
void cleanup();
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );
void checkErr(cl_int err, const char * name);
cl_platform_id platform = NULL; 
scoped_array<cl_device_id> devices; // num_devices elements
cl_context context = NULL;
scoped_array<cl_command_queue> queues; // num_devices elements
cl_program program = NULL;
scoped_array<cl_kernel> kernels; // num_devices elements
scoped_array<cl_event> kernel_events;
unsigned int num_devices=NUM_DEVICES, max_num_devices=MAX_NUM_DEVICES;
#endif
int forward_kernel_counter=0;
int reverse_kernel_counter=0;


// GACT first tile
int first_tile_size;
int first_tile_score_threshold;

//GACT extend 
int tile_size;
int tile_overlap;

//Multi-threading
int num_threads;

struct timeval start, end_time;

static std::string reference_string;
static std::string query_string;

uint32_t reference_length;
uint32_t query_length;

std::vector<long long int> reference_lengths;
std::vector<std::string> reference_seqs;

std::vector<long long int> reads_lengths;
std::vector<std::string> reads_seqs;
std::vector<std::string> rev_reads_seqs;

std::vector<std::vector<std::string> > reference_descrips;
std::vector<long long int> reference_fileposs;

std::vector<std::vector<std::string> > reads_descrips;
std::vector<long long int> reads_fileposs;

char* reference_char;
char** reads_char;
char** rev_reads_char;

std::map<int, uint32_t> chr_id_to_start_bin;
std::map<uint32_t, int> bin_to_chr_id;

SeedPosTable *sa;
DsoftHashTable *dsht;

std::mutex io_lock;

std::atomic<int> num_aligned(0);

std::string RevComp(std::string seq) {
    std::string rc = "";
    for (int i = seq.size()-1; i >= 0; i--) {
        if (seq[i] != 'a' && seq[i] != 'A' &&
                seq[i] != 'c' && seq[i] != 'C' &&
                seq[i] != 'g' && seq[i] != 'G' &&
                seq[i] != 't' && seq[i] != 'T' &&
                seq[i] != 'n' && seq[i] != 'N') {
            std::cerr<<"Bad Nt char: "<< seq[i] <<std::endl;
            exit(1);
        }
        else {
            switch (seq[i]) {
                case 'a': rc += 't';
                          break;
                case 'A': rc += 'T';
                          break;
                case 'c': rc += 'g';
                          break;
                case 'C': rc += 'G';
                          break;
                case 'g': rc += 'c';
                          break;
                case 'G': rc += 'C';
                          break;
                case 't': rc += 'a';
                          break;
                case 'T': rc += 'A';
                          break;
                case 'n': rc += 'n';
                          break;
                case 'N': rc += 'N';
                          break;
            }
        }
    }
    return rc;
}


bool CompareAlignments(Alignment a1, Alignment a2) {
    return (a1.score > a2.score);
}

int mode;
int l1l2enable;

void AlignRead (int start_read_num, int last_read_num, int num_threads) {

    uint32_t log_bin_size = (uint32_t) (log2(bin_size));
    int num_bins = 1 + (reference_length >> log_bin_size);
    uint64_t* candidate_hit_offset;

    uint32_t* nz_bins_array = new uint32_t[num_nz_bins];
    uint64_t* bin_count_offset_array = new uint64_t[num_bins];
    candidate_hit_offset = new uint64_t[max_candidates];

    for(int i=0; i < num_bins; i++) {
        bin_count_offset_array[i] = 0;
    }

    for (int k = start_read_num; k < last_read_num; k+=num_threads) {

        int len = reads_lengths[k];
        vector<Alignment> alignments;
        alignments.clear();
        int thread_id = k % num_threads;

        // Forward reads
        
        int num_candidates;
        if(l1l2enable)
           num_candidates = dsht->DSOFT(reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);
        else
           num_candidates = sa->DSOFT(reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);


        // TODO: Need to apply more filtering here to avoid duplicate and
        // unnecessary alignments
        for (int i = 0; i < num_candidates; i++) {
            uint32_t candidate_hit = (candidate_hit_offset[i] >> 32);
            uint32_t last_hit_offset = ((candidate_hit_offset[i] << 32) >> 32);

            int chr_id = bin_to_chr_id[candidate_hit/bin_size];
            std::string chrom = reference_descrips[chr_id][0];
            uint32_t start_bin = chr_id_to_start_bin[chr_id];

            uint32_t ref_pos = candidate_hit - (start_bin*bin_size);
            uint32_t query_pos = last_hit_offset;
            uint32_t ref_len = reference_lengths[chr_id];
            uint32_t query_len = reads_lengths[k];
            char* ref_start = reference_char + (start_bin*bin_size);
            char strand = '+';

            Alignment align = GACT(ref_start, reads_char[k], chrom, reads_descrips[k][0], gact_sub_mat, gap_open, gap_extend, tile_size, tile_overlap, ref_pos, query_pos, ref_len, query_len, strand, first_tile_score_threshold, mode, thread_id);
            alignments.push_back(align);
        }

        // Reverse complement reads
        if(l1l2enable)
             num_candidates = dsht->DSOFT(rev_reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);
        else
             num_candidates = sa->DSOFT(rev_reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);

        // TODO: Need to apply more filtering here to avoid duplicate and
        // unnecessary alignments
        for (int i = 0; i < num_candidates; i++) {
            uint32_t candidate_hit = (candidate_hit_offset[i] >> 32);
            uint32_t last_hit_offset = ((candidate_hit_offset[i] << 32) >> 32);

            int chr_id = bin_to_chr_id[candidate_hit/bin_size];
            std::string chrom = reference_descrips[chr_id][0];
            uint32_t start_bin = chr_id_to_start_bin[chr_id];

            uint32_t ref_pos = candidate_hit - (start_bin*bin_size);
            uint32_t query_pos = last_hit_offset;
            uint32_t ref_len = reference_lengths[chr_id];
            uint32_t query_len = reads_lengths[k];
            char* ref_start = reference_char + (start_bin*bin_size);
            char strand = '-';
            

            Alignment align = GACT(ref_start, rev_reads_char[k], chrom, reads_descrips[k][0], gact_sub_mat, gap_open, gap_extend, tile_size, tile_overlap, ref_pos, query_pos, ref_len, query_len, strand, first_tile_score_threshold, mode, thread_id);
            alignments.push_back(align);
        }

        std::stable_sort(alignments.begin(), alignments.end(), CompareAlignments);
        int num_alignments = alignments.size();
        int* flags = (int*) calloc(num_alignments, sizeof(int));

        for (int m=0; m < num_alignments; m++) {
            if (flags[m] < 0) {
                continue;
            }
            if (alignments[m].aligned_query_len == 0) {
                flags[m] = 0;
                continue;
            }
            uint32_t s1 = alignments[m].query_start; 
            uint32_t e1 = s1+alignments[m].aligned_query_len; 
            for (int n=m+1; n < num_alignments; n++) {
                if (flags[n] < 0) {
                    continue;
                }
                uint32_t s2 = alignments[n].query_start; 
                uint32_t e2 = s2+alignments[n].aligned_query_len; 
                uint32_t s = std::max(s1 , s2);
                uint32_t e = std::min(e1, e2);
                uint32_t overlap = 0;
                if (s < e) {
                    overlap = e-s;
                }
                if (2*overlap >= alignments[n].aligned_query_len) {
                    flags[n] = -1;
                }
            }
        }

        io_lock.lock();
        for (int m=0; m < num_alignments; m++) {
            if (flags[m] >= 0) {
                Alignment align = alignments[m];
                std::cout << "a score=" << align.score << std::endl;
                std::cout << "s\t" << align.ref_name << "\t" << 1+align.ref_start << "\t" << align.aligned_ref_len << "\t+\t" << align.ref_len << "\t" << align.aligned_ref_str << std::endl;
                std::cout << "s\t" << align.query_name << "\t" << 1+align.query_start << "\t" << align.aligned_query_len << "\t" << align.strand << "\t" << align.query_len << "\t" << align.aligned_query_str << std::endl;
                std::cout << std::endl;
            }
        }
        io_lock.unlock();
        int n = ++num_aligned;
        if (n % 100 == 0) {
            io_lock.lock();
            std::cerr << n << " reads aligned calling forward kernel " << forward_kernel_counter << " times and reverse kernel " << 
									  reverse_kernel_counter << " times\n";
            io_lock.unlock();
        }
    }

    delete[] bin_count_offset_array;
    delete[] nz_bins_array;
    delete[] candidate_hit_offset;
}


int main(int argc, char *argv[]) {

    if (argc < 5) {
        std::cerr << "Usage: ./darwin_ref_guided <REFERENCE>.fasta <READS>.fasta mode l1l2enable"<< endl;
	std::cerr << "(mode is 0: default darwin original, 1: cpu darwin-xl, 2:fpga darwin-xl)"<< endl;
	std::cerr << "(l1l2enable is 0: default darwin original, 1:l1l2enable)"<< endl;
        exit(1);
    }


    gettimeofday(&start, NULL);
    std::cerr << "\nReading configuration file ..." << std::endl;
    std::ifstream config_file("params.cfg");
    if (!config_file) {
        std::cerr << "Configuration file <params.cfg> not found! Exiting." << std::endl;
        exit(1);
    }
    ConfigFile cfg("params.cfg");

    // GACT scoring
    int sub_N = cfg.Value("GACT_scoring", "sub_N");
    for (int i = 0; i < 25; i++) {
        gact_sub_mat[i] = sub_N;
    }
    gact_sub_mat[0] = cfg.Value("GACT_scoring", "sub_AA");
    gact_sub_mat[1] = gact_sub_mat[5]  = cfg.Value("GACT_scoring", "sub_AC");
    gact_sub_mat[2] = gact_sub_mat[10] = cfg.Value("GACT_scoring", "sub_AG");
    gact_sub_mat[3] = gact_sub_mat[15] = cfg.Value("GACT_scoring", "sub_AT");
    gact_sub_mat[6] = cfg.Value("GACT_scoring", "sub_CC");
    gact_sub_mat[7] = gact_sub_mat[11] = cfg.Value("GACT_scoring", "sub_CG");
    gact_sub_mat[8] = gact_sub_mat[16] = cfg.Value("GACT_scoring", "sub_CT");
    gact_sub_mat[12] = cfg.Value("GACT_scoring", "sub_GG");
    gact_sub_mat[13] = gact_sub_mat[17] = cfg.Value("GACT_scoring", "sub_GT");
    gact_sub_mat[18] = cfg.Value("GACT_scoring", "sub_TT");
    gap_open        = cfg.Value("GACT_scoring", "gap_open");
    gap_extend      = cfg.Value("GACT_scoring", "gap_extend");

    // D-SOFT parameters
    seed_shape_str          = (std::string) cfg.Value("DSOFT_params", "seed_shape");
    bin_size                = cfg.Value("DSOFT_params", "bin_size");
    dsoft_threshold         = cfg.Value("DSOFT_params", "threshold");
    num_seeds               = cfg.Value("DSOFT_params", "num_seeds");
    seed_occurence_multiple = cfg.Value("DSOFT_params", "seed_occurence_multiple");
    max_candidates          = cfg.Value("DSOFT_params", "max_candidates");
    num_nz_bins             = cfg.Value("DSOFT_params", "num_nz_bins");

    // GACT first tile
    first_tile_size            = cfg.Value("GACT_first_tile", "first_tile_size");
    first_tile_score_threshold = cfg.Value("GACT_first_tile", "first_tile_score_threshold");
    std::cerr << "Running with configuration" << " bin size=" << bin_size << " threshold=" << dsoft_threshold << " ft_size=" << first_tile_size << " ft_threshold=" << first_tile_score_threshold << endl;

    // GACT extend
    tile_size    = cfg.Value("GACT_extend", "tile_size");
    tile_overlap = cfg.Value("GACT_extend", "tile_overlap");

    // Multi-threading
    num_threads = cfg.Value("Multithreading", "num_threads");

    seed_shape = seed_shape_str.c_str();

    std::string reference_filename(argv[1]);
    std::string reads_filename(argv[2]);
    mode = std::stoi(argv[3]);
    l1l2enable = std::stoi(argv[4]);

    if(mode == 1)
	std::cerr<< "Using Darwin-XL CPU implementation"<<std::endl;
    else if(mode == 2) {
	std::cerr<< "Using Darwin-XL FPGA implementation  supported on this"<<std::endl;
    }
    else
	std::cerr<< "Using Old Darwin CPU implementation"<<std::endl;

#ifdef FPGA
    if(!init()) {
        std::cerr << "Device initialization failed \n";
        exit(1);
    }
#endif

    gettimeofday(&end_time, NULL);
    long useconds = end_time.tv_usec - start.tv_usec;
    long seconds = end_time.tv_sec - start.tv_sec;
    long mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (reading configuration file): " << mseconds <<" msec" << std::endl;

    // LOAD REFERENCE
    std::cerr << "\nLoading reference genome ...\n";
    gettimeofday(&start, NULL);

    ParseFastaFile(reference_filename, reference_descrips, reference_seqs, reference_lengths, reference_fileposs);

    reference_string = "";

    int curr_bin = 0;

    for (size_t i=0; i < reference_seqs.size(); i++) {
        chr_id_to_start_bin[i] =  curr_bin;
        reference_string += reference_seqs[i];
        for (size_t j = 0; j < (reference_seqs[i].length() / bin_size); j++) {
            bin_to_chr_id[curr_bin++] = i;
        }
        if (reference_seqs[i].length() % bin_size > 0) {
            reference_string += std::string((bin_size - (reference_seqs[i].length() % bin_size)), 'N');
            bin_to_chr_id[curr_bin++] = i;
        }
    }

    reference_length = reference_string.length();
    std::cerr << "Reference length (after padding): " << (unsigned int) reference_length << std::endl;

    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (loading reference genome): " << mseconds <<" msec" << std::endl;

    // LOAD READS
    std::cerr << "\nLoading reads ...\n";
    gettimeofday(&start, NULL);

    ParseFastaFile(reads_filename, reads_descrips, reads_seqs, reads_lengths, reads_fileposs);
    int num_reads = reads_seqs.size();
    for (int i = 0; i < num_reads; i++) {
        std::string rev_read = RevComp(reads_seqs[i]);
        rev_reads_seqs.push_back(rev_read);
    }

    std::cerr << "Number of reads: " << num_reads << std::endl;

    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (loading reads): " << mseconds <<" msec" << std::endl;

    reference_char = (char*) reference_string.c_str();

    reads_char = new char*[num_reads];
    rev_reads_char = new char*[num_reads];
    for (int i =0; i < num_reads; i++) {
        reads_char[i] = (char*) reads_seqs[i].c_str();
        rev_reads_char[i] = (char*) rev_reads_seqs[i].c_str();
    }

    // CONSTRUCT SEED POSITION TABLE
    gettimeofday(&start, NULL);

    if(l1l2enable) {
       //dsht = new DsoftHashTable(reference_char, reference_length, seed_shape, seed_occurence_multiple, bin_size);
       std::cerr << "\nConstructing l1 l2 hash table ...\n";
       dsht = new DsoftHashTable(reference_char, reference_length, seed_shape, seed_occurence_multiple, bin_size);
    }
    else {
       std::cerr << "\nConstructing seed position table ...\n";
       sa = new SeedPosTable(reference_char, reference_length, seed_shape, seed_occurence_multiple, bin_size);
    }

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (seed position table construction): " << mseconds <<" msec" << std::endl;

    // RUN D-SOFT TO MAP READS
    gettimeofday(&start, NULL);

    std::cerr << "\nFinding candidate bin locations for each read: " << std::endl;

    std::vector<std::thread> align_threads;
    for (int k = 0; k < num_threads; k++) {
        align_threads.push_back(std::thread(AlignRead, k, num_reads, num_threads));
    }
    std::cerr << "Using " << align_threads.size() << " threads ...\n";
    for (auto& th : align_threads) th.join();
    std::cerr << "Synchronizing threads. " << num_aligned << " reads aligned\n";

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (seed table querying): " << mseconds <<" msec" << std::endl;

#ifdef FPGA
    cleanup();
#endif
    return 0;
}


#ifdef FPGA
/////// OPENCL HELPER FUNCTIONS ///////

bool init() {
  cl_int status;

  if(!setCwdToExeDir()) {
    return false;
  }

  // Get the OpenCL platform.
  platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform.\n");
    return false;
  }

  // Query the available OpenCL devices.
  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &max_num_devices));
  printf("Number of devices is %d\n", max_num_devices);

  // Create the context.
  //context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &status);
  context = clCreateContext(NULL, num_devices, &devices[0], NULL, NULL, &status);
  checkError(status, "Failed to create context");

  queues.reset(num_devices*num_threads);
  kernels.reset(num_devices*num_threads);
  kernel_events.reset(num_devices*num_threads*20);

  // Create command queue.
  for(int i = 0; i<num_devices ; i++) {
    for (int j = 0; j<num_threads ; j++){
          int id = (i*num_devices)+j;
          //queues[id] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &status);
          queues[id] = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
          checkError(status, "Failed to create command queue");
    }
  }


  // Create the program.
  std::string binary_file = getBoardBinaryFile(PRECOMPILED_BINARY, devices[0]);
  printf("Using AOCX: %s\n", binary_file.c_str());
  program = createProgramFromBinary(context, binary_file.c_str(), devices, num_devices);

  // Build the program that was just created.
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  for(int i= 0; i<num_devices ; i++) {
  for (int j = 0; j<num_threads ; j++){
        // Create the kernel - name passed in here must match kernel name in the
        // original CL file, that was compiled into an AOCX file using the AOC tool
        const char * kernel_name = "xl";  // Kernel name, as defined in the CL file
        kernels[(i*num_devices)+j] = clCreateKernel(program, kernel_name, &status);
        checkError(status, "Failed to create kernel");
  }
  }

  return true;
}


void cleanup() {
  for(int i = 0; i < num_devices; i++) {
    for(unsigned j = 0; j < num_threads; j++) {
      int id = (i*num_devices)+j;
      if(kernels[id]) {
        clReleaseKernel(kernels[id]);
      }
      if(queues[id]) {
        clReleaseCommandQueue(queues[id]);
      }
    }
  }
  if(program) {
    clReleaseProgram(program);
  }
  if(context) {
    clReleaseContext(context);
  }
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
   cl_ulong a;
   clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
   printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
   cl_uint a;
   clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
   printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
   cl_bool a;
   clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
   printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
   char a[STRING_BUFFER_LEN];
   clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
   printf("%-40s = %s\n", name, a);
}


// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

   printf("Querying device for info:\n");
   printf("========================\n");
   device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
   device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
   device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
   device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
   device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
   device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
   device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
   device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
   device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
   device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
   device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
   device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
   device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

   {
      cl_command_queue_properties ccp;
      clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
      printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
      printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
   }
}

/* Error checking */
void checkErr(cl_int err, const char * name)
{
        if (err != CL_SUCCESS) {
                printf("\n ERROR (%d): %s\n",err,name);
                exit(EXIT_FAILURE);
        }
}

#endif // end OpenCL helper

