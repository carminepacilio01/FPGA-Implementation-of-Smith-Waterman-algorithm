#include "host.h"
#include <random>
#include <time.h>

#include "common/xcl2.hpp"

#define MAX_SEQ_LEN 1000

#define NUM_KERNEL 1

#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
    BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
    BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
    BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
    BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
    BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
    BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
    BANK_NAME(30), BANK_NAME(31)};

void printConf(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement);
void fprintMatrix(int *P, int *D, int *Q, int lenA, int lenB);
int compute_golden(int lenA, char *seqA, int lenB, char *seqB, int wd, int ws, int gap_opening, int enlargement);
void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB);
int gen_rnd(int min, int max);

int main(int argc, char* argv[]){
	
	//TARGET_DEVICE macro needs to be passed from gcc command line
	// if(argc != 3) {
	// 	std::cout << "Usage: " << argv[0] <<" <xclbin> <dataset path>" << std::endl;
	// 	return EXIT_FAILURE;
	// }

    // FILE* f = fopen(argv[2], "r");

    if(argc < 2) {
		std::cout << "Usage: " << argv[0] <<" <xclbin>" << std::endl;
		return EXIT_FAILURE;
	}
    
    // double frequency = 350000000;

    srandom(2);
	// clock_t start, end;

	//	match score
	const int wd = 1;
	//	mismatch score
	const int ws = -1;

	const int gap_opening = -3;
	const int enlargement = -1;

	const int n = INPUT_SIZE;
	int lenA[INPUT_SIZE];
	int lenB[INPUT_SIZE];
	char seqA[INPUT_SIZE][MAX_DIM];
	char seqB[INPUT_SIZE][MAX_DIM];

	for(int i = 0; i < n; i++){

		//	generate random sequences
		//	length of the sequences
		lenA[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);
		lenB[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);

		lenA[i] += 1;
		lenB[i] += 1;

		//	sequences
		//seqA[i] = (char *) malloc(sizeof(char) * (lenA[i]));
		//seqB[i] = (char *) malloc(sizeof(char) * (lenB[i]));

		//	generate rand sequences
		seqA[i][0] = seqB[i][0] = '-';
		seqA[i][lenA[i]] = seqB[i][lenA[i]] = '\0';
		random_seq_gen(lenA[i], seqA[i], lenB[i], seqB[i]);

		//	Printing current configuration
	}

    int score[INPUT_SIZE];
	int golden_score[INPUT_SIZE];
	int n_ops = 0;

	double mean_golden_time = 0;
	double mean_golden_gcup = 0;

	for (int golden_rep = 0; golden_rep < n; golden_rep++) {

		auto start = clock();
		golden_score[golden_rep] = compute_golden(lenA[golden_rep], seqA[golden_rep], lenB[golden_rep], seqB[golden_rep], wd, ws, gap_opening, enlargement);
		auto end = clock();

		n_ops = lenA[golden_rep]*lenB[golden_rep];
		double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		double gcup = (double) (n_ops / cpu_time_used )  * 1e-9;

		mean_golden_time += cpu_time_used;
		mean_golden_gcup += gcup;
	}
	

/*
================================================================================================================================
	OPENCL STUFF
================================================================================================================================
*/

   	std::string binaryFile = argv[1]; // prendo il bitstream 
    auto fileBuf = xcl::read_binary_file(binaryFile); // leggi bitstream
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

    cl_int err;
    cl::Context context;
    cl::Kernel sw_maxi;
    cl::CommandQueue q;

    auto devices = xcl::get_xil_devices(); // lista di devices

    bool valid_device = false;

    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, sw_maxi= cl::Kernel(program, "decoderTop", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    cl_mem_ext_ptr_t lenT_buffer_ext;
    cl_mem_ext_ptr_t target_buffer_ext;
    cl_mem_ext_ptr_t lenD_buffer_ext;
    cl_mem_ext_ptr_t database_buffer_ext;
    cl_mem_ext_ptr_t score_buffer_ext;

    lenT_buffer_ext.param = 0;
    lenD_buffer_ext.param = 0;
    target_buffer_ext.param = 0;
    database_buffer_ext.param = 0;
    score_buffer_ext.param = 0;

    lenT_buffer_ext.flags = bank[0];
    target_buffer_ext.flags = bank[1];
    lenD_buffer_ext.flags = bank[2];
    database_buffer_ext.flags = bank[3];
    score_buffer_ext.flags = bank[4];

    lenT_buffer_ext.obj = lenA;
    target_buffer_ext.obj = seqA;
    lenD_buffer_ext.obj = lenB;
    database_buffer_ext.obj = seqB;
    score_buffer_ext.obj = score;

    cl::Buffer lenT_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, INPUT_SIZE * sizeof(int), &lenT_buffer_ext);
    cl::Buffer target_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, INPUT_SIZE * MAX_DIM * sizeof(char), &target_buffer_ext);
    cl::Buffer lenD_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, INPUT_SIZE * sizeof(int), &lenD_buffer_ext);
    cl::Buffer database_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, INPUT_SIZE * MAX_DIM * sizeof(char), &database_buffer_ext);
    cl::Buffer score_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, INPUT_SIZE * sizeof(int), &score_buffer_ext);

    q.finish();
/*
================================================================================================================================
	MY STUFF
================================================================================================================================
*/

	OCL_CHECK(err, err = sw_maxi.setArg(0, lenT_buffer));
	OCL_CHECK(err, err = sw_maxi.setArg(1, target_buffer));
    OCL_CHECK(err, err = sw_maxi.setArg(2, lenD_buffer));
    OCL_CHECK(err, err = sw_maxi.setArg(3, database_buffer));
    OCL_CHECK(err, err = sw_maxi.setArg(4, wd));
    OCL_CHECK(err, err = sw_maxi.setArg(5, ws));
    OCL_CHECK(err, err = sw_maxi.setArg(6, gap_opening));
    OCL_CHECK(err, err = sw_maxi.setArg(7, enlargement));
    OCL_CHECK(err, err = sw_maxi.setArg(8, score_buffer));

    // Data will be migrated to kernel space
    OCL_CHECK(err, q.enqueueMigrateMemObjects({lenT_buffer, target_buffer, lenD_buffer, database_buffer}, 0)); /*0 means from host*/
	OCL_CHECK(err, q.finish());

	auto start = std::chrono::high_resolution_clock::now();
	//Launch the Kernel
	q.enqueueTask(sw_maxi);
	OCL_CHECK(err, q.finish());
	
	auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    printf("Operation concluded in %f nanoseconds\n", (float)duration.count());
	
	
	//Data from Kernel to Host
	q.enqueueMigrateMemObjects({score_buffer}, CL_MIGRATE_MEM_OBJECT_HOST);
	q.finish();
	
	

	return 0;
}

//	Prints the current configuration
void printConf(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement) {
	std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
	std::cout << "+ Sequence A: [" << strlen(seqA) << "]: " << seqA << std::endl;
	std::cout << "+ Sequence B: [" << strlen(seqB) << "]: " << seqB << std::endl;
	std::cout << "+ Match Score: " << wd << std::endl;
	std::cout << "+ Mismatch Score: " << ws << std::endl;
	std::cout << "+ Gap Opening: " << gap_opening << std::endl;
	std::cout << "+ Enlargement: " << enlargement << std::endl;
	std::cout << "+++++++++++++++++++++" << std::endl;
}

int gen_rnd(int min, int max) {
     // Using random function to get random double value
    return (int) min + rand() % (max - min + 1);
}

void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB) {

	int i; 
	for(i = 1; i < lenA; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqA[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
			      (tmp_gen == 2) ? 'G' : 'T';
	}

	for(i = 1; i < lenB; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqB[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
				  (tmp_gen == 2) ? 'G' : 'T';
	}
}

int compute_golden(int lenA, char *seqA, int lenB, char *seqB, int wd, int ws, int gap_opening, int enlargement) {
	// Inizializza le matrici D, P, Q
	    std::vector< std::vector<int> > D(lenA, std::vector<int>(lenB, 0));
	    std::vector< std::vector<int> > P(lenA, std::vector<int>(lenB, std::numeric_limits<int>::min() / 2));
	    std::vector< std::vector<int> > Q(lenA, std::vector<int>(lenB, std::numeric_limits<int>::min() / 2));

	    int gap_penalty = gap_opening + enlargement * 1;
	    int max_score = 0;

	    for (int i = 1; i < lenA; ++i) {
	        for (int j = 1; j < lenB; ++j) {
	            // Calcola P[i][j]
	            P[i][j] = std::max(P[i-1][j] + enlargement, D[i-1][j] + gap_penalty);

	            // Calcola Q[i][j]
	            Q[i][j] = std::max(Q[i][j-1] + enlargement, D[i][j-1] + gap_penalty);

	            // Calcola D[i][j]
	            int match = (seqA[i] == seqB[j]) ? wd : ws;
	            D[i][j] = std::max(0, D[i-1][j-1] + match);
	            D[i][j] = std::max(D[i][j], P[i][j]);
	            D[i][j] = std::max(D[i][j], Q[i][j]);

	            // Aggiorna lo score massimo
	            max_score = std::max(max_score, D[i][j]);
	        }
	    }

	    return max_score;
}
