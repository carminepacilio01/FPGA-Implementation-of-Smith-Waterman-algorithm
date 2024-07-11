#include "host.h"

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

void printConf(char *target, char *database, int ws, int wd, int gap_opening, int enlargement);
int compute_golden(int lenT, char *target, int lenD, char *database, int wd, int ws, int gap_opening, int enlargement);
void random_seq_gen(int lenT, char *target, int lenD, char *database);
int gen_rnd(int min, int max);

int main(int argc, char* argv[]){

    if(argc < 2) {
		std::cout << "Usage: " << argv[0] <<" <xclbin>" << std::endl;
		return EXIT_FAILURE;
	}

    srandom(static_cast<unsigned>(time(0)));

	const int wd 			=  1; 	// match score
	const int ws 			= -1;	// mismatch score
	const int gap_opening	= -3;
	const int enlargement 	= -1;


    cl_int err;
    cl::Context context;
    cl::Kernel sw_maxi;

    cl::CommandQueue q;
	std::vector< int, aligned_allocator<int> > lenT(INPUT_SIZE);
    std::vector< char, aligned_allocator<char> > target(INPUT_SIZE * MAX_DIM);
    std::vector< int, aligned_allocator<int> > lenD(INPUT_SIZE);
	std::vector< char, aligned_allocator<char> > database(INPUT_SIZE * MAX_DIM);
    std::vector< int, aligned_allocator<int> > score(INPUT_SIZE);

/////////////////////////		DATASET GENERATION 		////////////////////////////////////

	printf("Generating %d random sequences pairs. \n", INPUT_SIZE);
	///////Generation of random sequences
	for(int i = 0; i < INPUT_SIZE; i++){

		//	generate random length of the sequences
		lenT[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);
		lenD[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);

		lenT[i] += 1;
		lenD[i] += 1;

		//	generate rand sequences
		target[i][0] = database[i][0] = '-';
		target[i][lenT[i]] = database[i][lenT[i]] = '\0';
		random_seq_gen(lenT[i], target[i], lenD[i], database[i]);
	}

/////////////////////////		OPENCL CONFIGURATION 		////////////////////////////////////

	printf("Programmin FPGA Device. \n");

   	std::string binaryFile = argv[1]; // prendo il bitstream 
    auto fileBuf = xcl::read_binary_file(binaryFile); // leggi bitstream
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

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
            OCL_CHECK(err, sw_maxi= cl::Kernel(program, "sw_maxi", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

	// Create device buffers
    std::vector<cl_mem_ext_ptr_t> lenT_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> target_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> lenD_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> database_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> score_buffer_ext(NUM_KERNEL);

    std::vector<cl::Buffer> lenT_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> target_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> lenD_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> database_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> score_buffer(NUM_KERNEL);

	for(int i = 0; i < NUM_KERNEL; i++) {

        lenT_buffer_ext[i].obj = lenT.data();
        lenT_buffer_ext[i].param = 0;
        lenT_buffer_ext[i].flags = bank[i*5];

        target_buffer_ext[i].obj = target.data();
        target_buffer_ext[i].param = 0;
        target_buffer_ext[i].flags = bank[i*5+1];
        
        lenD_buffer_ext[i].obj = lenD.data();
        lenD_buffer_ext[i].param = 0;
        lenD_buffer_ext[i].flags = bank[i*5+2];

        database_buffer_ext[i].obj = database.data();
        database_buffer_ext[i].param = 0;
        database_buffer_ext[i].flags = bank[i*5+3];

        score_buffer_ext[i].obj = score.data();
        score_buffer_ext[i].param = 0;
        score_buffer_ext[i].flags = bank[i*5+4];

    }

    for (int i = 0; i < NUM_KERNEL; i++) {
    	OCL_CHECK(err, lenT_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*INPUT_SIZE, &lenT_buffer_ext[i], &err));

        OCL_CHECK(err, target_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*INPUT_SIZE*MAX_DIM, &target_buffer_ext[i], &err));

        OCL_CHECK(err, lenD_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*INPUT_SIZE, &lenD_buffer_ext[i], &err));

        OCL_CHECK(err, database_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*INPUT_SIZE*MAX_DIM, &database_buffer_ext[i], &err));

        OCL_CHECK(err, score_buffer[i] = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*INPUT_SIZE, &score_buffer_ext[i], &err));
	}

    q.finish();

	printf("Copying input sequences on the FPGA. \n");
    // Data will be migrated to kernel space
	for(int i = 0; i < NUM_KERNEL; i++)
        err = commands.enqueueMigrateMemObjects({lenT_buffer[i], target_buffer[i], lenD_buffer[i], database_buffer[i]}, 0); /*0 means from host*/

    if (err != CL_SUCCESS) {
            printf("Error: Failed to write to device memory!\n");
            printf("Test failed\n");
            exit(1);
    }

	q.finish();

/////////////////////////		KERNEL EXCECUTION 		////////////////////////////////////
	
	// Set the arguments for kernel execution
	for (int i = 0; i < NUM_KERNEL; i++) {
		OCL_CHECK(err, err = sw_maxi.setArg(0, lenT_buffer[i]));
		OCL_CHECK(err, err = sw_maxi.setArg(1, target_buffer[i]));
		OCL_CHECK(err, err = sw_maxi.setArg(2, lenD_buffer[i]));
		OCL_CHECK(err, err = sw_maxi.setArg(3, database_buffer[i]));
		OCL_CHECK(err, err = sw_maxi.setArg(4, wd));
		OCL_CHECK(err, err = sw_maxi.setArg(5, ws));
		OCL_CHECK(err, err = sw_maxi.setArg(6, gap_opening));
		OCL_CHECK(err, err = sw_maxi.setArg(7, enlargement));
		OCL_CHECK(err, err = sw_maxi.setArg(8, score_buffer[i]));
		OCL_CHECK(err, err = sw_maxi.setArg(9, INPUT_SIZE));

		if (err != CL_SUCCESS) {
			printf("Error: Failed to set kernel arguments! %d\n", err);
			printf("Test failed\n");
			exit(1);
		}
	}

	q.finish();

	printf("Running FPGA accelerator. \n");
	auto start = std::chrono::high_resolution_clock::now();

	//Launch the Kernels
	for (int i = 0; i < NUM_KERNEL; ++i)
        err |= commands.enqueueTask(sw_maxi);


    if (err) {
        printf("Error: Failed to execute kernel! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

	q.finish();
	
	auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	float gcup = (double) (MAX_SEQ_LEN * MAX_SEQ_LEN / (float)duration.count() ) * 1e-9;
	
	//Data from Kernel to Host
	for (int i = 0; i < NUM_KERNEL; ++i) {
        err = commands.enqueueMigrateMemObjects({score_buffer[i]}, CL_MIGRATE_MEM_OBJECT_HOST);  
    }


    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

	q.finish();

	printf("Finished FPGA excecution. \n");
    printf("FPGA Kernel executed in %f ns \n", (float)duration.count());
	printf("GCUPS: %f \n",gcup);

/////////////////////////			TESTBENCH			////////////////////////////////////
	int golden_score[INPUT_SIZE];

	printf("Running Software version. \n");
	start = std::chrono::high_resolution_clock::now();

	for (int golden_rep = 0; golden_rep < INPUT_SIZE; golden_rep++) {
		golden_score[golden_rep] = compute_golden(lenT[golden_rep], target[golden_rep], lenD[golden_rep], database[golden_rep], wd, ws, gap_opening, enlargement);
	}

	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	
	printf("Software version executed in %f ns \n", (float)duration.count());

	////////test bench results
	printf("Comparing results. \n");
	for (int i=0; i<num; i++){
		if (scores[i]!=golden_scores[i]){
			printConf(target[i], database[i], ws, wd, gap_opening, enlargement);
            printf("HW: %d, SW: %d\n", scores[i], goleden_scores[i]);
            test_score=false;
        }
	}

	if (test_score) 
		printf("Test PASSED: Output matches reference.\n");
	else {
		printf("Test FAILED: Output does not match reference.\n");
	}
	

	return 0;
}

//	Prints the current configuration
void printConf(char *target, char *database, int ws, int wd, int gap_opening, int enlargement) {
	std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
	std::cout << "+ Sequence A: [" << strlen(target) << "]: " << target << std::endl;
	std::cout << "+ Sequence B: [" << strlen(database) << "]: " << database << std::endl;
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

void random_seq_gen(int lenT, char *target, int lenD, char *database) {

	int i; 
	for(i = 1; i < lenT; i++){
		int tmp_gen = gen_rnd(0, 3);
		target[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
			      (tmp_gen == 2) ? 'G' : 'T';
	}

	for(i = 1; i < lenD; i++){
		int tmp_gen = gen_rnd(0, 3);
		database[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
				  (tmp_gen == 2) ? 'G' : 'T';
	}
}

int compute_golden(int lenT, char *target, int lenD, char *database, int wd, int ws, int gap_opening, int enlargement) {
	// Inizializza le matrici D, P, Q
	    std::vector< std::vector<int> > D(lenT, std::vector<int>(lenD, 0));
	    std::vector< std::vector<int> > P(lenT, std::vector<int>(lenD, std::numeric_limits<int>::min() / 2));
	    std::vector< std::vector<int> > Q(lenT, std::vector<int>(lenD, std::numeric_limits<int>::min() / 2));

	    int gap_penalty = gap_opening + enlargement * 1;
	    int max_score = 0;

	    for (int i = 1; i < lenT; ++i) {
	        for (int j = 1; j < lenD; ++j) {
	            // Calcola P[i][j]
	            P[i][j] = std::max(P[i-1][j] + enlargement, D[i-1][j] + gap_penalty);

	            // Calcola Q[i][j]
	            Q[i][j] = std::max(Q[i][j-1] + enlargement, D[i][j-1] + gap_penalty);

	            // Calcola D[i][j]
	            int match = (target[i] == database[j]) ? wd : ws;
	            D[i][j] = std::max(0, D[i-1][j-1] + match);
	            D[i][j] = std::max(D[i][j], P[i][j]);
	            D[i][j] = std::max(D[i][j], Q[i][j]);

	            // Aggiorna lo score massimo
	            max_score = std::max(max_score, D[i][j]);
	        }
	    }

	    return max_score;
}
