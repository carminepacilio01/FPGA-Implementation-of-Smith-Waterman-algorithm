#include "host.h"

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

    cl_int 				err;
    cl::Context 		context;
    cl::Kernel 			krnl;
    cl::CommandQueue 	q;

	std::vector< int, aligned_allocator<int> > 		lenT(INPUT_SIZE);
    std::vector< int, aligned_allocator<int> > 		lenD(INPUT_SIZE);
    std::vector< int, aligned_allocator<int> > 		score(INPUT_SIZE);

	char target[INPUT_SIZE][MAX_DIM];
	char database[INPUT_SIZE][MAX_DIM];

	int cell_number;

/////////////////////////		DATASET GENERATION 		////////////////////////////////////

	std::cout << "Generating "<< INPUT_SIZE << " random sequence pairs." << std::endl;
	///////Generation of random sequences
	// Generation of random sequences
    for(int i = 0; i < INPUT_SIZE; i++){

		//	generate random sequences
		//	length of the sequences
		lenT[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);
		lenD[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);

		lenT[i] += 1;
		lenD[i] += 1;

		//	generate rand sequences
		target[i][0] = database[i][0] = '-';
		target[i][lenT[i]] = database[i][lenD[i]] = '\0';
		random_seq_gen(lenT[i], target[i], lenD[i], database[i]);
		
		cell_number += lenD[i] * lenT[i];
	}

	// reverse the string
		char t_rev[INPUT_SIZE][MAX_DIM];
		for(int i = 0; i < INPUT_SIZE; i++){
			for (int j = 0; j < lenT[i]; j++) {
					t_rev[i][j] = target[i][lenT[i] - j - 1];
				}
		}

/////////////////////////		OPENCL CONFIGURATION 		////////////////////////////////////

	std::cout << "Programmin FPGA Device. " << std::endl;;

   	std::string binaryFile = argv[1]; // prendo il bitstream 
    auto fileBuf = xcl::read_binary_file(binaryFile); // leggi bitstream
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

    auto devices = xcl::get_xil_devices(); // lista di devices

    bool valid_device = false;

    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, krnl= cl::Kernel(program, "sw_maxi", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

	// Create device buffers
    cl::Buffer lenT_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int)*INPUT_SIZE, lenT.data());
    cl::Buffer target_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(char)*INPUT_SIZE*MAX_DIM, t_rev);
    cl::Buffer lenD_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int)*INPUT_SIZE, lenD.data());
	cl::Buffer database_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(char)*INPUT_SIZE*MAX_DIM, database);
    cl::Buffer score_buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int)*INPUT_SIZE, score.data());

	std::cout <<("Copying input sequences on the FPGA. \n");

    // Data will be migrated to kernel space
    err = q.enqueueMigrateMemObjects({lenT_buffer, target_buffer, lenD_buffer, database_buffer}, 0); /*0 means from host*/
	q.finish();

/////////////////////////		KERNEL EXCECUTION 		////////////////////////////////////

	// Set the arguments for kernel execution
	OCL_CHECK(err, err = krnl.setArg(0, lenT_buffer));
	OCL_CHECK(err, err = krnl.setArg(1, target_buffer));
	OCL_CHECK(err, err = krnl.setArg(2, lenD_buffer));
	OCL_CHECK(err, err = krnl.setArg(3, database_buffer));
	OCL_CHECK(err, err = krnl.setArg(4, wd));
	OCL_CHECK(err, err = krnl.setArg(5, ws));
	OCL_CHECK(err, err = krnl.setArg(6, gap_opening));
	OCL_CHECK(err, err = krnl.setArg(7, enlargement));
	OCL_CHECK(err, err = krnl.setArg(8, score_buffer));

	if (err != CL_SUCCESS) {
		std::cout << "Error: Failed to set kernel arguments! " << err << std::endl;
		std::cout << "Test failed" << std::endl;
		return EXIT_FAILURE;;
	}

	std::cout <<("Running FPGA accelerator. \n");
	
	//Launch the Kernels
	auto start = std::chrono::high_resolution_clock::now();
	for(int ker = 0; ker < NUM_KER; ker++) {
		OCL_CHECK(err, err = krnl.setArg(9, ker * INPUT_SIZE / NUM_KER));
		OCL_CHECK(err, err = krnl.setArg(10, INPUT_SIZE / NUM_KER));
         q.enqueueTask(krnl);
     }
	q.finish();
	auto stop = std::chrono::high_resolution_clock::now();

	if (err != CL_SUCCESS) {
		std::cout << "Error: Failed to launch kernel(s)! " << err << std::endl;
		std::cout << "Test failed" << std::endl;
		return EXIT_FAILURE;;
	}
	
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	float gcup = (double) (cell_number / (float)duration.count());
	
	//Data from Kernel to Host
    err = q.enqueueMigrateMemObjects({score_buffer}, CL_MIGRATE_MEM_OBJECT_HOST);  
	if (err != CL_SUCCESS) {
		std::cout << "Error: Failed to retrive objects from kernel(s)!" << err << std::endl;
		std::cout << "Test failed" << std::endl;
		return EXIT_FAILURE;;
	}
	q.finish();

	std::cout << "Finished FPGA excecution." << std::endl;
    std::cout << "FPGA Kernel executed in " << (float)duration.count() * 1e-6 << "ms" << std::endl;
	std::cout << "GCUPS: " << gcup << std::endl;

/////////////////////////			TESTBENCH			////////////////////////////////////
	int golden_score[INPUT_SIZE];

	std::cout << "Running Software version." << std::endl;;
	start = std::chrono::high_resolution_clock::now();

	for (int golden_rep = 0; golden_rep < INPUT_SIZE; golden_rep++) {
		golden_score[golden_rep] = compute_golden(lenT[golden_rep], target[golden_rep], lenD[golden_rep], database[golden_rep], wd, ws, gap_opening, enlargement);
	}

	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	
	std::cout << "Software version executed in " <<  (float)duration.count() * 1e-6 << " ms " << std::endl;

	////////test bench results
	std::cout << "Comparing results. \n" << std::endl;
	bool test_score=true;
	for (int i=0; i < INPUT_SIZE; i++){
		if (score[i]!=golden_score[i]){
			printConf(target[i], database[i], ws, wd, gap_opening, enlargement);
            std::cout << "HW: "<< score[i] << ", SW: " << golden_score[i] << std::endl;
            test_score=false;
        }
	}

	if (test_score) 
		std::cout << "Test PASSED: Output matches reference." << std::endl;
	else {
		std::cout << "Test FAILED: Output does not match reference." << std::endl;
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

	char alphabet[4] = {'A', 'C', 'G', 'T'};
	int i;
	for(i = 1; i < lenT; i++){
		int tmp_gen = gen_rnd(0, 3);
		target[i] = alphabet[tmp_gen];
	}

	for(i = 1; i < lenD; i++){
		int tmp_gen = gen_rnd(0, 3);
		database[i] = alphabet[tmp_gen];
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
