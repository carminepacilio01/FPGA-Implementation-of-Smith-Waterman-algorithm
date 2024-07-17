# FPGA implementation of the Smith-Waterman-Gotoh algorithm for local DNA sequence alignment
In the rapidly evolving fields of bioinformatics and genomics, efficient sequence alignment is a crucial element for a wide range of molecular biology applications. As technologies continue to advance, the demand for faster and more effective algorithms to manage the escalating data volume is becoming increasingly evident. Current sequence alignment algorithms exhibit temporal and spatial complexities that are, at best, quadratic, rendering them inadequate for aligning large sequences. 
This repository showcases the source code for my project "Accelerated Smith-Waterman algorithm on FPGA" for PII at Politecnico di Milano. The project introduces an accelerated version of the affine gap-penalty Smith-Waterman algorithm for local sequence alignment. 
To accelerate this algorithm, I used Vitis HLS to write the top function and then used OpenCL and Vitis compiler to build the host and kernel for FPGA usage.

## Repository Structure
* `Build` Contains pre built binaries to replicate the experiment.
* `Design` Contains the source code for the c++ top function.
* `HostCode` Contains the host implementation.
* `TestBench` Contains the vitis HLS testbench.

## Running the application on FPGA
To run this implementation on FPGA, first clone the repository using:

```bash
git clone https://github.com/carminepacilio01/Accelerated-Smith-Waterman-algorithm-on-FPGA.git
```
Once cloned, source the necessary tools:
```bash
source /opt/xilinx/xrt/setup.sh
source /xilinx/software/Vitis/<Vitis_Version>/settings64.sh
```
Then, open the terminal in the `Build/<Build_Number>` folder and run the following:
```bash
./sw_maxi ./sw_maxi.xclbin
```
