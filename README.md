//Introduction
This is a tool for compressing collection of FASTA format genomes, 
which is based on reference sequence and uses the Suffix Array(SA) and 
the longest common prefix(LCP) to find the longest matched substrings（LMS）for the first-level matching 
and then uses a dynamic hash table for the second-level matching to complete the final compression. 
The program uses GPUs to parallelize the construction of SA, which is proposed in[1]. 
This package contains two programs, the TLMRGC is used for compression and the TLMRGCDE is used for decompression. 
The compression and decompression codes are stored in two run.cpp files, 
which are located in the directory TLMRGC(andTLMRGCDE)/src/run.cpp. 
The source code for building SA using GPUs is stored in the directory "TLMRGC/src/multi-gpu-suffix-array/src/". 
Meanwhile， the Verification code is also provided and stored in the directory “VERIFI/src/run.app”.
 
//Server Environment
The program is written in C++ and tested on Red Hat Enterprise Linux 7.9(64-bit), 
and with 2 RTX6000 GPUs with 24GB of RAM, and 2 * 2.6 GHz Intel Xeon Gold 6240 CPUs (18 cores) with 256GB RAM. 
This program requires Cmake 3.8 or greater but the version of GCC is no higher than 5.4.0. 
Meanwhile, this program also requires CUDA 10.0 (or 9.1) is installed in the root directory of the server.

Notice:
	The use of the original algorithm[1] for constructing SA using multiple GPUs needs to pay attention to the communication mode between GPUs, 
cuda10.0 is used when using NV-link to link GPUs for communication, and cuda9.1 is used when using PCIe to link GPUs for communication, 
and appropriate adjustments need to be made according to different communication modes during compiling the software. 

//Compile
During the compilation process, first locate the "Makefile" file, which is stored in the directory TLMRGC (or TLMRGCDE)/ ,
 and then use the command "make" to compile the source codes.
Compile command for compressor:
	cd TLMRGC
       	make
Compile command for decompressor:
	cd TLMRGCDE
       	make
Compile command for Verification:
cd VERIFI
make
//Compress , decompress and verification commands
Compress commands:
	(compress)   ./TLMRGC R-folder T-folder First-level folder R-T-file
R-folder is the reference genome;T-folder is the target genome; First-level folder is thte information sequence; R-T-file is the compressed result.
Example:    (./TLMRGC  hg17  Collection genomes  First-level  Compressed folder)
 Hg17 is the reference genome; Collection genomes is the set of the target genome; 
First-level is a folder store the information sequence which is the result of the First-level matching. 
Compressed folder is the compressed result.
 The default name of chromosomes is [chr1.fa, chr2.fa, chr3.fa, chr4.fa, chr5.fa, chr6.fa, chr7.fa, chr8.fa, chr9.fa, chr10.fa, chr11.fa, chr12.fa, 
chr13.fa, chr14.fa, chr15.fa, chr16.fa, chr17.fa, chr18.fa, chr19.fa, chr20.fa, chr21.fa, chr22.fa, chrX.fa, chrY.fa]

Decompress commands:
(decompress)  ./R-folder  Compressed-file  Second-level-folder  First-level-folder   De-file 
R-folder is the reference genome; Compressed-file is the compression result; 
Second-level-folder is the tempporary folder to store the decoded second-level matching information of the compression result; 
First-level-folderis the tempporary folder to store the decoded first-level matching information of the compression result; 
De-file is the decoded file;
Example:     (./TLMRGCDE   hg17   hg17-compressed  Second-level-folder First-level-folder   Decode)
hg17 is the reference genome; hg17--compressed is the compression result; 
Second-level-folder is the folder used to store the decoded information of the hg17-compressed; 
First-level-folder is the folder used to store the decoded information of First-level-folder; Decode is  the folder used to store the decoded result; 
Verification commands
(Verification)  ./ Decode genome    Verification
Example:     (./VERIFI   De-hg17   hg17)
De-hg17 is the result of the decoding, hg17 is the compressed genome;
The default name of the chromosome is the same as those used for compression.
