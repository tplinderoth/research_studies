/*
*
* replaceMaf.cpp
*
* Replaces parameter file 1 MAFs with parameter file 2 MAF values in output file.
* Output contains parameter file 1 likelihoods and inbreeding coefficients.
*
*/

#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<string>
#include<vector>
#include<sstream>

int main (int argc, char** argv) {
	if (argc < 6) {
		std::cerr << "\nreplaceMaf.cpp <nind1,nind2> <nsites> <par file 1> <par file 2> <outfile name>\n\n";
		return 0;
	}

	// parse input

	// number of individuals in each analysis
	std::stringstream s(argv[1]);
	std::string str;
	std::vector<unsigned int> nind;
	while(getline(s, str, ',')) { nind.push_back(stoi(str)); }

	// number sites - needs to be the same for both input parameter files
	unsigned int nsites = atoi(argv[2]);

	// open parameter files for reading
	const char* pf1 = argv[3];
	FILE* fpar1 = fopen(pf1, "rb");
	if (!fpar1) {
		std::cerr << "Unable to open input parameter file " << pf1 << "\n";
		return 1;
	}

	const char* pf2 = argv[4];
	FILE* fpar2 = fopen(pf2, "rb");
	if (!fpar2) {
		std::cerr << "Unable to open input parameter file " << pf2 << "\n";
		return 1;
	}

	// open parameter file stream for writing
	const char* of = argv[5];
	FILE* fparo = fopen(of, "wb");
	if (!fparo) {
		std::cerr << "Unable to open output file " << of << "\n";
		return 1;
	}
	
	// read parameter file 1
	size_t ll1sz = nind[0]+1; // includes global likelihood and all individual likelihoods
	double ll1[ll1sz];
	double f1[nind[0]];
	double maf1[nsites];
	
	fread(ll1, sizeof(double), ll1sz, fpar1);
	fread(f1, sizeof(double), nind[0], fpar1);
	fread(maf1, sizeof(double), nsites, fpar1);

	fclose(fpar1);
	
	// test print of parameter file 1
/*
	for (size_t i = 0; i < ll1sz; ++i) { 
		std::cout << ll1[i];
		if (i < ll1sz-1) std::cout << "\t"; else std::cout << "\n";
	}
	for (size_t i = 0; i < nind[0]; ++i) {
		std::cout << f1[i];
		if (i < nind[0]-1) std::cout << "\t"; else std::cout << "\n";
	}
	for (size_t i = 0; i < nsites; ++i) { std::cout << maf1[i] << "\n"; }
*/	

	// read paramter file 2

	size_t ll2sz = nind[1]+1; // includes global likeliood and all individual likelihoods
	double ll2[ll2sz];
	double f2[nind[1]];
	double maf2[nsites];

	fread(ll2, sizeof(double), ll2sz, fpar2);
	fread(f2, sizeof(double), nind[1], fpar2);
	fread(maf2, sizeof(double), nsites, fpar2);

	fclose(fpar2);

	// test print of parameter file 2
/*
	for (size_t i = 0; i < ll2sz; ++i) {
		std::cout << ll2[i];
		if (i < ll2sz-1) std::cout << "\t"; else std::cout << "\n";
	}
	for (size_t i = 0; i < nind[1]; ++i) {
		std::cout << f2[i];
		if (i < nind[1]-1) std::cout << "\t"; else std::cout << "\n";
	}
	for (size_t i = 0; i < nsites; ++i) { std::cout << maf2[i] << "\n"; }
*/

	// set par1 MAFs to par2 MAFs in output file
	size_t llcount = fwrite(ll1, sizeof(ll1[0]), ll1sz, fparo); // write par1 likelihoods
	size_t fcount = fwrite(f1, sizeof(f1[0]), nind[0], fparo); // write par1 F values
	size_t mafcount = fwrite(maf2, sizeof(maf2[0]), nsites, fparo); // write par2 MAFs
	fclose(fparo);

	std::cerr << "Wrote " << llcount << " likelihoods (includes global), " << fcount << " inbreeding coefficients, " << mafcount << "  MAF values\n";

	return 0;
}

