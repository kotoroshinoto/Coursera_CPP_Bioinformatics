#include "peptide.cpp"
#include <cstdint>
//implement main
int main() {
	try {
		InputData input;
		std::string aa_seq = input.next();
		Peptide prot = Peptide(aa_seq);
		Spectrum sizes=prot.linear_spectrum();
		std::cout << sizes.to_string() << std::endl << std::flush;
		return 0;
	} catch(std::exception& e){
		std::cerr << "\nException occurred: " << std::endl << std::flush;
		std::cerr << e.what() << std::endl << std::flush;
		std::exit(-1);
	} catch(...) {
		std::cerr << "\nUnexpected Exception" << std::endl << std::flush;
		std::exit(-1);
	}
}