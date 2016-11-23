#include "peptide.cpp"

//implement main
int main() {
	try {
		InputData input;
		std::string prot_seq = input.next();
		Peptide prot(prot_seq);
		std::cout << prot.number_of_sequences()<< std::endl << std::flush;
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