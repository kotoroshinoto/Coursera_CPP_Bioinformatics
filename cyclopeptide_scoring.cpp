#include "peptide.cpp"

int main() {
	try {
		InputData input;
		std::string aa_seq = input.next();
		Peptide prot = Peptide(aa_seq);
//		std::cerr << prot.to_string() <<std::endl<<std::flush;
		std::vector<std::size_t> spectrum;
		std::string whitespace("\\s");
		input.next_into_int_type_vector<std::size_t>(spectrum, whitespace);
//		std::cerr << join(spectrum, " ") <<std::endl<<std::flush;
		std::size_t score = prot.score_spectrum(spectrum, true);
		std::cout << score <<std::endl<<std::flush;
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