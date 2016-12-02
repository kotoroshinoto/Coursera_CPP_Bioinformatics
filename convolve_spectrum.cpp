#include "peptide.cpp"

//implement main
int main() {
	try {
		std::string whitespace("\\s");
		InputData input;
		std::vector<std::size_t> spectrumv;
		input.next_into_int_type_vector<std::size_t>(spectrumv, whitespace);
		Spectrum spectrum(spectrumv);
		Spectrum conv=spectrum.convolve();
//		std::cout<<spectrum.to_string()<<std::endl;
		std::cout<<conv.to_string()<<std::endl;
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
