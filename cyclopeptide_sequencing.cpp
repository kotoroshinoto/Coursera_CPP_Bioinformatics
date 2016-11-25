#include "peptide.cpp"

//implement main
int main() {
	try {
		std::string whitespace("\\s");
		InputData input;
		std::vector<std::size_t> spectrum;
		input.next_into_int_type_vector<std::size_t>(spectrum, whitespace);
		std::vector<Peptide> seqs = Peptide::cyclopeptide_sequencing(spectrum);
		for(size_t i =0; i < seqs.size(); i++){
			if(i != 0){
				std::cout<<" ";
			}
			std::cout<<seqs[i].to_mass_string();
		}
		std::cout<<std::endl;
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