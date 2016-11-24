#include "peptide.cpp"

//implement main
int main() {
	try {
		std::string whitespace("\\s");
		InputData input;
		std::vector<std::size_t> spectrum;
		input.next_into_size_t_vector(spectrum, whitespace);
		std::vector<Peptide> seqs = Peptide::cyclopeptide_sequencing(spectrum);
		for(size_t i =0; i < seqs.size(); i++){
			std::vector<size_t> masses;
			for(size_t j=0; j<seqs[i].size();j++){
				masses.push_back(seqs[i].aa_mass_at(i));
			}
			if(i != 0){
				std::cout<<" ";
			}
			std::cout<<join(masses,"-");
		}
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