#include "peptide.cpp"

//implement main
int main() {
	try {
		std::string whitespace("\\s");
		InputData input;
		std::size_t N = input.next_as_int_type<std::size_t>();
		std::vector<std::size_t> spectrumv;
		input.next_into_int_type_vector<std::size_t>(spectrumv, whitespace);
		Spectrum spectrum(spectrumv);
		std::cout<<"# of masses in extended list: "<<Peptide::get_aa_mass_list_extended().size()<<std::endl;
		std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
		std::vector<Peptide> leaders = Peptide::leaderboard_cyclopeptide_sequencing_all(spectrum, N, true);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::chrono::duration<double> score_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		std::cout<<"duration: "<<score_time.count()<<std::endl;
		std::cout<<Peptide::pep_container_to_mass_string(leaders.begin(), leaders.end(), " ")<<std::endl;
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
