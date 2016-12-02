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
		std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
		Peptide leader = Peptide::leaderboard_cyclopeptide_sequencing(spectrum, N);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::chrono::duration<double> score_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
		std::cout<<"duration: "<<score_time.count()<<std::endl;
		std::cout<<leader.to_mass_string()<<std::endl;
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