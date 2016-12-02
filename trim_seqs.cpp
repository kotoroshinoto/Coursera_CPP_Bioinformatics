#include "peptide.cpp"

//implement main
int main() {
	try {
		std::string whitespace("\\s");
		InputData input;
		std::vector<std::string> prots;
		input.next_into_str_vector(prots, whitespace);
		std::vector<std::size_t> spectrum;
		input.next_into_int_type_vector<std::size_t>(spectrum, whitespace);
		std::size_t N = input.next_as_int_type<std::size_t>();
		std::vector<Peptide> leaders;
		std::map<Peptide, std::size_t> scores;
		for(std::size_t i=0; i<prots.size();i++){
			Peptide p = Peptide(prots[i]);
			std::size_t score = p.score_spectrum(spectrum);
			scores[p] = score;
			leaders.push_back(p);
		}
		trim_leaderboard(leaders,spectrum,scores,N);
		std::cout<<Peptide::pep_container_to_string(leaders.begin(), leaders.end(), " ")<<std::endl;
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