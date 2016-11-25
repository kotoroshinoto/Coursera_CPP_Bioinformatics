#include "peptide.cpp"

std::vector<std::size_t> masses(Peptide::get_aa_mass_list());

void test_print(std::size_t n){
	std::cout<<"test_print: "<<n<<std::endl<<std::flush;
}
std::size_t get_count_for_compomer(std::vector<std::size_t> &c){
	
}
std::size_t get_peps_with_mass_recur(std::size_t i, std::size_t m, std::vector<std::size_t> &c){


	return 0;
}

std::size_t count_peps_with_mass(std::size_t m){
	std::sort(masses.begin(),masses.end());
	std::vector<std::size_t> c;
	for(std::size_t i=0;i<masses.size();i++){
		c.push_back(0);
	}
	return get_peps_with_mass_recur(masses.size(),m,c);
}

//implement main
int main() {
	try {
		InputData input;
		std::size_t N = input.next_as_int_type<std::size_t>();
		std::cout<<count_peps_with_mass(N)<<std::endl;
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