#include "peptide.cpp"

//implement main
int main() {
	try {
		InputData input;
		std::size_t N = input.next_as_int_type<std::size_t>();
		std::cout << ((N * (N+1))/2)+1 <<std::endl;
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