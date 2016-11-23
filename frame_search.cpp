#include "peptide.cpp"

//implement main
int main() {
	try {
		InputData input;
		std::vector<std::string> dna_seq_v;
		input.rest_into_str_vector(dna_seq_v);
		std::string aa_seq = dna_seq_v[dna_seq_v.size()-1];
		dna_seq_v.pop_back();
		std::string dna_seq = join(dna_seq_v,"");
		Peptide prot = Peptide(aa_seq);
		std::cout<<"size of sequence: "<<dna_seq.size()<<std::endl;
		std::cout<<"search_sequence: "<<prot.to_string()<<std::endl<<std::flush;
		std::vector<std::string> subseqs=prot.subseqs_encoding(dna_seq);
		std::cout << join(subseqs,"\n")<< std::endl << std::flush;
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