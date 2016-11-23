#include <iostream>
#include <queue>
#include <stack>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <regex>
#include <iterator>
#include <type_traits>
#include <cctype>

// define static functions
std::size_t compute_original_index(std::size_t aa_pos, std::size_t orig_len, std::size_t offset, std::size_t search_aa_len, bool is_reverse);

// define Objects
class intVectorFiller {
    std::vector<int>& v;
public:
    intVectorFiller(std::vector<int>& v);
    void operator()(std::string& item);
};

class InputData {
private:
    // regexes
    static std::regex whitespace_trimmer_regex_lead_and_trail;
    static std::regex whitespace_trimmer_regex_lead_only;
    static std::regex whitespace_trimmer_regex_trail_only;
    static std::regex whitespace_re;
    std::queue<std::string> data;
public:
    InputData();
    int next_as_int();
    void next_into_int_vector(std::vector<int> &intv, std::string &separator);
    void next_into_str_vector(std::vector<std::string> &strv, std::string &separator);
    void rest_into_int_vector(std::vector<int> &intv);
    void rest_into_str_vector(std::vector<std::string> &strv);
    std::string next();
    std::size_t size();
    std::string StripWhitespace(std::string& text);
    std::string StripLeadingWhitespace(std::string& text);
    std::string StripTrailingWhitespace(std::string& text);
	static int stoi(std::string &text);
};

//instantiate static variables
std::regex InputData::whitespace_trimmer_regex_lead_and_trail("\\s+$|^\\s+");
std::regex InputData::whitespace_trimmer_regex_lead_only("^\\s*");
std::regex InputData::whitespace_trimmer_regex_trail_only("\\s*$");
std::regex InputData::whitespace_re("\\s+");

//implement object functions

intVectorFiller::intVectorFiller(std::vector<int>& v): v(v){}

void intVectorFiller::operator()(std::string& item) {
	v.push_back(InputData::stoi(item));
}

int InputData::stoi(std::string &text){
	if(text.find(',') >= 0) {
		std::string t(text);
		std::cout << "text: " << t <<  std::flush;
		t.erase(std::remove(t.begin(), t.end(), ','), t.end());
		std::cout << "\tt: " << t << std::flush;
		int i = std::stoi(t);
		std::cout << "\tint: " << i << std::endl << std::flush;
		return i;
	} else {
		return std::stoi(text);
	}
}

InputData::InputData() {
    std::string inputline;
    while(std::getline(std::cin, inputline)){
        inputline = StripWhitespace(inputline);
        this->data.push(inputline);
    }
}

int InputData::next_as_int(){
    int next = InputData::stoi(this->data.front());
    this->data.pop();
    return next;
}

std::string InputData::next(){
    std::string next(this->data.front().c_str());
    this->data.pop();
    return next;
}

std::size_t InputData::size(){
    return this->data.size();
}

void InputData::next_into_int_vector(std::vector<int> &intv, std::string &separator) {
    std::vector<std::string> str_form;
    this->next_into_str_vector(str_form, separator);
    std::for_each(str_form.begin(),str_form.end(),intVectorFiller(intv));
}

void InputData::next_into_str_vector(std::vector<std::string> &strv, std::string &separator) {
    std::regex re(separator);
    std::string &next = this->data.front();
    std::sregex_token_iterator first{next.begin(), next.end(), re, -1}, last;
    std::copy(first, last, std::back_inserter<std::vector<std::string> >(strv));
    this->data.pop();
}

void InputData::rest_into_int_vector(std::vector<int> &intv) {
    std::vector<std::string> str_form;
    this->rest_into_str_vector(str_form);
    std::for_each(str_form.begin(),str_form.end(),intVectorFiller(intv));
}

void InputData::rest_into_str_vector(std::vector<std::string> &strv) {
    while(this->data.size() > 0){
        strv.push_back(this->next());
    }
}

std::string InputData::StripWhitespace(std::string& text){
    return std::regex_replace(text, whitespace_trimmer_regex_lead_and_trail, "", std::regex_constants::match_any | std::regex_constants::format_sed);
}

std::string InputData::StripLeadingWhitespace(std::string& text){
    return std::regex_replace(text, whitespace_trimmer_regex_lead_only, "", std::regex_constants::match_any | std::regex_constants::format_sed);
}

std::string InputData::StripTrailingWhitespace(std::string& text){
    return std::regex_replace(text, whitespace_trimmer_regex_trail_only, "", std::regex_constants::match_any | std::regex_constants::format_sed);
}
//template <typename T, template<typename...> class C>
template <typename T, template<typename...> class Container>
std::string join(Container<T> &c, std::string s){
	if(c.size() == 0){
		return "";
	}
	std::stringstream ss;
	std::copy(c.begin(), c.end()-1, std::ostream_iterator<T>(ss, s.c_str()));
	ss << c[c.size() -1];
	return ss.str();
}

class Peptide {
public:
	enum PeptideCode {
		STOP,
		HISTIDINE,
		GLUTAMINE,
		PROLINE,
		ARGININE,
		LEUCINE,
		ASPARTIC_ACID,
		GLUTAMIC_ACID,
		ALANINE,
		GLYCINE,
		VALINE,
		TYROSINE,
		SERINE,
		CYSTEINE,
		TRYPTOPHAN,
		PHENYLALANINE,
		ASPARAGINE,
		LYSINE,
		THREONINE,
		ISOLEUCINE,
		METHIONINE
	};
	static Peptide translateDNA(std::string& dnaseq);
	static std::string transcribeDNA(std::string& dnaseq);
	static Peptide translateRNA(std::string& rnaseq);
	static std::string PeptideCode_to_string(PeptideCode aa);
	static std::string PeptideCode_to_abbrev(PeptideCode aa);
	static char PeptideCode_to_symbol(PeptideCode aa);
	static PeptideCode codon_to_peptide(char cod1, char cod2, char cod3);
	static int number_of_codons(PeptideCode aa);
	static PeptideCode symbol_to_peptide(char symb);
	static std::string revcomp(std::string& dna, bool rna_mode=false);
	static char complement(char nuc, bool rna_mode=false);
	Peptide(std::string aastr);
	Peptide();
	Peptide(PeptideCode& aa);
	Peptide(std::vector<PeptideCode>& aa);
	void addAminoAcid(PeptideCode aa);
	std::string to_string(bool include_stops=false);
	std::string to_abbrev_string(bool include_stops=false);
	std::string to_fullword_string(bool include_stops=false);
	std::size_t number_of_sequences(bool include_stops=false);
	std::vector<std::string> subseqs_encoding(std::string& dna);
	std::size_t size(bool include_stops=false);
	std::size_t nuc_size(bool include_stops=false);
	bool operator ==(const Peptide &prot2);
	bool isSameAs(const Peptide &prot2, bool include_stops=false);
	Peptide subseq(std::size_t i, std::size_t len);
//	Peptide addAminoAcid_new(PeptideCode aa);
//	std::vector<Peptide> branch();
//	bool bind();
private:
//	int totalmass;
//	std::vector<int> masses;
	std::vector<PeptideCode> sequence;
};

Peptide::Peptide():sequence(){}

Peptide::Peptide(PeptideCode& aa):sequence(){
	this->addAminoAcid(aa);
}

Peptide::Peptide(std::vector<PeptideCode>& aaseq):sequence(){
	std::for_each(aaseq.begin(), aaseq.end(),std::bind1st(std::mem_fun(&Peptide::addAminoAcid), this));
}


std::string Peptide::to_string(bool include_stops){
	std::vector<char> aastr;
	for(std::size_t i = 0; i < this->sequence.size() ;i++){
		Peptide::PeptideCode current = this->sequence[i];
		if((! include_stops) && current == Peptide::PeptideCode::STOP){
			break;
		}
		aastr.push_back(Peptide::PeptideCode_to_symbol(current));
	}
	return join(aastr, "");
}

std::string Peptide::to_abbrev_string(bool include_stops){
	std::vector<std::string> abbrs;
	for(std::size_t i = 0; i < this->sequence.size() ;i++){
		Peptide::PeptideCode current = this->sequence[i];
		if((! include_stops) && current == Peptide::PeptideCode::STOP){
			break;
		}
		abbrs.push_back(Peptide::PeptideCode_to_abbrev(current));
	}
	return join(abbrs, "");
}

std::string Peptide::to_fullword_string(bool include_stops){
	std::vector<std::string> words;
	for(std::size_t i = 0; i < this->sequence.size() ;i++){
		Peptide::PeptideCode current = this->sequence[i];
		if((! include_stops) && current == Peptide::PeptideCode::STOP){
			break;
		}
		words.push_back(Peptide::PeptideCode_to_string(current));
	}
	return join(words, "");
}

Peptide::Peptide(std::string aastr):sequence(){
	for(std::size_t i =0;i < aastr.size();i++){
		this->addAminoAcid(Peptide::symbol_to_peptide(aastr[i]));
	}
//	std::cout<<"number of peptide entries: "<<this->sequence.size()<<std::endl<<std::flush;
}

std::size_t Peptide::number_of_sequences(bool include_stops){
	std::size_t total = 1;
	int num_codons;
	for(std::size_t i =0;i < this->sequence.size();i++){
		Peptide::PeptideCode current = this->sequence[i];
		if((! include_stops) && current == Peptide::PeptideCode::STOP){
			break;
		}
		num_codons = Peptide::number_of_codons(current);
//		std::cout<< Peptide::PeptideCode_to_string(current) << "\t" <<Peptide::PeptideCode_to_symbol(current) << "\t" << Peptide::number_of_codons(current)<<std::endl<<std::flush;
		total *= num_codons;
	}
	return total;
}

void Peptide::addAminoAcid(PeptideCode aa){
	this->sequence.push_back(aa);
}


Peptide Peptide::translateDNA(std::string& dnaseq){
//	std::cout<<"[Peptide::translateDNA] before: '"<<dnaseq<<"'"<<std::endl<<std::flush;
	std::string rnaseq(Peptide::transcribeDNA(dnaseq));
//	std::cout<<"[Peptide::translateDNA] after: '"<<rnaseq<<"'"<<std::endl<<std::flush;
	return Peptide::translateRNA(rnaseq);
}

std::string Peptide::transcribeDNA(std::string& dnaseq){
	std::string rnaseq(dnaseq);
	std::replace(rnaseq.begin(),rnaseq.end(),'T','U');
	return rnaseq;
}

Peptide Peptide::translateRNA(std::string& rnaseq){
	std::vector<PeptideCode> aaseq;
	for(std::size_t i = 0; i+2 < rnaseq.size() ; i+=3){
		aaseq.push_back(Peptide::codon_to_peptide(rnaseq[i],rnaseq[i+1],rnaseq[i+2]));
	}
	Peptide prot(aaseq);
	return prot;
}

bool Peptide::isSameAs(const Peptide &prot2, bool include_stops){
	if(this->sequence.size() != prot2.sequence.size()){
		return false;
	}
	for(std::size_t i = 0; i < this->sequence.size() ;i++){
		if(this->sequence[i] != prot2.sequence[i]){
			return false;
		}
		if((!include_stops) && this->sequence[i] == Peptide::PeptideCode::STOP){
			break;
		}
	}
	return true;
}

bool Peptide::operator ==(const Peptide &prot2){
	return this->isSameAs(prot2, true);
}

Peptide Peptide::subseq(std::size_t i, std::size_t len){
	Peptide newpep;
	for(std::size_t j = i; j - i < len && j < this->sequence.size(); j++){
		newpep.addAminoAcid(this->sequence[j]);
	}
	return newpep;
}

std::size_t compute_original_index(std::size_t aa_pos, std::size_t orig_len, std::size_t offset, std::size_t search_aa_len, bool is_reverse){
	//MA len 2 A is pos[1]
	//ATGGCC GCC is pos[3]
//	std::cout<<"\taa pos: "<<aa_pos;
//	std::cout<<"\toriginal len: "<<orig_len;
	std::size_t nuc_pos = aa_pos * 3;
//	std::cout<<"\tnuc_pos: "<<nuc_pos;
	if(is_reverse){

		std::size_t rev_index = nuc_pos + offset; // position in original reversed string before bases were truncated away
		std::size_t rev_index_end = rev_index + (search_aa_len*3) - 1; //0-3 in 4 lenstr 0 + 4 -1
//		std::cout<<"\trev_index: "<<rev_index;
//		std::cout<<"\trev_index_end: "<<rev_index_end;
		// len 6
		// ACTCCG
		//     ^(4)
		// CGGAGT
		//  ^(1) 6- 4 -1  len - oldindex - 1 = revindex
		//  GGAGT
		//  ^(0) new index is len - revindex - 1 - offset

		//len - oldindex - 1 = revindex
		//len - oldindex = revindex + 1
		//len - revindex =  + 1 + oldindex
		//len - revindex - 1 =  oldindex
//		std::cout<<std::endl;
		return orig_len - rev_index_end - 1; //originally 6 long, position 4 in original is 1 in reverse, and 0 in offset of 1
	} else{
		//ATCGGA T is [1]
		// TCGGA T is [0]
//		std::cout<<std::endl;
		return nuc_pos + offset;
	}
}

std::vector<std::string> Peptide::subseqs_encoding(std::string& dna){
//	std::cout<<"[Peptide::subseqs_encoding] start"<<std::endl<<std::flush;
	std::vector<std::string> strv;
	if(this->sequence.size() * 3 > dna.size()){
//		std::cout<<"sequence too short "<<(this->sequence.size()*3)<<" > "<<dna.size()<<std::endl<<std::flush;
		return strv;
	}
	std::string rdna = Peptide::revcomp(dna);
	if(this->sequence.size() * 3 == dna.size()){
		//only one frame is possible
		Peptide aastr = Peptide::translateDNA(dna);
		Peptide aarstr = Peptide::translateDNA(rdna);
		if(this->operator==(aastr)){
			strv.push_back(dna);
		}
		if(this->operator==(aarstr)){
			strv.push_back(rdna);
		}
	}
	std::vector<std::string> frames;
	std::vector<Peptide> transframes;
	std::vector<std::size_t> coords;
	std::vector<bool> is_reverse;
	//extra: number that have to be removed from end to translate
	//only add frames that CAN contain the codon sequence
	for(std::size_t i =0; i<3 && this->sequence.size() * 3 <= dna.size();i++){
//		std::cout<<"[Peptide::subseqs_encoding] building frames i: "<<i<<std::endl<<std::flush;
		std::string subdna = dna.substr(i,dna.size()-i); //ATCG -> TGC
		std::string subrdna = rdna.substr(i,dna.size()-i); // ATGC -> CGTA -> GCAT -> CAT
//		std::cout<<std::string(i,' ')<<subdna<<std::endl;
//		std::cout<<std::string(i,' ')<<subrdna<<std::endl<<std::flush;
		frames.push_back(subdna);
		frames.push_back(subrdna);

		transframes.push_back(Peptide::translateDNA(subdna));
		coords.push_back(i);
		is_reverse.push_back(false);
		transframes.push_back(Peptide::translateDNA(subrdna));
		coords.push_back(i);
		is_reverse.push_back(true);
	}
//	std::cout<<"[Peptide::subseqs_encoding] frames.size(): "<<frames.size()<<"\ttransframes.size(): "<<transframes.size()<<std::endl<<std::flush;

	for(std::size_t i = 0; i < transframes.size() ; i++){
//		std::stringstream ss2;
//		{
//			std::stringstream ss1;
//			ss1<<"\toffset: "<<coords[i]<<"\tis_reverse: "<<is_reverse[i];
//		 	std::cout<<transframes[i].to_string(true)<< ss1.str()<<" len: "<<transframes[i].size(true)<<std::endl;
//		}
		std::string pointers(dna.size(), ' ');
		std::size_t offset = coords[i];
		for(std::size_t j = 0;this->size(true) + j <= transframes[i].size(true) && j < transframes[i].size(true) ; j++) {
			//make subseqs of same size as self
			Peptide subseq = transframes[i].subseq(j,this->sequence.size());
			std::size_t substr_pos = j * 3;
			std::size_t orig_index = compute_original_index(j,dna.size(),offset, this->size(true),is_reverse[i]);

			if(this->isSameAs(subseq)){
//				ss2<<"j: "<<j;
//				ss2<<" offset: "<<offset<<"\tsubstr_pos:"<<substr_pos<<"\torig_index: "<<orig_index<<std::endl<<std::flush;
				strv.push_back(dna.substr(orig_index,this->nuc_size(true)));
				pointers[orig_index] = '^';
			}
		}
//		std::cout<<dna<<std::endl<<pointers<<std::endl<<ss2.str()<<std::endl<<std::flush;
	}
//	std::cout<<"[Peptide::subseqs_encoding] end"<<std::endl<<std::flush;
	return strv;
}

std::size_t Peptide::size(bool include_stops){
	if(include_stops){
		return this->sequence.size();

	}else {
		return this->to_string(include_stops).size();
	}
}

std::size_t Peptide::nuc_size(bool include_stops){
	if(include_stops){
		return this->sequence.size() * 3;

	}else {
		return this->to_string(include_stops).size() * 3;
	}
}
char Peptide::complement(char nuc, bool rna_mode){
	switch(nuc){
		default:
			throw std::runtime_error("[Peptide::complement] invalid nucleotide");
		case 'A':
			if(rna_mode){
				return 'U';
			}
			return 'T';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'T':
			if(rna_mode){
				throw std::runtime_error("[Peptide::complement] DNA character in RNA: ");
			}
			return 'A';
		case 'U':
			if(!rna_mode){
				throw std::runtime_error("[Peptide::complement] RNA character in DNA: ");
			}
			return 'A';
	}
}

std::string Peptide::revcomp(std::string& dna, bool rna_mode){
	std::string rev(dna);
	std::reverse(rev.begin(),rev.end());
	std::vector<char> revcomp;
	for(std::size_t i = 0; i<rev.size() ; i++){
		revcomp.push_back(Peptide::complement(rev[i],rna_mode));
	}
	return join(revcomp, "");
}

Peptide::PeptideCode Peptide::symbol_to_peptide(char symb){
//	int chari(symb);
//	std::cout<<"\naa symbol ascii value: "<<chari<<std::endl<<std::flush;
	switch ((char)std::toupper(symb)){
		default:
			throw std::runtime_error("[Peptide::symbol_to_peptide] invalid amino acid symbol");
		case '*':
			return PeptideCode::STOP;
		case 'H':
			return PeptideCode::HISTIDINE;
		case 'Q':
			return PeptideCode::GLUTAMINE;
		case 'P':
			return PeptideCode::PROLINE;
		case 'R':
			return PeptideCode::ARGININE;
		case 'L':
			return PeptideCode::LEUCINE;
		case 'D':
			return PeptideCode::ASPARTIC_ACID;
		case 'E':
			return PeptideCode::GLUTAMIC_ACID;
		case 'A':
			return PeptideCode::ALANINE;
		case 'G':
			return PeptideCode::GLYCINE;
		case 'V':
			return PeptideCode::VALINE;
		case 'Y':
			return PeptideCode::TYROSINE;
		case 'S':
			return PeptideCode::SERINE;
		case 'C':
			return PeptideCode::CYSTEINE;
		case 'W':
			return PeptideCode::TRYPTOPHAN;
		case 'F':
			return PeptideCode::PHENYLALANINE;
		case 'N':
			return PeptideCode::ASPARAGINE;
		case 'K':
			return PeptideCode::LYSINE;
		case 'T':
			return PeptideCode::THREONINE;
		case 'I':
			return PeptideCode::ISOLEUCINE;
		case 'M':
			return PeptideCode::METHIONINE;
	}
}

int Peptide::number_of_codons(Peptide::PeptideCode aa){
	switch (aa){
		case ISOLEUCINE:
		case STOP:
			return 3;
		case HISTIDINE:
		case GLUTAMINE:
		case ASPARTIC_ACID:
		case GLUTAMIC_ACID:
		case TYROSINE:
		case CYSTEINE:
		case PHENYLALANINE:
		case ASPARAGINE:
		case LYSINE:
			return 2;
		case ALANINE:
		case PROLINE:
		case GLYCINE:
		case VALINE:
		case THREONINE:
			return 4;
		case ARGININE:
		case LEUCINE:
		case SERINE:
			return 6;
		case TRYPTOPHAN:
		case METHIONINE:
			return 1;
	}
}

std::string Peptide::PeptideCode_to_string(Peptide::PeptideCode aa){
	switch (aa){
		case STOP:
			return "STOP";
		case HISTIDINE:
			return "Histidine";
		case GLUTAMINE:
			return "Glutamine";
		case PROLINE:
			return "Proline";
		case ARGININE:
			return "Arginine";
		case LEUCINE:
			return "Leucine";
		case ASPARTIC_ACID:
			return "Aspartic acid";
		case GLUTAMIC_ACID:
			return "Glutamic acid";
		case ALANINE:
			return "Alanine";
		case GLYCINE:
			return "Glycine";
		case VALINE:
			return "Valine";
		case TYROSINE:
			return "Tyrosine";
		case SERINE:
			return "Serine";
		case CYSTEINE:
			return "Cysteine";
		case TRYPTOPHAN:
			return "Tryptophan";
		case PHENYLALANINE:
			return "Phenylalanine";
		case ASPARAGINE:
			return "Asparagine";
		case LYSINE:
			return "Lysine";
		case THREONINE:
			return "Threonine";
		case ISOLEUCINE:
			return "Isoleucine";
		case METHIONINE:
			return "Methionine";
	}
}
std::string Peptide::PeptideCode_to_abbrev(PeptideCode aa){
	switch (aa){
		case STOP:
			return "STP";
		case HISTIDINE:
			return "His";
		case GLUTAMINE:
			return "Glu";
		case PROLINE:
			return "Pro";
		case ARGININE:
			return "Arg";
		case LEUCINE:
			return "Leu";
		case ASPARTIC_ACID:
			return "Asp";
		case GLUTAMIC_ACID:
			return "Glu";
		case ALANINE:
			return "Ala";
		case GLYCINE:
			return "Gly";
		case VALINE:
			return "Val";
		case TYROSINE:
			return "Tyr";
		case SERINE:
			return "Ser";
		case CYSTEINE:
			return "Cys";
		case TRYPTOPHAN:
			return "Trp";
		case PHENYLALANINE:
			return "Phe";
		case ASPARAGINE:
			return "Asn";
		case LYSINE:
			return "Lys";
		case THREONINE:
			return "Thr";
		case ISOLEUCINE:
			return "Ile";
		case METHIONINE:
			return "Met";
	}
}
char Peptide::PeptideCode_to_symbol(PeptideCode aa){
	switch (aa){
		case STOP:
			return '*';
		case HISTIDINE:
			return 'H';
		case GLUTAMINE:
			return 'Q';
		case PROLINE:
			return 'P';
		case ARGININE:
			return 'R';
		case LEUCINE:
			return 'L';
		case ASPARTIC_ACID:
			return 'D';
		case GLUTAMIC_ACID:
			return 'E';
		case ALANINE:
			return 'A';
		case GLYCINE:
			return 'G';
		case VALINE:
			return 'V';
		case TYROSINE:
			return 'Y';
		case SERINE:
			return 'S';
		case CYSTEINE:
			return 'C';
		case TRYPTOPHAN:
			return 'W';
		case PHENYLALANINE:
			return 'F';
		case ASPARAGINE:
			return 'N';
		case LYSINE:
			return 'K';
		case THREONINE:
			return 'T';
		case ISOLEUCINE:
			return 'I';
		case METHIONINE:
			return 'M';
	}
}

Peptide::PeptideCode Peptide::codon_to_peptide(char cod1, char cod2, char cod3) {
	char N1 = (char) std::toupper(cod1);
	char N2 = (char) std::toupper(cod2);
	char N3 = (char) std::toupper(cod3);
//	std::cout<<N1<<N2<<N3<<std::endl<<std::flush;
	switch(std::toupper(N1)) {
		case 'A': // A--
			switch(N2) {
				case 'A': //AA-
					switch(N3) {
						case 'C': //AAC
						case 'U': //AAU
							return PeptideCode::ASPARAGINE;
						case 'A': //AAA
						case 'G': //AAG
							return PeptideCode::LYSINE;
						default:
							break;
					}
				case 'U': //AU-
					switch(N3) {
						case 'G': // AUG
							return PeptideCode::METHIONINE;
						case 'A': // AUA
						case 'U': // AUU
						case 'C': // AUC
							return PeptideCode::ISOLEUCINE;
						default:
							break;
					}
				case 'G': // AG-
					switch(N3) {
						case 'A': // AGA
						case 'G': // AGG
							return PeptideCode::ARGININE;
						case 'U': // AGU
						case 'C': // AGC
							return PeptideCode::SERINE;
						default:
							break;
					}
				case 'C': // AC-
					return PeptideCode::THREONINE;
				default:
					break;
			}
		case 'U': //U--
			switch(N2) {
				case 'A': //UA-
					switch(N3) {
						case 'A':
						case 'G':
							return PeptideCode::STOP;
						case 'U':
						case 'C':
							return PeptideCode::TYROSINE;
						default:
							break;
					}
				case 'U': //UU-
					switch(N3) {
						case 'A':
						case 'G':
							return PeptideCode::LEUCINE;
						case 'U':
						case 'C':
							return PeptideCode::PHENYLALANINE;
						default:
							break;
					}
				case 'G'://UG-
					switch(N3) {
						case 'A': // UGA
							return PeptideCode::STOP;
						case 'G': // UGG
							return PeptideCode::TRYPTOPHAN;
						case 'U': // UGU
						case 'C': // UGC
							return PeptideCode::CYSTEINE;
						default:
							break;
					}
				case 'C'://UC-
					return PeptideCode::SERINE;
				default:
					break;
			}
		case 'G'://G--
			switch(N2) {
				case 'A'://GA-
					switch(N3) {
						case 'A'://GAA
						case 'G'://GAG
							return PeptideCode::GLUTAMIC_ACID;
						case 'U'://GAU
						case 'C'://GAC
							return PeptideCode::ASPARTIC_ACID;
						default:
							break;
					}
				case 'U'://GU-
					return PeptideCode::VALINE;
				case 'G'://GG-
					return PeptideCode::GLYCINE;
				case 'C'://GC-
					return PeptideCode::ALANINE;
				default:
					break;
			}
		case 'C':
			switch(N2) {
				case 'A':
					switch(N3) {
						case 'A':
						case 'G':
							return PeptideCode::GLUTAMINE;
						case 'U':
						case 'C':
							return PeptideCode::HISTIDINE;
						default:
							break;
					}
				case 'U':
					return PeptideCode::LEUCINE;
				case 'G':
					return PeptideCode::ARGININE;
				case 'C':
					return PeptideCode::PROLINE;
				default:
					break;
			}
		default:
			break;
	}
	throw std::runtime_error("codon_to_peptide: INVALID RNA CHARACTER");
}


//implement static functions
