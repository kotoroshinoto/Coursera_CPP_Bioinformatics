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
#include <map>
#include <chrono>
#include <numeric>
#include <cstdint>
#include <random>
// define static functions
std::size_t compute_original_index(std::size_t aa_pos, std::size_t orig_len, std::size_t offset, std::size_t search_aa_len, bool is_reverse);

template<typename int_type>
int_type str_to_int_type(std::string &text);

// define Objects
template<typename int_type>
class intVectorFiller {
    std::vector<int_type>& v;
public:
    intVectorFiller(std::vector<int_type>& v);
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
	template<typename int_type>
	int_type next_as_int_type();
	template<typename int_type>
	void next_into_int_type_vector(std::vector<int_type> &intv, std::string &separator);
	template<typename int_type>
	void rest_into_int_type_vector(std::vector<int_type> &intv);
    void next_into_str_vector(std::vector<std::string> &strv, std::string &separator);
    void rest_into_str_vector(std::vector<std::string> &strv);
    std::string next();
    std::size_t size();
    std::string StripWhitespace(std::string& text);
    std::string StripLeadingWhitespace(std::string& text);
    std::string StripTrailingWhitespace(std::string& text);
};

//instantiate static variables
std::regex InputData::whitespace_trimmer_regex_lead_and_trail("\\s+$|^\\s+");
std::regex InputData::whitespace_trimmer_regex_lead_only("^\\s*");
std::regex InputData::whitespace_trimmer_regex_trail_only("\\s*$");
std::regex InputData::whitespace_re("\\s+");

//implement object functions

InputData::InputData() {
	std::string inputline;
	while(std::getline(std::cin, inputline)){
		inputline = StripWhitespace(inputline);
		this->data.push(inputline);
	}
}

template<typename int_type>
int_type str_to_int_type(std::string &text){
	std::string t(text);
	if(t.find(',') >= 0) {
		t.erase(std::remove(t.begin(), t.end(), ','), t.end());
	}
	std::stringstream ss(t);
	int_type i;
	ss>>i;
	return i;
}

template<typename int_type>
intVectorFiller<int_type>::intVectorFiller(std::vector<int_type>& v): v(v){}

template<typename int_type>
void intVectorFiller<int_type>::operator()(std::string& item) {
	v.push_back(str_to_int_type<int_type>(item));
}

template<typename int_type>
int_type InputData::next_as_int_type(){
	int_type next = str_to_int_type<int_type>(this->data.front());
	this->data.pop();
	return next;
}

template<typename int_type>
void InputData::next_into_int_type_vector(std::vector<int_type> &intv, std::string &separator){
	std::vector<std::string> str_form;
	this->next_into_str_vector(str_form, separator);
	std::for_each(str_form.begin(),str_form.end(),intVectorFiller<int_type>(intv));
}

template<typename int_type>
void InputData::rest_into_int_type_vector(std::vector<int_type> &intv){
	std::vector<std::string> str_form;
	this->rest_into_str_vector(str_form);
	std::for_each(str_form.begin(),str_form.end(),intVectorFiller<int_type>(intv));
}

std::string InputData::next(){
    std::string next(this->data.front().c_str());
    this->data.pop();
    return next;
}

std::size_t InputData::size(){
    return this->data.size();
}

void InputData::rest_into_str_vector(std::vector<std::string> &strv) {
	while(this->data.size() > 0){
		strv.push_back(this->next());
	}
}

void InputData::next_into_str_vector(std::vector<std::string> &strv, std::string &separator) {
    std::regex re(separator);
    std::string &next = this->data.front();
    std::sregex_token_iterator first{next.begin(), next.end(), re, -1}, last;
    std::copy(first, last, std::back_inserter<std::vector<std::string> >(strv));
    this->data.pop();
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
struct Spectrum{
	std::vector<std::size_t> masses;
	std::map<std::size_t, std::size_t> mass_count;
	Spectrum(std::vector<std::size_t>& spectrum);
	std::vector<std::size_t> to_vector();
	std::string to_string();
	std::string to_table_string();
	std::size_t parentmass();
	bool operator==(Spectrum other);
	Spectrum convolve();
};

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
		METHIONINE,
		NONPROTEINOGENIC
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
	static std::size_t getMassFor(PeptideCode aa);
	static Peptide::PeptideCode getPeptideCodeFor(std::size_t mass);
	static std::vector<Peptide> cyclopeptide_sequencing(Spectrum& spectrum);
	static Peptide leaderboard_cyclopeptide_sequencing(Spectrum& spectrum, std::size_t N, bool extended=false);
	static std::vector<Peptide> leaderboard_cyclopeptide_sequencing_all(Spectrum& spectrum, std::size_t N, bool extended=false);
	template <typename Iter>
	static std::string pep_container_to_string(Iter it, Iter end, std::string separator="\n", bool include_stops=false);
	bool is_consistent_with_spectrum(Spectrum& spectrum);\
	Spectrum * linspec_ptr();
	Spectrum * cycspec_ptr();
//	static std::size_t get_peps_with_mass(std::size_t mass);
	Peptide();
	Peptide(std::string aastr);
	Peptide(PeptideCode& aa);
	Peptide(std::vector<PeptideCode>& aa);
	Peptide(std::vector<std::size_t>& mass_seq);
	Peptide(const Peptide& other);
	Peptide& operator=(const Peptide & other);
	~ Peptide ();
	void addAminoAcid(PeptideCode aa);
	void addAminoAcidMass(std::size_t mass);
	std::string to_mass_string(bool include_stops=false);
	template <typename Iter>
	static std::string pep_container_to_mass_string(Iter it, Iter end, std::string separator="\n", bool include_stops=false);
	std::string to_string(bool include_stops=false);
	std::string to_abbrev_string(bool include_stops=false);
	std::string to_fullword_string(bool include_stops=false);
	std::size_t number_of_sequences(bool include_stops=false);
	std::vector<std::string> subseqs_encoding(std::string& dna);
	std::size_t size(bool include_stops=false);
	std::size_t nuc_size(bool include_stops=false);
	bool operator ==(const Peptide prot2) const;
	bool operator !=(const Peptide prot2) const;
	bool operator <(const Peptide prot2) const;
	bool operator <=(const Peptide prot2) const;
	bool operator >(const Peptide prot2) const;
	bool operator >=(const Peptide prot2) const;
	bool isSameAs(const Peptide &prot2, bool include_stops=false);
	Peptide subseq(std::size_t i, std::size_t len);
	Spectrum &linear_spectrum();
	Spectrum &cyclic_spectrum();
	std::size_t get_mass();
	Peptide prefix(std::size_t len);
	std::vector<Peptide> expand(bool extended=false);
	bool has_spectrum(Spectrum& spectrum, bool cyclic=false);
	std::size_t score_spectrum(Spectrum& spectrum, bool cyclic=false);
	PeptideCode aa_at(std::size_t index);
	std::size_t aa_mass_at(std::size_t index);
	static std::vector<std::size_t> get_aa_mass_list();
	static std::vector<std::size_t> get_aa_mass_list_extended();
//	bool bind();
    static void set_extended(Spectrum s,std::size_t M);
private:
	std::size_t totalmass;
	std::vector<PeptideCode> sequence;
	std::vector<std::size_t> mass_sequence;
	Spectrum * linspectrum;
	Spectrum * cycspectrum;
	static std::vector<PeptideCode> aa_list;
	static std::vector<std::size_t> aa_mass_list;
	static std::vector<std::size_t> aa_mass_list_extended;
	void invalidate_spectra();
	void generate_linear_spectrum();
	void generate_cyclic_spectrum();

};

Peptide::Peptide(const Peptide& other):totalmass(other.totalmass), sequence(other.sequence), mass_sequence(other.mass_sequence),linspectrum(NULL),cycspectrum(NULL){
//	std::cerr<<"[copy constructor] executed"<<std::endl;
	if(other.linspectrum != NULL){
//		std::cerr<<"[copy constructor] linspectrum was not NULL"<<std::endl;
		this->linspectrum = new Spectrum(*(other.linspectrum));
//	} else {
//		std::cerr<<"[copy constructor] linspectrum was NULL"<<std::endl;
	}
	if(other.cycspectrum != NULL){
//		std::cerr<<"[copy constructor] cycspectrum was not NULL"<<std::endl;
		this->cycspectrum = new Spectrum(*(other.cycspectrum));
//	} else {
//		std::cerr<<"[copy constructor] cycspectrum was NULL"<<std::endl;
	}
}

Peptide& Peptide::operator=(const Peptide & other){
//	std::cerr<<"[operator =] executed"<<std::endl;
	this->invalidate_spectra();
	this->totalmass = other.totalmass;
	this->sequence = other.sequence;
	this->mass_sequence = other.mass_sequence;
	if(other.linspectrum != NULL){
//		std::cerr<<"[operator =] linspectrum was not NULL"<<std::endl;
		this->linspectrum = new Spectrum(*(other.linspectrum));
//	} else {
//		std::cerr<<"[operator =] linspectrum was NULL"<<std::endl;
	}
	if(other.cycspectrum != NULL){
//		std::cerr<<"[operator =] cycspectrum was not NULL"<<std::endl;
		this->cycspectrum = new Spectrum(*(other.cycspectrum));
//	} else {
//		std::cerr<<"[operator =] cycspectrum was NULL"<<std::endl;
	}
}

Spectrum * Peptide::linspec_ptr(){ return this->linspectrum;}
Spectrum * Peptide::cycspec_ptr(){return this->cycspectrum;}

Spectrum Spectrum::convolve(){
	std::vector<std::size_t> spec = this->to_vector();
	std::vector<std::size_t> conv_v;
	for(std::size_t i=0;i<spec.size()-1;i++){
		std::size_t smallmass = spec[i];
		for(std::size_t j=i+1;j<spec.size();j++){
			std::size_t bigmass = spec[j];
			if(smallmass != bigmass) {
//				std::cout << "smallmass: " << smallmass << ", bigmass:" << bigmass << std::endl;
				conv_v.push_back(bigmass - smallmass);
			}
		}
	}
	return Spectrum(conv_v);
}

bool Spectrum::operator==(Spectrum other){
	return this->mass_count == other.mass_count;
}

Spectrum::Spectrum(std::vector<std::size_t>& spectrum){
	for(std::size_t i =0;i< spectrum.size(); i++ ){
		std::size_t mass = spectrum[i];
		if(!this->mass_count.count(mass)){
			this->masses.push_back(mass);
			this->mass_count[mass] = std::count(spectrum.begin(),spectrum.end(), mass);
		}
	}
	std::sort(this->masses.begin(),this->masses.end());
}

std::size_t Spectrum::parentmass(){
	return this->masses[this->masses.size()-1];
}

std::vector<std::size_t> Spectrum::to_vector(){
	std::vector<std::size_t> spect;
	for(std::size_t i =0 ;i<this->masses.size();i++){
		std::size_t mass = this->masses[i];
		std::size_t count = this->mass_count[mass];
		for(std::size_t j = 0; j<count;j++){
			spect.push_back(mass);
		}
	}
	return spect;
}

std::string Spectrum::to_string(){
	std::stringstream ss;
	for(std::size_t i=0;i<this->masses.size();i++){
		for(std::size_t j = 0; j<this->mass_count[this->masses[i]]; j++) {
			if(!(i == 0 && j == 0)){
				ss <<" ";
			}
			ss << this->masses[i];
		}
	}
	return ss.str();
}

std::string Spectrum::to_table_string(){
	std::stringstream ss;
	for(std::size_t i=0;i<this->masses.size();i++){
		ss<<this->masses[i]<<"\t"<<this->mass_count[this->masses[i]]<<std::endl;
	}
	return ss.str();
}

std::vector<std::size_t> Peptide::get_aa_mass_list(){
	return Peptide::aa_mass_list;
}

std::vector<std::size_t> Peptide::get_aa_mass_list_extended(){
	return Peptide::aa_mass_list_extended;
}

Peptide::PeptideCode Peptide::aa_at(std::size_t index){
	return this->sequence[index];
}
std::size_t Peptide::aa_mass_at(std::size_t index){
	return Peptide::getMassFor(this->sequence[index]);
}

std::size_t Peptide::get_mass(){
	return this->totalmass;
}

void Peptide::invalidate_spectra(){
	if(this->cycspectrum != NULL){
		delete this->cycspectrum;
		this->cycspectrum = NULL;
	}
	if(this->linspectrum != NULL){
		delete this->linspectrum;
		this->linspectrum = NULL;
	}
}

void Peptide::addAminoAcid(PeptideCode aa){
	this->sequence.push_back(aa);
	this->mass_sequence.push_back(Peptide::getMassFor(aa));
	this->totalmass += this->mass_sequence.back();
	this->invalidate_spectra();
}

void Peptide::addAminoAcidMass(std::size_t mass){
	this->mass_sequence.push_back(mass);
	this->sequence.push_back(Peptide::getPeptideCodeFor(mass));
	this->totalmass += mass;
	this->invalidate_spectra();
}

//Peptide::Peptide(Peptide& other):totalmass(other.totalmass), sequence(other.sequence){}
Peptide::~Peptide (){
	this->invalidate_spectra();
}
Peptide::Peptide():totalmass(0), sequence(), mass_sequence(), linspectrum(NULL), cycspectrum(NULL){}

Peptide::Peptide(PeptideCode& aa):totalmass(0), sequence(), mass_sequence(), linspectrum(NULL), cycspectrum(NULL){
	this->addAminoAcid(aa);
}

Peptide::Peptide(std::vector<PeptideCode>& aaseq):totalmass(0), sequence(), mass_sequence(), linspectrum(NULL), cycspectrum(NULL){
	std::for_each(aaseq.begin(), aaseq.end(),std::bind1st(std::mem_fun(&Peptide::addAminoAcid), this));
}

Peptide::Peptide(std::vector<std::size_t>& mass_seq):totalmass(0), sequence(), mass_sequence(), linspectrum(NULL), cycspectrum(NULL){
	std::for_each(mass_seq.begin(), mass_seq.end(),std::bind1st(std::mem_fun(&Peptide::addAminoAcidMass), this));
}

Peptide::Peptide(std::string aastr):totalmass(0), sequence(), mass_sequence(), linspectrum(NULL), cycspectrum(NULL){
	for(std::size_t i =0;i < aastr.size();i++){
		this->addAminoAcid(Peptide::symbol_to_peptide(aastr[i]));
	}
//	std::cout<<"number of peptide entries: "<<this->sequence.size()<<std::endl<<std::flush;
}

Peptide Peptide::prefix(std::size_t len){
	return this->subseq(0,len);
}

template <typename Iter>
std::string Peptide::pep_container_to_string(Iter it, Iter end, std::string separator, bool include_stops){
	std::vector<std::string> pstr;
	for(; it < end ; it++){
		Peptide current = *it;
		pstr.push_back(current.to_string(include_stops));
	}
	return join(pstr, separator);
}

template <typename Iter>
std::string Peptide::pep_container_to_mass_string(Iter it, Iter end, std::string separator, bool include_stops){
	std::vector<std::string> pstr;
	for(; it < end ; it++){
		Peptide current = *it;
		pstr.push_back(current.to_mass_string(include_stops));
	}
	return join(pstr, separator);
}

std::string Peptide::to_mass_string(bool include_stops){
	if(!include_stops) {
		std::vector<std::size_t> sizes;
		for (std::size_t i = 0; i < this->mass_sequence.size(); i++) {
			std::size_t current = this->mass_sequence[i];
			if (current == 0) {
				break;
			}
			sizes.push_back(current);
		}
		return join(sizes, "-");
	} else {
		return join(this->mass_sequence, "-");
	}

}

std::string Peptide::to_string(bool include_stops){
	std::vector<char> aastr;
	for(std::vector<Peptide::PeptideCode>::iterator it = this->sequence.begin(); it < this->sequence.end(); it++){
		Peptide::PeptideCode current = *it;
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
	if(this->mass_sequence.size() != prot2.mass_sequence.size()){
		return false;
	}
	//then by mass
	if(this->totalmass != prot2.totalmass){
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
	for(std::size_t i = 0; i < this->mass_sequence.size() ;i++){
		if(this->mass_sequence[i] != prot2.mass_sequence[i]){
			return false;
		}
		if((!include_stops) && this->mass_sequence[i] == 0){
			break;
		}
	}
	return true;
}

bool Peptide::operator ==(const Peptide prot2) const{
	//first by length
	if(this->sequence.size() != prot2.sequence.size()){
		return false;
	}

	if(this->mass_sequence.size() != prot2.mass_sequence.size()){
		return false;
	}

	//then by mass
	if(this->totalmass != prot2.totalmass){
		return false;
	}

//	std::string mystr(this->to_string()), prot2str(prot2.to_string());
	for(std::size_t i = 0; i< this->sequence.size();i++){
		if(this->sequence[i] != prot2.sequence[i]){
			return false;
		}
	}
	for(std::size_t i = 0; i< this->mass_sequence.size();i++){
		if(this->mass_sequence[i] != prot2.mass_sequence[i]){
			return false;
		}
	}
	return true;
}


bool Peptide::operator !=(const Peptide prot2) const{
	return ! (this->operator==(prot2));
}


bool Peptide::operator <(const Peptide prot2) const{
	//first by length
	if(this->sequence.size() < prot2.sequence.size()){
		return true;
	} else if (this->sequence.size() > prot2.sequence.size()){
		return false;
	}
	//then by mass
	if(this->totalmass < prot2.totalmass){
		return true;
	} else if (this->totalmass > prot2.totalmass){
		return false;
	}

	for(std::size_t i = 0; i< this->mass_sequence.size();i++){
		if(this->mass_sequence[i] < prot2.mass_sequence[i]){
			return true;
		} else if(this->mass_sequence[i] > prot2.mass_sequence[i]){
			return false;
		}
	}
	return false;
}

bool Peptide::operator >(const Peptide prot2) const{

	//first by length
	if(this->sequence.size() > prot2.sequence.size()){
		return true;
	} else if (this->sequence.size() < prot2.sequence.size()){
		return false;
	}
	//then by mass
	if(this->totalmass > prot2.totalmass){
		return true;
	} else if (this->totalmass < prot2.totalmass){
		return false;
	}

	for(std::size_t i = 0; i< this->mass_sequence.size();i++){
		if(this->mass_sequence[i] > prot2.mass_sequence[i]){
			return true;
		} else if(this->mass_sequence[i] < prot2.mass_sequence[i]){
			return false;
		}
	}
	return false;
}
bool Peptide::operator <=(const Peptide prot2) const{

	//first by length
	if(this->sequence.size() < prot2.sequence.size()){
		return true;
	} else if (this->sequence.size() > prot2.sequence.size()){
		return false;
	}
	//then by mass
	if(this->totalmass < prot2.totalmass){
		return true;
	} else if (this->totalmass > prot2.totalmass){
		return false;
	}
//	std::string mystr(this->to_string()), prot2str(prot2.to_string());
	for(std::size_t i = 0; i< this->mass_sequence.size();i++){
		if(this->mass_sequence[i] < prot2.mass_sequence[i]){
			return true;
		} else if(this->mass_sequence[i] > prot2.mass_sequence[i]){
			return false;
		}
	}
	return true;
}
bool Peptide::operator >=(const Peptide prot2) const{


	//first by length
	if(this->sequence.size() > prot2.sequence.size()){
		return true;
	} else if (this->sequence.size() < prot2.sequence.size()){
		return false;
	}
	//then by mass
	if(this->totalmass > prot2.totalmass){
		return true;
	} else if (this->totalmass < prot2.totalmass){
		return false;
	}

	for(std::size_t i = 0; i< this->mass_sequence.size();i++){
		if(this->mass_sequence[i] > prot2.mass_sequence[i]){
			return true;
		} else if(this->mass_sequence[i] < prot2.mass_sequence[i]){
			return false;
		}
	}
	return true;
}

Peptide Peptide::subseq(std::size_t i, std::size_t len){
	Peptide newpep;
	for(std::size_t j = i; j - i < len && j < this->sequence.size(); j++){
		if(this->sequence[j] == PeptideCode::NONPROTEINOGENIC){
			newpep.addAminoAcidMass(this->mass_sequence[j]);
		} else {
			newpep.addAminoAcid(this->sequence[j]);
		}
	}
	return newpep;
}

std::size_t compute_original_index(std::size_t aa_pos, std::size_t orig_len, std::size_t offset, std::size_t search_aa_len, bool is_reverse){
	std::size_t nuc_pos = aa_pos * 3;
	if(is_reverse){
		std::size_t rev_index = nuc_pos + offset; // position in original reversed string before bases were truncated away
		std::size_t rev_index_end = rev_index + (search_aa_len*3) - 1; //0-3 in 4 lenstr 0 + 4 -1
		return orig_len - rev_index_end - 1; //originally 6 long, position 4 in original is 1 in reverse, and 0 in offset of 1
	} else{
		return nuc_pos + offset;
	}
}

std::vector<std::string> Peptide::subseqs_encoding(std::string& dna){
	std::vector<std::string> strv;
	if(this->sequence.size() * 3 > dna.size()){
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
	for(std::size_t i =0; i<3 && this->sequence.size() * 3 <= dna.size();i++){
		std::string subdna = dna.substr(i,dna.size()-i); //ATCG -> TGC
		std::string subrdna = rdna.substr(i,dna.size()-i); // ATGC -> CGTA -> GCAT -> CAT
		frames.push_back(subdna);
		frames.push_back(subrdna);
		transframes.push_back(Peptide::translateDNA(subdna));
		coords.push_back(i);
		is_reverse.push_back(false);
		transframes.push_back(Peptide::translateDNA(subrdna));
		coords.push_back(i);
		is_reverse.push_back(true);
	}
	for(std::size_t i = 0; i < transframes.size() ; i++){
		std::string pointers(dna.size(), ' ');
		std::size_t offset = coords[i];
		for(std::size_t j = 0;this->size(true) + j <= transframes[i].size(true) && j < transframes[i].size(true) ; j++) {
			//make subseqs of same size as self
			Peptide subseq = transframes[i].subseq(j,this->sequence.size());
			std::size_t substr_pos = j * 3;
			std::size_t orig_index = compute_original_index(j,dna.size(),offset, this->size(true),is_reverse[i]);
			if(this->isSameAs(subseq)){
				strv.push_back(dna.substr(orig_index,this->nuc_size(true)));
				pointers[orig_index] = '^';
			}
		}
	}
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
		case '@':
			return PeptideCode::NONPROTEINOGENIC;
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
		case NONPROTEINOGENIC:
			return 0;
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
		case NONPROTEINOGENIC:
			return "Non_Proteinogenic";
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
		case NONPROTEINOGENIC:
			return "@@@";
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
		case NONPROTEINOGENIC:
			return '@';
	}
}

Peptide::PeptideCode Peptide::getPeptideCodeFor(std::size_t mass){
	switch(mass){
		case 0:
			return PeptideCode::STOP;
		case 137:
			return PeptideCode::HISTIDINE;
		case 128:
			return PeptideCode::LYSINE;
//			return PeptideCode::GLUTAMINE;
		case 97:
			return PeptideCode::PROLINE;
		case 156:
			return PeptideCode::ARGININE;
		case 113:
			return PeptideCode::ISOLEUCINE;
//			return PeptideCode::LEUCINE;
		case 115:
			return PeptideCode::ASPARTIC_ACID;
		case 129:
			return PeptideCode::GLUTAMIC_ACID;
		case 71:
			return PeptideCode::ALANINE;
		case 57:
			return PeptideCode::GLYCINE;
		case 99:
			return PeptideCode::VALINE;
		case 163:
			return PeptideCode::TYROSINE;
		case 87:
			return PeptideCode::SERINE;
		case 103:
			return PeptideCode::CYSTEINE;
		case 186:
			return PeptideCode::TRYPTOPHAN;
		case 147:
			return PeptideCode::PHENYLALANINE;
		case 114:
			return PeptideCode::ASPARAGINE;
		case 101:
			return PeptideCode::THREONINE;
		case 131:
			return PeptideCode::METHIONINE;
		default:
			return PeptideCode::NONPROTEINOGENIC;
		case 1 ... 56:
		case 201 ... SIZE_MAX:
			std::stringstream ss;
			ss<< "[Peptide::getPeptideCodeFor] no peptidecode exists for mass: "<<mass;
			throw std::runtime_error(ss.str());
	}
}

std::size_t Peptide::getMassFor(PeptideCode aa){
	switch (aa){
		case STOP:
			return 0;
		case HISTIDINE:
			return 137;
		case GLUTAMINE:
		case LYSINE:
			return 128;
		case PROLINE:
			return 97;
		case ARGININE:
			return 156;
		case LEUCINE:
		case ISOLEUCINE:
			return 113;
		case ASPARTIC_ACID:
			return 115;
		case GLUTAMIC_ACID:
			return 129;
		case ALANINE:
			return 71;
		case GLYCINE:
			return 57;
		case VALINE:
			return 99;
		case TYROSINE:
			return 163;
		case SERINE:
			return 87;
		case CYSTEINE:
			return 103;
		case TRYPTOPHAN:
			return 186;
		case PHENYLALANINE:
			return 147;
		case ASPARAGINE:
			return 114;
		case THREONINE:
			return 101;
		case METHIONINE:
			return 131;
		case NONPROTEINOGENIC:
			throw std::runtime_error("[Peptide::getMassFor] cannot return mass for NONPROTEINOGENIC");
	}
}

Peptide::PeptideCode Peptide::codon_to_peptide(char cod1, char cod2, char cod3) {
	char N1 = (char) std::toupper(cod1);
	char N2 = (char) std::toupper(cod2);
	char N3 = (char) std::toupper(cod3);
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
					break;
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
					break;
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
					break;
				case 'C': // AC-
					return PeptideCode::THREONINE;
				default:
					break;
			}
			break;
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
					break;
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
					break;
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
					break;
				case 'C'://UC-
					return PeptideCode::SERINE;
				default:
					break;
			}
			break;
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
					break;
				case 'U'://GU-
					return PeptideCode::VALINE;
				case 'G'://GG-
					return PeptideCode::GLYCINE;
				case 'C'://GC-
					return PeptideCode::ALANINE;
				default:
					break;
			}
			break;
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
					break;
				case 'U':
					return PeptideCode::LEUCINE;
				case 'G':
					return PeptideCode::ARGININE;
				case 'C':
					return PeptideCode::PROLINE;
				default:
					break;
			}
			break;
		default:
			break;
	}
	throw std::runtime_error("codon_to_peptide: INVALID RNA CHARACTER");
}

Spectrum &Peptide::linear_spectrum(){
	this->generate_linear_spectrum();
	return *(this->linspectrum);
}

Spectrum &Peptide::cyclic_spectrum(){
	this->generate_cyclic_spectrum();
	return *(this->cycspectrum);
}


void Peptide::generate_linear_spectrum(){
	if(this->linspectrum != NULL){
		return;// time saver, no need to recount and recompile this if its been done
	}
	std::vector< std::size_t > prefix_mass;
	prefix_mass.push_back(0);
	for(std::size_t i=0;i<this->sequence.size();i++){
		prefix_mass.push_back(prefix_mass[i] + this->mass_sequence[i]);
	}
	std::vector< std::size_t > linspec;
	linspec.push_back(0);
	for(std::size_t i = 0;i < this->sequence.size();i++){
		for(std::size_t j = i+1;j <= this->sequence.size();j++){
			linspec.push_back(prefix_mass[j] - prefix_mass[i]);
		}
	}
	std::sort(linspec.begin(),linspec.end());
	this->linspectrum = new Spectrum(linspec);
}

void Peptide::generate_cyclic_spectrum(){
	if(this->cycspectrum != NULL){
		return; // time saver, no need to recount and recompile this if its been done
	}
	std::vector< std::size_t > prefix_mass;
	prefix_mass.push_back(0);
	for(std::size_t i=0;i<this->sequence.size();i++){
		prefix_mass.push_back(prefix_mass[i] + this->mass_sequence[i]);
	}
	std::size_t peptide_mass = prefix_mass[this->sequence.size()];
	std::vector< std::size_t > cycspec;
	cycspec.push_back(0);
	for(std::size_t i = 0;i < this->sequence.size();i++){
		for(std::size_t j = i+1;j <= this->sequence.size();j++){
			cycspec.push_back(prefix_mass[j] - prefix_mass[i]);
			if(i > 0 && j < this->sequence.size()){
				cycspec.push_back(peptide_mass- (prefix_mass[j] - prefix_mass[i]));
			}
		}
	}
	std::sort(cycspec.begin(),cycspec.end());
	this->cycspectrum = new Spectrum(cycspec);
}

std::vector<std::size_t> Peptide::aa_mass_list= {57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186};
std::vector<std::size_t> extended_list_initializer(){
	std::vector<std::size_t> aa_masses(200-75+1);
	std::iota(aa_masses.begin(), aa_masses.end(), 57);
	return aa_masses;
}
std::vector<std::size_t> Peptide::aa_mass_list_extended(extended_list_initializer());
void Peptide::set_extended(Spectrum s, std::size_t M){
    Spectrum conv(s.convolve());
    std::vector<std::pair<std::size_t,std::size_t>> m_highest;
    m_highest.reserve(conv.masses.size());
    for (size_t i = 0; i < conv.masses.size(); i++) {
        std::size_t mass = conv.masses[i];
//        std::cout<<"mass: "<<mass<<std::endl;
        if(mass >=57 && mass <=200){
            m_highest.push_back(std::pair<std::size_t, std::size_t>(conv.mass_count[mass],mass));
        }
    }
    std::sort(m_highest.begin(), m_highest.end(),std::greater<std::pair<std::size_t,std::size_t> >());
    std::vector<std::size_t> new_aa_mass_list;
    std::pair<std::size_t,std::size_t> & mth_mass= m_highest[M-1];
//    std::cout<<"mth mult: "<<mth_mass.first <<"; mass: "<<mth_mass.second<<std::endl;
//    std::cout<<"m_highest.size(): "<<m_highest.size()<<std::endl;
    for(std::size_t i=M;i< m_highest.size();i++){
        std::pair<std::size_t,std::size_t> & current= m_highest[i];
//        std::cout<<"ith mult: "<<current.first <<"; ith mass: "<<current.second<<std::endl;
        if(current.first < mth_mass.first){
            m_highest.erase(m_highest.begin()+i,m_highest.end());
            break;
        }
    }
    new_aa_mass_list.reserve(m_highest.size());
    for(std::size_t i=0;i<m_highest.size();i++){
        std::pair<std::size_t,std::size_t> & current = m_highest[i];
        new_aa_mass_list.push_back(current.second);
    }
    std::sort(new_aa_mass_list.begin(),new_aa_mass_list.end());
    Peptide::aa_mass_list_extended.swap(new_aa_mass_list);
}

std::vector<Peptide::PeptideCode> Peptide::aa_list={
		Peptide::PeptideCode::GLYCINE,  //G
		Peptide::PeptideCode::ALANINE,  //A
		Peptide::PeptideCode::SERINE,  //S
		Peptide::PeptideCode::PROLINE,  //P
		Peptide::PeptideCode::VALINE,  //V
		Peptide::PeptideCode::THREONINE,  //T
		Peptide::PeptideCode::CYSTEINE,  //C
		Peptide::PeptideCode::ISOLEUCINE,  //I
//		Peptide::PeptideCode::LEUCINE,  //L mass is redundant
		Peptide::PeptideCode::ASPARAGINE,  //N
		Peptide::PeptideCode::ASPARTIC_ACID,  //D
		Peptide::PeptideCode::LYSINE,  //K
//		Peptide::PeptideCode::GLUTAMINE,  //Q mass is redundant
		Peptide::PeptideCode::GLUTAMIC_ACID,  //E
		Peptide::PeptideCode::METHIONINE,  //M
		Peptide::PeptideCode::HISTIDINE,  //H
		Peptide::PeptideCode::PHENYLALANINE,  //F
		Peptide::PeptideCode::ARGININE,  //R
		Peptide::PeptideCode::TYROSINE,  //Y
		Peptide::PeptideCode::TRYPTOPHAN,  //W
};

std::vector<Peptide> Peptide::expand(bool extended){
	std::vector<Peptide> new_peps;
	if(extended){
		for (std::size_t i = 0; i < Peptide::aa_mass_list_extended.size(); i++) {
			Peptide newpep = (*this);
			newpep.addAminoAcidMass(Peptide::aa_mass_list_extended[i]);
			new_peps.push_back(newpep);
		}
	} else {
		for (std::size_t i = 0; i < Peptide::aa_list.size(); i++) {
			Peptide newpep = (*this);
			newpep.addAminoAcid(Peptide::aa_list[i]);
			new_peps.push_back(newpep);
		}
	}
	return new_peps;
}

std::size_t Peptide::score_spectrum(Spectrum& spectrum, bool cyclic){
	Spectrum &myspectrum(cyclic ? this->cyclic_spectrum() : this->linear_spectrum());
	std::size_t score = 0;
	// we have x count in our spectrum
	// other spectrum has y count for mass
	// x > y, y contribution to count
	// x < y, x contribution to the count
	// min-value contributes
	if(myspectrum.masses.size() <= spectrum.masses.size()){
		for(std::size_t i = 0; i < myspectrum.masses.size(); i++){
			std::size_t mass = myspectrum.masses[i];
			if(spectrum.mass_count.count(mass) == 1){
				score += std::min(myspectrum.mass_count[mass],spectrum.mass_count[mass]);
			}
		}
	} else {
		for(std::size_t i = 0; i < spectrum.masses.size(); i++){
			std::size_t mass = spectrum.masses[i];
			if(myspectrum.mass_count.count(mass) == 1) {
				score += std::min(myspectrum.mass_count[mass],spectrum.mass_count[mass]);
			}
		}
	}
	return score;
}

bool Peptide::has_spectrum(Spectrum& spectrum, bool cyclic){
	Spectrum &myspectrum(cyclic ? this->cyclic_spectrum() : this->linear_spectrum());
	return myspectrum == spectrum;
}

bool Peptide::is_consistent_with_spectrum(Spectrum& spectrum){
	//spectrum is allowed to have values that we don't.
	//we are not allowed to have values that are not in spectrum;
	//we cannot have more counted values for any given mass than spectrum
	this->generate_linear_spectrum();
	if(this->linspectrum->masses.size() > spectrum.masses.size()){
		//we have more masses, might as well quit now
		return false;
	}
	for(std::size_t j =0;j< this->linspectrum->masses.size(); j++ ){ // for each unique value found in our spectrum
		std::size_t mass = this->linspectrum->masses[j];
		if(spectrum.mass_count.count(mass) != 1){ // if our mass is not in spectrum, we're inconsistent
			return false;
		} else if(spectrum.mass_count[mass] < this->linspectrum->mass_count[mass]){
			//if there are more of ANY value in my spectrum, I am not a subspectrum and not consistent
			return false;
		}
	}
	return true;
}

std::vector<Peptide> Peptide::cyclopeptide_sequencing(Spectrum& spectrum){
	std::size_t spectrum_parentmass = spectrum.parentmass();
	std::queue<Peptide> peps;
	peps.push(Peptide());
	std::vector<Peptide> output;
	while(peps.size() > 0){
		Peptide currpep = peps.front();
		peps.pop();
		if(currpep.get_mass() > spectrum_parentmass){
			continue;
		}
		else if (currpep.get_mass() == spectrum_parentmass){ //if mass of peptide is consistent with given spectrum's ParentMass
			if(currpep.has_spectrum(spectrum, true)) {//if peptide's cyclic spectrum is consistent with given spectrum
				output.push_back(currpep);
			}
			//remove peptide from peptide collection i.e. do not calculate its expansion
			//given my change to the algorithm, doing nothing accomplishes this
		} else if(currpep.is_consistent_with_spectrum(spectrum)){
			std::vector<Peptide> expanded = currpep.expand();
			for(std::size_t i = 0; i < expanded.size();i++){
				peps.push(expanded[i]);
			}
		}
		//else if Peptide's spectrum is inconsistent with given spectrum
		//remove peptide from peptide collection i.e. do not calculate its expansion
		//^ was rephrased as if consistent DO calculate its expansion
	}
	return output;
}

void set_score(std::map<Peptide, std::size_t> &scoremap, Peptide &pep, Spectrum &spec, bool cyclic=false){
	if(scoremap.count(pep) != 1){
		scoremap[pep] = pep.score_spectrum(spec, cyclic);
	}
}
std::size_t get_score(std::map<Peptide, std::size_t> &scoremap, Peptide &pep, Spectrum &spec, bool cyclic=false){
	set_score(scoremap,pep,spec,cyclic);
	return scoremap[pep];
}

void remove_score(std::map<Peptide, std::size_t> &scoremap, Peptide &pep){
	if(scoremap.count(pep) == 1) {
		scoremap.erase(pep);
	}
}

class ScorePeptideSorter{
public:
	ScorePeptideSorter(std::map<Peptide, std::size_t> &scores);
	bool operator()(const Peptide p1,const Peptide p2)const ;
private:
	std::map<Peptide, std::size_t> &scores;
};

ScorePeptideSorter::ScorePeptideSorter(std::map<Peptide, std::size_t> &_scores): scores(_scores){}

bool ScorePeptideSorter::operator()(const Peptide p1,const Peptide p2)const {
	return scores[p1] < scores[p2];
}

void trim_leaderboard(std::vector<Peptide> &leaderboard, std::map<Peptide, std::size_t> &linscores, std::map<Peptide, std::size_t> &cycscores, std::size_t N){
	std::cout<<"[trim_leaderboard] begin"<<std::endl<<std::flush;
	if(leaderboard.size() > N){
		std::vector<Peptide> removed;
		std::cout<<"pre-trim leaderboard size: "<<leaderboard.size()<<std::endl;
	//	std::cout<< Peptide::pep_container_to_string(leaderboard.begin(), leaderboard.end(), " ")<<std::endl<<std::endl;
		std::sort(leaderboard.rbegin(), leaderboard.rend(), ScorePeptideSorter(linscores));
	//	std::cout<< Peptide::pep_container_to_string(leaderboard.begin(), leaderboard.end(), " ")<<std::endl;
		std::size_t nth_score = linscores[leaderboard[N-1]];
//		std::cout<<"nth score: "<<nth_score<<std::endl;
		for(std::size_t i=N;i<leaderboard.size();i++){
			std::size_t score = linscores[leaderboard[i]];
//			std::cout<<"ith score: "<<score<<std::endl;
			if( score < nth_score){
				removed = std::vector<Peptide>(leaderboard.begin()+i,leaderboard.end());
				leaderboard.erase(leaderboard.begin()+i,leaderboard.end());
//				leaderboard.resize(i);
//				std::cout<<"leaderboard truncated using i: "<<i<<std::endl;
				break;
			}
		}
		for(std::size_t i=0;i<removed.size();i++){
			//remove scores from scoremaps
			remove_score(linscores,removed[i]);
			remove_score(cycscores,removed[i]);
		}
		std::cout<<"post-trim leaderboard size: "<<leaderboard.size()<<std::endl;
	} else {
		std::cout<<"trim not required: "<<std::endl;
	}
	std::cout<<"[trim_leaderboard] end"<<std::endl<<std::flush;
}

void expand_leaderboard(std::vector<Peptide> &leaderboard, std::size_t parentmass, std::map<Peptide, std::size_t> &linscores, std::map<Peptide, std::size_t> &cycscores, Spectrum& spectrum, bool extended=false){
	std::cout<<"[expand_leaderboard] begin"<<std::endl<<std::flush;
	//pull old leaderboard into new vector
	std::vector<Peptide> _leaderboard;
	_leaderboard.reserve(leaderboard.size() * 20);
	//expand all peptides in leaderboard
//	std::cout<<"[expand]size of old leaderboard: "<<leaderboard.size()<<std::endl;

	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	std::chrono::duration<double> score_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::duration<double>::zero());
	for(std::size_t i= 0; i < leaderboard.size(); i++){
		Peptide &currpep = leaderboard[i];
		std::vector<Peptide> expandv = currpep.expand(extended);
//		std::cout<<"number of expanded peps: "<<expandv.size()<<std::endl;
		for(std::size_t j = 0; j < expandv.size();j++){
			Peptide &newpep = expandv[j];
//			std::cout<<"expanded pep: "<<newpep.to_string()<<std::endl;
			//avoid adding peptides that are too heavy
			if(newpep.get_mass() <= parentmass){
//				std::cout<<"keeping this pep"<<std::endl;
				_leaderboard.push_back(newpep);
				std::chrono::steady_clock::time_point score_start = std::chrono::steady_clock::now();
				set_score(linscores,newpep,spectrum,false);
				set_score(cycscores,newpep,spectrum,true);
				std::chrono::steady_clock::time_point score_end = std::chrono::steady_clock::now();
				score_time += std::chrono::duration_cast<std::chrono::duration<double> >(score_end - score_start);
			}
		}
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::chrono::duration<double> total_time = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
	leaderboard.swap(_leaderboard);
//	std::cout<<"[expand]size of new leaderboard: "<<leaderboard.size()<<std::endl;
	std::cout<<"[expand_leaderboard] end {total time elapsed: "<<total_time.count()<<"; time spent on scoring: "<<score_time.count()<<"; rest: "<<(total_time-score_time).count()<<"}"<<std::endl<<std::flush;
}

void update_leader(std::vector<Peptide> &leaderboard, std::map<Peptide, std::size_t> &cycscores, Peptide &leader_pep, std::size_t spectrum_parentmass){
	std::cout<<"[update_leader] begin"<<std::endl<<std::flush;
	for(std::size_t i = 0 ; i < leaderboard.size(); i++){
		Peptide &currpep = leaderboard[i];
		if (currpep.get_mass() == spectrum_parentmass) {
			//if mass of peptide is consistent with given spectrum's ParentMass
			if (cycscores[currpep] > cycscores[leader_pep]) {
				//if peptide's cyclic spectrum has leading score
				leader_pep = currpep;
			}
		}
	}
	std::cout<<"[update_leader] end"<<std::endl<<std::flush;
}

void update_leaders(std::vector<Peptide> &leaderboard, std::map<Peptide, std::size_t> &cycscores, std::vector<Peptide> &leader_peps, std::size_t &leading_score, std::size_t spectrum_parentmass){
	std::cout<<"[update_leaders] begin"<<std::endl<<std::flush;
	bool score_changed = false;
	for(std::size_t i = 0 ; i < leaderboard.size(); i++){
		Peptide &currpep = leaderboard[i];
		if (currpep.get_mass() == spectrum_parentmass) {
			//if mass of peptide is consistent with given spectrum's ParentMass
			if (cycscores[currpep] > leading_score) {
				//if peptide's cyclic spectrum has leading score
				leading_score = cycscores[currpep];
				score_changed = true;
			}
		}
	}
	if(score_changed) {
		//replace leader peptides with the new, higher scoring peptides
		std::vector<Peptide> new_leaders;
		for (std::size_t i = 0; i < leaderboard.size(); i++) {
			Peptide &currpep = leaderboard[i];
			if (currpep.get_mass() == spectrum_parentmass) {
				//if mass of peptide is consistent with given spectrum's ParentMass
				if (cycscores[currpep] == leading_score) {
					//if peptide's cyclic spectrum has leading score
					new_leaders.push_back(currpep);
				}
			}
		}
		leader_peps.swap(new_leaders);
	} else {
		//check leaderboard for peptides with scores that match existing leader score
		for (std::size_t i = 0; i < leaderboard.size(); i++) {
			Peptide &currpep = leaderboard[i];
			if (currpep.get_mass() == spectrum_parentmass) {
				//if mass of peptide is consistent with given spectrum's ParentMass
				if (cycscores[currpep] == leading_score) {
					//if peptide's cyclic spectrum has leading score
					leader_peps.push_back(currpep);
				}
			}
		}
	}
	std::cout<<"[update_leaders] end"<<std::endl<<std::flush;
}

std::vector<Peptide> Peptide::leaderboard_cyclopeptide_sequencing_all(Spectrum& spectrum, std::size_t N, bool extended){
	std::cout<<"[leaderboard_cyclopeptide_sequencing_all] begin "<<std::endl<<std::flush;
	std::cout<<"N: "<<N<<std::endl<<std::flush;
	std::size_t spectrum_parentmass = spectrum.parentmass();
	std::vector<Peptide> leaderboard;
	std::map<Peptide, std::size_t> linscores;
	std::map<Peptide, std::size_t> cycscores;
	std::vector<Peptide> leader_peps;
	std::size_t leading_score(0);
	leader_peps.push_back(Peptide());
	leaderboard.push_back(leader_peps[0]);
	while(leaderboard.size() > 0){
		expand_leaderboard(leaderboard,spectrum_parentmass, linscores, cycscores, spectrum, extended);
		update_leaders(leaderboard, cycscores, leader_peps, leading_score, spectrum_parentmass);
		trim_leaderboard(leaderboard, linscores, cycscores, N);
	}
	std::cout<<"[leaderboard_cyclopeptide_sequencing_all] completed; Number of leaders: "<<leader_peps.size()<< "; final leading score: "<<leading_score<<std::endl<<std::flush;
	return leader_peps;
}

Peptide Peptide::leaderboard_cyclopeptide_sequencing(Spectrum& spectrum, std::size_t N, bool extended){
	std::vector<Peptide> leader_peps = Peptide::leaderboard_cyclopeptide_sequencing_all(spectrum,N,extended);
	if(leader_peps.size() == 0){
		return Peptide();
	}
	if (leader_peps.size() == 1){
		return leader_peps[0];
	}
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<std::size_t> dis(0,leader_peps.size()-1);
	return leader_peps[dis(gen)];
}

//implement static functions

//std::size_t Peptide::get_peps_with_mass(std::size_t mass){
//	std::queue<Peptide> peps;
//	peps.push(Peptide());
//	std::vector<Peptide> output;
//	while(peps.size() > 0){
//		Peptide currpep = peps.front();
//		std::size_t currmass = currpep.get_mass();
//		peps.pop();
//		if(currmass > mass){
//			continue;
//		}
//		else if (currmass == mass){
//			output.push_back(currpep);
//		} else { // currmass < mass
//			std::vector<Peptide> expanded = currpep.expand();
//			for(std::size_t i = 0; i < expanded.size();i++){
//				peps.push(expanded[i]);
//			}
//		}
//	}
//	return output.size();
//}