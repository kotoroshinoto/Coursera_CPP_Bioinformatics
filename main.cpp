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

// define static functions


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
template <class T>
std::string join(std::vector<T> &v, std::string s){
	std::stringstream ss;
	std::copy(v.begin(), v.end()-1, std::ostream_iterator<T>(ss, s.c_str()));
	ss << v[v.size() -1];
	return ss.str();
}

//implement static functions
//implement main
int main() {
	std::string separator("\\s+");
    std::queue<std::string> tests;
    InputData input;
    std::vector<int> int_v, int_v_comma;
	std::vector<std::string> str_v, str_cols;
	input.next_into_int_vector(int_v, separator);
	std::cout << "int_v length: " << int_v.size()<<std::endl<<std::flush;
//	std::copy(int_v.begin(), int_v.end(), std::ostream_iterator<int>(std::cout, ";"));
	std::cout<<join(int_v, ";")<< std::endl << std::flush;
	input.next_into_str_vector(str_cols, separator);
	std::cout << "str_cols length: " << str_cols.size()<<std::endl<<std::flush;
	std::cout<<join(str_cols, ";")<< std::endl << std::flush;
	input.next_into_int_vector(int_v_comma, separator);
	std::cout << "int_v_comma length: " << int_v_comma.size()<<std::endl<<std::flush;
	std::cout<<join(int_v_comma, ";")<< std::endl << std::flush;
    input.rest_into_str_vector(str_v);
	std::cout << "str_v length: " << str_v.size()<<std::endl<<std::flush;
	std::cout<<join(str_v, ";")<< std::endl << std::flush;
    return 0;
}