#include "peptide.cpp"
#include <array>
//Σ = {a1,...,ak}
//k = |Σ|
//µ: Σ -> N
//s = s1 ... sn over Σ
//µ(s):= summation from i=1 to n of µ(s[i])
//
/*

*/
std::vector<std::size_t> masses(Peptide::get_aa_mass_list());

std::vector<std::vector<std::size_t> > ERT; // N[i][r] i : k and r from 1 ... a1-1
std::vector<std::vector<bool> > ERT_INF; //
std::size_t gcd(std::size_t a, std::size_t b){
    std::size_t cd;
    for(int i=1;i<=a&&i<=b;i++){
        if(a%i==0 && b%i == 0 ){

            cd=i;

        }
    }
    return cd;
}

void generate_ERT() {
    std::size_t i,p,j,d,q,r,n;
    //columns correspond to a[i]
    //rows correspond to a[1] from 1 ... a1-1
    ERT.reserve(masses.size()+1);
    ERT_INF.reserve(masses.size()+1);
    for(i = 0;i<=masses.size();i++){
        ERT.emplace_back();
        ERT_INF.emplace_back();
        ERT[i].reserve(masses[0]);
        ERT_INF[i].reserve(masses[0]);
        for(r = 0;r<masses[0];r++){
            if(r > 0 && i == 0){
                ERT[i].push_back(0);
                ERT_INF[i].push_back(true);
            }else{
                ERT[i].push_back(0);
                ERT_INF[i].push_back(false);
            }
        }
    }
    for(i = 1; i < masses.size();i++){
        d = gcd(masses[0], masses[i]);
        for(p=0;p<d-1;p++){
            //find n = min{N[q]}
            if(ERT_INF[][]){
                for(j = 1;  j <masses[0];j++){
                    n = n + masses[i];
                    r = n % masses[i];
                    n = std::min(n,N[r]);
                    N[r] = n;
                }
            }
        }
    }


}


std::size_t get_count_for_compomer(std::vector<std::size_t> &c){
    std::size_t k = masses.size();
    std::size_t N = 0;
    std::vector<bool> denoms;
    denoms.reserve(k);
    for(std::size_t i=0;i<c.size();i++){
        N+=c[i];
        denoms[i] = factorial(c[i]);
    }
    double numer = factorial(N);
    for(std::size_t i=0;i<denoms.size();i++){
        numer /= denoms[i];
    }
    return (std::size_t )numer;
}

std::size_t get_peps_with_mass_recur(std::size_t i, std::size_t m, std::vector<std::size_t> &c){
    std::cout<<"get_peps_with_mass_recur start i: "<<i<<" m: "<<m<<" c: "<<join(c, ",")<<std::endl<<std::flush;
    if(i == 0){
        return get_count_for_compomer(c);
    }
    std::size_t count =0;
    if(B(i,m)){
        count += get_peps_with_mass_recur(i-1,m,c);
    }
    std::size_t ai = masses[i];
    if(ai > m){
        return 0;
    }
    if(B(i,m-ai)){
        std::vector<std::size_t> recur_c(c);
        recur_c[i]++;
        count += get_peps_with_mass_recur(i,m-ai,recur_c);
    }
    std::cout<<"get_peps_with_mass_recur end"<<std::endl<<std::flush;
    return count;
}

std::size_t count_peps_with_mass(std::size_t M){
    std::sort(masses.begin(),masses.end());
    generate_C(M);
    std::size_t k=masses.size();
    std::vector<std::size_t> c;
    c.reserve(k);
    for(std::size_t i=0;i<k;i++){
        c.push_back(0);
    }
    return get_peps_with_mass_recur(k,M,c);
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