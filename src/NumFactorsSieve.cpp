#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector NumFactorsSieve(int n) {
    std::vector<int> numFacs(n, 1);
    int i, j;
    
    for (i = 2; i <= n; i++) {
        for (j = i; j <= n; j+=i) {numFacs[j - 1]++;}
    }
    
    return wrap(numFacs);
}
