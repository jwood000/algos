/*! Factoring with Multiple Polynomial Quadratic Sieve.\
 *  
 *  From https://codegolf.stackexchange.com/ (Credit to user primo for answer)
 *      P., & Chowdhury, S. (2012, October 7). Fastest semiprime factorization. Retrieved October 06, 2017, from 
 *          https://codegolf.stackexchange.com/questions/8629/fastest-semiprime-factorization/9088#9088
 *          
 *  Paper 1
 *      URL: http://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866119-8/S0025-5718-1987-0866119-8.pdf
 *      citation: Silverman, R. D. (1987). The Multiple Polynomial Quadratic Sieve. 
 *                  Mathematics of Computation, 48(177), 329-339. doi:10.2307/2007894
 *
 *  Paper 2
 *      URL: http://www.cs.virginia.edu/crab/QFS_Simple.pdf
 *      author: Eric Landquist
 *      date: December 14, 2001
 *      title: The Quadratic Sieve Factoring Algorithm
 *
 *  Paper 3
 *      URL: https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/
 *      author: MSDN Archive
 *      date: June 19, 2006
 *      title: Factoring large numbers with quadratic sieve
 *
 *  Paper 4
 *      URL: http://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf
 *
 *  Paper 5
 *      URL: http://library.msri.org/books/Book44/files/03carl.pdf
 *      citation: Pomerance, C. (2008). Smooth numbers and the quadratic sieve. In Algorithmic
 *       Number Theory Lattices, Number Fields, Curves and Cryptography (pp. 69-81). 
 *       Cambridge: Cambridge University Press.
 *
 *  Paper 6
 *      URL: http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf
 *      author: Chad Seibert
 *      date: March 16, 2011
 *      title: Integer Factorization using the Quadratic Sieve
 * 
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
// This is needed as cinttypes is C++11
#include <inttypes.h>
#include <math.h>
#include "Rgmp.h"
#include "quadraticsieve.h"

typedef std::vector<signed long int> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

static inline v1d outersect (v1d x, v1d y) {
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    unsigned long int lenX = x.size(), lenY = y.size();
    v1d v(lenX + lenY);
    std::vector<signed long int>::iterator it, it2;
    it = std::set_difference(x.begin(), x.end(), y.begin(), y.end(), v.begin());
    it2 = std::set_difference(y.begin(), y.end(), x.begin(), x.end(), it);
    v.resize(it2 - v.begin());
    return v;
}

static bigvec getFacs (bigvec v1, bigvec v2, mpz_t n) {
    mpz_t x, y;
    bigvec myAns;

    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(y, 1);

    signed long int sv1 = v1.size();
    signed long int sv2 = v2.size();
    signed long int i;

    for (i = 0; i < sv1; i++) {mpz_mul(x, x, v1[i].value.getValue());}
    mpz_mod(x, x, n);

    for (i = 0; i < sv2; i++) {mpz_mul(y, y, v2[i].value.getValue());}
    mpz_mod(y, y, n);

    myAns.push_back(x);
    myAns.push_back(y);
    std::sort(myAns.value.begin(), myAns.value.end());

    mpz_clear(x);
    mpz_clear(y);

    return myAns;
}

static void reduceMatrix (unsigned long int n1,
                          unsigned long int n2, v2d & nullMat, v1d & myCols) {
    unsigned long int i, j, myMin, temp, myMax = 0;
    std::vector<signed long int>::iterator it;
    v1d myOnes;

    for (i = 0; i < n1; i++) {
        myOnes.reserve(n2);
        for (j = myMax; j < n2; j++) {
            if (nullMat[j][i] == 1) {
                myOnes.push_back(j);
            }
        }
        if (myOnes.size() > 0) {
            myMin = myOnes[0];
            if (myMin != myMax) {
                for (j = 0; j < n1; j++) {
                    temp = nullMat[myMin][j];
                    nullMat[myMin][j] = nullMat[myMax][j];
                    nullMat[myMax][j] = temp;
                }
            }
            myOnes.erase(myOnes.begin());
            for (it = myOnes.begin(); it < myOnes.end(); it++) {
                for (j = 0; j < n1; j++) {
                    nullMat[*it][j] = (nullMat[*it][j] + nullMat[myMax][j]) % 2;
                }
            }
            myMax++;
        }
        myOnes.clear();
    }

    unsigned long int newLen;

    if (myMax < n2 && myMax != 0) {
        for (j = myMax; j < n2; j++) {
            nullMat.erase(nullMat.begin() + j);
        }
        newLen = nullMat.size();
    } else {
        newLen = 0;
    }

    bool allZero;

    if (newLen > 0) {
        i = 0;
        while (i < newLen) {
            allZero = true;
            for (j = 0; j < n1; j++) {
                if (nullMat[i][j] != 0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                nullMat.erase(nullMat.begin() + i);
                newLen--;
                continue;
            }
            if (nullMat[i][i] != 1) {
                for (j = 0; j < n1; j++) {
                    if (nullMat[i][j] == 1) {
                        myMin = j;
                        break;
                    }
                }
                for (j = 0; j < newLen; j++) {
                    temp = nullMat[j][i];
                    nullMat[j][i] = nullMat[j][myMin];
                    nullMat[j][myMin] = temp;
                }
                temp = myCols[i];
                myCols[i] = myCols[myMin];
                myCols[myMin] = temp;
            }
            i++;
        }
    }
}

static bool solutionSearch (v2d mat, mpz_t M2, mpz_t n, v1d FB) {
    unsigned long int nrow = mat.size(), ncol = mat[0].size();
    signed long int i, j, k, r = 0;
    v2d nullMat;

    for (j = 0; j < ncol; j++) {
        i = 0;
        while (mat[i][j] == 0 && i < nrow) {i++;}
        if (i < nrow) {
            nullMat.push_back(v1d(nrow, 0));
            for (k = 0; k < nrow; k++) {nullMat[r][k] = mat[k][j] % 2;}
            r++;
        }
    }

    v1d myCols(ncol, 0);
    for (i = 0; i < myCols.size(); i++) {myCols[i] = i;}
    reduceMatrix (ncol, nrow, nullMat, myCols);

    unsigned long int tLen, newNrow = nullMat.size();
    std::vector<signed long int>::iterator it, it2;
    v2d myList(ncol, v1d());
    v1d freeVariables, temp;
    freeVariables.reserve(ncol);
    
    if (ncol > newNrow) {
        for (i = newNrow + 1; i < ncol; i++) {freeVariables.push_back(myCols[i]);}
        for (it = freeVariables.begin(); it < freeVariables.end(); it++) {
            myList[*it].push_back(*it);
        }
    }
    
    bool allBigNewNrow;
    unsigned long int lenI;
    v1d newVec;
    
    if (newNrow > 0) {
        for (i = newNrow; i > 0; i--) {
            temp.reserve(ncol);
            for (j = 0; j < ncol; j++) {
                if (nullMat[i][j] == 1) {
                    temp.push_back(j);
                }
            }
            if (temp.size() == 1) {
                for (j = 0; j < nrow; j++) {nullMat[j][temp[0]] = 0;}
                myList[myCols[i]].clear();
                myList[myCols[i]].push_back(0);
            } else {
                temp.clear();
                temp.reserve(ncol);
                allBigNewNrow = true;
                for (j = i+1; j < ncol; j++) {
                    if (nullMat[i][j] == 1) {
                        temp.push_back(j);
                        if (allBigNewNrow) {
                            if (j < newNrow) {  // possibly need to change to <=
                                allBigNewNrow = false;
                            }
                        }
                    }
                }
                
                if (allBigNewNrow) {
                    myList[myCols[i]].clear();
                    for (j = 0; j < temp.size(); j++) {
                        tLen = myList[myCols[temp[j]]].size();
                        if (tLen > 0) {
                            for (k = 0; k < tLen; j++) {
                                myList[myCols[i]].push_back(myList[myCols[temp[k]]][j]);
                            }
                        }
                    }
                } else {
                    lenI = myList[myCols[i]].size();
                    for (it = temp.begin(); it < temp.end(); it++) {
                        Rprintf("\n%d", lenI);
                        newVec = outersect(myList[myCols[i]], myList[myCols[*it]]);
                        myList[myCols[i]].assign(newVec.begin(), newVec.end());
                    }
                }
            }
        }
    }
    
    unsigned long int lenFree = freeVariables.size();
    
    if (lenFree > 0) {
        
    }
    return(true);
}

// ## still in SOLUTION SEARCH
//     listSol <- nullMat[[1]]; freeVar <- nullMat[[2]]; LF <- length(freeVar)
//     if (LF > 0L) {
//         for (i in 2:min(10^8,(2^LF + 1L))) {
//             PosAns <- MyIntToBit(i, LF)
//             posVec <- sapply(listSol, function(x) {
//                 t <- which(freeVar %in% x)
//                 if (length(t)==0L) {
//                     0
//                 } else {
//                     sum(PosAns[t])%%2L
//                 }
//             })
//             ansVec <- which(posVec==1L)
//             if (length(ansVec)>0) {
//
//                 if (length(ansVec) > 1L) {
//                     myY <- apply(mymat[ansVec,],2,sum)
//                 } else {
//                     myY <- mymat[ansVec,]
//                 }
//
//                 if (sum(myY %% 2) < 1) {
//                     myY <- as.integer(myY/2)
//                     myY <- pow.bigz(FB,myY[-1])
//                     temp <- GetFacs(M2[ansVec], myY, n)
//                     if (!(1==temp[1]) && !(1==temp[2])) {
//                         return(temp)
//                     }
//                 }
//             }
//         }
//     }
// }

static v1d getPrimesQuadRes (mpz_t myN, double n) {
    std::vector<bool> primes(n+1, true);
    v1d myps;
    myps.reserve(floor(2*n/log(n)));

    signed long int lastP = 3;
    signed long int fsqr = floor(sqrt(n));
    signed long int k, ind, j;

    for (j = 4; j <= n; j += 2) {primes[j] = false;}

    while (lastP <= fsqr) {
        for (j = lastP*lastP; j <= n; j += 2*lastP) {primes[j] = false;}
        k = lastP + 2;
        ind = 2;
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        lastP += ind;
    }

    mpz_t test, jmpz, temp;
    mpz_init(test);
    mpz_init(jmpz);
    mpz_init(temp);

    myps.push_back(2);

    for (j = 3; j <= n; j += 2) {
        if (primes[j]) {
            mpz_set_ui(jmpz, j);
            mpz_set_ui(temp, j);
            mpz_sub_ui(temp, jmpz, 1);
            mpz_div_2exp(temp, temp, 1);
            mpz_powm(test,myN,temp,jmpz);
            if (mpz_cmp_ui(test, 1) == 0) {
                myps.push_back(j);
            }
        }
    }

    mpz_clear(jmpz); mpz_clear(temp); mpz_clear(test);
    return myps;
}

static inline unsigned long int positive_modulo (signed long int i,
                                                  signed long int n) 
{
    if (i < 0) {
        return (i % n + n);
    } else {
        return (i % n);
    }
}

static v3d SieveLists (signed long int facLim,
                       std::vector<signed long int> FBase,
                       signed long int vecLen,
                       signed long int myLow,
                       v2d sieveD)
{
    v3d outList(facLim, v2d(2, v1d()));
    unsigned long int tLen1, tLen2;
    signed long int i, j, modLow, vStrt1, vStrt2;

    for (i = 2; i < facLim; i++) {
        modLow = positive_modulo(myLow, FBase[i]);
        vStrt1 = ((sieveD[i][0] - myLow) % FBase[i]) + 1;
        tLen1 = (vecLen - vStrt1)/FBase[i];
        outList[i][0].reserve(tLen1);
        for (j = vStrt1; j <= vecLen; j += FBase[i]) {outList[i][0].push_back(j-1);}

        vStrt2 = ((sieveD[i][1] - myLow) % FBase[i]) + 1;
        tLen2 = (vecLen - vStrt2)/FBase[i];
        outList[i][1].reserve(tLen2);
        for (j = vStrt2; j <= vecLen; j += FBase[i]) {outList[i][1].push_back(j-1);}
    }

    return outList;
}

void quadraticSieve (mpz_t myNum, double fudge1,
                     double fudge2,
                     unsigned long int LenB) 
{
    unsigned long int digCount = mpz_sizeinbase(myNum, 10);
    unsigned long int bits = mpz_sizeinbase(myNum, 2);
    unsigned long int myMalloc = bits + 5;

    // The two values below are the slope and y-intercept
    // of the linear model obtained by running the
    // following in R:
    // DigSize <- c(4,10,15,20,23)
    // f_Pos <- c(0.5,0.25,0.15,0.1,0.05)
    // MSize <- c(5000,7000,10000,12500,15000)
    // LM2 <- lm(MSize ~ DigSize)
    // m2 <- summary(LM2)$coefficients[2,1]
    // b2 <- summary(LM2)$coefficients[1,1]
    double m2 = 524.0137221, b2 = 2354.20240137;

    // LM1 <- lm(f_Pos ~ DigSize)
    // m1 <- summary(LM1)$coefficients[2,1]
    // b1 <- summary(LM1)$coefficients[1,1]
    double m1 = -0.022384219554, b1 = 0.532332761578;

    unsigned long int myTarget, facSize;
    double sqrLogLog, LimB, lognum = bits/log2(exp(1));
    sqrLogLog = sqrt(lognum*log(lognum));
    mpz_t currP, nextP, resTest, CP1;
    mpz_init2(currP, myMalloc);
    mpz_init2(nextP, myMalloc);
    mpz_init2(CP1, myMalloc);
    mpz_init2(resTest, myMalloc);
    v1d facBase;

    if (digCount < 24) {
        if (fabs(fudge1) < 0.0001) {fudge1 = digCount*m1 + b1;}
        if (LenB == 0) {LenB = ceil(digCount*m2 + b2);}
        LimB = exp((.5+fudge1)*sqrLogLog);
        facBase = getPrimesQuadRes(myNum, LimB);
    } else if (digCount < 67) {
        // These values were obtained from "The Multiple Polynomial
        // Quadratic Sieve" by Robert D. Silverman
        // DigSize <- c(24,30,36,42,48,54,60,66)
        // FBSize <- c(100,200,400,900,1200,2000,3000,4500)
        // MSize <- c(5,25,25,50,100,250,350,500)
        //
        // rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
        // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
        // rawCoef
        //    intercept          x^1          x^2          x^3          x^4
        // 3637.0670996 -391.8275012   15.1541456   -0.2475566    0.0016806

        if (fudge1 == 0) {
            fudge1 = -0.4;
            LimB = exp((.5+fudge1)*sqrLogLog);

            myTarget = ceil(-391.8275012*digCount + 15.1541456*pow(digCount,2) -
                0.2475566*pow(digCount,3) + 0.0016806*pow(digCount,4) + 3637.0671);

            while (LimB < myTarget) {
                LimB = exp((.5+fudge1)*sqrLogLog);
                fudge1 += 0.001;
            }

            facBase = getPrimesQuadRes(myNum, LimB);
            facSize = facBase.size();

            while (facSize < myTarget) {
                fudge1 += 0.005;
                LimB = exp((.5+fudge1)*sqrLogLog);
                mpz_set_ui(currP, facBase.back());
                mpz_nextprime(nextP, currP);
                while (mpz_cmp_ui(nextP, LimB) < 0) {
                    mpz_set(currP, nextP);
                    mpz_nextprime(nextP, currP);
                    mpz_sub_ui(CP1, currP, 1);
                    mpz_div_2exp(CP1, CP1, 1);
                    mpz_powm(resTest,myNum,CP1,currP);
                    if (mpz_cmp_ui(resTest, 1) == 0) {
                        facBase.push_back(mpz_get_ui(currP));
                    }
                }
                facSize = facBase.size();
            }
        } else {
            LimB = exp((.5+fudge1)*sqrLogLog);
            facBase = getPrimesQuadRes(myNum, LimB);
        }
        // rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
        // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
        // rawCoef
        //    intercept           x^1           x^2           x^3           x^4
        // -1650.8252165   176.9043861    -6.7567603     0.1085362    -0.0005955

        if (LenB==0L) {
            LenB = 1000*ceil(176.9043861*digCount - 6.7567603*pow(digCount,2) +
                0.1085362*pow(digCount,3) - 0.0005955*pow(digCount,4) - 1650.8252165);
        }
    }

    facSize = facBase.size();

    mpz_t sqrtInt;
    mpz_init(sqrtInt);
    mpz_sqrt(sqrtInt, myNum);
    v1d myInterval;

    signed long int i, j;
    signed long int Lower = -1*LenB, Upper = LenB, LenB2 = 2*LenB+1;
    myInterval.reserve(LenB2);
    for (i = Lower; i <= Upper; i++) {myInterval.push_back(i);}

    v2d SieveDist(facSize, v1d(2));
    SieveDist[0][0] = 1;
    SieveDist[0][1] = 1;

    mpz_t TS[13];
    for (i = 0; i < 13; i++) {mpz_init2(TS[i], myMalloc);}
    unsigned long int pow2;
    signed long int iter1, iter2;
    mpz_set_ui(TS[12], 2);

    for (i = 1; i < facSize; i++) {
        mpz_set_ui(TS[5], facBase[i]);
        mpz_set_ui(TS[0], facBase[i]);
        mpz_sub_ui(TS[0], TS[0], 1);
        mpz_set(TS[1], TS[0]);
        mpz_set_ui(TS[7], 2);
        pow2 = mpz_scan1 (TS[1], 0);
        mpz_div_2exp (TS[1], TS[1], pow2);

        if (pow2 == 1) {
            mpz_add_ui (TS[4], TS[5], 1);
            mpz_div_2exp (TS[4], TS[4], 2);
            mpz_powm (TS[2], myNum, TS[4], TS[5]);
            mpz_neg (TS[4], TS[2]);
            mpz_mod (TS[3], TS[4], TS[5]);
        } else {
            mpz_div_2exp (TS[4], TS[0], 1);
            mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            while (mpz_cmp_ui(TS[6], 1) == 0) {
                mpz_add_ui(TS[7], TS[7], 1);
                mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            }

            mpz_add_ui(TS[4], TS[1], 1);
            mpz_div_2exp(TS[4], TS[4], 1);
            mpz_powm(TS[10], myNum, TS[4], TS[5]);
            mpz_powm(TS[8], myNum, TS[1], TS[5]);
            mpz_powm(TS[9], TS[7], TS[1], TS[5]);

            iter1 = pow2;
            iter2 = 1;
            mpz_mod(TS[11], TS[8], TS[5]);

            while ((mpz_cmp_ui(TS[11], 1) != 0) && (iter2 != 0)) {
                iter2 = 0;
                mpz_mod(TS[11], TS[8], TS[5]);
                while (mpz_cmp_ui(TS[11], 1) != 0) {
                    iter2++;
                    mpz_pow_ui(TS[4], TS[12], iter2);
                    mpz_powm(TS[11], TS[8], TS[4], TS[5]);
                }
                if (iter2 != 0) {
                    mpz_pow_ui(TS[4], TS[12], iter1-iter2-1);
                    mpz_powm(TS[4], TS[9], TS[4], TS[5]);
                    mpz_mul(TS[4], TS[4], TS[10]);
                    mpz_mod(TS[10], TS[4], TS[5]);

                    mpz_pow_ui(TS[4], TS[12], iter1-iter2);
                    mpz_powm(TS[9], TS[9], TS[4], TS[5]);

                    mpz_mul(TS[4], TS[8], TS[9]);
                    mpz_mod(TS[8], TS[4], TS[5]);
                    iter1 = iter2;
                }
                mpz_set_ui(TS[11], 0);
            }
            mpz_set(TS[2], TS[10]);
            mpz_sub(TS[4], TS[5], TS[10]);
            mpz_mod(TS[3], TS[4], TS[5]);
        }
        SieveDist[i][0] = mpz_get_si(TS[2]);
        SieveDist[i][1] = mpz_get_si(TS[3]);
    }

    for (i = 0; i < 13; i++) {mpz_clear(TS[i]);}

    v3d CoolList;
    CoolList = SieveLists(facSize, facBase, LenB2, Lower, SieveDist);

    mpz_t largeInterval[LenB2];
    mpz_t sqrDiff[LenB2];

    for (i = 0; i < LenB2; i++) {
        mpz_init2(largeInterval[i], myMalloc);
        mpz_init2(sqrDiff[i], myMalloc);
    }

    mpz_sub_ui(largeInterval[0], sqrtInt, Upper);
    for (i = 1; i < LenB2; i++) {mpz_add_ui(largeInterval[i], largeInterval[i-1], 1);}

    mpz_t temp;
    mpz_init2(temp, myMalloc);

    Rprintf("\nLenB = %d", LenB);

    for (i = 0; i < LenB2; i++) {
        mpz_pow_ui(temp, largeInterval[i], 2);
        mpz_sub(sqrDiff[i], temp, myNum);
    }

    mpz_mul_ui(temp, myNum, 2);
    mpz_sqrt(temp, temp);
    mpz_mul_ui(temp, temp, Upper);

    if (fudge2 == 0) {
        if (digCount < 8) {
            fudge2 = 0.1;
        } else if (digCount < 12) {
            fudge2 = 0.4;
        } else if (digCount < 30) {
            fudge2 = 0.7;
        } if (digCount < 50) {
            fudge2 = 0.9;
        } else {
            fudge2 = 1.5;
        }
    }

    double theCut = fudge2 * mpz_sizeinbase(temp, 10);
    mpz_t myPrimes[facSize];
    for (i = 0; i < facSize; i++) {mpz_init2(myPrimes[i], myMalloc);}
    std::vector<double> LnFB, myLogs(LenB2, 0);
    LnFB.reserve(facSize);
    std::vector<signed long int>::iterator it;

    j = 0;
    for (it = facBase.begin(); it < facBase.end(); it++) {
        mpz_set_ui(myPrimes[j], *it);
        LnFB.push_back(log(*it));
        j++;
    }

    unsigned long int tempSize, myIndex, minPrime;
    mpz_mul_ui(temp, sqrtInt, Upper);
    minPrime = mpz_sizeinbase(temp, 10) * 2;
    Rprintf("\nminPrimeforreal = %d", minPrime);
    minPrime = 2;

    for (i = 2; i < facSize; i++) {
        if (facBase[i] > minPrime) {
            tempSize = CoolList[i][0].size();
            for (j = 0; j < tempSize; j++) {myLogs[CoolList[i][0][j]] = myLogs[CoolList[i][0][j]] + LnFB[i];}
            tempSize = CoolList[i][1].size();
            for (j = 0; j < tempSize; j++) {myLogs[CoolList[i][1][j]] = myLogs[CoolList[i][1][j]] + LnFB[i];}
        }
    }

    v1d largeLogs;

    for (i = 0; i < LenB2; i++) {
        if (myLogs[i] > 23.91848) {
            largeLogs.push_back(i);
        }
    }

    unsigned long int largeLogsSize = largeLogs.size();
    mpz_t testInterval[largeLogsSize], newSqrDiff[largeLogsSize];

    for (i = 0; i < largeLogsSize; i++) {
        mpz_init2(testInterval[i], myMalloc);
        mpz_init2(newSqrDiff[i], myMalloc);
        mpz_set(testInterval[i], largeInterval[largeLogs[i]]);
        mpz_set(newSqrDiff[i], sqrDiff[largeLogs[i]]);
    }
    Rprintf("\nminP = %d", minPrime);
    Rprintf("\ntheCut = %f\n", theCut);
    bool GoForIt = false;
    v2d myMat(largeLogsSize, v1d(facSize+1, 0));
    i = 0;

    for (i = 0; i < largeLogsSize; i++) {
        if (mpz_cmp_ui(newSqrDiff[i], 0) < 0) {myMat[i][0] = 1;}
    }

    mpz_t rem, quot;
    mpz_init2(rem, myMalloc);
    mpz_init2(quot, myMalloc);
    v1d sFacs;

    if (largeLogsSize > 0) {
        GoForIt = true;
        bool divides = true;

        for (j = 0; j < largeLogsSize; j++) {
            for (i = 0; i < facSize; i++) {
                while (divides) {
                    mpz_fdiv_qr_ui(quot, rem, newSqrDiff[j], facBase[i]);
                    divides = (mpz_cmp_ui(rem, 0) == 0);
                    if (divides) {
                        mpz_set(newSqrDiff[j], quot);
                        myMat[j][i+1]++;
                    }
                }
                divides = true;
            }
            if (mpz_cmp_ui(newSqrDiff[j], 1) == 0) {
                sFacs.push_back(j);
                gmp_printf ("%Zd\n", sqrDiff[largeLogs[j]]);
            }
        }
    }

    unsigned long int lenM = sFacs.size();
    bool bSolution = false;

    if (lenM > 0) {
        mpz_t newTestInt[lenM];
        v2d newMat(lenM, v1d(facSize+1, 0));
        for (i = 0; i < lenM; i++) {
            mpz_init2(newTestInt[i], myMalloc);
            mpz_set(newTestInt[i], testInterval[sFacs[i]]);
            for (j = 0; j <= facSize; j++) {
                newMat[i][j] = myMat[sFacs[i]][j];
            }
        }
        bSolution = solutionSearch(newMat, myNum, myNum, facBase);
    }

    for (i = 0; i < largeLogsSize; i++) {
        mpz_clear(testInterval[i]);
        mpz_clear(newSqrDiff[i]);
    }

    for (i = 0; i < facSize; i++) {mpz_clear(myPrimes[i]);}
    for (i = 0; i < LenB2; i++) {
        mpz_clear(largeInterval[i]);
        mpz_clear(sqrDiff[i]);
    }

    mpz_clear(rem);
    mpz_clear(quot);
    mpz_clear(temp);
    mpz_clear(sqrtInt);
}

// testNum <- gmp::prod.bigz(gmp::nextprime(gmp::urand.bigz(2,35)))

// // 185081
