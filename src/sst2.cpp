#include <ctime>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


/* name spaces */
using namespace Rcpp;
using namespace std;


// Helpers
// Thanks: http://stackoverflow.com/questions/4003232/how-to-code-a-modulo-operator-in-c-c-obj-c-that-handles-negative-numbers
int mod (int a, int b) {
  if(b < 0) //you can check for b == 0 separately and do what you want
    return mod(-a, -b);
  int ret = a % b;
  if(ret < 0)
    ret += b;
  return ret;
}


/*
 * Find next neighbor in a wrapped world
 *
 * i: current cell
 * w: number of cells
 * dir: direction of movement
 */

int nei (int i, int w, int dir) {
  if (dir == 0) // stay
    return i;
  if (dir == 1) // east
    return mod((i % w) + 1, w) + w * (i / w);
  if (dir == 2) // west
    return mod((i % w) - 1, w) + w * (i / w);
  if (dir == 3) // north
    return (i % w) + w * mod(((i / w) + 1), w);
  if (dir == 4) // south
    return (i % w) + w * mod(((i / w) - 1), w);
  return -1;
}

/*
 * Find next neighbor, stop at boundary
 *
 * i: current cell
 * w: number of cells
 * dir: direction of movement
 */

int nei2 (int i, int w, int dir) {
  if (dir == 0) // stay
    return i;
  else {
    int x = i % w;
    int y = i / w;
    if (dir == 1) x++; // east
    if (dir == 2) x--; // west
    if (dir == 3) y++; // north
    if (dir == 4) y--; // south
    if (x >= 0 & y >= 0 & x < w & y < w) {
      return x + w * y;
    } else {
     // Rcout << "outside the box " << x << std::endl;
      return -1;
    }
  }
}

// http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
double frand(double upper) {
  double f = (double)rand() / RAND_MAX;
  return f * upper;
}

// http://stackoverflow.com/questions/1761626/weighted-random-numbers
int rsamp(vector<int>& nums, vector<double>& probs) {

  double sum_of_weight = 0, rnd;
  vector<int>::iterator it_num;
  vector<double>::iterator it_prob;

  // sum weights
  for (it_prob = probs.begin(); it_prob < probs.end(); it_prob++) {
    sum_of_weight += *it_prob;
  }

  rnd = frand(sum_of_weight);

  for(it_prob = probs.begin(), it_num = nums.begin(); it_prob < probs.end(); it_prob++, it_num++) {
    if(rnd < *it_prob)
      return *it_num;
    rnd -= *it_prob;
  }
}

int rsamp2(IntegerVector nums, NumericVector probs) {

  double sum_of_weight = sum(probs), rnd;
  rnd = frand(sum_of_weight);
  int n = probs.size();

  for (int i = 0; i < n; i++) {
    if (rnd < probs(i))
      return nums(i);
    rnd -= probs(i);
  }
}

// Transistion probabilities

//[[Rcpp::export]]
NumericMatrix tpm_func(double alpha, NumericVector omegas, NumericMatrix resources, NumericVector d, int nc) {
  int i, ii, j, k, n = nc * nc, l = omegas.length();
  double one;
  NumericMatrix tp(n, 5), num(5);
  IntegerVector p(5);

  for (i = 0; i < n; i++) {
    /*
     * get neighbour pixels (north = up)
     * 0 = stay
     * 1 = east
     * 2 = west
     * 3 = north
     * 4 = south
     */
    for (ii = 0; ii < 5; ii++)
      p(ii) = nei(i, nc, ii);

    for (j = 0; j < 5; j++) {
      one = -alpha * d[j];
      for (int k = 0; k < l; k++) {
        one += resources(p[j], k) * omegas(k);
      }
      num[j] = exp(one);
    }
    tp(i, _) = num / sum(num);
  }
  return tp;
}

// Wrapped world
//[[Rcpp::export]]
List walk(NumericMatrix tpm, int n, NumericVector xy0, int nc, NumericVector dp, int init_dir = -1, int boundary = 1, int max_try = 100) {

  NumericVector dp_stay = NumericVector::create(1, 1, 1, 1, 1);

  NumericVector dp_e = NumericVector::create(dp[1], dp[0], dp[2], dp[1], dp[1]);
  NumericVector dp_w = NumericVector::create(dp[1], dp[2], dp[0], dp[1], dp[2]);

  NumericVector dp_s = NumericVector::create(dp[1], dp[1], dp[1], dp[2], dp[0]);
  NumericVector dp_n = NumericVector::create(dp[1], dp[1], dp[1], dp[0], dp[2]);


  IntegerVector c(n), w, ch = Range(0, 4); // choices for the new yx
  int i, k, t;
  c(0) = xy0(0) + xy0(1) * nc;

  for (i = 1; i < n; i++) {
    k = c(i-1);

    // find next cell
    if (i == 1) {
      if (init_dir == -1) {
        w = rsamp2(ch, tpm(k, _)); // weights of the neighbours -> multiply by directional persistence
      } else {
        w = NumericVector::create(init_dir);
      }
    }

    int wlast = w(0);
    if (wlast == 0) {
      w = rsamp2(ch, tpm(k, _) * dp_stay); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 1) {
      w = rsamp2(ch, tpm(k, _) * dp_e); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 2) {
      w = rsamp2(ch, tpm(k, _) * dp_w); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 3) {
      w = rsamp2(ch, tpm(k, _) * dp_n); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 4) {
      w = rsamp2(ch, tpm(k, _) * dp_s); // weights of the neighbours -> multiply by directional persistence
    }
    /*
     * Different boundary conditions:
     * 1: Wrapped world
     * 2: Stop if animal steps out of boundary
     * 3: Reflective boundaries
     */

    if (boundary == 1) {
      // Wrapped
      c(i) = nei(k, nc, w(0));
    } else if (boundary == 2) {
      // Stop if animal hits boundary

      int next = nei2(k, nc, w(0));
      if (next == -1) { // nei2 return -1 if a step is outside the boundary
        return List::create(Named("steps") = c,
                            Named("last_step") = i,
                            Named("outcome") = -1);
      } else {
        c(i) = next;
      }
    } else if (boundary == 3) {
      // Reflective boundary
      int next = nei2(k, nc, w(0));
      while (next == -1 & t < max_try) {// nei2 return -1 if a step is outside the boundary
        w = rsamp2(ch, tpm(k, _)); // weights of the neighbours
        next = nei2(k, nc, w(0));
        t++;
      }
      c(i) = next;
      t = 0;
    }
  }
  return List::create(Named("steps") = c,
                      Named("outcome") = 1);
}


//[[Rcpp::export]]
List ud_func(NumericMatrix tpm, int n, NumericVector xy0, int nc, int burnin, NumericVector dp, int boundary = 1, int init_dir = -1,  int max_try = 100) {

  IntegerVector w, ch = Range(0, 4); // choices for the new yx

  IntegerVector ud(tpm.nrow(), 0);

  // first cell
  int k = xy0(0) + (nc - xy0(1)) * nc; // index for neighbouring cell
  int next, t;
  ud(k)++;

  // directional persistence
  NumericVector dp_stay = NumericVector::create(1, 1, 1, 1, 1);

  NumericVector dp_e = NumericVector::create(dp[1], dp[0], dp[2], dp[1], dp[1]);
  NumericVector dp_w = NumericVector::create(dp[1], dp[2], dp[0], dp[1], dp[2]);

  NumericVector dp_s = NumericVector::create(dp[1], dp[1], dp[1], dp[2], dp[0]);
  NumericVector dp_n = NumericVector::create(dp[1], dp[1], dp[1], dp[0], dp[2]);

  for (int i = 1; i < (n + burnin); i++) {
    // w = rsamp2(ch, tpm(k, _)); // weights of the neighbours
    // 2017-04-01: Removed normalization in favor for speed
    // w = RcppArmadillo::sample(ch, 1, true, (NumericVector)(tpm(k, _) / sum(tpm(k, _)))); // weights of the neighbours


    // find next cell
    if (i == 1) {
      if (init_dir == -1) {
        w = rsamp2(ch, tpm(k, _)); // weights of the neighbours -> multiply by directional persistence
      } else {
        w = NumericVector::create(init_dir);
      }
    }

    int wlast = w(0);
    if (wlast == 0) {
      w = rsamp2(ch, tpm(k, _) * dp_stay); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 1) {
      w = rsamp2(ch, tpm(k, _) * dp_e); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 2) {
      w = rsamp2(ch, tpm(k, _) * dp_w); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 3) {
      w = rsamp2(ch, tpm(k, _) * dp_n); // weights of the neighbours -> multiply by directional persistence
    } else if (wlast == 4) {
      w = rsamp2(ch, tpm(k, _) * dp_s); // weights of the neighbours -> multiply by directional persistence
    }
    /*
     * Different boundary conditions:
     * 1: Wrapped world
     * 2: Stop if animal steps out of boundary
     * 3: Reflective boundaries
     */

    if (boundary == 1) {
      // Wrapped
      next = nei(k, nc, w(0));
    } else if (boundary == 2) {
      // Stop if animal hits boundary

      next = nei2(k, nc, w(0));
      if (next == -1) { // nei2 return -1 if a step is outside the boundary
        return List::create(Named("steps") = ud,
                            Named("last_step") = i,
                            Named("outcome") = -1);
      }
    } else if (boundary == 3) {
      // Reflective boundary
      next = nei2(k, nc, w(0));
      while (next == -1 & t < max_try) {// nei2 return -1 if a step is outside the boundary
        w = rsamp2(ch, tpm(k, _)); // weights of the neighbours
        next = nei2(k, nc, w(0));
        t++;
      }
      t = 0;
    }

    if (i > burnin) {
      ud(next)++;
    }

    k = next;
  }

  return List::create(
    Named("ud") = ud,
    Named("lastPos") = k);
}


