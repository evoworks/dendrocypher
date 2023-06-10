
#ifndef TOOLS_H_
#define TOOLS_H_

#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <ctime>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
using namespace std;


namespace tools
{
   // Math stuff ..
   int gcd(int a, int b);
   int factorial(int n);
   double sterling_factorial(int n);                                        // Approximate the factorial function with Sterling's formulas; fairly good approximation for larger values of "n"
   double double_rand(double min, double max);                              // Returns a pseudorandom double precision number within the specified range
   double rand_gaussian(double mu=0.0, double sigma=1.0);                   // Returns a pseudorandom double precision number that follows a normal distribution
   void int_rand_series(vector<int>* vec, int size, int min, int max);      // Sets size of vector and fills it with random integers sampled WITHOUT from the interval (min,max)

   string intersection(string a, string b);
   string my_union(string a, string b);
   string unique_dna(string str);

   double compute_mean(vector<double>* vec);   // Computes mean of an vector of doubles
   double compute_sd(vector<double>* vec);     // Computes standard deviation of values in a vector of double

   // probability distributions:
   //    - these are probability density functions
   //    - take a random variable as input and return a probability
   double binomial(int n, double p, int k);                               // IMP: Add some notes ...
   //double normal_apprx_binomial( int k);                                // k will be converted to k-.5 < x < k+.5
   double poisson_apprx_binomial(int n, double p, int k);                 // Typically use this when n is large (>50) and p is small (such that n*p < 10)
   double normal_PDF(double x, double mean=0.0, double sigma=1.0);        // Returns the *height* of the probability distribution at the value of x
   double p_uniform(double x, double lower, double upper);                // Returns probability of a uniform random variable x over the interval lower to upper

   //misc string manipulations and copies
   void print_string_vector(vector<string> *my_vector);               // just print the contents to the screen
   void str_tok(string my_str, vector<string> *my_vector);            // breaks text within "my_str" into a stream of tokens (space delimited) which are then stored as separate elements within a vector
   void find_replace(string &str, const string &find_this, const string &replacement); // a method for "find and replace" within a string
   void find_erase(string &str, const string &find_this);             // a method for finding all instances of a substring within a string and erasing it
   void remove_whitespace(string &input);                             // remove all white space from a string
   void copy(vector<double>* copy , const vector<double>* primary);   // generic vector copy method; use to copy site likelihoods from likelihood engine to local vector
   void copy(vector<bool>* copy , const vector<bool>* primary);       // generic boolean vector copy method; use to copy optimize_this flags
   void copy(vector<string>* copy , const vector<string>* primary);   // generic string vector copy method
   string string_tolower(string str);                                 // "for loop" method for converting a string to all LOWERcase characters; see WARNING in .cpp code
   string STL_string_tolower(string str);                             // STL algorithm method for converting a string to all lower case characters; same warnings
   string string_toupper(string str);                                 // "for loop" method for converting a string to all UPPERcase characters; see WARNING in .cpp code
   string get_current_date();                                         // returns a string holding the current date in format: day/month/year

   double str_to_double(string str);                                  // returns double; if conversion fails it returns 0
   int str_to_int(string str);                                        // returns int; if conversion fails it returns 0
   string int_to_str(int number);                                     // returns string; if conversion fails it returns "fail"; NOTE: C++11 now has std::to_string

   string parse_keyword(string input);                                // parse keyword from keyword(#,#,#)
   void parse_args(string input, vector<string> *my_vector);          // parse arguments from between () and load them into args vector
   void expand_series(vector<string> *a_vec);                         // given a vector with  a "-" in the string at index=0, it expands series; e.g., start with 1-4 in a single element of the vector and end with 1, 2, 3, 4 in separate elements of the array
   int string_to_bool(string str);                                    // if string=="false", return 0; if string=="true", return 1; any other string, return -1

   //big stuff
   double log_factorial(int n);                                       // To avoid factorial you can work in the log domain; hence this function return the log(n!) rather than n!
}

#endif /* TOOLS_H_ */
