#include <cmath>
#include <cassert>
#include <map>
#include "statistics.h"

//	Euler's number
const double e = exp(1.0);

// Map containing z-scores and corresponding percentiles
const std::map<double, double> z_percentiles = { {-3.09, 0.001}, {-3.08, 0.001}, {-3.07, 0.0011}, {-3.06, 0.0011}, {-3.05, 0.0011},
	{-3.04, 0.0012}, {-3.03, 0.0012}, {-3.02, 0.0013}, {-3.01, 0.0013}, {-3.00, 0.0013},
	{-2.99, 0.0014}, {-2.98, 0.0014}, {-2.97, 0.0015}, {-2.96, 0.0015}, {-2.95, 0.0016},
	{-2.94, 0.0016}, {-2.93, 0.0017}, {-2.92, 0.0018}, {-2.91, 0.0018}, {-2.90, 0.0019},
	{-2.89, 0.0019}
};

/**
*	Calculates the factorial of an integer.
* 
*	If the integer is greater than 20, return -1, indicating that the number is too large to store its factorial.
* 
*	@param number The integer whose factorial is to be calculated.
*	@return The factorial of the integer, or -1 if it is too large.
*/
unsigned long long factorial(unsigned int number)
{
	if (number == 0 || number == 1)
		return 1;

	if (number > 20)
	{
		// number too large to compute factorial
		return -1;
	}

	unsigned long long result = 1;
	for (int i = 2; i <= number; i++)
		result *= i;
	return result;
}

/**
*	Calculates the binomial coefficient, i.e the number of n ways to choose unordered collections from k elements.
* 
*	@param n Number of ways to choose.
*	@param k Total number of elements.
*	@return The binomial coefficient, or -1 if either n or k is larger than 20, as we cannot compute that factorial.
*/
unsigned long long choose(int n, int k)
{
	if (n > 20 || k > 20)
		return -1;	// number too large to compute binomial coefficient
	return factorial(n) / (factorial(k) * factorial(n - k));
}






/**
*	Calculates the standard score (z-score) of a data point.
*	
*	@param datapoint The observation.
*	@param mean The population mean.
*	@param StdDev The population standard deviation.
*	@return The standard score, i.e, the number of standard deviations the observed data point deviates from the mean by.
*/
inline double z(double datapoint, double mean, double StdDev)
{
	return ((datapoint - mean) / StdDev);
}

/**
*	Calculates the CDF of a Gaussian distribution.
* 
*	@param z The z-score.
*	@return The cumulative distribution function.
*/
inline double standard_normal_cdf(double z)
{
	return 0.5 * (1 + erf(z / sqrt(2)));
}

/**
*	Calculates the nth percentile of a set of normally distributed data.
* 
*	@param container The pointer to the vector containing the set.
*	@param percent The percentile the caller is interested in.
*	@return The data point that is the nth percentile.
*/
inline double nth_percentile(std::vector<int>& container, double percent)
{
	// unimplemented

}

/**
*	Calculates the percentile of a specific data point in a standard normal distribution.
* 
*	@param container Pointer to the vector containing the set of data.
*	@return The percentile the data point lies in.
*/
inline double z_to_percentile(std::vector<int>& container, double z)
{
	return standard_normal_cdf(z) * 100;
}

/**
*	Calculates the mean of a set.
*
*	@param container Container that contains the set.
*	@return The mean.
*/
inline double mean(std::vector<int>& container)
{
	long sum = 0;
	for (int i : container)
		sum += i;
	return sum / container.size();
}

/**
*	Calculates the standard deviation of a set.
* 
*	@param container Container that contains the set.
*	@return The standard deviation.
*/
inline double StdDev(std::vector<int>& container)
{
	double average = mean(container);
	long sum = 0;
	for (int i : container)
	{
		sum += pow(i - average, 2);
	}
	return sqrt(sum / container.size() - 1);
}

/**
*	Calculates the variance of a set.
*
*	@param container Container that contains the set.
*	@return The variance.
*/
inline double variance(std::vector<int>& container)
{
	return pow(StdDev(container), 2);
}



/**
*	Calculates the PMF of a uniform discrete random variable.
*
*	@param a The minimum value in the range of values.
*	@param b The maximum value in the range of values.
*	@return The probability of selecting any given value in the range of values.
*/
inline double uni_pmf(int a, int b, int x)
{
	if (x < a || x > b)
		return 0.0;
	else
		return 1.0 / (b - a + 1.0);
}

/**
*	Calculates the expected value of a uniform discrete random variable.
* 
*	@param a The minimum value in the range of values.
*  	@param b The maximum value in the range of values.
*	@return The expected value.
*/
inline double uni_mean(int a, int b)
{
	return 0.5 * (a + b);
}

/**
*	Calculates the variance of a uniform discrete random variable.
* 
*	@param a The minimum value in the range of values.
*  	@param b The maximum value in the range of values.
*	@return The variance.
*/
inline double uni_variance(int a, int b)
{
	return ((pow(b - a + 1, 2) - 1) / 12);
}

/**
*	Calculates the PMF of a Bernoulli random variable.
* 
*	This function only functions correctly when parameter k is 0 or 1.
* 
*	@param p The probability of success on a single trial.
*	@param k The event we are testing for. 0 represents failure, 1 represents success.
*	@return The probability of parameter k, whether failure or success.
*/
inline double bern_pmf(double p, unsigned int k)
{
	assert(k == 0 || k == 1);
	if (k == 0)
		return 1 - p;
	else
		return p;
}

/**
*	Calculates the mean of a Bernoulli random variable.
* 
*	@param p The probability of success on a single trial.
*	@return The expected value (aka mean).
*/
inline double bern_mean(double p)
{
	return p;
}

/**
*	Calculates the variance of a Bernoulli random variable.
*	@param p The probability of success on a single trial.
*	@return The variance.
*/
inline double bern_variance(double p)
{
	return p * (1 - p);
}

/**
*	Calculates the PMF of a binomial random variable.
* 
*	If either n or k is very large or p is very small, approximate this RV with a Poisson RV.
* 
*	@param p The probability of success on a single trial.
*	@param n The number of trials.
*	@param k The number of successes being tested for.
*	@return The probability of observing k successes in n trials.
*/
inline double bin_pmf(unsigned int n, double p, unsigned int k)
{
	if (n * p < 5 || n > 20 || k > 20)
	{
		return (choose(n, k) * pow(p, k) * pow(1 - p, n - k));
	}
	else {
		return pois_pmf(n * p, k);
	}
}

/**
*	Calculates the mean/expected value of a binomial random variable.
* 
*	@param p The probability of success on a single trial.
*	@param n The number of trials.
*	@return The mean/expected value.
*/
inline double bin_mean(unsigned int n, double p)
{
	return n * p;
}

/**
*	Calculates the variance of a binomial random variable.
* 
*	@param p The probability of success on a single trial.
*	@param n The number of trials.
*	@return The variance.
*/
inline double bin_variance(unsigned int n, double p)
{
	return n * p * (1 - p);
}

/**
*	Calculates the PMF of a geometric random variable.
* 
*	@param p The probability of success on one trial.
*	@param k Number of trials being tested for.
*	@return The probability of needing k trials to observe the first success.
*/
inline double geo_pmf(double p, unsigned int k)
{
	return pow(1 - p, k - 1);
}

/**
*	Calculates the mean/expected value of a geometric random variable.
* 
*	@param p Probability of success on a single trial.
*	@return The expected value of number of successes needed to achieve first success.
*/
inline double geo_mean(double p)
{
	return 1 / p;
}

/**
*	Calculates the variance of a geometric random variable.
* 
*	@param p Probability of success on a single trial.
*	@return The variance.
*/
inline double geo_variance(double p)
{
	return ((1 - p) / pow(p, 2));
}

/**
*	Calculates the PMF of a negative binomial/pascal random variable.
* 
*	Uses Poisson distribution to approximate it if either r or k are too large ( > 21).
* 
*	@param r Number of successes we are testing for.
*	@param k Number of trials.
*	@param p Probability of success on a single trial.
*	@return Probability of needing k trials to achieve r successes.
*/
inline double nb_pmf(int r, int k, double p)
{
	if (r > 20 || k > 20)
	{
		//	approximate using Poisson distribution, as calculating the PMF directly is too computationally expensive
		return pois_pmf(r * (1 - p) / p, k);
	}
	else {
		return (choose(k - 1, r - 1) * pow(p, r) * pow(1 - p, k - r));
	}
}

/**
*	Calculates the mean of a negative binomial random variable.
* 
*	@param r Number of successes.
*	@param p Probability of success on a single trial.
*	@return Average number of trials needed to achieve r successes.
*/
inline double nb_mean(int r, double p)
{
	return r / p;
}

/**
*	Calculates the variance of a negative binomial random variable.
* 
*	@param r Number of successes.
*	@param p Probability of success on a single trial.
*	@return The variance.
*/
inline double nb_variance(int r, double p)
{
	return (r * (1 - p) / pow(p, 2));
}

/**
*	Calculates the PMF of a hypergeometric random variable.
*	
*	Uses binomial distribution to approximate if the sample size is large enough.
* 
*	@param N Population size.
*	@param M Number of successes in population.
*	@param n Sample size.
*	@param k Number of successes whose PMF is to be calculated.
*	@return Probability of drawing k successes with n total trials from a population with N size and M successes.
*/
inline double hg_pmf(int N, int M, int n, int k)
{
	assert(n <= N && k <= n);
	if (n / N <= 0.05)
	{
		// Use binomial distribution to approximate.
		return M / N;
	}
	else {
		return (choose(M, k) * choose(N - M, n - k)) / choose(N, n);
	}
}

/**
*	Calculates the expected value/mean of a hypergeometric random variable.
* 
*	@param N Population size.
*	@param M Number of successes in population.
*	@param n Sample size.
*	@return The average number of successes in the sample.
*/
inline double hg_mean(int N, int M, int n)
{
	return (n * M / N);
}

/**
*	Calculates the variance of a hypergeometric random variable.
* 
*	@param N Population size.
*	@param M Number of successes in population.
*	@param n Sample size.
*	@return The variance of the number of successes in the sample.
*/
inline double hg_variance(int N, int M, int n)
{
	return ((n * M / N) * (1 - M / N) * ((N - n) / (N - 1)));
}

/**
*	Calculates the PMF of a Poisson random variable.
* 
*	This function may also be used to approximate a Binomial distribution under certain conditions. 
*	See documentation of bin_pmf function.
* 
*	@param mean The mean number of occurrences in a time frame.
*	@param successes Number of successes whose PMF is being calculated.
*	@return Probability of observing the specified number of successes.
*/
inline double pois_pmf(double mean, int k)
{
	return (pow(mean, k) * pow(e, -mean)) / factorial(k);
}






//	Constructor for Uni class. Sets minimum and maximum values.
Uni::Uni(int a, int b)
{
	this->a = a;
	this->b = b;
}

//	PMF, mean, and variance of this uniform discrete RV.
double Uni::pmf(int k)
{
	return uni_pmf(this->a, this->b, k);
}

double Uni::ExpectedValue()
{
	return uni_mean(this->a, this->b);
}

double Uni::variance()
{
	return uni_variance(this->a, this->b);
}







//	Default constructor for Bern class. Sets probability to -1.
Bern::Bern() 
{
	this->p = -1.0;
}

// Constructor for Bern class. Sets the probability.
Bern::Bern(double p)
{
	assert(p >= 0);
	this->p = p;
}

//	Sets the probability of this Bernoulli random variable.
void Bern::setProbability(double p)
{
	this->p = p;
}

//	Calculates the PMF, expected value, and variance for this Bernoulli RV.
double Bern::pmf(int k)
{
	assert(k >= 0);
	return bern_pmf(this->p, k);
}

double Bern::ExpectedValue()
{
	return bern_mean(this->p);
}

double Bern::variance()
{
	return bern_variance(this->p);
}





// Constructor for Bin class. Sets the probability.
Bin::Bin(int n, double p)
{
	this->n = n;
}

// Computes PMF of this binomial RV.
double Bin::pmf(int k)
{
	return bin_pmf(n, this->p, k);
}

// Calculates the mean of this binomial RV.
double Bin::ExpectedValue()
{
	return bin_mean(this->n, this->p);
}

// Calculates the variance of this binomial RV.
double Bin::variance()
{
	return bin_variance(this->n, this->p);
}





// Constructor for Geo class. Sets the probability.
Geo::Geo(double p)
{
	this->p = p;
}

// Computes PMF of geometric RV.
double Geo::pmf(int k)
{
	return geo_pmf(this->p, k);
}

// Computes mean of this geometric RV.
double Geo::ExpectedValue()
{
	return geo_mean(this->p);
}

// Computes the variance of this geometric RV.
double Geo::variance()
{
	return geo_variance(this->p);
}






// Constructor for NB class. Sets the number of successes and probability.
NB::NB(int r, double p)
{
	this->r = r;
	this->p = p;
}

// Constructor for NB class with only success parameter.
NB::NB(int r)
{
	this->r = r;
}

// Calculates the PMF for this negative binomial RV.
double NB::pmf(int k)
{
	return nb_pmf(this->r, k, this->p);
}

// Calculates the mean for this negative binomial RV.
double NB::ExpectedValue()
{
	return nb_mean(this->r, this->p);
}

// Calculates the variance for this negative binomial RV.
double NB::variance()
{
	return nb_variance(this->r, this->p);
}






// Constructor for HG class. Sets the population size, number of successes, and the sample size.
HG::HG(int N, int M, int n)
{
	this->N = n;
	this->M = M;
	this->n = n;
}

// Calculates the PMF of this hypergeometric RV.
double HG::pmf(int k)
{
	return hg_pmf(this->N, this->M, this->n, k);
}

// Calculates the mean of this hypergeometric RV.
double HG::ExpectedValue()
{
	return hg_mean(this->N, this->M, this->n);
}

// Calculates the variance of this hypergeometric RV.
double HG::variance()
{
	return hg_variance(this->N, this->M, this->n);
}






//	Constructor for Pois class. Sets the mean.
Pois::Pois(double mean)
{
	assert(mean >= 0.0);
	this->mean = mean;
}

//	Overloading constructor for Pois class. Takes time frame and average number of occurrences to calculate the mean.
Pois::Pois(unsigned int occurrences, double unitTimeFrame)
{
	this->mean = occurrences * unitTimeFrame;
}

// Calculates the PMF for this Poisson RV.
double Pois::pmf(int k)
{
	assert(k >= 0);
	return pois_pmf(this->mean, k);
}

double Pois::ExpectedValue()
{
	return this->mean;
}

double Pois::variance()
{
	return this->mean;
}
