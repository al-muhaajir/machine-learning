#pragma once
#include <vector>
#ifndef STATISTICS_H
#define STATISTICS_H

// gaussian statistics
inline double z(double datapoint, double mean, double StdDev);
inline double standard_normal_cdf(double z);

inline double nth_percentile(std::vector<int>& container, double percent);
inline double z_to_percentile(std::vector<int>& container, double z);
inline double mean(std::vector<int>& container);
inline double median(std::vector<int>& container);
inline double mode(std::vector<int>& container);
inline double StdDev(std::vector<int>& container);
inline double variance(std::vector<int>& container);

// normal random variables
inline double n_pmf(double x, double mean, double variance);

// discrete random variables
inline double uni_pmf(int a, int b, int x);
inline double uni_mean(int a, int b);
inline double uni_variance(int a, int b);

inline double bern_pmf(double p, unsigned int k);
inline double bern_mean(double p);
inline double bern_variance(double p);

inline double bin_pmf(unsigned int n, double p, unsigned int k);
inline double bin_mean(unsigned int n, double p);
inline double bin_variance(unsigned int n, double p);

inline double geo_pmf(double p, unsigned int k);
inline double geo_mean(double p);
inline double geo_variance(double p);

inline double nb_pmf(int r, int k, double p);
inline double nb_mean(int r, double p);
inline double nb_variance(int r, double p);

inline double hg_pmf(int N, int M, int n, int k);
inline double hg_mean(int N, int M, int n);
inline double hg_variance(int N, int M, int n);

inline double pois_pmf(double mean, int k);



/*
*	Models a normal random variable.
*/
class N
{
private:
	double mean;
	double variance;
public:
	N(double mean, double variance);
	double pmf(int x);
};

/*
*	Abstract class that models a discrete random variable.
*/
class DiscreteRandomVariable
{
public:
	virtual double pmf(int k) = 0;
	virtual double ExpectedValue() = 0;
	virtual double variance() = 0;
};

/*
*	Models a uniform discrete random variable.
*/
class Uni : public DiscreteRandomVariable
{
private:
	int a;
	int b;
public:
	Uni(int a, int b);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a Bernoulli random variable.
*/
class Bern : public DiscreteRandomVariable
{
protected:
	double p;
public:
	Bern();
	Bern(double p);
	void setProbability(double p);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a binomial random variable.
*/
class Bin : public Bern
{
private:
	int n;
public:
	Bin(int n, double p);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a geometric random variable.
*/
class Geo : public Bern
{
public:
	Geo(double p);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a negative binomial random variable.
*/
class NB : public Bern
{
private:
	int r;
public:
	NB(int r, double p);
	NB(int r);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a geometric random variable.
*/
class HG : public DiscreteRandomVariable
{
private:
	int N;
	int M;
	int n;
public:
	HG(int N, int M, int n);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Models a Poisson random variable.
*/
class Pois : public DiscreteRandomVariable
{
private:
	double mean;
public:
	Pois(double mean);
	Pois(unsigned int occurrences, double unitTimeFrame);
	double pmf(int k) override;
	double ExpectedValue() override;
	double variance() override;
};

/*
*	Abstract class that models a continuous random variable.
*/
class ContinuousRandomVariable
{

};

#endif
