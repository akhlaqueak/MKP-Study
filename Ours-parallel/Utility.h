#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <tuple>
#include <queue>
#include <set>
#include <cassert>

#define TIME_NOW std::chrono::steady_clock::now()
#define NDEBUG
using namespace std;
//#define NO_TRUSS_PRUNE // This controls truss pruning in preprocessing, not in BBSearch

using ui = unsigned int; // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges
// typedef std::vector<ui> vecui ;
using vecui = std::vector<ui>;
#define pb push_back
#define mp make_pair

#define mmax(a,b) ((a)>(b)?(a):(b))
#define mmin(a,b) ((a)<(b)?(a):(b))

class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<ui> sequence;
		if(number == 0) sequence.pb(0);
		while(number > 0) {
			sequence.pb(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(ui i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}
};

class MBitSet
{
private:
public:
	ui n;
	ui cap;
	ui *buf;

	MBitSet()
	{
		buf = nullptr;
		cap = n = 0;
	}
	MBitSet(ui _cap)
	{
		cap = _cap;
		n = (cap >> 5) + 1;
		buf = new ui[n];
		fill(buf, buf + n, 0);
		// for (ui i = 0; i < n; ++i)
		// 	buf[i] = 0;
	}
	void init(ui _cap)
	{
		cap = _cap;
		n = (cap >> 5) + 1;
		buf = new ui[n];
		fill(buf, buf + n, 0);
	}

	// MBitSet(const MBitSet& src){
	// 	cap = src.cap;
	// 	n = src.n;
	// 	buf = src.buf;
	// 	fill(buf, buf + n, 0);
	// }
	void nullify(){

		buf=nullptr;
	}
	~MBitSet()
	{
		// todo doing double free, see what causing it and fix

		if (buf != nullptr){
			delete[] buf;
			buf = nullptr;
		}
	}
	void reset()
	{
		fill(buf, buf + n, 0);
	}
	void reset(ui n)
	{
		fill(buf, buf + (n>>5)+1, 0);
	}
	void setAll()
	{
		fill(buf, buf + n, 0xffffffff);
	}
	// FLIP all the bits
	void flip()
	{
		for (ui i = 0; i < n; ++i)
			buf[i] = ~buf[i];
	}
	void set(ui x)
	{
		// assert(x < cap);
		buf[x >> 5] |= (ui)1 << (x & 31);
	}

	bool test(ui x)
	{
		// cout << x << " " << n << " " << cap << endl;
		return buf[x >> 5] >> (x & 31) & 1;
	}

	bool empty()
	{
		for (ui i = 0; i < n; ++i)
			if (buf[i])
				return false;
		return true;
	}
	void setup(vecui* adj, ui m)
	{
		m = (m >> 5) + 1;
		fill(buf, buf + m, 0);
		
		for (ui u : *adj)
			set(u);	
	}
};
// class MBitSet
// {
// private:
// public:
// 	ui n;
// 	ui cap;
// 	char *buf;

// 	MBitSet()
// 	{
// 		buf = nullptr;
// 		cap = n = 0;
// 	}
// 	MBitSet(ui _cap)
// 	{
// 		cap = _cap;
// 		n = cap;
// 		buf = new char[n];
// 		fill(buf, buf + n, 0);
// 		// for (ui i = 0; i < n; ++i)
// 		// 	buf[i] = 0;
// 	}
// 	void init(ui _cap)
// 	{
// 		cap = _cap;
// 		n = cap;
// 		buf = new char[n];
// 		fill(buf, buf + n, 0);
// 	}
// 	~MBitSet()
// 	{
// 		// todo doing double free, see what causing it and fix

// 		if (buf != nullptr)
// 			// delete[] buf;
// 			buf = nullptr;
// 	}
// 	void reset()
// 	{
// 		fill(buf, buf + n, 0);
// 	}
// 	void reset(ui n)
// 	{
// 		fill(buf, buf + n, 0);
// 	}
// 	void setAll()
// 	{
// 		fill(buf, buf + n, 0xff);
// 	}
// 	// FLIP all the bits
// 	void flip()
// 	{
// 		for (ui i = 0; i < n; ++i)
// 			buf[i] = ~buf[i];
// 	}
// 	void set(ui x)
// 	{
// 		// assert(x < cap);
// 		buf[x] = 1;
// 	}

// 	bool test(ui x)
// 	{
// 		// cout << x << " " << n << " " << cap << endl;
// 		return buf[x];
// 	}

// 	bool empty()
// 	{
// 		for (ui i = 0; i < n; ++i)
// 			if (buf[i])
// 				return false;
// 		return true;
// 	}
// 	void setup(vecui &adj, ui m)
// 	{
// 		fill(buf, buf + m, 0);
		
// 		for (ui u : adj)
// 			set(u);	
// 	}
// };
#endif
