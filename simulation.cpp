
#ifdef _MSC_VER
#	include "StdAfx.h"
#endif

#include <boost/pool/pool_alloc.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

#include <set>
#include <vector>
#include <map>

#include <algorithm>

#include <cmath>
#include <ctime>

template <class T, class ForwardIterator>
T Average(ForwardIterator begin, ForwardIterator end)
{
	T value = T();
	int len = 0;
	
	for (; begin != end; ++begin)
		value += *begin, ++len;
	
	return value / len;
}

template <class T, class ForwardIterator>
T Variance(ForwardIterator begin, ForwardIterator end, T average)
{
	T value = T();
	int len = 0;
	
	for (; begin != end; ++begin)
		value += (*begin - average)*(*begin - average), ++len;
	
	return value / len;
}

struct RandomGenerator {
	int operator () (int max) {
		return rand() % max;
	}
} RandomGenerator;

class Node;

typedef std::set< Node *, std::less<Node *>, boost::fast_pool_allocator<Node *> > NodeSet;


typedef enum Action {
	Increase,
	Decrease
};

class Node
{
public:
	double capacity, link_capacity;
	int N, window_size, steps, last_update, last_wakeup, last_changed;
	int window_size_index, new_window_size, max_window_size;
	Action last_action;
	
	NodeSet upload_peers, download_peers;
	
	double total;
	double avg, last_avg;
	
	Node(double capacity, int initN, int sim_window_size) : capacity(capacity), link_capacity(capacity), N(initN), window_size(sim_window_size), steps(0), last_update(0), last_wakeup(0), last_changed(0), last_action(Decrease), total(0.0), avg(0.0), last_avg(0.0) 
		{
			window_size_index = 0;
   			 new_window_size = sim_window_size;
   			 max_window_size = sim_window_size;
		}
};

struct LinkGreater {
	bool operator () (Node *lhs, Node *rhs) {
		return lhs->link_capacity > rhs->link_capacity;
	}
} LinkGreater;

template <typename ContainerT>
void make_union(const NodeSet &a, const NodeSet &b, ContainerT &result) {
	std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.begin()));
}

template <typename ContainerT>
void make_intersection(const NodeSet &a, const NodeSet &b, ContainerT &result) {
	std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.begin()));
}

template <typename ContainerT>
void make_difference(const NodeSet &a, const NodeSet &b, ContainerT &result) {
	std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.begin()));
}

double download_rate(const NodeSet &peers) {
	double rate = 0.0;
	
	for (NodeSet::iterator it = peers.begin(); it != peers.end(); ++it)
		rate += (*it)->link_capacity;
	
	return rate;
}

typedef int (*init_func_t)(double c, int n);

class Simulation
{
public:
	std::vector<Node> nodes;
	NodeSet nodeset;
	
	int _fx, _gx;
	double alpha;

	NodeSet upload, drop, add, extra, changed;
	
	Simulation(std::vector<double> &capacities, int fx, int gx, init_func_t func, double alpha, int window_size) : _fx(fx), _gx(gx), alpha(alpha) {
		nodes.reserve(capacities.size());
		
		for (size_t i = 0; i < capacities.size(); i++)
			nodes.push_back( Node(capacities[i], (*func)(capacities[i], capacities.size()), window_size) );
		
		for (size_t i = 0; i < nodes.size(); i++)
			nodeset.insert( &nodes[i] );
		
		for (size_t i = 0; i < nodes.size(); i++)
			update_connections(nodes[i], 0);
	}
	
	void fx(Node &node) {
		node.N = std::min(node.N + _fx, std::min<int>(nodes.size()-2,  (int)(node.capacity/5.0)));
	}
	
	void gx(Node &node) {
		node.N = std::max(node.N - _gx, 1);
	}
	
	void update(int step) {
		Node &node = nodes[rand() % nodes.size()];
		
		dynamic_wakeup(node, step);
	}

	void dynamic_wakeup(Node &node, int step) {
		NodeSet trading_peers;
		make_intersection(node.download_peers, node.upload_peers, trading_peers);
		
		update_total(node, step);
		
		int elapsed = step - node.last_wakeup;
			double avg = node.total / elapsed;
			node.total = 0;
			
		if (node.steps == 0) {
		  node.avg = avg;
		} else {
		  node.avg = alpha * node.avg + (1.0 - alpha) * avg;
		}
		
		node.last_wakeup = step;
		
		if (++node.window_size_index == node.new_window_size) {
		  update_N(node, step);
		  node.window_size_index = 0;
		  int new_size = node.max_window_size + 5 - node.upload_peers.size();
		  if (new_size > node.max_window_size) {
			new_size = node.max_window_size;
		  }
		  node.new_window_size = new_size;
		}
		
		update_connections(node, step);
	}
	
	void wakeup(Node &node, int step) {
		NodeSet trading_peers;
		make_intersection(node.download_peers, node.upload_peers, trading_peers);
		
		update_total(node, step);
		
		int elapsed = step - node.last_wakeup;
		double avg = node.total / elapsed;
		node.total = 0;
		
		if (node.steps == 0)
			node.avg = avg;
		else
			node.avg = alpha * node.avg + (1.0 - alpha) * avg;
		
		node.last_wakeup = step;
		
		if (++node.steps >= node.window_size) {
			update_N(node, step);
			node.steps = 0;
		}
		
		update_connections(node, step);
	}
	
	void relative_wakeup(Node &node, int step) {
		node.window_size = std::max<int>((int) sqrt(node.capacity), node.window_size);
		wakeup(node, step);
	}
	
	void update_connections(Node &node, int step) {
		upload.clear();
		drop.clear();
		add.clear();
		extra.clear();
		changed.clear();
		
		//		std::cout << "calculating best contributors..." << std::endl;
		
		std::vector< Node * > temp(node.download_peers.begin(), node.download_peers.end());
		std::sort(temp.begin(), temp.end(), LinkGreater);
		if (temp.size() > node.N)
			temp.resize(node.N);
		
		upload.insert(temp.begin(), temp.end());
		
		//		std::cout << "calculating drop..." << std::endl;
		
		make_difference(node.upload_peers, upload, drop);
		
		//		std::cout << "calculating add..." << std::endl;
		
		make_difference(upload, node.upload_peers, add);
		
		int diff = node.N + 1 - upload.size();
		
		//		std::cout << "finding " << diff << " unchokes..." << std::endl;
		
		get_extra_peers(node, upload, diff, extra);
		
		make_union(upload, extra, upload);
		make_union(add, extra, add);
		
		//		std::cout << "updating changed peers..." << std::endl;
		
		make_union(add, drop, changed);
		for (NodeSet::iterator it = changed.begin(); it != changed.end(); ++it) {
			update_total(**it, step);
		}
		
		//		std::cout << "dropping..." << std::endl;
		
		for (NodeSet::iterator it = drop.begin(); it != drop.end(); ++it) {
			Node *peer = *it;
			peer->download_peers.erase(&node);
			peer->link_capacity = peer->capacity / std::max<int>(peer->upload_peers.size(), 1);
		}
		
		//		std::cout << "adding..." << std::endl;
		
		for (NodeSet::iterator it = add.begin(); it != add.end(); ++it) {
			Node *peer = *it;
			peer->download_peers.insert(&node);
			peer->link_capacity = peer->capacity / std::max<int>(peer->upload_peers.size(), 1);
		}
		
		//		std::cout << "setting upload peers..." << std::endl;
		
		node.upload_peers.swap(upload);
		node.link_capacity = node.capacity / std::max<int>(node.upload_peers.size(), 1);
	}
	
	void get_extra_peers(Node &node, NodeSet upload, int count, NodeSet &add) {
		NodeSet used;
		used.insert(&node);
		
		std::vector<Node *> vpool;
		make_difference(nodeset, used, vpool);
		
		while (upload.size() < node.N+1) {
			int i = rand() % vpool.size();
			if (upload.insert(vpool[i]).second)
				add.insert(vpool[i]);
		}
	}
	
	void update_total(Node &node, int step) {
		int elapsed = step - node.last_changed;
		NodeSet trading_peers;
		make_intersection(node.download_peers, node.upload_peers, trading_peers);
		double add = elapsed * download_rate(trading_peers);
		node.total += add;
		node.last_changed = step;
	}
	
	void update_N(Node &node, int step) {
		for (NodeSet::iterator it = node.upload_peers.begin(); it != node.upload_peers.end(); ++it) {
			update_total(**it, step);
		}
		
		double avg = node.avg;
		
		if (avg > 0.0) {
			double f = avg / node.last_avg;
			
			if (f < 1.0 - 0.01) {
				if (node.last_action == Decrease) {
					fx(node);
					node.last_action = Increase;
				} else {
					gx(node);
					node.last_action = Decrease;
				}
			} else if (f > 1.0 - 0.01) {
				if (node.last_action == Decrease)
					gx(node);
				else
					fx(node);
			}
		}
		
		node.last_avg = avg;
		node.avg = 0.0;
		node.last_update = step;
	}
};

class Output
{
	Simulation &sim;
	
	std::ofstream variance, download_rates;
public:
	Output(Simulation &simulation) : sim(simulation), variance("variance.txt"), download_rates("download_rates.txt") {}
	
	void dump_variance() {
	}
	
	void dump_download_speeds() {
	}
	
	void dump_adjacency() {
	}
	
	void output(int step) {
		bool disabled_adjacency = true;
		bool disabled_variance = false;
		bool disabled_download_rates = true;
		
		int download_connections = 0;
		std::vector<double> downloads;
		
		std::ostringstream str;
					
		if (disabled_download_rates == false) {			
			download_rates << step;
		}
		
		for (int i = 0; i < sim.nodes.size(); i++) {
			Node &node = sim.nodes[i];
			
			download_connections += node.download_peers.size();
			double D = download_rate(node.download_peers);
			downloads.push_back(D);
			
			NodeSet peers;
			make_union(node.upload_peers, node.download_peers, peers);
			NodeSet trading;
			make_intersection(node.upload_peers, node.download_peers, trading);
		
			if (disabled_adjacency == false) {
				str << "adjacency" << step << ".txt";
				std::ofstream adjacency(str.str().c_str());
				for (NodeSet::iterator it = peers.begin(); it != peers.end(); ++it) {
					Node *peer = *it;
					adjacency << node.capacity << '\t';
					if (node.upload_peers.find(peer) != node.upload_peers.end())
						adjacency << peer->capacity;
					adjacency << '\t';
					if (node.download_peers.find(peer) != node.download_peers.end())
						adjacency << peer->capacity;
					adjacency << '\t';
					if (trading.find(peer) != trading.end())
						adjacency << peer->capacity;
					adjacency << '\n';
				}
			}
			
			if (disabled_download_rates == false) {
				download_rates << '\t' << D;
			}
		}
		
		if (disabled_download_rates == false) {
			download_rates << std::endl;
		}
		
		double download_average = Average<double>(downloads.begin(), downloads.end());
		double download_variance = Variance<double>(downloads.begin(), downloads.end(), download_average);
		
		if (disabled_variance == false) {
			variance << step << '\t' << download_connections << '\t' << download_variance << std::endl;
		}
	}
};

int init_fixed(double c, int n) {
	return std::min(5, (int)(c/5.0));
}

int init_sqrt(double c, int n) {
	return std::min( std::min((int)sqrt(c), (int)(c/5.0)), n-2);
}

int init_random(double c, int n) {
	return std::min( rand() % n-2, (int)(c/5.0) );
}

int init_full(double c, int n) {
	return std::min( n-2, (int)(c/5.0) );
}

void read_capacities(const std::string &file, std::vector<double> &capacities)
{
	std::ifstream f(file.c_str());
	
	while (!f.eof()) {
		double d;
		f >> d;
		capacities.push_back(int(d));
	}
}

void run(const std::string &capacities_file, int count, int fx, int gx, int seed, init_func_t init, double alpha, int window_size, int limit)
{
	std::vector<double> capacities;
	read_capacities(capacities_file, capacities);
	
	srand(1);
	
	std::random_shuffle(capacities.begin(), capacities.end(), RandomGenerator );
	if (capacities.size() > count)
		capacities.resize(count);
	std::sort(capacities.begin(), capacities.end(), std::greater<double>());
	
	if (true)
	{
		std::cout << "[";
		
		for (int i = 0; i < capacities.size()-1; i++)
			std::cout << capacities[i] << ", ";
		std::cout << capacities.back();
		
		std::cout << "]" << std::endl;
	}
	
	srand(seed);
	
	Simulation *sim = new Simulation(capacities, fx, gx, init, alpha, window_size);
	
	Output output(*sim);
	
	for (int i = 0; i <= limit; i++) {
		sim->update(i);
		
		if (i % (10000) == 0)
			output.output(i);
	}
	
	delete sim;
}

int main(int argc, char * const argv[])
{
	if (argc < 10) {
		std::cout << "usage: " << argv[0] << " <capacities-file> <int count> <int fx> <int gx> <int seed> <fixed|sqrt|random|full> <double alpha> <int window-size> <int limit>" << std::endl;
		return 0;
	}
	
	std::map<std::string, init_func_t> init_funcs;
	init_funcs[std::string("fixed")] = &init_fixed;
	init_funcs[std::string("sqrt")] = &init_sqrt;
	init_funcs[std::string("random")] = &init_random;
	init_funcs[std::string("full")] = &init_full;
	
	time_t start, end;
	
	time(&start);
	
	run(argv[1], strtol(argv[2], NULL, 0), strtol(argv[3], NULL, 0), strtol(argv[4], NULL, 0), strtol(argv[5], NULL, 0), init_funcs[argv[6]], strtod(argv[7], NULL), strtol(argv[8], NULL, 0), strtol(argv[9], NULL, 0));
	
	time(&end);
	double elapsed = difftime(end, start);
	
	std::cout << "\nsimulation completed in " << elapsed << " seconds.\n" << std::endl;
}
