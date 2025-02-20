#include "Graph.h"
#include <atomic>
#include <omp.h>
double threshold=1e9;
Timer thresh, branchings, bounding;
std::atomic<ui> best_solution_size(0);
std::vector<ui> kplex;

#include "KPlex_BB_matrix.h"
#include "KPlex_BB.h"
#include "CTPrune.h"
using namespace std;


Graph::Graph(const char *_dir, const int _K) {
	dir = string(_dir);
	K = _K;

	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	edgelist_pointer = nullptr;

	kplex.clear();
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pend_buf != nullptr) {
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(edgelist_pointer != nullptr) {
		delete[] edgelist_pointer;
		edgelist_pointer = nullptr;
	}
}

void Graph::read_graph_binary() {
	    FILE *f = Utility::open_file(dir.c_str(), "rb");
        ui tt;
        fread(&tt, sizeof(ui), 1, f);
        if (tt != sizeof(ui)) {
            printf("sizeof unsigned ui is different: file %u, machine %lu\n", tt, sizeof(ui));
        }
        fread(&n, sizeof(ui), 1, f);	// the number of vertices
        fread(&m, sizeof(ui), 1, f); // the number of edges (twice the acutal number).
 
        ui *degree = new ui[n];
        fread(degree, sizeof(ui), n, f);
        if (pstart != nullptr) delete[] pstart;
        pstart = new ui[n + 1];
        if (edges != nullptr) delete[] edges;
        edges = new ui[m];

        pstart[0] = 0;
        for (ui i = 0; i < n; i++) {
            if (degree[i] > 0){
                fread(edges + pstart[i], sizeof(ui), degree[i], f);
                //std::sort(edges+pstart[i], edges + pstart[i] + degree[i]);
            }
            pstart[i + 1] = pstart[i] + degree[i];
        }
        fclose(f);        
        delete[] degree; 
// 	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
// 	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

// 	ui tt;
// 	fread(&tt, sizeof(int), 1, f);
// 	if(tt != sizeof(int)) {
// 		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
// 		return ;
// 	}
// 	fread(&n, sizeof(int), 1, f);
// 	fread(&m, sizeof(int), 1, f);

// 	printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

// 	ui *degree = new ui[n];
// 	fread(degree, sizeof(int), n, f);

// #ifndef NDEBUG
// 	long long sum = 0;
// 	for(ui i = 0;i < n;i ++) sum += degree[i];
// 	if(sum != m) printf("m not equal sum of degrees\n");
// #endif

// 	fclose(f);

// 	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

// 	if(pstart == nullptr) pstart = new ept[n+1];
// 	if(edges == nullptr) edges = new ui[m];

// 	pstart[0] = 0;
// 	for(ui i = 0;i < n;i ++) {
// 		if(degree[i] > 0) {
// 			fread(edges+pstart[i], sizeof(int), degree[i], f);

// 			// remove self loops and parallel edges
// 			ui *buff = edges+pstart[i];
// 			sort(buff, buff+degree[i]);
// 			ui idx = 0;
// 			for(ui j = 0;j < degree[i];j ++) {
// 				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
// 				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
// 				buff[idx ++] = buff[j];
// 			}
// 			degree[i] = idx;
// 		}

// 		pstart[i+1] = pstart[i] + degree[i];
// 	}

// 	fclose(f);

// 	delete[] degree;
}

void Graph::read_graph() {
	printf("# Start reading graph from an edgelist file, Require files \"edges.txt\"\n");
	printf("# Note that this function is not optimized. Reading from a binary file will be faster\n");
	FILE *f = Utility::open_file((dir + string("/edges.txt")).c_str(), "r");

	fscanf(f, "%u%u", &n, &m);
	m *= 2;
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	vector<pair<ui,ui> > vp;
	for(ui i = 0;i < m/2;i ++) {
		ui a, b;
		fscanf(f, "%u%u", &a, &b);
		if(a >= n || b >= n) {
			printf("!!! Vertex IDs must be between 0 and n-1. Exit !!!\n");
			return ;
		}
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ept[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pstart[i+1] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edges[pstart[i+1] ++] = vp[idx ++].second;
	}

	fclose(f);

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::output_one_kplex() {
	FILE *fout = Utility::open_file("kplexes.txt", "w");
	fprintf(fout, "%lu\n", kplex.size());
	sort(kplex.begin(), kplex.end());
	for(ui i = 0;i < kplex.size();i ++) fprintf(fout, " %u", kplex[i]);
	fprintf(fout, "\n");
	fclose(fout);
}

void Graph::verify_kplex() {
	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	FILE *fin = Utility::open_file("kplexes.txt", "r");

	ui kplex_size = n, kplex_n, idx = 0;
	char ok = 1;
	while(fscanf(fin, "%u", &kplex_n) == 1) {
		++ idx;
		if(kplex_size == n) {
			kplex_size = kplex_n;
			printf("k-plex sizes: %u\n", kplex_size);
		}
		if(kplex_n != kplex_size) printf("!!! WA k-plex size: %u!\n", kplex_n);
		vector<ui> kplex;
		for (ui i = 0; i < kplex_n; i++) {
			ui tmp;
			fscanf(fin, "%u", &tmp);
			kplex.pb(tmp);
		}

		for (ui i = 0; i < kplex.size(); i++) {
			if (vis[kplex[i]]) {
				printf("WA k-plex! Duplicate vertex: %u\n", idx);
				ok = 0;
				break;
			}
			vis[kplex[i]] = 1;
		}
		for(ui i = 0;i < kplex.size();i ++) {
			ui d = 0;
			for(ui j = pstart[kplex[i]];j < pstart[kplex[i]+1];j ++) if(vis[edges[j]]) ++ d;
			if(d + K < kplex.size()) {
				ok = 0;
				printf("WA k-plex! Not enough neighbors!\n");
			}
		}
		for(ui i = 0;i < kplex.size();i ++) vis[kplex[i]] = 0;
	}
	if(ok) printf("Correct k-plexes!\n");
	fclose(fin);

	delete[] vis;
}

void Graph::kPlex_degen() {
	Timer t;


	assert(K > 0&&K < n);

	kplex.clear();
	heuristic_kplex_max_degree(10);

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;

	if(kplex.size() < UB) printf("\tHeuristic k-plex Size: %lu, UB: %u, Total Time: %s (microseconds)\n", kplex.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
	else printf("\tMaximum k-plex Size: %lu, Total Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::kPlex_exact(int mode) {
	Timer t;

	auto nn=n;
	auto mm=m;

	assert(K > 0);
	if(K <= 1) {
		printf("\tFor k <= 1, please invoke clique computation algorithms\n");
		return ;
	}
	if(K >= n) {
	printf(">>%s \tMaxKPlex_Size: %lu t_Total: %f t_Bounding: %f\n", dir.substr(dir.find_last_of("/")+1).c_str(), kplex.size(), t.elapsed()/1e6, bounding.ticktock());
	return ;
	}

	kplex.resize(2*K-2); // trivial cases
	heuristic_kplex_max_degree(10);

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);
	assert(kplex.size() >= K);
	ui search_time=0;
	if(kplex.size() < UB) {
		ui old_size = kplex.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];

		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);


		ego_degen(n, m, peel_sequence, pstart, edges, degree, rid, vis, heap, true);

		if(kplex.size() > old_size) {
			old_size = kplex.size();
			for(ui i = 0;i < kplex.size();i ++) {
				assert(kplex[i] < n);
				kplex[i] = out_mapping[kplex[i]];
			}
		}
		ui presize=kplex.size();
		if(kplex.size()+1 > 2*K) CTPrune::core_truss_copruning(n, m, kplex.size()+1-K, kplex.size()+1-2*K, peel_sequence, out_mapping, rid, pstart, edges, degree, true);
		else core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);



			ui *peel_sequence_rid = core;
			for(ui i= 0;i < n;i ++) peel_sequence_rid[peel_sequence[i]] = i;

			memset(vis, 0, sizeof(char)*n);

			if(pend == nullptr) pend = new ept[n+1];
			reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
			best_solution_size.store(kplex.size());
			cout<<"Best solution size: "<<best_solution_size.load()<<endl;
			thread_local ui ts_time=0;
#pragma omp parallel
		{

			thread_local Timer tt;
			vector<ui> ids, kplex_local=kplex;
			vector<pair<ui,ui> > vp;


			char* exists = new char[n];
			ui* degree = new ui[n];
			ui* rid = new ui[n];

			std::fill(exists, exists+n, 0);
			// std::fill(degree, degree+n, 0);
			// std::fill(rid, rid+n, 0);


			KPLEX_BB_MATRIX *kplex_solver_m = new KPLEX_BB_MATRIX();
			kplex_solver_m->allocateMemory(n);

			ui search_cnt = 0;
			double min_density = 1, total_density = 0;
#pragma omp for schedule(dynamic)
			for(ui i = 0;i < n;i ++) {
				ui best_sz = best_solution_size.load();
				if(best_sz >= UB) continue;
				ui u = peel_sequence[i];

				if(pend[u]-pstart[u]+K <= best_sz||n-i < best_sz) continue;

				fflush(stdout);
				
				if(best_sz >= 2*K-1) 
				extract_subgraph_with_prune(u, best_sz+1-K, best_sz+1-2*K, best_sz+3-2*K, peel_sequence_rid, degree, ids, rid, vp, exists, pstart, pend, edges);
				else extract_subgraph_wo_prune(u, peel_sequence_rid, ids, rid, vp, exists, pstart, pend, edges);

				if(ids.empty()||ids.size() <= best_sz) continue;

				double density = vp.size()*2/(double)ids.size()/(ids.size()-1);
				++ search_cnt;
				total_density += density;
				if(density < min_density) min_density = density;
				ui t_old_size = kplex_local.size();
					kplex_solver_m->load_graph(ids, vp);
				#pragma omp taskgroup
				{
					kplex_solver_m->kPlex(K, UB, kplex_local, true);
				}
			}
			// delete kplex_solver_m;

			ts_time = tt.elapsed();
			if(search_cnt == 0) printf("search_cnt: 0, ave_density: 1, min_density: 1\n");
			// else printf("search_cnt: %u, ave_density: %.5lf, min_density: %.5lf\n", search_cnt, total_density/search_cnt, min_density);
			// printf("*** Search time: %s\n", Utility::integer_to_string(tt.elapsed()).c_str());
		}
	        #pragma omp parallel reduction(max : search_time)
			{
			 search_time=ts_time;
			}
		if(kplex.size()>presize) for(ui i=0;i<kplex.size();i++)kplex[i]=out_mapping[kplex[i]];
		delete[] out_mapping;
		delete[] rid;
	}

	delete heap;
	delete[] core;
	delete[] peel_sequence;
	delete[] vis;
	delete[] degree;
	printf(">>%s \tMaxKPlex_Size: %lu t_Total: %f t_search: %f\n", dir.substr(dir.find_last_of("/")+1).c_str(), kplex.size(), t.elapsed()/1e6, search_time/1e6);


	// printf("\tMaximum kPlex Size: %lu, Total Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		for(ui j = pstart[u];j < pend[u]&&rid[edges[j]] > rid[u];j ++) {
			ui v = edges[j];
			edges[pend[v] ++] = u;
			assert(pend[v] <= pstart[v+1]);
		}
	}
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend[i] == pstart[i+1]);
#endif
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		while(end < pstart[i+1]&&rid[edges[end]] > rid[i]) ++ end;
	}
}

// each of u's neighbors must have at least triangle_threshold common neighbors with u
// each of u's non-neighbor must have at least cn_threshold common neighbors with u
// after pruning, u must have at least degree_threshold neighbors
ui Graph::extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!exists[i]);
#endif

	ids.clear(); vp.clear();
	ids.push_back(u); exists[u] = 1;
	for(ept i = pstart[u];i < pend[u];i ++) {
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]); exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u+1]||p_rid[edges[pend[u]]] < p_rid[u]);

	// Utility::print_array("ids1", ids.data(), 0, ids.size(), 0);

	ui *Q = rid;
	ui Q_n = 0;
	for(ui i = 1;i < ids.size();i ++) {
		ui v = ids[i];
		degree[v] = 0;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] > p_rid[u];j ++) {
			if(exists[edges[j]]) ++ degree[v];
		}
		if(degree[v] < triangle_threshold) Q[Q_n++] = v;
	}
	for(ui i = 0;i < Q_n;i ++) {
		ui v = Q[i];
		exists[v] = 3;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] > p_rid[u];j ++) if(exists[edges[j]] == 2) {
			if(degree[edges[j]] == triangle_threshold) Q[Q_n++] = edges[j];
			-- degree[edges[j]];
		}
	}
	assert(Q_n < ids.size());
	if(ids.size() - Q_n  - 1 < degree_threshold) {
		for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
		ids.clear();
		return 0;
	}

	ui old_size = ids.size();
	for(ui i = 1;i < old_size;i ++) if(exists[ids[i]] == 2) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] > p_rid[u];j ++) {
			if(!exists[edges[j]]) {
				ids.push_back(edges[j]);
				exists[edges[j]] = 1;
				degree[edges[j]] = 1;
			}
			else ++ degree[edges[j]];
		}
	}

	ui new_size = 1;
	for(ui i = 1;i < old_size;i ++) {
		if(exists[ids[i]] == 3) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);
	for(ui i = old_size;i < ids.size();i ++) {
		if(degree[ids[i]] < cn_threshold) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for(ui i = 0;i < ids.size();i ++) rid[ids[i]] = i;
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) {
			assert(rid[v] < ids.size()&&rid[edges[j]] < ids.size());
			vp.push_back(make_pair(rid[v], rid[edges[j]]));
		}
	}
	for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
	return old_size;
}

ui Graph::extract_subgraph_wo_prune(ui u, const ui *p_rid, vector<ui> &ids, ui *rid, vector<pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges) {
	ids.clear(); vp.clear();
	ids.push_back(u); exists[u] = 1; rid[u] = 0;
	for(ept i = pstart[u];i < pend[u];i ++) {
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]); exists[edges[i]] = 1; rid[edges[i]] = ids.size()-1;
	}
	assert(pend[u] >= pstart[u+1]||p_rid[edges[pend[u]]] < p_rid[u]);
	ui old_size = ids.size();
	for(ui i = 1;i < old_size;i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] > p_rid[u];j ++) {
			ui w = edges[j];
			if(exists[w]) continue;
			ids.push_back(w); exists[w] = 1; rid[w] = ids.size()-1;
		}
	}
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}
	for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
	return old_size;
}

void Graph::write_subgraph(ui n, const vector<pair<int,int> > &edge_list) {
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for(ui i = 0;i < edge_list.size();i ++) fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer) {
	ids_n = 0; vp.clear();
	ids[ids_n++] = u; exists[u] = 1; rid[u] = 0;
	ui u_n = pstart[u];
	for(ept i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		ui v = edges[i];
		rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
	}
	pend[u] = u_n;
	ui old_size = ids_n;
	for(ui i = 1;i < old_size;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(ept j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			ui v = edges[j];
			if(exists[v]) continue;
			rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
		}
		pend[u] = u_n;
	}
	for(ui i = 0;i < old_size;i ++) {
		u = ids[i];
		for(ept j = pstart[u];j < pend[u];j ++) if(edges[j] > u) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(ui i = old_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(ept j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
}

void Graph::extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, vector<pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer) {
	vp.clear();
	for(ui i = 0;i < ids_n;i ++) {
		rid[ids[i]] = i;
		exists[ids[i]] = 1;
	}
	for(ui i = 0;i < ids_n;i ++) {
		ui u = ids[i];
		for(ept j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]&&edges[j] > u&&exists[edges[j]]) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui,ui> > &vp, ui *Q, ui* degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer) {
	vp.clear();
	ids_n = 0; ids[ids_n++] = u; exists[u] = 1;
	ui u_n = pstart[u];
	for(ept i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		ui v = edges[i];
		ids[ids_n++] = v; exists[v] = 2;
	}
	pend[u] = u_n;
	
	ui Q_n = 0;
	for(ui i = 1;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for(ept j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(exists[edges[j]] == 2) ++ degree[u];
		}
		pend[u] = u_n;
		if(degree[u]+2*K <= kplex.size()) Q[Q_n++] = u;
	}
	for(ui i = 0;i < Q_n;i ++) {
		u = Q[i];
		exists[u] = 10;
		for(ept j = pstart[u];j < pend[u];j ++) if(exists[edges[j]] == 2) {
			if( (degree[edges[j]]--) + 2*K == kplex.size()+1) {
				assert(Q_n < m/2);
				Q[Q_n++] = edges[j];
			}
		}
	}
	assert(Q_n <= ids_n);
	if(ids_n - 1 - Q_n + K <= kplex.size()) {
		for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
		ids_n = 0;
		return ;
	}
	
	ui nr_size = ids_n;
	for(ui i = 1;i < nr_size;i ++) if(exists[ids[i]] == 2) {
		u = ids[i];
		for(ept j = pstart[u];j < pend[u];j ++) {
			if(!exists[edges[j]]) {
				ids[ids_n++] = edges[j];
				exists[edges[j]] = 3;
				degree[edges[j]] = 1;
			}
			else if(exists[edges[j]] == 3) ++ degree[edges[j]];
		}
	}

#ifndef NDEBUG
	//printf("Entire list: ");
	//for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
	//printf("\n");
#endif

	ui new_size = 1;
	for(ui i = 1;i < nr_size;i ++) {
		if(exists[ids[i]] == 10) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if(new_size + Q_n != nr_size) {
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for(ui i = 0;i < new_size;i ++) printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for(ui i = 0;i < Q_n;i ++) printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	ui old_nr_size = nr_size;
	nr_size = new_size;
	for(ui i = old_nr_size;i < ids_n;i ++) {
		if(degree[ids[i]] + 2*K <= kplex.size()+2) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for(ui i = 1;i < nr_size;i ++) assert(exists[ids[i]] == 2);
	for(ui i = nr_size;i < ids_n;i ++) assert(exists[ids[i]] == 3);
#endif

	//for(ui i = 0;i < ids_n;i ++) printf(" %u", ids[i]);
	//printf("\n");

	for(ui i = 0;i < ids_n;i ++) {
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for(ui i = 0;i < nr_size;i ++) {
		u = ids[i];
		for(ept j = pstart[u];j < pend[u];j ++) if(exists[edges[j]]&&edges[j] > u) {
			assert(!deleted[edgelist_pointer[j]]);
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(ui i = nr_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(ept j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(exists[i] == 0);
#endif
}

// max-degree-based heuristic k-plex computation
void Graph::heuristic_kplex_max_degree(ui processed_threshold) {
	Timer t;
	// assert(kplex.empty());
	ui *head = new ui[n];
	ui *next = new ui[n];
	ui *degree = new ui[n];

	ui *vis = new ui[n];
	memset(vis, 0, sizeof(ui)*n);

	int max_degree = 0;
	for(ui i = 0;i < n;i ++) head[i] = n;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1]-pstart[i];
		if(degree[i] > max_degree) max_degree = degree[i];
		next[i] = head[degree[i]];
		head[degree[i]] = i;
	} // continued... 

	for(ui processed_vertices = 0;max_degree + K >= kplex.size()&&processed_vertices < processed_threshold;processed_vertices ++) {
		ui u = n;
		while(max_degree >= 0&&max_degree + K >= kplex.size()&&u == n) {
			for(ui v = head[max_degree];v != n;) {
				ui tmp = next[v];
				if(degree[v] == max_degree) {
					u = v;
					head[max_degree] = tmp;
					break;
				}
				else if(degree[v] + K >= kplex.size()) {
					next[v] = head[degree[v]];
					head[degree[v]] = v;
				}
				v = tmp;
			}
			if(u == n) {
				head[max_degree] = n;
				-- max_degree;
			}
		}
		if(u == n) break;

		vis[u] = 1;
		for(ui k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k]]) -- degree[edges[k]];

		vector<ui> vs;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vs.pb(edges[j]);

		vector<ui> vs_deg(vs.size());
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 2;
		for(ui j = 0;j < vs.size();j ++) {
			ui v = vs[j], d = 0;
			for(ui k = pstart[v];k < pstart[v+1];k ++) {
				if(vis[edges[k]] == 2) ++ d;
			}
			vs_deg[j] = d;
		}
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 0;

		vector<ui> res; res.pb(u);
		ui vs_size = vs.size();
		while(vs_size > 0&&res.size() + vs_size + K - 1 > kplex.size()) {
		// while(vs_size > 0) {
			ui idx = 0;
			for(ui j = 1;j < vs_size;j ++) {
				if(vs_deg[j] > vs_deg[idx]) idx = j;
				else if(vs_deg[j] == vs_deg[idx]&&degree[vs[j]] > degree[vs[idx]]) idx = j;
			}
			u = vs[idx];

			ui new_size = 0;
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vis[edges[j]] = 2;
			for(ui j = 0;j < vs_size;j ++) if(vis[vs[j]]) {
				if(j != new_size) swap(vs[new_size], vs[j]);
				vs_deg[new_size] = vs_deg[j];
				++ new_size;
			}
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 2) vis[edges[j]] = 0;

			res.pb(u);
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = k+2;
			for(ui j = new_size;j < vs_size;j ++) {
				ui v = vs[j];
				for(ui k = pstart[v];k < pstart[v+1];k ++) {
					if(vis[edges[k]] >= 2) -- vs_deg[vis[edges[k]]-2];
				}
			}
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = 0;

			vs_size = new_size;
		}

		// TO DO: extend res to be a maximal k-plex

		if(res.size() > kplex.size()) kplex = res;
	}

	delete[] vis;
	delete[] head;
	delete[] next;
	delete[] degree;

	printf("*** Heuristic kplex size: %lu, time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;

	ui threshold = (kplex.size()+1 > K? kplex.size()+1-K: 0);

	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] < threshold) peel_sequence[queue_n ++] = i;
	for(ui i = 0;i < queue_n;i ++) {
		ui u = peel_sequence[i]; degree[u] = 0;
		for(ept j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) peel_sequence[queue_n ++] = edges[j];
		}
	}
	ui UB = n;
	if(queue_n == n) UB = kplex.size();

	memset(vis, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		if(degree[i] >= threshold) peel_sequence[queue_n + (new_size ++)] = i;
		else {
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	if(new_size != 0) {
		heap->init(new_size, new_size-1, peel_sequence+queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		UB = 0;
		for(ui i = 0;i < new_size;i ++) {
			ui u, key;
			heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + K;
			if(new_size - i < t_UB) t_UB = new_size - i;
			if(t_UB > UB) UB = t_UB;

			if(idx == n&&key + K >= new_size - i) idx = i;
			vis[u] = 1;

			for(ept j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) {
				heap->decrement(edges[j], 1);
			}
		}

		if(output) printf("*** Degeneracy k-plex size: %u, max_core: %u, UB: %u, Time: %s (microseconds)\n", new_size-idx, max_core, UB, Utility::integer_to_string(t.elapsed()).c_str());

		if(new_size - idx > kplex.size()) {
			kplex.clear();
			for(ui i = idx;i < new_size;i ++) kplex.pb(peel_sequence[queue_n + i]);
			if(!output) printf("Find a k-plex of size: %u\n", new_size - idx);
		}
	}

	return UB;
}

void Graph::ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;
	if(pend == nullptr) pend = new ept[n+1];
	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);

	if(pend_buf == nullptr) pend_buf = new ept[n+1];
	if(edgelist_pointer == nullptr) edgelist_pointer = new ui[m];
	ui *pstart_s = pend_buf;
	ui *pend_s = rid;
	ui *edges_s = edgelist_pointer;

	vector<ui> Q;
	vector<ui> vs;
	memset(vis, 0, sizeof(char)*n);
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		if(pend[u] - pstart[u] < kplex.size()) continue;

		vs.clear();
		for(ui j = pstart[u];j < pend[u];j ++) {
			vs.push_back(edges[j]);
			vis[edges[j]] = 1;
			degree[edges[j]] = 0;
		}
		for(ui j = 0;j < vs.size();j ++) for(ui k = pstart[vs[j]];k < pend[vs[j]];k ++) if(vis[edges[k]]) {
			++ degree[vs[j]]; ++ degree[edges[k]];
		}
		pend_s[vs[0]] = pstart_s[vs[0]] = 0;
		for(ui j = 1;j < vs.size();j ++) pend_s[vs[j]] = pstart_s[vs[j]] = pstart_s[vs[j-1]] + degree[vs[j-1]];
		for(ui j = 0;j < vs.size();j ++) for(ui k = pstart[vs[j]];k < pend[vs[j]];k ++) if(vis[edges[k]]) {
			edges_s[pend_s[vs[j]]++] = edges[k];
			edges_s[pend_s[edges[k]]++] = vs[j];
		}

		ui threshold = (kplex.size() > K? kplex.size()-K: 0); // all vertices with degree < threshold can be pruned
		Q.clear();
		for(ui j = 0;j < vs.size();j ++) if(degree[vs[j]] < threshold) {
			Q.push_back(vs[j]);
			vis[vs[j]] = 0;
		}
		for(ui j = 0;j < Q.size();j ++) for(ui k = pstart_s[Q[j]];k < pend_s[Q[j]];k ++) if(vis[edges_s[k]]) {
			if( (degree[edges_s[k]]--) == threshold) {
				Q.push_back(edges_s[k]);
				vis[edges_s[k]] = 0;
			}
		}
		ui cnt = 0;
		for(ui j = 0;j < vs.size();j ++) if(vis[vs[j]]) vs[cnt++] = vs[j];
		assert(cnt + Q.size() == vs.size());
		vs.resize(cnt);
		if(cnt == 0) continue;

		heap->init(vs.size(), vs.size()-1, vs.data(), degree);
		bool found = false;
		for(ui ii = 0;ii < vs.size();ii ++) {
			ui v, key;
			heap->pop_min(v, key);
			if(found) {
				kplex.push_back(v);
				continue;
			}

			if(vs.size()-ii+1 <= kplex.size()) break;

			if(key + K >= vs.size() - ii) {
				kplex.clear();
				kplex.push_back(u);
				kplex.push_back(v);
				found = true;
				continue;
			}

			vis[v] = 0;
			for(ept j = pstart_s[v];j < pend_s[v];j ++) if(vis[edges_s[j]]) heap->decrement(edges_s[j], 1);
		}
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 0;
	}
	for(ui i = 0;i < n;i ++) pend_buf[i] = pend[i];
	for(ui i = 0;i < n;i ++) for(ept j = pstart[i];j < pend[i];j ++) edges[pend_buf[edges[j]]++] = i;
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend_buf[i] == pstart[i+1]);
#endif

	if(output) printf("*** EGo-Degen kPlex size: %lu, Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}


// in_mapping and out_mapping can be the same array
void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output) {
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] + K > kplex.size()) {
		rid[i] = cnt;
		if(in_mapping == nullptr) out_mapping[cnt] = i;
		else out_mapping[cnt] = in_mapping[i];
		++ cnt;
	}

	if(cnt != n) {
		cnt = 0;
		ept pos = 0;
		for(ui i = 0;i < n;i ++) if(core[i] + K > kplex.size()) {
			ept t_start = pstart[i];
			pstart[cnt] = pos;
			for(ept j = t_start;j < pstart[i+1];j ++) if(core[edges[j]] + K > kplex.size()) {
				edges[pos ++] = rid[edges[j]];
			}
			++ cnt;
		}
		pstart[cnt] = pos;

		//printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n-cnt-1]] == 0||core[peel_sequence[n-cnt-1]] + K <= kplex.size());
		assert(cnt == 0||core[peel_sequence[n-cnt]] + K > kplex.size());
		for(ui i = 0;i < cnt;i ++) {
			peel_sequence[i] = rid[peel_sequence[n-cnt+i]];
			core[i] = core[out_mapping[i]];
		}

		if(pos > 0&&pos < m/2) {
			ept *pstart_new = new ept[cnt+1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept)*(cnt+1));
			memcpy(edges_new, edges, sizeof(ui)*pos);
			delete[] pstart; pstart = pstart_new;
			delete[] edges; edges = edges_new;
		}

		n = cnt;
		m = pos;
	}

	if(output) printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());
}

// orient graph
void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(ui i = 0;i < n;i ++) {
		ept &end = pend[i] = pstart[i];
		for(ept j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for(int i = 0;i < n;i ++) sum += pend[i] - pstart[i];
	assert(sum*2 == m);
#endif
}

// oriented triangle counting
void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj) {
	memset(adj, 0, sizeof(ui)*n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui)*m);
	for(ui u = 0;u < n;u ++) {
		for(ept j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(ept j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];
			for(ept k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
				++ tri_cnt[j];
				++ tri_cnt[k];
				++ tri_cnt[adj[edges[k]]-1];
				++ cnt;
			}
		}

		for(ept j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

#ifndef NDEBUG
	//printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf) {
	for(ui i = 0;i < n;i ++) pend2[i] = pend[i];
	ept pos = 0;
	for(ui i = 0;i < n;i ++) {
		for(ept j = pstart[i];j < pend[i];j ++) {
			tri_cnt[pos>>1] = edgelist_pointer[j]; edge_list[pos++] = i; edge_list[pos++] = edges[j];

			ept &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos>>1)-1;
			edges[k ++] = i;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend2[i] == pstart[i+1]);
#endif

	for(ui i = 0;i < n;i ++) {
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for(ui i = 0;i < n;i ++) {
		for(ept j = pend2[i];j < pstart[i+1];j ++) {
			ept &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k ++] = i;
		}
	}

	ept *ids = pend2;
	for(ui i = 0;i < n;i ++) {
		if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while(j < pend[i]&&k < pstart[i+1]) {
			if(edges[j] < edges[k]) {
				ids[pos] = edges[j];
				buf[pos ++] = edgelist_pointer[j ++];
			}
			else {
				ids[pos] = edges[k];
				buf[pos ++] = edgelist_pointer[k ++];
			}
		}
		while(j < pend[i]) {
			ids[pos] = edges[j];
			buf[pos ++] = edgelist_pointer[j ++];
		}
		while(k < pstart[i+1]) {
			ids[pos] = edges[k];
			buf[pos ++] = edgelist_pointer[k ++];
		}
		for(ept j = 0;j < pos;j ++) {
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}

char Graph::find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges) {
	if(b >= e) return 0;

	while(b+1 < e) {
		idx = b + (e-b)/2;
		if(edges[idx] > w) e = idx;
		else b = idx;
	}

	if(edges[b] == w) {
		idx = edgelist_pointer[b];
		if(!deleted[idx]) return 1;
	}

	return 0;
}

// return the number of peeled edges
ept Graph::peeling(ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists) {
	ept Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	if(initialize_Qe) {
		ept active_edgelist_newn = 0;
		for(ept j = 0;j < active_edgelist_n;j ++) if(!deleted[active_edgelist[j]]) {
			if(tri_cnt[active_edgelist[j]] < t_threshold) Qe[Qe_n++] = active_edgelist[j];
			else active_edgelist[active_edgelist_newn ++] = active_edgelist[j];
		}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	//printf("%lu\n", Qe_n);

	ept deleted_edges_n = 0;
	ui Qv_idx = 0;
	while(Qv_idx < Qv_n || Qe_n) {
		if(Qe_n == 0) {
			//printf("hit\n");
			ui u = Qv[Qv_idx ++]; // delete u from the graph due to have a degree < d_threshold
			ept u_n = pstart[u];
			for(ept k = pstart[u];k < pend[u];k ++) if(!deleted[edgelist_pointer[k]]) {
				edges[u_n] = edges[k]; edgelist_pointer[u_n++] = edgelist_pointer[k];
				exists[edges[k]] = 1;
			}
			pend[u] = u_n;

			for(ept k = pstart[u];k < pend[u];k ++) deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			linear_heap->decrement(u, degree[u]); degree[u] = 0;
			//printf("Removed %u\n", u);

			for(ept k= pstart[u];k < pend[u];k ++) {
				ui v = edges[k];
#ifndef NO_TRUSS_PRUNE
				ept v_n = pstart[v];
				for(ept x = pstart[v];x < pend[v];x ++) if(!deleted[edgelist_pointer[x]]) {
					edges[v_n] = edges[x]; edgelist_pointer[v_n++] = edgelist_pointer[x];
					if(edges[x] > v&&exists[edges[x]]) {
						if( (tri_cnt[edgelist_pointer[x]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[x];
					}
				}
				pend[v] = v_n;
#endif

				linear_heap->decrement(v, 1);
				if( (degree[v]--) == d_threshold) Qv[Qv_n++] = v;
			}

			for(ept k = pstart[u];k < pend[u];k ++) exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for(ept j = 0;j < Qe_n;j ++) {
			ept idx = Qe[j];
			ui u = edge_list[idx<<1], v = edge_list[(idx<<1)+1];
			ui tri_n = tri_cnt[idx];
			//printf("remove %u %u\n", u, v);
			deleted[idx] = 1;
			linear_heap->decrement(u, 1);
			linear_heap->decrement(v, 1);
			if( (degree[u] --) == d_threshold) Qv[Qv_n++] = u;
			if( (degree[v] --) == d_threshold) Qv[Qv_n++] = v;
			deleted_edges_n ++;
			
			if(degree[u] < degree[v]) swap(u,v);
			//printf("here\n");

			if(degree[u] > degree[v]*2) { // binary search
			//if(false) {
				ept v_n = pstart[v], start = pstart[u];
				for(ept k = pstart[v];k < pend[v];k ++) if(!deleted[edgelist_pointer[k]]) {
					edges[v_n] = edges[k]; edgelist_pointer[v_n++] = edgelist_pointer[k];

					if(tri_n&&find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges)) {
						-- tri_n;
						if( (tri_cnt[idx]--) == t_threshold) Qe[Qe_n++] = idx;
						if( (tri_cnt[edgelist_pointer[k]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[k];
					}
				}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else { // sorted_merge
				ept ii = pstart[u], jj = pstart[v];
				ept u_n = pstart[u], v_n = pstart[v];

				while(true) {
					while(ii < pend[u]&&deleted[edgelist_pointer[ii]]) ++ ii;
					while(jj < pend[v]&&deleted[edgelist_pointer[jj]]) ++ jj;
					if(ii >= pend[u]||jj >= pend[v]) break;

					if(edges[ii] == edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if( (tri_cnt[edgelist_pointer[ii]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[ii];
						if( (tri_cnt[edgelist_pointer[jj]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[jj];

						++ ii;
						++ jj;
					}
					else if(edges[ii] < edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ ii;
					}
					else {
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++ jj;
					}
				}
				while(ii < pend[u]) {
					if(!deleted[edgelist_pointer[ii]]) {
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ ii;
				}
				while(jj < pend[v]) {
					if(!deleted[edgelist_pointer[jj]]) {
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++ jj;
				}
				pend[u] = u_n; pend[v] = v_n;
			}
			//printf("finish %u %u\n", u, v);
		}
		Qe_n = 0;
	}
#ifndef NDEBUG
	printf("*** Truss removed %s undirected edges\n", Utility::integer_to_string(deleted_edges_n).c_str());
#endif
	return deleted_edges_n;
}


#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage() {
	printf("Example usage: ./kPlexT -g path_to_graph -k 3 -o -b\n");
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** kPlexT (Debug) build at %s %s ***\n", __TIME__, __DATE__);
	printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
	printf("**** kPlexT (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif
#define LEN_LIMIT (1<<10)
char filename[LEN_LIMIT];
	printf("\n-----------------------------------------------------------------------------------------\n");
	if (argc == 3) {
		strncpy(filename, argv[1], LEN_LIMIT);
		int k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read_graph_binary();
		graph->kPlex_exact(false);
		graph->output_one_kplex();
		delete graph;
		// delete graph; // there are some bugs in releasing memory
	}
	else if (argc == 4) {
		strncpy(filename, argv[1], LEN_LIMIT);
		int k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read_graph_binary();
		graph->kPlex_exact(false);
		// graph->kPlex_dense();
		graph->output_one_kplex();
		delete graph;
	}
	else printf("[usage]: exe file k [dense]\n");
	printf("-----------------------------------------------------------------------------------------\n\n");
	return 0;
}
