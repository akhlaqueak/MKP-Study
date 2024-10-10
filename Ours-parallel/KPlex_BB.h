#ifndef _KPLEX_BB_
#define _KPLEX_BB_

#include "Utility.h"
#include "Timer.h"

class KPLEX_BB {
private:
	ui n;
	ept *pstart;
	ept *pend;
	ui *edges;
	ept edges_cap;

	ui *degree;
	ui *degree_in_S;

	ui K;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;

	ui *neighbors;
	ui *nonneighbors;

	ui *S2;
	ui *SR; // union of S and R, where S is at the front
	ui *SR_rid; // reverse ID for SR
	std::queue<ui> Qv;
	ui *level_id;

	ui *buf;
	ui *buf1;
	ui *buf2;
	char *vis;

	std::vector<std::pair<ui,ui> > vp;

public:
	KPLEX_BB() {
		n = 0;
		edges_cap = 0;
		pstart = NULL;
		pend = NULL;
		edges = NULL;
		degree = degree_in_S = NULL;

		best_solution = NULL;
		K = best_solution_size = _UB_ = 0;

		S2 = SR = SR_rid = NULL;
		level_id = NULL;

		neighbors = nonneighbors = NULL;
		buf = buf1 = buf2 = NULL;
		vis = NULL;
	}

	~KPLEX_BB() {
		if(pstart != NULL) {
			delete[] pstart;
			pstart = NULL;
		}
		if(pend != NULL) {
			delete[] pend;
			pend = NULL;
		}
		if(edges != NULL) {
			delete[] edges;
			edges = NULL;
		}
		if(degree != NULL) {
			delete[] degree;
			degree = NULL;
		}
		if(degree_in_S != NULL) {
			delete[] degree_in_S;
			degree_in_S = NULL;
		}
		if(best_solution != NULL) {
			delete[] best_solution;
			best_solution = NULL;
		}
		if(S2 != NULL) {
			delete[] S2;
			S2 = NULL;
		}
		if(SR != NULL) {
			delete[] SR;
			SR = NULL;
		}
		if(SR_rid != NULL) {
			delete[] SR_rid;
			SR_rid = NULL;
		}
		if(level_id != NULL) {
			delete[] level_id;
			level_id = NULL;
		}
		if(neighbors != NULL) {
			delete[] neighbors;
			neighbors = NULL;
		}
		if(nonneighbors != NULL) {
			delete[] nonneighbors;
			nonneighbors = NULL;
		}
		if(buf != NULL) {
			delete[] buf;
			buf = NULL;
		}
		if(buf1 != NULL) {
			delete[] buf1;
			buf1 = NULL;
		}
		if(buf2 != NULL) {
			delete[] buf2;
			buf2 = NULL;
		}
		if(vis != NULL) {
			delete[] vis;
			vis = NULL;
		}
	}

	void allocateMemory(ui n) {
		if(n <= 0) return ;
		pstart = new ept[n+1];
		pend = new ept[n+1];
		edges = new ui[1];
		edges_cap = 1;

		degree = new ui[n+1];
		degree_in_S = new ui[n+1];
		best_solution = new ui[n+1];
		S2 = new ui[n+1];
		SR = new ui[n+1];
		SR_rid = new ui[n+1];
		neighbors = new ui[n+1];
		nonneighbors = new ui[n+1];
		level_id = new ui[n+1];

		buf = new ui[n+1];
		buf1 = new ui[n+1];
		buf2 = new ui[n+1];
		vis = new char[n+1];
	}

	void load_graph(ui _n, const std::vector<std::pair<ui,ui> > &vp) {
		n = _n;
		if(vp.size()*2 > edges_cap) {
			do {
				edges_cap *= 2;
			} while(vp.size()*2 > edges_cap);
			delete[] edges; edges = new ui[edges_cap];
		}
		for(ui i = 0;i < n;i ++) degree[i] = 0;
		for(ept i = 0;i < vp.size();i ++) {
			assert(vp[i].first < n&&vp[i].second < n);
			++ degree[vp[i].first];
			++ degree[vp[i].second];
		}
		for(ui i = 1;i < n;i ++) {
			degree[i] += degree[i-1];
			pstart[i] = degree[i-1];
		}
		pstart[0] = 0;
		for(ept i = 0;i < vp.size();i ++) {
			ui a = vp[i].first, b = vp[i].second;
			edges[pstart[a]++] = b;
			edges[pstart[b]++] = a;
		}
		for(ui i = n;i > 0;i --) pstart[i] = pstart[i-1];
		pstart[0] = 0;
	}

	void load_graph(ui n_, ui *_pstart, ui *_pend, ui *_edges) {
		ept m = 0;
		n = n_;
		for(ui i = 0;i < n;i ++) m += _pend[i]-_pstart[i];
		if(m > edges_cap) {
			do {
				edges_cap *= 2;
			} while(m > edges_cap);
			delete[] edges; edges = new ui[edges_cap];
		}
		m = 0;
		for(ui i = 0;i < n;i ++) {
			pstart[i] = m;
			for(ept j = _pstart[i];j < _pend[i];j ++) {
				assert(_edges[j] < n);
				edges[m ++] = _edges[j];
			}
		}
		pstart[n] = m;
		printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m/2, double(m)/n/(n-1));
	}

	void kPlex(ui K_, ui UB_, std::vector<ui> &kplex, bool must_include_0) {
		K = K_;
		_UB_ = UB_;
		if(K <= 1) {
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return ;
		}
		best_solution_size = kplex.size();
		ui R_end;
		initialization(R_end);
		if(R_end&&best_solution_size < _UB_) BB_search(0, R_end, 1, must_include_0);
		if(best_solution_size > kplex.size()) {
			kplex.clear();
			for(int i = 0;i < best_solution_size;i ++) kplex.push_back(best_solution[i]);
		}
	}

	int main(int argc, char *argv[]) {
		if(argc < 3) {
			printf("Usage: [1]exe [2]dir [3]k [4 option] lb_of_max_kplex_size\n");
			return 0;
		}
		readGraph_binary(argv[1]);
		printf("Finish reading graph\n");
		K = atoi(argv[2]);
		if(K == 1) {
			printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
			return 0;
		}
		best_solution_size = 1;
		_UB_ = n;
		if(argc >= 4) {
			best_solution_size = atoi(argv[3]);
			printf("initial lb: %u\n", best_solution_size);
		}
		ui R_end;
		Timer t;
		initialization(R_end);
		if(R_end) BB_search(0, R_end, 1, 0);
		printf("Maximum %u-plex size: %u, time excluding reading: %s (micro seconds)\n", K, best_solution_size, Utility::integer_to_string(t.elapsed()).c_str());
		return 0;
	}

private:
	void readGraph_binary(char* dir) {
		FILE *f = Utility::open_file( (std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");

		int tt;
		fread(&tt, sizeof(int), 1, f);
		if(tt != sizeof(int)) {
			printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
			return ;
		}
		ui m;
		fread(&n, sizeof(int), 1, f);
		fread(&m, sizeof(int), 1, f);
		printf("n = %u, m = %u\n", n, m/2);

		allocateMemory(n);
		if(m > edges_cap) {
			do {
				edges_cap *= 2;
			} while(m > edges_cap);
			delete[] edges; edges = new ui[edges_cap];
		}

		fread(degree, sizeof(ui), n, f);
		fclose(f);

		f = Utility::open_file( (std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");
		pstart[0] = 0;
		for(ui i = 0;i < n;i ++) {
			pstart[i+1] = pstart[i] + degree[i];
			if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		}
		fclose(f);
	}

	void initialization(ui &R_end) {
		memset(vis, 0, sizeof(char)*(n+1));
		memset(degree_in_S, 0, sizeof(ui)*(n+1));
		for(ui i = 0;i < n;i ++) level_id[i] = n;
		for(ui i = 0;i < n;i ++) SR[i] = SR_rid[i] = i;
		for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];
		while(!Qv.empty()) Qv.pop();

		for(ui i = 0;i < n;i ++) {
			degree[i] = pstart[i+1] - pstart[i];
			if(degree[i] + K <= best_solution_size) {
				level_id[i] = 0;
				Qv.push(i);
			}
		}

		R_end = n;
		if(!remove_vertices_with_prune(0, R_end, 0)) R_end = 0;
	}

	void reorganize_edges(ui S_end, ui R_end, ui level) {
		assert(level > 0);
		for(ui i = 0;i < R_end;i ++) {
			ui u = SR[i];
			ui non_neighbors_n = 0, end = pstart[u];
			for(ui j = pstart[u];j < pend[u]&&level_id[edges[j]] >= level-1;j ++) {
				assert(level_id[edges[j]] == level-1||level_id[edges[j]] == n);
				if(level_id[edges[j]] >= level) edges[end ++] = edges[j];
				else nonneighbors[non_neighbors_n ++] = edges[j];
			}
			assert(degree[u] == end-pstart[u]);
			for(ui j = 0;j < non_neighbors_n;j ++) edges[end ++] = nonneighbors[j];
			assert((end < pend[u]&&level_id[edges[end]] < level-1)||end == pend[u]);
#ifndef NDEBUG
			for(ui j = end;j < pend[u];j ++) {
				if(level_id[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", level_id[edges[j]], level);
				assert(level_id[edges[j]] < level);
			}
#endif
		}
	}

	void compute_a_heuristic_solution_and_prune(ui &R_end, ui level) {
		// the following computes the degeneracy ordering and a heuristic solution
#ifndef NDEBUG
		for(ui i = 0;i < R_end;i ++) {
			ui u = SR[i];
			assert(degree[u] + K > best_solution_size);
			ui end = pstart[u];
			while(end < pend[u]&&level_id[edges[end]] >= level) ++ end;
#ifndef NDEBUG
			for(ui j = end;j < pend[u];j ++) {
				if(level_id[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", level_id[edges[j]], level);
				assert(level_id[edges[j]] < level);
			}
#endif
			if(degree[u] != end - pstart[u]) printf("degree[u]: %u, %u\n", degree[u], end-pstart[u]);
			assert(degree[u] == end - pstart[u]);
		}
#endif
		ui *core = neighbors;
		ui *rid = nonneighbors;
		ui *id = buf;
		ui *t_degree = buf1;
		for(ui i = 0;i < R_end;i ++) {
			id[i] = 0;
			t_degree[SR[i]] = degree[SR[i]];
		}
		for(ui i = 0;i < R_end;i ++) ++ id[t_degree[SR[i]]];
		for(ui i = 1;i < R_end;i ++) id[i] += id[i-1];

		for(ui i = 0;i < R_end;i ++) rid[SR[i]] = -- id[t_degree[SR[i]]];
		for(ui i = 0;i < R_end;i ++) id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for(ui i = 0, j = 0;i <= R_end;i ++) {
			while(j < R_end&&t_degree[id[j]] < i) ++ j;
			degree_start[i] = j;
		}

		ui max_core = 0, pre_solution_size = best_solution_size;
		for(ui i = 0;i < R_end;i ++) {
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if(t_degree[u] > max_core) max_core = t_degree[u];
			core[u] = max_core;

			if(t_degree[u] + K >= R_end - i&&R_end - i > best_solution_size) {
				best_solution_size = R_end - i;
				for(ui j = i;j < R_end;j ++) best_solution[j-i] = id[j];
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			++ degree_start[t_degree[u]];
			if(t_degree[u] == 0) continue;

			degree_start[t_degree[u]-1] = degree_start[t_degree[u]];
			for(ui j = pstart[u];j < pend[u]&&level_id[edges[j]] >= level;j ++) if(rid[edges[j]] > i) {
				ui v = edges[j];
				ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
				std::swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++ degree_start[t_degree[v]];
				-- t_degree[v];
			}
		}

		assert(Qv.empty());
		for(ui i = 0;i < R_end;i ++) if(core[SR[i]] + K <= best_solution_size) {
			assert(level_id[SR[i]] > level);
			level_id[SR[i]] = level;
			Qv.push(SR[i]);
		}
		remove_vertices_with_prune(0, R_end, level);
	}

	void store_solution(ui size) {
		if(size <= best_solution_size) {
			printf("!!! the solution to store is no larger than the current best solution!");
			return ;
		}
		best_solution_size = size;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
	}

	bool is_kplex(ui R_end) {
		for(ui i = 0;i < R_end;i ++) if(degree[SR[i]] + K < R_end) return false;
		return true;
	}

	void BB_search(ui S_end, ui R_end, ui level, bool choose_zero) {
		if(S_end > best_solution_size) store_solution(S_end);
		if(R_end > best_solution_size&&is_kplex(R_end)) store_solution(R_end);
		if(R_end <= best_solution_size+1 || best_solution_size >= _UB_) return ;

#ifndef NDEBUG
		for(ui i = 0;i < S_end;i ++) {
			assert(degree[SR[i]] + K > best_solution_size);
			assert(degree_in_S[SR[i]] + K >= S_end);
		}
#endif

		ui old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());

		if(level > 1) reorganize_edges(S_end, R_end, level);
		if(S_end == 0) compute_a_heuristic_solution_and_prune(R_end, level);

		if(choose_zero&&SR_rid[0] < R_end&&!move_u_to_S_with_prune(0, S_end, R_end, level)) {
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return ;
		}

		// apply reduction rules
		ui S2_n = 0;
		for(ui i = 0;i < S_end;i ++) if(R_end - degree[SR[i]] > K) S2[S2_n++] = SR[i];
		if(S2_n >= 2) {
			collect_removable_vertices_based_on_total_edges(S2_n, S_end, R_end, level);
			if(!remove_vertices_with_prune(S_end, R_end, level)) {
				restore_SR(S_end, R_end, old_S_end, old_R_end, level);
				return ;
			}
		}

		// greedily add vertices to S
		if(!greedily_add_vertices_to_S(S_end, R_end, level)) {
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return ;
		}

		if(S_end > best_solution_size) store_solution(S_end);
		if(R_end > best_solution_size&&is_kplex(R_end)) store_solution(R_end);
		if(R_end <= best_solution_size+1 || best_solution_size >= _UB_) {
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return ;
		}

#ifndef NDEBUG
		for(ui i = 0;i < R_end;i ++) {
			ui u = SR[i], cnt = 0, cnt2 = 0;
			for(ui j = pstart[u];j < pend[u]&&level_id[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				++ cnt;
				if(SR_rid[edges[j]] < S_end) ++ cnt2;
			}
			assert(degree[u] == cnt);
			assert(degree_in_S[u] == cnt2);
		}
#endif

		ui u = choose_branch_vertex(S_end, R_end, level);
		assert(degree[u] + K > best_solution_size&&degree[u] + K > S_end);

		// the first branch includes u into S
		ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
		if(move_u_to_S_with_prune(u, S_end, R_end, level)) {
			BB_search(S_end, R_end, level+1, false);
		}
		if(best_solution_size >= _UB_) return ;
		assert(S_end == t_old_S_end + 1&&SR[S_end-1] == u);
		restore_SR(S_end, R_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
		for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif

		// the second branch exclude u from S
		assert(Qv.empty());
		ui v = n, candidates_n = 0; // the unique non-neighbor of u, and moreover v has exactly k non-neighbors
		ui *candidates = S2;
		if(degree[u]+2 == R_end) {
			ui neighbors_n = 0, non_neighbors_n = 0;
			get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, non_neighbors_n);
			assert(non_neighbors_n == 1);
			v = nonneighbors[0];
			if(degree[v]+K+1 != R_end) v = n;
			else {
				if(SR_rid[v] >= S_end) candidates[candidates_n++] = v;
				ui neighbors_n = 0, non_neighbors_n = 0;
				get_neighbors_and_non_neighbors(v, 0, R_end, level, neighbors_n, non_neighbors_n);
				for(ui i = 0;i < non_neighbors_n;i ++) if(nonneighbors[i] != u&&SR_rid[nonneighbors[i]] >= S_end) candidates[candidates_n++] = nonneighbors[i];
			}
		}
		bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
		if(succeed&&best_solution_size > pre_best_solution_size) succeed = collect_removable_vertices(S_end, R_end, level);
		if(succeed) succeed = remove_vertices_with_prune(S_end, R_end, level);
		if(succeed&&(v == n||greedily_add_nonneighbors(candidates, candidates_n, S_end, R_end, level))) {
			BB_search(S_end, R_end, level+1, false);
		}
		if(best_solution_size >= _UB_) return ;
		assert(S_end >= old_S_end&&R_end <= old_R_end);
		restore_SR(S_end, R_end, old_S_end, old_R_end, level);
	}

	void collect_removable_vertices_based_on_total_edges(ui S2_n, ui S_end, ui R_end, ui level) {
		vp.resize(R_end - S_end);
		ui max_nn = 0;
		for(ui i = S_end;i < R_end;i ++) {
			ui nn = 0;
			if(S2_n != S_end) {
				for(ept j = pstart[SR[i]];j < pend[SR[i]]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 1;
				for(ui j = 0;j < S2_n;j ++) if(!vis[S2[j]]) ++ nn;
				for(ept j = pstart[SR[i]];j < pend[SR[i]]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 0;
			}
			else nn = S_end - degree_in_S[SR[i]];
			if(nn > max_nn) max_nn = nn;
			vp[i-S_end].second = nn;
		}
		ui *cnt = neighbors;
		for(ui i = 0;i <= max_nn;i ++) cnt[i] = 0;
		for(ui i = 0;i < vp.size();i ++) ++ cnt[vp[i].second];
		for(ui i = 0;i < max_nn;i ++) cnt[i+1] += cnt[i];
		for(ui i = max_nn;i > 0;i --) cnt[i] = cnt[i-1];
		cnt[0] = 0;
		ui *ids = nonneighbors;
		for(ui i = 0;i < vp.size();i ++) ids[cnt[vp[i].second]++] = i;

		ui total_support = 0;
		for(ui i = 0;i < S2_n;i ++) total_support += K-S_end+degree_in_S[S2[i]];

		ui new_n = 0;
		while(!Qv.empty()) Qv.pop();
		for(ui i = 0;i < vp.size();i ++) {
			ui idx = ids[i], v = SR[S_end+ids[i]];
			ui t_support = total_support - vp[idx].second;
			for(ept j = pstart[v];j < pend[v]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 1;
			ui j = 0, v_support = K-1-S_end+degree_in_S[v], ub = S_end+1;
			while(true) {
				if(j == new_n) j = i+1;
				if(j >= vp.size()||ub > best_solution_size||ub + vp.size() - j <= best_solution_size) break;
				ui u = SR[S_end+ids[j]], nn = vp[ids[j]].second;
				if(t_support < nn || (nn != 0&&ub + (t_support/nn) <= best_solution_size)) break;
				if(vis[u]) {
					t_support -= nn;
					++ ub;
				}
				else if(v_support > 0) {
					-- v_support;
					t_support -= nn;
					++ ub;
				}
				++ j;
			}
			for(ept j = pstart[v];j < pend[v]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 0;
			if(ub <= best_solution_size) {
				level_id[v] = level;
				Qv.push(v);
			}
			else ids[new_n++] = ids[i];
		}
	}

	bool greedily_add_vertices_to_S(ui &S_end, ui &R_end, ui level) {
		while(true) {
			ui *candidates = S2;
			ui candidates_n = 0;
			for(ui i = S_end;i < R_end;i ++) {
				ui u = SR[i];
				if(R_end - degree[u] > K) continue;

				ui neighbors_n = 0, non_neighbors_n = 0;
				get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, non_neighbors_n);
				assert(non_neighbors_n < K);
				bool OK = true;
				for(ui j = 0;j < non_neighbors_n;j ++) if(R_end - degree[nonneighbors[j]] > K) {
					OK = false;
					break;
				}
				if(OK) candidates[candidates_n ++] = u;
			}

			if(!candidates_n) break;

			while(candidates_n) {
				ui u = candidates[--candidates_n];
				assert(SR_rid[u] >= S_end);
				if(SR_rid[u] >= R_end) return false;

				if(!move_u_to_S_with_prune(u, S_end, R_end, level)) return false;
			}
		}
		return true;
	}

	bool greedily_add_nonneighbors(ui *candidates, ui candidates_n, ui &S_end, ui &R_end, ui level) {
		while(candidates_n) {
			ui u = candidates[--candidates_n];
			assert(SR_rid[u] >= S_end);
			if(SR_rid[u] >= R_end||!move_u_to_S_with_prune(u, S_end, R_end, level)) return false;
		}
		return true;
	}



	void get_neighbors_and_non_neighbors(ui u, ui idx_start, ui idx_end, ui level, ui &neighbors_n, ui &non_neighbors_n) {
		neighbors_n = non_neighbors_n = 0;
		for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < idx_end) vis[edges[i]] = 1;
		for(ui i = idx_start;i < idx_end;i ++) if(SR[i] != u) {
			if(vis[SR[i]]) neighbors[neighbors_n++] = SR[i];
			else nonneighbors[non_neighbors_n++] = SR[i];
		}
		for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) vis[edges[i]] = 0;
	}

	bool move_u_to_S_with_prune(ui u, ui &S_end, ui &R_end, ui level) {
		assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
		assert(degree_in_S[u] + K > S_end);
#ifndef NDEBUG
		for(ui i = 0;i < S_end;i ++) assert(degree_in_S[SR[i]] + K >= S_end);
		for(ui i = S_end;i < R_end;i ++) assert(degree_in_S[SR[i]] + K > S_end);
#endif
		if(SR_rid[u] != S_end) swap_pos(S_end, SR_rid[u]);
		++ S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end-1);
		for(ui i = 0;i < neighbors_n;i ++) ++ degree_in_S[neighbors[i]];

		while(!Qv.empty()) Qv.pop();
		// reduction rules based on the fact that each vertex can have at most k-1 nonneighbors
		for(ui i = 0;i < nonneighbors_n;i ++) {
			ui v = nonneighbors[i];
			if(SR_rid[v] >= S_end) {
				if(level_id[v] == level) continue;
				if(S_end - degree_in_S[v] >= K||S_end - degree_in_S[u] == K) {
					level_id[v] = level;
					Qv.push(v);
				}
			}
			else if(S_end - degree_in_S[v] == K) {
				for(ept j = pstart[v];j < pend[v]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 1;
				for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!vis[SR[j]]) {
					level_id[SR[j]] = level;
					Qv.push(SR[j]);
				}
				for(ept j = pstart[v];j < pend[v]&&level_id[edges[j]] >= level;j ++) vis[edges[j]] = 0;
			}
		}
		return remove_vertices_with_prune(S_end, R_end, level);
	}

	bool remove_vertices_with_prune(ui S_end, ui &R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			bool terminate = false;
			ui neighbors_n = 0;
			for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
				ui w = edges[i];
				neighbors[neighbors_n++] = w;
				-- degree[w];
				if(degree[w] + K <= best_solution_size) {
					if(SR_rid[w] < S_end) terminate = true; // UB1
					else if(level_id[w] > level) { // RR3
						level_id[w] = level;
						Qv.push(w);
					}
				}
			}
			if(terminate) {
				for(ui i = 0;i < neighbors_n;i ++) ++ degree[neighbors[i]];
				level_id[u] = n;
				++ R_end;
				return false;
			}
		}

		return true;
	}

	void restore_SR(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

		for(;R_end < old_R_end;R_end ++) { // insert u back into R
			ui u = SR[R_end];
			assert(level_id[u] == level&&SR_rid[u] == R_end);
			level_id[u] = n;

			degree[u] = degree_in_S[u] = 0;
			for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) {
				ui w = edges[i];
				assert(SR[SR_rid[w]] == w);
				if(SR_rid[w] < R_end) {
					++ degree[w];
					++ degree[u];
				}
				if(SR_rid[w] < S_end) ++ degree_in_S[u];
			}
		}

		for(;S_end > old_S_end;S_end --) { // move u from S to R
			ui u = SR[S_end-1];
			assert(SR_rid[u] == S_end-1);

			for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) {
				ui w = edges[i];
				if(SR_rid[w] < R_end) -- degree_in_S[w];
			}
		}
	}

	bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
		assert(S_end);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool terminate = false;
		for(ept i = pstart[u];i < pend[u]&&level_id[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
			ui v = edges[i];
			-- degree_in_S[v];
			-- degree[v];
			if(degree[v] + K <= best_solution_size) {
				if(SR_rid[v] < S_end) terminate = true;
				else {
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if(terminate) return false;
		return true;
	}

	bool collect_removable_vertices(ui S_end, ui R_end, ui level) {
		for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K <= best_solution_size) return false;

		for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level){
			ui v = SR[i];
			if(S_end - degree_in_S[v] >= K||degree[v] + K <= best_solution_size) {
				level_id[v] = level;
				Qv.push(v);
				continue;
			}
		}

		return true;
	}

	void swap_pos(ui i, ui j) {
		std::swap(SR[i], SR[j]);
		SR_rid[SR[i]] = i;
		SR_rid[SR[j]] = j;
	}

	ui choose_branch_vertex(ui S_end, ui R_end, ui level) {
		ui *D = buf;
		ui D_n = 0;
		for(ui i = 0;i < R_end;i ++) if(R_end - degree[SR[i]] > K) D[D_n++] = SR[i];
		assert(D_n != 0);

		ui min_degree_in_S = n;
		for(ui i = 0;i < D_n;i ++) if(degree_in_S[D[i]] < min_degree_in_S) min_degree_in_S = degree_in_S[D[i]];
		ui u = n, min_degree =n;
		for(ui i = 0;i < D_n;i ++) if(degree_in_S[D[i]] == min_degree_in_S&&degree[D[i]] < min_degree) {
			min_degree = degree[D[i]];
			u = D[i];
		}
		assert(u != n);

		if(SR_rid[u] < S_end) {
			ui max_degree = 0, b = n;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(u, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for(ui i = 0;i < nonneighbors_n;i ++) if(degree[nonneighbors[i]] > max_degree) {
				max_degree = degree[nonneighbors[i]];
				b = nonneighbors[i];
			}
			return b;
		}
		else if(degree_in_S[u] < S_end||R_end - degree[u] > K+1) return u;
		else {
			ui max_degree = 0, w = n;
			for(ui i = S_end;i < R_end;i ++) if(degree[SR[i]] > max_degree) {
				max_degree = degree[SR[i]];
				w = SR[i];
			}
			if(degree[w]+1 >= R_end) {
				printf("!!! WA degree[w]: %u, R_end: %u\n", degree[w], R_end);
			}
			assert(degree[w]+1 < R_end);
			if(R_end-degree[w] == 2) return w;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(w, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for(ui i = 0;i < nonneighbors_n;i ++) if(R_end - degree[nonneighbors[i]] == K+1) return nonneighbors[i];
		}

		printf("!!! WA in choose_branch_vertex\n");
		return n;
	}
};
#endif
