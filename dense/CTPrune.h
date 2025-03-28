/*
 * CTPrune.h
 *
 *  Created on: 31 May 2023
 *      Author: ljchang@outlook.com
 */

#ifndef CTPRUNE_H_
#define CTPRUNE_H_

#include <cstring>
#include <algorithm>

using ui = unsigned int; // vertex type
// using ui = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges

namespace CTPrune {
	// orient graph
	void orient_graph(ui n, ui m, ui *peel_sequence, ui *pstart, ui *pend, ui *edges, ui *rid) {
		for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
		for(ui i = 0;i < n;i ++) {
			ui &end = pend[i] = pstart[i];
			for(ui j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
		}
	}

	// oriented triangle counting
	void oriented_triangle_counting(ui n, ui m, ui *pstart, ui *pend, ui *edges, ui *tri_cnt, ui *adj) {
		memset(adj, 0, sizeof(ui)*n);
		// long long cnt = 0;
		memset(tri_cnt, 0, sizeof(ui)*m);
		for(ui u = 0;u < n;u ++) {
			for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

			for(ui j = pstart[u];j < pend[u];j ++) {
				ui v = edges[j];
				for(ui k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
					++ tri_cnt[j];
					++ tri_cnt[k];
					++ tri_cnt[adj[edges[k]]-1];
					// ++ cnt;
				}
			}

			for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
		}
	}

	bool remove_and_shrink_oriented_tri(ui &n, ui &m, ui triangle_threshold, ui *out_mapping, ui *peel_sequence, ui *pstart, ui *pend, ui *edges, ui *tri_cnt, ui *rid, ui *degree) {
		for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];
		ui removed_edges = 0;
		for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pend[i];j ++) if(tri_cnt[j] < triangle_threshold) {
			-- degree[i]; -- degree[edges[j]];
			++ removed_edges;
		}

		if(removed_edges <= m/4) return false;

		ui cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
		ui t_cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[t_cnt++] = rid[peel_sequence[i]];

		ui pos = 0; cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			ui start = pstart[i];
			pstart[cnt] = pos;
			for(ui j = start;j < pend[i];j ++) if(tri_cnt[j] >= triangle_threshold) edges[pos ++] = rid[edges[j]];
			pend[cnt] = pos;
			pos = degree[i] + pstart[cnt];
			++ cnt;
		}
		pstart[cnt] = m = pos;
		n = cnt;

		return true;
	}

	// reorganize the adjacency lists
	// and sort each adjacency list to be in increasing order
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *pstart, ui *pend, ui *pend2, ui *edges, ui *edges_pointer, ui *buf) {
		for(ui i = 0;i < n;i ++) pend2[i] = pend[i];
		for(ui i = 0;i < n;i ++) {
			for(ui j = pstart[i];j < pend[i];j ++) {
				ui &k = pend2[edges[j]];
				edges[k] = i; tri_cnt[k] = tri_cnt[j];
				++ k;
			}
		}

		for(ui i = 0;i < n;i ++) {
			pend2[i] = pend[i];
			pend[i] = pstart[i];
		}
		for(ui i = 0;i < n;i ++) {
			for(ui j = pend2[i];j < pstart[i+1];j ++) {
				ui &k = pend[edges[j]];
				edges[k] = i; tri_cnt[k] = tri_cnt[j];
				edges_pointer[k] = j; edges_pointer[j] = k;
				++ k;
			}
		}

		ui *pointer = pend2;
		for(ui i = 0;i < n;i ++) {
			if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
			ui j = pstart[i], k = pend[i], pos = 0;
			while(j < pend[i]||k < pstart[i+1]) {
				if(k >= pstart[i+1]||(j < pend[i]&&edges[j] < edges[k])) {
					buf[pos] = edges[j];
					pointer[pos ++] = edges_pointer[j ++];
				}
				else {
					buf[pos] = edges[k];
					pointer[pos ++] = edges_pointer[k ++];
				}
			}
			for(ui j = 0;j < pos;j ++) {
				ui idx = pstart[i]+j, k = pointer[j];
				edges[idx] = buf[j];
				tri_cnt[idx] = tri_cnt[k];
				edges_pointer[idx] = k; edges_pointer[k] = idx;
			}
		}
	}

	bool find(ui w, ui b, ui e, char *deleted, ui &idx, ui *edges) {
		if(b >= e) return false;

		while(b+1 < e) {
			idx = b + (e-b)/2;
			if(edges[idx] > w) e = idx;
			else b = idx;
		}

		idx = b;
		if(edges[idx] == w&&!deleted[idx]) return true;
		return false;
	}

	void compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted, ui *pstart, ui *pend, ui *edges) {
		ui end = pstart[u];
		for(ui i = pstart[u];i < pend[u];i ++) if(!deleted[i]) {
			edges[end] = edges[i];
			tri_cnt[end] = tri_cnt[i];
			edges_pointer[end] = edges_pointer[i];
			edges_pointer[edges_pointer[end]] = end;
			deleted[end] = 0;
			++ end;
		}
		pend[u] = end;
	}

	void truss_peeling(ui degree_threshold, ui *Qv, ui &Qv_n, ui triangle_threshold, ui *Qe, ui Qe_n, ui *tri_cnt, ui *edges_pointer, char *deleted, ui *degree, ui *pstart, ui *pend, ui *edges) {
		for(ui j = 0;j < Qe_n;j += 2) {
			ui u = Qe[j], v = Qe[j+1], idx;
			find(v, pstart[u], pend[u], deleted, idx, edges);

			ui tri_n = tri_cnt[idx];
			deleted[idx] = deleted[edges_pointer[idx]] = 1;
			if(Qv != nullptr) {
				if(degree[u] == degree_threshold) Qv[Qv_n++] = u;
				if(degree[v] == degree_threshold) Qv[Qv_n++] = v;
			}
			-- degree[u]; -- degree[v];
			if(pend[u]-pstart[u] > degree[u]*2) compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			if(pend[v]-pstart[v] > degree[v]*2) compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			if(pend[u]-pstart[u] < pend[v]-pstart[v]) std::swap(u,v);

			if(pend[u]-pstart[u] > (pend[v]-pstart[v])*2) { // binary search
				for(ui k = pstart[v];k < pend[v];k ++) if(!deleted[k]) {
					if(tri_n&&find(edges[k], pstart[u], pend[u], deleted, idx, edges)) {
						-- tri_n;
						-- tri_cnt[edges_pointer[idx]];
						if( (tri_cnt[idx]--) == triangle_threshold) {
							Qe[Qe_n++] = u; Qe[Qe_n++] = edges[idx];
						}
						-- tri_cnt[edges_pointer[k]];
						if( (tri_cnt[k]--) == triangle_threshold) {
							Qe[Qe_n++] = v; Qe[Qe_n++] = edges[k];
						}
					}
				}
			}
			else { // sorted_merge
				ui ii = pstart[u], jj = pstart[v];
				while(ii < pend[u]&&jj < pend[v]) {
					if(edges[ii] == edges[jj]) {
						if(!deleted[ii]&&!deleted[jj]) {
							-- tri_n;
							-- tri_cnt[edges_pointer[ii]];
							if( (tri_cnt[ii]--) == triangle_threshold) {
								Qe[Qe_n++] = u; Qe[Qe_n++] = edges[ii];
							}
							-- tri_cnt[edges_pointer[jj]];
							if( (tri_cnt[jj]--) == triangle_threshold) {
								Qe[Qe_n++] = v; Qe[Qe_n++] = edges[jj];
							}
						}

						++ ii;
						++ jj;
					}
					else if(edges[ii] < edges[jj]) ++ ii;
					else ++ jj;
				}
			}
		}
	}

	void core_truss_copeeling(ui n, ui degree_threshold, ui triangle_threshold, ui *Qv, ui *Qe, ui *tri_cnt, ui *edges_pointer, char *deleted, char *exists, ui *degree, ui *pstart, ui *pend, ui *edges) {
		ui Qe_n = 0;
		for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pend[i];j ++) if(tri_cnt[j] < triangle_threshold&&edges[j] > i) {
			Qe[Qe_n++] = i; Qe[Qe_n++] = edges[j];
		}
		ui Qv_n = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] < degree_threshold) Qv[Qv_n++] = i;
		while(Qe_n || Qv_n) {
			while(Qe_n == 0&&Qv_n != 0) {
				ui u = Qv[-- Qv_n]; // delete u from the graph due to have a degree < degree_threshold
				if(degree[u] == 0) continue;

				if(pend[u]-pstart[u] != degree[u]) compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
				assert(pend[u] - pstart[u] == degree[u]);

				for(ui i = pstart[u];i < pend[u];i ++) deleted[i] = deleted[edges_pointer[i]] = exists[edges[i]] = 1;
				degree[u] = 0;

				for(ui i = pstart[u];i < pend[u];i ++) {
					ui v = edges[i];
					if( (degree[v]--) == degree_threshold) Qv[Qv_n++] = v;
					if(pend[v]-pstart[v] > degree[v]*2) compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges);

					for(ui j = pstart[v];j < pend[v];j ++) if(!deleted[j]&&edges[j] > v&&exists[edges[j]]) {
						-- tri_cnt[edges_pointer[j]];
						if( (tri_cnt[j]--) == triangle_threshold) {
							Qe[Qe_n++] = v, Qe[Qe_n++] = edges[j];
						}
					}
				}

				for(ui i = pstart[u];i < pend[u];i ++) exists[edges[i]] = 0;
			}
			truss_peeling(degree_threshold, Qv, Qv_n, triangle_threshold, Qe, Qe_n, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);
			Qe_n = 0;
		}
	}

	// reduce the graph to its maximal subgraph with and minimum edge triangle count at least triangle_threshold
	void truss_pruning(ui &n, ui &m, ui triangle_threshold, ui *peel_sequence, ui *out_mapping, ui *rid, ui *pstart, ui *edges, ui *degree, bool output) {
		if(triangle_threshold == 0) {
			printf("!!! Triangle_threshold is 0\n");
			return ;
		}

		ui *pend = new ui[n+1];
		orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);
		ui *tri_cnt = new ui[m];
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

		while(n&&remove_and_shrink_oriented_tri(n, m, triangle_threshold, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree)) {
			oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
		}

		ui *pend_buf = new ui[n+1];
		ui *edges_pointer = new ui[m];
		reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid);
		delete[] pend_buf; pend_buf = nullptr;

		for(ui i = 0;i < n;i ++) {
			pend[i] = pstart[i+1];
			degree[i] = pstart[i+1] - pstart[i];
		}
		ui *Qe = new ui[m];
		char *deleted = new char[m];
		memset(deleted, 0, sizeof(char)*m);
		ui Qe_n = 0;
		for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pend[i];j ++) if(tri_cnt[j] < triangle_threshold&&edges[j] > i) {
			Qe[Qe_n++] = i; Qe[Qe_n++] = edges[j];
		}
		truss_peeling(0, nullptr, Qe_n, triangle_threshold, Qe, Qe_n, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);

		ui cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
		ui t_cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[t_cnt++] = rid[peel_sequence[i]];
		ui pos = 0; cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			ui start = pstart[i];
			pstart[cnt] = pos;
			for(ui j = start;j < pend[i];j ++) if(!deleted[j]) edges[pos++] = rid[edges[j]];
			++ cnt;
		}
		pstart[cnt] = m = pos;
		n = cnt;

		delete[] Qe;
		delete[] deleted;
		delete[] edges_pointer;
		delete[] pend;
		delete[] tri_cnt;

		if(output) printf("*** After truss_pruning: n = %u, m = %lu (undirected)\n", n, m/2);
	}

	// reduce the graph to its maximal subgraph with minimum degree at least degree_threshold and minimum edge triangle count at least triangle_threshold
	void core_truss_copruning(ui &n, ui &m, ui degree_threshold, ui triangle_threshold, ui *peel_sequence, ui *out_mapping, ui *rid, ui *pstart, ui *edges, ui *degree, bool output) {
		if(triangle_threshold == 0) {
			printf("!!! Triangle_threshold is 0\n");
			return ;
		}
		if(degree_threshold <= triangle_threshold+1) {
			printf("!!! Degree_threshold <= triangle_threshold + 1, please invoke truss_pruning\n");
			return ;
		}

		ui *pend = new ui[n+1];
		orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);
		ui *tri_cnt = new ui[m];
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

		while(n&&remove_and_shrink_oriented_tri(n, m, triangle_threshold, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree)) {
			oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
		}

		ui *pend_buf = new ui[n+1];
		ui *edges_pointer = new ui[m];
		reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid);
		delete[] pend_buf; pend_buf = nullptr;

		for(ui i = 0;i < n;i ++) {
			pend[i] = pstart[i+1];
			degree[i] = pstart[i+1] - pstart[i];
		}
		ui *Qe = new ui[m];
		ui *Qv = new ui[n];
		char *deleted = new char[m];
		memset(deleted, 0, sizeof(char)*m);
		char *exists = new char[n];
		memset(exists, 0, sizeof(char)*n);
		core_truss_copeeling(n, degree_threshold, triangle_threshold, Qv, Qe, tri_cnt, edges_pointer, deleted, exists, degree, pstart, pend, edges);

		ui cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
		ui t_cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[t_cnt++] = rid[peel_sequence[i]];
		ui pos = 0; cnt = 0;
		for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
			ui start = pstart[i];
			pstart[cnt] = pos;
			for(ui j = start;j < pend[i];j ++) if(!deleted[j]) edges[pos++] = rid[edges[j]];
			++ cnt;
		}
		pstart[cnt] = m = pos;
		n = cnt;

		delete[] Qv;
		delete[] Qe;
		delete[] deleted;
		delete[] edges_pointer;
		delete[] pend;
		delete[] tri_cnt;

		if(output) printf("*** After core_truss_copruning: n = %u, m = %lu (undirected)\n", n, m/2);
	}

	void check_core_pruning(ui n, ui m, ui degree_threshold, ui* pstart, ui *edges) {
		for(ui i = 0;i < n;i ++) {
			if(pstart[i+1]-pstart[i] < degree_threshold) printf("!!! WA degree in check_core_pruning\n");
			for(ui j = pstart[i];j < pstart[i+1];j ++) if(edges[j] >= n) printf("!!! WA edge in check_core_pruning\n");
		}
	}

	void check_truss_pruning(ui n, ui m, ui triangle_threshold, ui* pstart, ui *edges, char *exists) {
		memset(exists, 0, sizeof(char)*n);
		for(ui i = 0;i < n;i ++) {
			for(ui j = pstart[i];j < pstart[i+1];j ++) exists[edges[j]] = 1;
			for(ui j = pstart[i];j < pstart[i+1];j ++) {
				ui v = edges[j], cn = 0;
				for(ui k = pstart[v];k < pstart[v+1];k ++) if(exists[edges[k]]) ++ cn;
				if(cn < triangle_threshold) printf("!!! WA triangle in check_truss_pruning\n");
			}
			for(ui j = pstart[i];j < pstart[i+1];j ++) exists[edges[j]] = 0;
		}
	}
}

#endif /* CTPRUNE_H_ */
