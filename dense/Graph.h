#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph
{
private:
	std::string dir; // input graph directory
	ui n;			 // number of nodes of the graph
	ept m;			 // number of edges of the graph
	ui K;			 // the value of k in k-plex

	ept *pstart; // offset of neighbors of nodes
	ept *pend;	 // used in search
	ept *pend_buf;
	ui *edges; // adjacent ids of edges

	std::vector<ui> kplex;

	ui *s_degree;
	ept *s_pstart;
	ept *s_pend;
	ui *s_edges;
	ui *s_peel_sequence;
	ui *s_core;
	char *s_vis;
	ListLinearHeap *s_heap;

	ui *s_edgelist_pointer;
	ui *s_tri_cnt;
	ui *s_edge_list;
	ui *s_active_edgelist;
	char *s_deleted;
	bool dense_search = false;

public:
	ui best_n_edges=0;
	Graph(const char *_dir, const int _K);
	~Graph();

	void read();
	void read_csv();
	void read_graph();

	void write();
	void verify_kplex();

	void search();
	void search_dense();

private:
	void write_subgraph(ui n, const std::vector<std::pair<ui, ui > > &edge_list);
	void heuristic_kplex_max_degree(ui processed_threshold);
	void extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui > > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);
	void extract_graph(ui n, ui m, ui *deg, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui > > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui > > &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer);

	ui degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output);
	void shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *pstart, ui *edges, bool output);
	void oriented_triangle_counting(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj);
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf);
	ept peeling(ui critical_vertex, ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists);
	char find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges);

	// functions for subgraph processing
	void load_graph_from_edgelist(ui _n, const std::vector<std::pair<ui, ui > > &edge_list, ui &n, ept &m, ui *degree, ept *pstart, ui *edges);
	void subgraph_prune(ui *ids, ui &_n, std::vector<std::pair<ui, ui > > &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists);
	void ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, ui *edgelist_pointer, bool output);
	void write(vector<ui> &kplex, ui *pstart, ui* pend, ui *edges, bool append = false)
	{

		FILE *fout = append ? Utility::open_file("dense_kplexes.txt", "a") : Utility::open_file("dense_kplexes.txt", "w");
		ui ne = 0;
		for (ui u : kplex)
			fprintf(fout, "%u ", u);
		fprintf(fout, "\n");

		for (ui u : kplex)
		{
			fprintf(fout, "%u ", u);
			ui *st = edges + pstart[u];
			ui *en = edges + pend[u];
			for (ui v : kplex)
			{
				if (binary_search(st, en, v))
				{
					fprintf(fout, "%u ", v);
					ne++;
				}
			}
			fprintf(fout, "\n");
		}
		fprintf(fout, "total no. of edges: %u\n", ne);
		fprintf(fout, "\n");
		fclose(fout);
	}
};
#endif
