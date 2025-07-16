//connected dominating set interdiction
//2024/6/24 version
//process blocking node prior to not blocking node
//lower bound in bote associated wite block and non block
//UB uses Stoer_Wagner's algoritem
//no infeasibility check when block an edge
//worst bound first search
//infeasibility check with updated ST identification

#include <iostream>

#include <vector>
using namespace std;
//To be able to use vector
#include<queue>
//To be able to use queue
#include <list>
//To be able to use forwardlist
#include <list>
//To be able to use list
#include <utility> 
//to be able to use pair
#include <fstream>
//To be able to use fin
#include <sstream>
//to be able to read tee elements of an string (istringstream comamnd)
#include <time.h>

#include <algorithm>    // std::max
//used in main searce tree

class edge
{
public:
	int head;
	int tail;
	double blocking_cost;
	edge(int ed, int tl, double bc)
	{
		head = ed;
		tail = tl;
		blocking_cost = bc;
	}

	edge()
	{
		head = 0;
		tail = 0;
		blocking_cost = 0;
	}
};

class CDSI_instance
{
public:
	vector<double> vw;//weiget of vertices
	vector <vector<int>> two_d;//store edge id correspond to one_d
	vector<vector<int>> neigebor;//store tee neigebor of eace vertex
	//==-1 means no edge
	vector <edge> one_d;
	double desired_ST_weiget;
	double gamma;

	CDSI_instance()
	{
		vw = vector<double>();
		two_d = vector <vector<int>>();
		one_d = vector<edge>();
		desired_ST_weiget = -1;
		gamma = 0;
	}
};

class node_EST
{
public:
	vector<int> solution;//format of oned edge
	vector<int> candidate; //format of oned edge
	vector<int> degree;
	vector<char> visitedV;//'T' if vertex is visited
	int weight;
	node_EST(vector<int>& bs, vector<int>& cs, vector<int>& m, int w)
	{
		solution = bs;
		candidate = cs;
		degree = m;
		weight = w;
	}
	node_EST(int s)
	{
		visitedV.resize(s, 'F');
		weight = 0;
	}
	node_EST()
	{
		weight = 0;
	}
};

class graph
{
public:
	vector <vector<double>> two_d;//store edge costs of combined nodes
	vector<vector<int>> vset;//store tee vertices (ID of instance) in eace combined nodes

	graph()
	{
		two_d = vector<vector<double>>();
		vset = vector <vector<int>>();
	}
};

void d_sort_ratio(vector<pair<int, double>>& candidate);

class node
{
public:
	vector<int> status;//0 block,1 not block,2 candidate,format of oned
	vector<int> candidate;//id of oned
	vector<int> Splusq;
	int e;
	double solution_cost;
	bool b;//true if it's associated wite block e
	double LB;
	//Parametrized Constructor
	node(vector<int>& ub, vector<int>& cc, int ve, double sc, bool vb, double lb)
	{
		status = ub;
		candidate = cc;
		e = ve;
		solution_cost = sc;
		b = vb;
		LB = lb;
	}
	node(CDSI_instance& instance)
	{//create root node
		double temp_c = 1e7;
		vector<pair<int, double>> tke(instance.one_d.size());

		for (int i = 0; i < instance.one_d.size(); i++) {
			tke[i] = { i, instance.one_d[i].blocking_cost };
			temp_c = min(temp_c, instance.one_d[i].blocking_cost);
		}

		d_sort_ratio(tke);

		candidate.resize(instance.one_d.size());
		transform(tke.begin(), tke.end(), candidate.begin(), [](const pair<int, double>& p) { return p.first; });

		LB = temp_c;

		status.resize(instance.one_d.size(), 2);
		e = -1;
		solution_cost = 0;
		b = true;
	}
	node()
	{
		e = -1;
		solution_cost = 0;
	}
};

class EST
{//store an EST
public:
	vector<int> EST1;//id of oned in instance
	double cost;//cost of minimum cost edge
	vector<char> EST2;//T if element is in the EST
	int me;//minimum cost edge
	//Parametrized Constructor
	EST(vector<int> egc1, int sc, vector<char> egc2, int e)
	{
		me = e;
		cost = sc;
		EST1 = egc1;
		EST2 = egc2;
	}
	EST()
	{
		me = -1;
		cost = 0;
	}
};

void readdata(CDSI_instance& instance, string filename);

void runmodel(string);

double MINIMUMCUT(CDSI_instance& instance, vector<int>& edge_MC);
void MINIMUMCUTPHASE(CDSI_instance& instance, graph& G, vector<int>& edge_MC, double& MC, int a);

int c_brance(CDSI_instance& instance, int t, vector<node_EST>& tree);

void branch(CDSI_instance& instance, list<node>& tree, int& node_count);
//return true if tee candidate set of current node is feasibl
char ContainEST(CDSI_instance& instance, node& tq, vector<EST>& estlist, clock_t start,double);
bool isconnectCon(CDSI_instance& instance, vector<int>& he);
bool find_EST_nb(vector<EST>& estlist, vector<int>& nb, int eq);
bool find_EST_CBB(CDSI_instance& instance, vector<int>& Splusq, vector<int>&, int eq, clock_t start,double);
bool check_EST(CDSI_instance& instance, vector<int>& edges, vector<vector<char>>&);

double lower_bound(CDSI_instance& instance, node q, vector<EST>& estlist, clock_t start);
int find_EST_union(vector<EST>& estlist, vector<int>& he, int start);
int insert_EST(CDSI_instance& instance, vector<EST>& estlist, vector<int>& ste);
bool isconnectLB(CDSI_instance& instance, vector<int>& he, int nhe);
vector<int> find_ST(CDSI_instance& instance, vector<int>& he, vector<int>& heud);
bool check_STW(CDSI_instance& instance, vector<int>& STE);

double LB_unprocessed_nodes(CDSI_instance& instance, list<node>& tree);

int main()
{
	ofstream fout;
	fout.open("CDSI_CBB_v2.csv");
	fout << "File Name,|V|,density,Eta,number_of_EST,number_of_cycle,Time Limit (s),Run Time (s),Gap,Best Sol,mincut,Alternative solution, Best Bound, Nodes #,processed_node,fathomed by bound,fathomed by infeasibility, fathomed by feasibility,update solution" << endl;
	fout.close();

	ifstream fin;
	fin.open("filename-real.txt");
	string str3;
	int number_of_files;
	fin >> number_of_files;
	//finish reading the current line
	getline(fin, str3);
	for (int i = 0; i < number_of_files; i++)
	{
		getline(fin, str3);
		runmodel(str3);
	}
	fin.open("filename_random.txt");
	fin >> number_of_files;
	//finish reading the current line
	getline(fin, str3);
	for (int i = 0; i < number_of_files; i++)
	{
		getline(fin, str3);
		runmodel(str3);
	}
	return 0;
}

void runmodel(string filename)
{
	double total_time_limit = 3600.0;
	CDSI_instance instance;
	readdata(instance, filename);

	cout << filename << " start" << endl;
	clock_t pre_process_start = clock();
	list<node> tree;
	tree.push_back(node(instance));

	vector<EST> EST_list;
	vector<int> incumbent_solution;
	vector<int> UB_solution;
	double mincut = MINIMUMCUT(instance, UB_solution);
	double incumbent_cost = mincut;
	vector<int> min_cut_S;
	min_cut_S.resize(instance.one_d.size(), 1);
	for (int i = 0; i < UB_solution.size(); i++)
	{
		min_cut_S[i] = 0;
	}
	clock_t now = clock();

	int node_count = 1;//the root node
	int fathomed_by_bound = 0;
	int fathomed_by_infeasibility = 0;
	int fathomed_by_feasibility = 0;
	int update_solution = 0;
	//current_node = branch(instance, current_node, tree, node_count);
	int processed_node = 1;
	int aos = 0;
	double rest_LB = 0;
	bool solved = true;
	while (tree.size() != 0)
	{
		//check time limit
		clock_t now = clock();
		if (static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC > total_time_limit)
		{
			double time = static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC;
			cout << filename << ' ' << time << "s, TL, start to calculate gap." << endl;
			//calculate lower bound for all unprocessed nodes
			rest_LB = LB_unprocessed_nodes(instance, tree);

			solved = false;
			break;
		}
		char infeasibility = 'N';

		if (!tree.front().b)
		{//associated wite not blocking
			infeasibility = ContainEST(instance, tree.front(), EST_list, pre_process_start,total_time_limit);
			if (infeasibility == 'L')
			{//time limit
				cout << filename << ' ' << static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC << "s, TL, start to calculate gap." << endl;
				//calculate lower bound for all unprocessed nodes
				rest_LB = LB_unprocessed_nodes(instance, tree);
				solved = false;
				break;
			}
		}
		if (infeasibility == 'N')
		{//contain EST
			if (tree.front().candidate.size() == 0)
			{//fateom by feasibility
				fathomed_by_feasibility++;
				if (tree.front().solution_cost == incumbent_cost)
				{
					bool equal = true;
					for (int i = 0; i < min_cut_S.size(); i++)
					{
						if (min_cut_S[i] != tree.front().status[i])
						{
							equal = false;
							break;
						}
					}
					if (!equal)
					{
						aos++;
					}
				}
				if (tree.front().solution_cost < incumbent_cost)
				{
					incumbent_cost = tree.front().solution_cost;
					incumbent_solution = tree.front().status;
					update_solution++;
				}
				tree.pop_front();
			}
			else
			{
				double newLB = lower_bound(instance, tree.front(), EST_list, pre_process_start);
				if (newLB > tree.front().LB)
				{
					tree.front().LB = newLB;
				}
				if (tree.front().LB < incumbent_cost)
				{
					branch(instance, tree, node_count);
				}
				else//fateom by bound
				{
					fathomed_by_bound++;
					tree.pop_front();
				}
			}
		}
		else
		{
			fathomed_by_infeasibility++;
			tree.pop_front();
		}
		processed_node++;
	}
	clock_t end = clock();
	cout << filename << " end" << endl;
	ofstream fout;
	fout.open("CDSI_CBB_v2.csv", std::ios_base::app);
	fout << filename << ',' << instance.vw.size();
	fout << ',' << static_cast<double>(instance.one_d.size()) / static_cast<double>((instance.vw.size() * (instance.vw.size() - 1)) / 2);
	fout << ',' << instance.desired_ST_weiget;
	fout << ',' << EST_list.size() << ",0";
	fout << ',' << total_time_limit << ",";

	if (solved)
	{
		fout << static_cast<double>(end - pre_process_start) / CLOCKS_PER_SEC;
		fout << ",0,";
		fout << incumbent_cost << ',';
		fout << mincut << ",";
		fout << aos << ',';
		fout << incumbent_cost << ",";
	}
	else
	{
		fout << "TL" << ',';
		fout << static_cast<double> (incumbent_cost - rest_LB) / (incumbent_cost) << ',';
		fout << incumbent_cost << ',';
		fout << mincut << ",";
		fout << aos << ',';
		fout << rest_LB << ",";
	}
	fout << node_count << ',';
	fout << processed_node << ',';
	fout << fathomed_by_bound << ',';
	fout << fathomed_by_infeasibility << ',';
	fout << fathomed_by_feasibility << ',';
	fout << update_solution << ',';
	if (update_solution > 0)
	{
		for (int i = 0; i < incumbent_solution.size(); i++)
		{
			if (incumbent_solution[i] == 0)
			{
				fout << '(' << instance.one_d[i].head + 1 << ';';
				fout << instance.one_d[i].tail + 1 << ')';
			}
		}
	}
	else
	{
		for (int i = 0; i < UB_solution.size(); i++)
		{
			fout << '(' << instance.one_d[UB_solution[i]].head + 1 << ';';
			fout << instance.one_d[UB_solution[i]].tail + 1 << ')';
		}
	}
	fout << ',' << static_cast<double>(end - pre_process_start) / CLOCKS_PER_SEC;
	fout << endl;
	fout.close();
}

double MINIMUMCUT(CDSI_instance& instance, vector<int>& edge_MC)
{
	graph G;
	G.vset.resize(instance.vw.size());
	vector<vector<double>> edge_weiget(instance.two_d.size());
	for (int i = 0; i < instance.two_d.size(); i++)
	{
		G.vset[i].assign(1, i);
		edge_weiget[i].resize(instance.two_d[i].size(), 0);
	}
	for (int i = 0; i < instance.two_d[i].size() - 1; i++)
	{
		for (int j = i + 1; j < instance.two_d[i].size(); j++)
		{
			if (instance.two_d[i][j] != -1)
			{
				edge_weiget[i][j] = instance.one_d[instance.two_d[i][j]].blocking_cost;
				edge_weiget[j][i] = edge_weiget[i][j];
			}
		}
	}
	G.two_d = edge_weiget;

	double minimumcut = 100 * instance.one_d.size();
	while (G.vset.size() > 1)
	{
		MINIMUMCUTPHASE(instance, G, edge_MC, minimumcut, 1);
	}
	return minimumcut;
}
void MINIMUMCUTPHASE(CDSI_instance& instance, graph& G, vector<int>& edge_MC, double& MC, int a)
{//G store edge weiget, a is starting vertex
	vector<int> A;
	list<int> VminusA;
	A.push_back(a);
	for (int i = 0; i < G.vset.size(); i++)
	{
		if (i != a)
		{
			VminusA.push_front(i);
		}
	}
	while (A.size() < G.vset.size() - 1)
	{//ween a size equal to g size-1,we need to process tee last one
		double MaxWeiget = 0;
		list<int>::iterator MaxV;//point to v before max
		for (list<int>::iterator it = VminusA.begin(); it != VminusA.end(); it++)
		{
			double temp_w = 0;
			for (int i = 0; i < A.size(); i++)
			{
				if (G.two_d[A[i]][*it] > 0)
				{
					temp_w += G.two_d[A[i]][*it];
				}
			}
			if (temp_w > MaxWeiget)
			{
				MaxWeiget = temp_w;
				MaxV = it;
			}
		}
		//add to A tee most tigetly connected vertex
		A.push_back(*MaxV);
		VminusA.erase(MaxV);
	}
	//store tee cut - of - tee - pease and serink G by merging tee two vertices added last
	double cutofthephase = 0;
	int s = A.back();//the vertex before last
	int t = VminusA.front();//last vertex in V
	//process cut between s and t
	cutofthephase += G.two_d[s][t];
	for (int i = 0; i < A.size(); i++)
	{
		if (G.two_d[A[i]][t] > 0 && A[i] != s)
		{
			cutofthephase += G.two_d[A[i]][t];
			G.two_d[A[i]][s] += G.two_d[A[i]][t];
			G.two_d[s][A[i]] = G.two_d[A[i]][s];
		}
	}

	if (cutofthephase < MC)
	{
		vector<int> edgeofcut;//edges in the cut wite id of oned in instance
		MC = cutofthephase;
		for (int i = 0; i < G.vset[t].size(); i++)
		{
			int ti = G.vset[t][i];
			for (int j = 0; j < A.size(); j++)
			{
				if (G.two_d[A[j]][t] > 0)
				{
					for (int k = 0; k < G.vset[A[j]].size(); k++)
					{
						if (instance.two_d[ti][G.vset[A[j]][k]] != -1)
						{
							edgeofcut.push_back(instance.two_d[ti][G.vset[A[j]][k]]);
						}
					}
				}
			}
		}
		edge_MC = edgeofcut;
	}
	G.vset[s].insert(G.vset[s].end(), G.vset[t].begin(), G.vset[t].end());
	G.vset[t] = G.vset.back();
	G.vset.pop_back();
	for (int i = 0; i < G.two_d.size(); i++)
	{
		G.two_d[i][t] = G.two_d[i].back();
		G.two_d[i].pop_back();
	}
	G.two_d[t] = G.two_d.back();
	G.two_d.pop_back();
}

void readdata(CDSI_instance& instance, string filename)
{//return number of edges
	ifstream fin;
	fin.open(filename);
	string str3;
	int number_of_points, number_of_edges;
	fin >> number_of_points;
	fin >> number_of_edges;
	fin >> instance.desired_ST_weiget;
	instance.two_d.resize(number_of_points);
	instance.neigebor.resize(number_of_points);
	instance.vw.resize(number_of_points);
	for (int i = 0; i < number_of_points; i++)
	{
		instance.two_d[i].resize(number_of_points, -1);
	}
	getline(fin, str3);//finise current line
	for (int i = 0; i < number_of_points; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		double temp_w;
		ss0 >> temp_w;
		instance.vw[i] = temp_w;
	}
	int one_code = 0;
	for (int i = 0; i < number_of_edges; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		int head;
		ss0 >> head;
		head--;
		int tail;
		ss0 >> tail;
		tail--;
		double temp_cost;
		ss0 >> temp_cost;

		if (instance.two_d[head][tail] == -1)
		{
			instance.one_d.push_back(edge(head, tail, temp_cost));
			instance.two_d[head][tail] = one_code;
			instance.two_d[tail][head] = one_code;
			instance.neigebor[head].push_back(tail);
			instance.neigebor[tail].push_back(head);
			one_code++;
		}
	}
}

void branch(CDSI_instance& instance, list<node>& tree, int& node_count)
{//return position of tree node
	int e = tree.front().candidate.back();
	tree.front().candidate.pop_back();
	tree.front().e = e;
	//create tq associated wite not blocking e
	node tqnode = tree.front();
	tqnode.Splusq.push_back(e);
	tqnode.status[e] = 1;
	tqnode.b = false;
	//create tp associated with blocking e
	node tpnode = tree.front();
	tpnode.solution_cost += instance.one_d[e].blocking_cost;
	tpnode.status[e] = 0;
	tpnode.b = true;
	tree.pop_front();
	node_count += 2;//block and not block
	for (auto it = tree.begin(); it != tree.end(); it++)
	{
		if ((*it).LB >= tpnode.LB)
		{
			tree.insert(it, tqnode);
			tree.insert(it, tpnode);
			return;
		}
	}
	tree.push_back(tpnode);
	tree.push_back(tqnode);
}

char ContainEST(CDSI_instance& instance, node& tq, vector<EST>& estlist, clock_t start,double ttl)
{//return L if reach time limit, C if S+q contain est, N if does not contain
	if (tq.Splusq.size() < instance.vw.size() - 1)
	{
		return 'N';
	}
	if (!isconnectCon(instance, tq.status))
	{
		return 'N';
	}
	if (check_STW(instance, tq.Splusq))
	{
		return 'C';
	}
	if (find_EST_nb(estlist, tq.status, tq.e))
	{
		return 'C';
	}

	vector<int> e;
	bool TL = find_EST_CBB(instance, tq.Splusq, e, tq.e, start,ttl);
	if (TL)
	{//time limit reach
		return 'L';
	}
	else if (e.size() == 0)
	{
		return 'N';
	}
	else
	{
		insert_EST(instance, estlist, e);
		return 'C';
	}
}
bool isconnectCon(CDSI_instance& instance, vector<int>& he)
{//return true if connected
	vector<char> visitedV;//bool vector
	visitedV.resize(instance.vw.size(), 'F');
	queue<int> unvisitedV;
	unvisitedV.push(0);
	int countv = 0;
	while (!unvisitedV.empty())
	{
		int cv = unvisitedV.front();
		unvisitedV.pop();
		if (visitedV[cv] == 'F')
		{
			visitedV[cv] = 'T';
			countv++;
			for (int i = 0; i < instance.neigebor[cv].size(); i++)
			{
				if (visitedV[instance.neigebor[cv][i]] == 'F' && he[instance.two_d[cv][instance.neigebor[cv][i]]] == 1)
				{
					unvisitedV.push(instance.neigebor[cv][i]);
				}
			}
		}
	}
	if (countv == instance.vw.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool find_EST_nb(vector<EST>& estlist, vector<int>& nb, int eq)
{//check if there exist EST in S+q containing eq, return true if find EST
	for (int i = 0; i < estlist.size(); i++)
	{
		if (estlist[i].EST2[eq] == 'T')
		{
			bool contain = true;
			for (int j = 0; j < estlist[i].EST1.size(); j++)
			{
				if (nb[estlist[i].EST1[j]] != 1)
				{
					contain = false;
					break;
				}
			}
			if (contain)
			{
				return true;
			}
		}
	}
	return false;
}
bool find_EST_CBB(CDSI_instance& instance, vector<int>& Splusq, vector<int>& solution, int eq, clock_t start,double ttl)
{//return true if reach time limit
	vector<node_EST> tree;
	//initialize the root node
	vector<int> edges = Splusq;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i] == eq)
		{
			edges[i] = edges.back();
			edges.pop_back();
			break;
		}
	}
	tree.push_back(node_EST());
	tree[0].candidate = edges;
	tree[0].solution.push_back(eq);
	tree[0].visitedV.resize(instance.vw.size(), 'F');
	tree[0].visitedV[instance.one_d[eq].head] = 'T';
	tree[0].visitedV[instance.one_d[eq].tail] = 'T';
	tree[0].degree.resize(instance.vw.size(), 0);
	tree[0].degree[instance.one_d[eq].head] = 1;
	tree[0].degree[instance.one_d[eq].tail] = 1;
	//searce tree
	int current_node = 0;
	while (current_node != -1)
	{
		clock_t now = clock();
		if (static_cast<double>(now - start) / CLOCKS_PER_SEC > ttl)
		{
			return true;
		}
		if (tree[current_node].solution.size() == instance.vw.size() - 1)
		{
			solution = tree[current_node].solution;
			return false;
		}
		else
		{
			current_node = c_brance(instance, current_node, tree);
		}
	}
	return false;
}
int c_brance(CDSI_instance& instance, int t, vector<node_EST>& tree)
{
	if (tree[t].candidate.size() == 0)
	{
		return t - 1;
	}
	if (tree[t].solution.size() + tree[t].candidate.size() < instance.vw.size() - 1)
	{//not enough to make a spanning tree;
		return t - 1;
	}
	double min_weight = 10000;
	int min_e = -1;
	int min_v = -1;
	int min_nv = -1;
	for (int i = 0; i < tree[t].candidate.size(); i++)
	{
		edge temp_e = instance.one_d[tree[t].candidate[i]];
		if (tree[t].visitedV[temp_e.head] != tree[t].visitedV[temp_e.tail])
		{
			int notvisit = (tree[t].visitedV[temp_e.head] == 'T') ? temp_e.tail : temp_e.head;
			int visit = (tree[t].visitedV[temp_e.head] == 'T') ? temp_e.head : temp_e.tail;
			if (instance.vw[notvisit] < min_weight)
			{
				min_weight = instance.vw[notvisit];
				min_v = visit;
				min_nv = notvisit;
				min_e = i;
			}
			else if (instance.vw[notvisit] == min_weight)
			{
				if (tree[t].degree[visit] > tree[t].degree[min_v])
				{
					min_v = visit;
					min_nv = notvisit;
					min_e = i;
				}
			}
		}
		else if (tree[t].visitedV[temp_e.head] == 'T')
		{
			tree[t].candidate[i] = tree[t].candidate.back();
			tree[t].candidate.pop_back();
			i--;
			if (tree[t].candidate.empty()) {
				return t - 1;
			}
		}
	}
	if (min_e == -1)
	{
		return t - 1;
	}
	tree.resize(t + 2);
	tree[t + 1] = tree[t];
	tree[t].candidate[min_e] = tree[t].candidate.back();
	tree[t].candidate.pop_back();
	tree[t + 1].degree[min_v]++;
	if (tree[t + 1].degree[min_v] == 2)
	{
		tree[t + 1].weight += instance.vw[min_v];
		if (tree[t + 1].weight >= instance.desired_ST_weiget)
		{
			return t;
		}
	}
	tree[t + 1].solution.push_back(tree[t + 1].candidate[min_e]);
	tree[t + 1].candidate[min_e] = tree[t + 1].candidate.back();
	tree[t + 1].candidate.pop_back();
	tree[t + 1].degree[min_nv]++;
	tree[t + 1].visitedV[min_nv] = 'T';
	return t + 1;
}
bool check_EST(CDSI_instance& instance, vector<int>& edges, vector<vector<char>>& edge_m)
{//return true if weight less than eta
	double weight = 0;
	vector<int> degree(instance.vw.size());
	int h, t;
	for (int i = 0; i < edges.size(); i++)
	{
		h = instance.one_d[edges[i]].head;
		t = instance.one_d[edges[i]].tail;
		degree[h]++;
		degree[t]++;
	}
	//check weight
	for (int i = 0; i < degree.size(); i++)
	{
		if (degree[i] > 1)
		{
			weight += instance.vw[i];
		}
	}
	if (weight >= instance.desired_ST_weiget)
	{
		return false;
	}
	else
	{
		return true;
	}
}

double lower_bound(CDSI_instance& instance, node q, vector<EST>& estlist, clock_t start)
{
	vector<int> he = q.status;//1 if edge in H S+, 2 if in H Uq,format of oned
	double c_t = 0;
	vector<int> heud = q.candidate;//id of edge in H that is in undecided set
	int nhe = q.Splusq.size() + q.candidate.size();//number of edges in H

	int ESTid = find_EST_union(estlist, he, 0);
	while (ESTid < estlist.size())
	{
		bool egceud = false;
		if (he[estlist[ESTid].me] == 2)
		{//is in Ud
			c_t += estlist[ESTid].cost;
			for (int i = 0; i < heud.size(); i++)
			{
				if (estlist[ESTid].EST2[heud[i]])
				{//edge in EST and edge is in ud
					he[heud[i]] = 0;//remove from H
					nhe--;
					heud[i] = heud.back();
					heud.pop_back();
					i--;
				}
			}
		}
		else
		{
			double c_min = 100 * instance.one_d.size();
			for (int i = 0; i < heud.size(); i++)
			{
				if (estlist[ESTid].EST2[heud[i]])
				{//edge in EST and edge is in ud
					he[heud[i]] = 0;//remove from H
					nhe--;
					if (instance.one_d[heud[i]].blocking_cost < c_min)
					{
						c_min = instance.one_d[heud[i]].blocking_cost;
					}
					heud[i] = heud.back();
					heud.pop_back();
					i--;
				}
			}
			c_t += c_min;
		}
		ESTid = find_EST_union(estlist, he, ESTid + 1);
	}

	while (heud.size() > 0 && isconnectLB(instance, he, nhe))
	{
		vector<int> STE = find_ST(instance, he, heud);
		if (check_STW(instance, STE))
		{
			int estpos = insert_EST(instance, estlist, STE);
			if (q.status[estlist[estpos].me] == 2)
			{//is Ud
				c_t += estlist[estpos].cost;
				for (int i = 0; i < heud.size(); i++)
				{
					if (estlist[estpos].EST2[heud[i]])
					{//edge in EST and edge is in ud
						he[heud[i]] = 0;//remove from H
						nhe--;
						heud[i] = heud.back();
						heud.pop_back();
						i--;
					}
				}
			}
			else
			{
				double c_min = 100 * instance.one_d.size();
				for (int i = 0; i < heud.size(); i++)
				{
					if (estlist[estpos].EST2[heud[i]])
					{//edge in EST and edge is in ud
						he[heud[i]] = 0;//remove from H
						nhe--;
						if (instance.one_d[heud[i]].blocking_cost < c_min)
						{
							c_min = instance.one_d[heud[i]].blocking_cost;
						}
						heud[i] = heud.back();
						heud.pop_back();
						i--;
					}
				}
				c_t += c_min;
			}
		}
		else
		{
			for (int i = 0; i < STE.size(); i++)
			{
				if (he[STE[i]] == 2)
				{//edge in EST and edge is in ud
					he[STE[i]] = 0;//remove from H
					nhe--;
					int ii;
					for (int j = 0; j < heud.size(); j++)
					{
						if (heud[j] == STE[i])
						{
							ii = j;
							break;
						}
					}
					heud[ii] = heud.back();
					heud.pop_back();
				}
			}
		}
	}
	return c_t + q.solution_cost;
}
int find_EST_union(vector<EST>& estlist, vector<int>& he, int start)
{//check if there exist EST in given edges, from start in the list
	for (int i = start; i < estlist.size(); i++)
	{
		bool contain = true;
		for (int j = 0; j < estlist[i].EST1.size(); j++)
		{
			if (he[estlist[i].EST1[j]] == 0)
			{//element in EST is not contained
				contain = false;
				break;
			}
		}
		if (contain)
		{
			return i;//contain this EST
		}
	}
	return estlist.size();
}
bool isconnectLB(CDSI_instance& instance, vector<int>& he, int nhe)
{//return true if connected
	if (nhe < instance.vw.size() - 1)
	{
		return false;
	}
	vector<char> visitedV;//bool vector
	visitedV.resize(instance.vw.size(), 'F');
	queue<int> unvisitedV;
	unvisitedV.push(0);
	int countv = 0;
	while (!unvisitedV.empty())
	{
		int cv = unvisitedV.front();
		unvisitedV.pop();
		if (visitedV[cv] == 'F')
		{
			visitedV[cv] = 'T';
			countv++;
			for (int i = 0; i < instance.neigebor[cv].size(); i++)
			{
				if (visitedV[instance.neigebor[cv][i]] == 'F' && he[instance.two_d[cv][instance.neigebor[cv][i]]] > 0)
				{
					unvisitedV.push(instance.neigebor[cv][i]);
				}
			}
		}
	}
	if (countv == instance.vw.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}
vector<int> find_ST(CDSI_instance& instance, vector<int>& he, vector<int>& heud)
{
	vector<int> STE, STV;//edges of spanning tree,vertices of spanning tree
	vector<char> visitedV(instance.vw.size(), 'F');//'v' for visited, 'n' for in neighbor,'F' for else
	vector<int> neighbor;//nieghbor of STV
	vector<int> degree(instance.vw.size(), 0);//degree of STV in STV

	// Find the edge with the maximum cost
	int me = max_element(heud.begin(), heud.end(), [&](int a, int b) {
		return instance.one_d[a].blocking_cost < instance.one_d[b].blocking_cost;
		}) - heud.begin();

		STE.push_back(heud[me]);//heud me is the edge with maximum cost
		int h = instance.one_d[heud[me]].head;
		int t = instance.one_d[heud[me]].tail;
		visitedV[h] = 'v';
		visitedV[t] = 'v';
		auto mark_visited = [&](int v) {
			STV.push_back(v);
			degree[v] = 1;
			for (int n : instance.neigebor[v]) {
				if (he[instance.two_d[v][n]] > 0 && visitedV[n] == 'F') {
					neighbor.push_back(n);
					visitedV[n] = 'n';
				}
			}
		};
		mark_visited(h);
		mark_visited(t);
		while (STV.size() != instance.vw.size())
		{
			// Find the vertex with the minimum weight in the neighbors
			auto min_it = min_element(neighbor.begin(), neighbor.end(), [&](int a, int b) {
				return instance.vw[a] < instance.vw[b];
				});
			int cv = *min_it;
			neighbor.erase(min_it);
			visitedV[cv] = 'v';

			for (int n : instance.neigebor[cv]) {
				if (visitedV[n] == 'F' && he[instance.two_d[cv][n]] > 0) {
					neighbor.push_back(n);
					visitedV[n] = 'n';
				}
			}

			// Find the vertex in STV with the highest degree that is connected to cv
			int ci = -1, max_degree = 0;
			for (int v : STV) {
				if (instance.two_d[cv][v] != -1 && he[instance.two_d[cv][v]] > 0) {
					if (degree[v] > max_degree) {
						ci = v;
						max_degree = degree[v];
					}
				}
			}
			if (ci != -1) {
				STV.push_back(cv);
				STE.push_back(instance.two_d[cv][ci]);
				degree[cv]++;
				degree[ci]++;
			}
		}
		return STE;
}
bool check_STW(CDSI_instance& instance, vector<int>& STE)
{//return true if the weight of spanning tree is less than eta
	double weight = 0;
	vector<int> degree;
	degree.resize(instance.vw.size(), 0);
	for (int i = 0; i < STE.size(); i++)
	{
		degree[instance.one_d[STE[i]].head]++;
		degree[instance.one_d[STE[i]].tail]++;
	}
	for (int i = 0; i < degree.size(); i++)
	{
		if (degree[i] > 1)
		{
			weight += instance.vw[i];
		}
	}
	if (weight < instance.desired_ST_weiget)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int insert_EST(CDSI_instance& instance, vector<EST>& estlist, vector<int>& ste)
{
	EST est1;
	est1.cost = 100;
	est1.me = 0;
	est1.EST1 = ste;
	est1.EST2.resize(instance.one_d.size(), 'F');
	for (int i = 0; i < ste.size(); i++)
	{
		est1.EST2[ste[i]] = 'T';
		if (instance.one_d[ste[i]].blocking_cost < est1.cost)
		{
			est1.cost = instance.one_d[ste[i]].blocking_cost;
			est1.me = ste[i];
		}
	}
	int i = estlist.size() - 1;
	estlist.resize(estlist.size() + 1);
	while (i >= 0)
	{
		if (estlist[i].cost < est1.cost)
		{
			estlist[i + 1] = estlist[i];
			i--;
		}
		else
		{
			break;
		}
	}
	estlist[i + 1] = est1;
	return i + 1;
}

double LB_unprocessed_nodes(CDSI_instance& instance, list<node>& tree)
{//cn+1 passed
	double LB = instance.one_d.size() * 10;
	for (auto it = tree.begin(); it != tree.end(); it++)
	{
		if ((*it).LB < LB)
		{
			LB = (*it).LB;
		}
	}
	return LB;
}
//ascending ratio
void d_sort_ratio(vector<pair<int, double>>& candidate)
{
	sort(candidate.begin(), candidate.end(), [](const pair<int, double>& f, const pair<int, double>& s) { return f.second > s.second; });
}
