//add a lazy cut of an EST each time
//use the a,b seperator based IP formulation to detect EST
//use min cut at the beginning
#include <vector>
using namespace std;
#include <list>
//To be able to use list
//To be able to use vector
#include <queue>
#include "gurobi_c++.h"
//To be able to use Gurobi
#include <fstream>
//To be able to use fin
#include <sstream>
//to be able to read the elements of an string (istringstream comamnd)
#include <time.h>

#include <algorithm>//sort
#include <regex>

using namespace std;

class vertex
{
public:
	double vertex_weight;

	vertex()
	{
		vertex_weight = 0;
	}
};

class edge
{
public:
	int head;
	int tail;
	double blocking_costs;
	edge(int hd, int tl, double bc)
	{
		head = hd;
		tail = tl;
		blocking_costs = bc;
	}

	edge()
	{
		head = 0;
		tail = 0;
		blocking_costs = 0;
	}
};

class CDSI_instance
{
public:
	vector <vertex> vertices;
	vector <vector<int>> two_d;//store edge id correspond to one_d
	//==-1 means no edge
	vector <vector<int>> neighbor;//store neighbors of every vertices
	vector <edge> one_d;
	double desired_clique_weight;
	int total_edge;//total possible edge size,n*(n-1)/2
	double total_weight;//total weight of all vertices
	CDSI_instance()
	{
		vertices = vector<vertex>();
		two_d = vector <vector<int>>();
		neighbor = vector <vector<int>>();
		one_d = vector<edge>();
		desired_clique_weight = -1;
		total_weight = 0;
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

class lazy_status
{
public:
	bool not_TL;
	double solution;
	double bound;
	double total_time_limit;
	clock_t start;
	int number_of_lazy;
	int gurobi_lazy;
	bool take_start;
	int aos = 0;//alternative optimal solution other than min-cut
	double root_ub;
	double root_lb;

	lazy_status(double tt, clock_t st)
	{
		not_TL = true;
		solution = 0;
		bound = 0;
		total_time_limit = tt;
		start = st;
		number_of_lazy = 0;
		take_start = true;
		aos = 0;
	}
	lazy_status()
	{
		not_TL = true;
		solution = 0;
		bound = 0;
		total_time_limit = 0;
		start = clock();
		number_of_lazy = 0;
		take_start = true;
		aos = 0;
	}
};

bool get_EST(vector<int>& remain_edges, vector<int>& eta_EST, CDSI_instance& instance, char* eta_ST, clock_t start, clock_t total_time_limit,string);
int get_ST(CDSI_instance& instance, vector<vector<int>>& neighbor, int D, vector<char>& D2, vector<int>& STE);
bool check_full_cycle(CDSI_instance& instance, vector<int>& remain_edges, vector<char>& D2, int D_size, int headD, int tailD, vector<int> EST);
char check_arc_cycle(CDSI_instance& instance, vector<int>& remain_edges, vector<char>& D2, int D_size, int headD, int tailD, vector<int> EST);

class mycallback : public GRBCallback
{
public:
	int numx;//number of variable x
	GRBVar* varsx;
	CDSI_instance instance;
	int* tree;
	int* arc;
	int* cycle;
	lazy_status* lazy_s;
	double minS;
	vector<int> min_cut;
	string filename;
	mycallback(int xnumvars, GRBVar* xvarsx, CDSI_instance& xinstance, int* t, int* a, int* c, lazy_status* ls, double mS, vector<int>& mcut,string fn)
	{
		numx = xnumvars;
		varsx = xvarsx;
		instance = xinstance;
		lazy_s = ls;
		tree = t;
		arc = a;
		cycle = c;
		minS = mS;
		min_cut = mcut;
		filename = fn;
	}
protected:
	void callback() {
		try {
			if (where == GRB_CB_MIPSOL) {
				// MIP solution callback
				if ((*lazy_s).not_TL) {
					(*lazy_s).bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
					(*lazy_s).solution = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
				}
				double* x = getSolution(varsx, numx);
				vector<int> remain_edges;
				int edge_size = 0;
				remain_edges.resize(instance.one_d.size(), 0);
				for (int i = 0; i < numx; i++)
				{
					if (x[i] < 0.5)
					{//not block
						remain_edges[i] = 1;
						edge_size++;
					}
				}
				clock_t now = clock();

				vector<int> eta_est;
				char eta_ST = 0;
				bool reach_time_limit = get_EST(remain_edges, eta_est, instance, &eta_ST, (*lazy_s).start, (*lazy_s).total_time_limit,filename);
				if (eta_est.size() > 0)
				{
					GRBLinExpr lhs = NULL;
					for (int i = 0; i < eta_est.size(); i++)
					{
						lhs += varsx[eta_est[i]];
					}
					if (eta_ST == 'C')
					{
						(*cycle)++;
						addLazy(lhs >= 2);
					}
					else
					{
						addLazy(lhs >= 1);
						if (eta_ST == 'A')
						{
							(*arc)++;
						}
						else
						{
							(*tree)++;
						}
					}
					((*lazy_s).number_of_lazy)++;
				}
				else if (reach_time_limit)
				{
					(*lazy_s).bound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
					(*lazy_s).solution = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
					(*lazy_s).not_TL = false;
					ofstream fout;
					fout.open(filename, std::ios_base::app);
					fout << "Manual abort called now" << endl;
					fout.close();
					abort();
				}
				else
				{
					(*lazy_s).solution = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
					int temps = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
					if (temps == static_cast<int>(minS))
					{
						bool equal=true;
						vector<int> temp_set;
						for (int i = 0; i < numx; i++)
						{
							if (x[i] > 0.5)
							{//not block
								temp_set.push_back(i);
							}
						}
						if (temp_set.size() != min_cut.size())
						{
							equal = false;
						}
						else
						{
							sort(temp_set.begin(), temp_set.end());
							for (int i = 0; i < min_cut.size(); i++)
							{
								if (min_cut[i] != temp_set[i])
								{
									equal = false;
									break;
								}
							}
						}
						if (!equal)
						{
							(*lazy_s).aos++;
						}
					}
					else if (temps < static_cast<int>(minS))
					{
						(*lazy_s).aos++;
					}
				}
				delete[] x;
			}
			else if (where == GRB_CB_MIP) {
				// MIP callback
				if ((*lazy_s).not_TL) {
					(*lazy_s).solution = getDoubleInfo(GRB_CB_MIP_OBJBST);
					(*lazy_s).bound = getDoubleInfo(GRB_CB_MIP_OBJBND);
					if (getDoubleInfo(GRB_CB_MIP_NODCNT) == 0) {
						(*lazy_s).root_lb = getDoubleInfo(GRB_CB_MIP_OBJBND);
						(*lazy_s).root_ub = getDoubleInfo(GRB_CB_MIP_OBJBST);
					}
				}
			}
			else if (where == GRB_CB_MIPNODE) {
			// MIP node callback
				if ((*lazy_s).not_TL) {
					(*lazy_s).solution = getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
					(*lazy_s).bound = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
				}
			}
			else if (where == GRB_CB_MESSAGE) {
				// Message callback
				string msg = getStringInfo(GRB_CB_MSG_STRING);
				if (msg == "User MIP start did not produce a new incumbent solution\n")
				{
					(*lazy_s).take_start = false;
				}
				regex pattern(R"(Lazy constraints:\s*(\d+))");
				smatch match;
				if (regex_search(msg, match, pattern)) {
					(*lazy_s).gurobi_lazy = stoi(match[1].str());
				}
				ofstream fout;
				fout.open(filename, std::ios_base::app);
				fout << msg << endl;
				fout.close();

			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};

int readdata(CDSI_instance& instance, string filename);

void runmodel(string);

double MINIMUMCUT(CDSI_instance& instance, vector<int>& edge_MC);
void MINIMUMCUTPHASE(CDSI_instance& instance, graph& G, vector<int>& edge_MC, double& MC, int a);

int main()
{
	ofstream fout;
	fout.open("CDSI_lazy.csv");
	fout << "File Name,|V|,density,Eta,number_of_lazy,gurobi_#lazy,#tree,#arc,#cycle,Time Limit (s),Run Time (s),Gurobi gap,Gurobi ub,Gurobi lb,manual gap, manual ub, manual lb,root gap,root ub,root lb,take_start, Nodes #,min_cut,mincut not optimal, alternative solution" << endl;
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
	fin.close();
	fin.open("filename-random.txt");
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
	ofstream fout;
	string fn = filename + "_log.txt";
	fout.open(fn);
	fout.close();
	double total_time_limit = 3600.0;
	CDSI_instance instance;
	int number_of_edges = readdata(instance, filename);

	clock_t pre_process_start = clock();

	//Defining Gurobi environment
	GRBEnv* env = 0;
	env = new GRBEnv();

	// Model
	GRBModel model = GRBModel(*env);
	model.set(GRB_StringAttr_ModelName, filename);

	// Must set LazyConstraints parameter when using lazy constraints
	model.set(GRB_IntParam_LazyConstraints, 1);
	// Decision variables
	GRBVar* x = 0;//edges
	x = model.addVars(instance.one_d.size(), GRB_BINARY);
	for (int i = 0; i < instance.one_d.size(); i++)
	{
		x[i].set(GRB_DoubleAttr_Obj, instance.one_d[i].blocking_costs);
	}
	// The objective is minimization
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	model.update();

	model.set(GRB_IntParam_Threads, 1);
	model.update();

	//Set initial solution
	vector<int> UB_solution;
	double incumbent_cost = MINIMUMCUT(instance, UB_solution);
	sort(UB_solution.begin(), UB_solution.end());
	int ubi = 0;
	for (int i = 0; i < instance.one_d.size(); i++)
	{
		if (ubi < UB_solution.size() && i == UB_solution[ubi])
		{
			x[i].set(GRB_DoubleAttr_Start, 1);
			ubi++;
		}
		else
		{
			x[i].set(GRB_DoubleAttr_Start, 0);
		}
	}

	// Create a callback object and associate it with the model
	int ntree = 0, narc = 0, ncycle = 0;
	lazy_status ls = lazy_status(total_time_limit, pre_process_start);
	mycallback cb = mycallback(instance.one_d.size(), x, instance, &ntree, &narc, &ncycle, &ls, incumbent_cost, UB_solution,fn);

	model.setCallback(&cb);
	clock_t now = clock();
	double limit = total_time_limit - static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC;
	model.set(GRB_DoubleParam_TimeLimit, limit);
	model.optimize();
	clock_t end_time = clock();

	fout.open("CDSI_lazy.csv", std::ios_base::app);
	fout << filename << ',' << instance.vertices.size();
	double density = number_of_edges / (static_cast<double>(instance.vertices.size() * (instance.vertices.size() - 1)) / 2);
	fout << ',' << density << ',' << instance.desired_clique_weight;
	fout << ',' << ls.number_of_lazy;
	fout<<','<<ls.gurobi_lazy << "," << ntree << "," << narc << "," << ncycle;
	fout << ',' << total_time_limit << ",";

	if (model.get(GRB_IntAttr_Status) == GRB_INTERRUPTED)
	{//time limit reached in call back
		fout << "LTL" << ',';
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		fout << static_cast<double>(end_time - pre_process_start) / CLOCKS_PER_SEC << ",";
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
	{
		fout << "TL" << ',';
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
	{
		fout << "IN" << ',';
	}
	else
	{
		fout << "ELSE" << ',';
	}
	fout << model.get(GRB_DoubleAttr_MIPGap) << ',';
	fout << model.get(GRB_DoubleAttr_ObjVal) << ',';
	fout << model.get(GRB_DoubleAttr_ObjBoundC) << ',';
	fout << (ls.solution - ls.bound) / ls.solution << ',';
	fout << ls.solution << ',';
	fout << ls.bound << ',';
	fout << (ls.root_ub - ls.root_lb) / ls.root_ub << ',';
	fout << ls.root_ub << ',';
	fout << ls.root_lb << ',';
	fout << ls.take_start << ',';
	fout << model.get(GRB_DoubleAttr_NodeCount) << ',';
	fout << incumbent_cost << ',';
	if (ls.solution < incumbent_cost)
	{
		fout << "Yes" << ',';
	}
	else
	{
		fout << "No" << ',';
	}
	fout << ls.aos << ',';
	fout << endl;
	fout.close();
	delete[] x;
	delete env;

}

double MINIMUMCUT(CDSI_instance& instance, vector<int>& edge_MC)
{
	graph G;
	G.vset.resize(instance.vertices.size());
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
				edge_weiget[i][j] = instance.one_d[instance.two_d[i][j]].blocking_costs;
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

bool get_EST(vector<int>& remain_edges, vector<int>& eta_EST, CDSI_instance& instance, char* eta_ST, clock_t start, clock_t total_time_limit,string fn)
{//return true if reach time limit, else return false, eta_ST is 'A' for arc,'C' for cycle,'T' for tree
	//Defining Gurobi environment
	GRBEnv* sub_env = new GRBEnv();

	// Model
	GRBModel sub_model = GRBModel(*sub_env);

	// Decision variables
	GRBVar* y = NULL;
	y = sub_model.addVars(instance.vertices.size(), GRB_BINARY);
	GRBVar**** z = NULL;
	int ab_size = instance.total_edge - instance.one_d.size();
	z = new GRBVar ***[instance.vertices.size()];
	vector<vector<int>> neighbor;
	vector<vector<int>> not_adjacent;
	neighbor.resize(instance.vertices.size());
	not_adjacent.resize(instance.vertices.size() - 1);
	for (int i = 0; i < instance.vertices.size(); i++)
	{
		for (int j = i + 1; j < instance.vertices.size(); j++)
		{
			if (instance.two_d[i][j] == -1)
			{
				not_adjacent[i].push_back(j);
			}
			else if (remain_edges[instance.two_d[i][j]] == 0)
			{
				not_adjacent[i].push_back(j);
			}
			else
			{
				neighbor[i].push_back(j);
				neighbor[j].push_back(i);
			}
		}
	}
	for (int a = 0; a < not_adjacent.size(); a++)
	{
		z[a] = new GRBVar * *[not_adjacent[a].size()];
		for (int b = 0; b < not_adjacent[a].size(); b++)
		{
			z[a][b] = new GRBVar * [instance.vertices.size()];
			for (int i = 0; i < instance.vertices.size(); i++)
			{
				z[a][b][i] = sub_model.addVars(neighbor[i].size());
				for (int j = 0; j < neighbor[i].size(); j++)
				{
					z[a][b][i][j].set(GRB_DoubleAttr_LB, 0);
					z[a][b][i][j].set(GRB_DoubleAttr_UB, 1);
				}
			}
		}
	}
	clock_t now = clock();
	double passing_time = static_cast<double>(now - start) / CLOCKS_PER_SEC;
	double optimization_time_limit = total_time_limit - passing_time;
	if (optimization_time_limit <= 0)
	{//reach time limit
		ofstream fout;
		fout.open(fn, std::ios_base::app);
		fout << "Time limit in callback model building." << endl;
		return true;
	}
	
	GRBLinExpr lh3a = 0, lh3c = 0, lh3d = 0, lh3e = 0;
	GRBLinExpr lh3b = 0;
	for (int i = 0; i < instance.vertices.size(); i++)
	{
		lh3a = y[i];
		for (int j = 0; j < neighbor[i].size(); j++)
		{
			lh3a += y[neighbor[i][j]];
		}
		sub_model.addConstr(lh3a, GRB_GREATER_EQUAL, 1);
		lh3b += y[i] * instance.vertices[i].vertex_weight;
	}
	sub_model.addConstr(lh3b, GRB_LESS_EQUAL, instance.desired_clique_weight - 0.00001);
	for (int a = 0; a < not_adjacent.size(); a++)
	{
		for (int b = 0; b < not_adjacent[a].size(); b++)
		{
			int nab = not_adjacent[a][b];
			for (int i = 0; i < instance.vertices.size(); i++)
			{
				lh3c = -y[i];
				for (int j = 0; j < neighbor[i].size(); j++)
				{
					lh3c += z[a][b][i][j];
				}
				sub_model.addConstr(lh3c, GRB_LESS_EQUAL, 0);
				if (i != a && i != nab)
				{
					lh3e = 0;
					for (int j = 0; j < neighbor[i].size(); j++)
					{
						lh3e += z[a][b][i][j];
						int nij = neighbor[i][j];
						int niji;
						for (int k = 0; k < neighbor[nij].size(); k++)
						{
							if (neighbor[nij][k] == i)
							{
								niji = k;
								break;
							}
						}
						lh3e -= z[a][b][nij][niji];
					}
					sub_model.addConstr(lh3e, GRB_EQUAL, 0);
				}
			}
			lh3d = y[a] + y[not_adjacent[a][b]];
			for (int i = 0; i < neighbor[a].size(); i++)
			{
				int nai = neighbor[a][i];
				int naia;
				for (int j = 0; j < neighbor[nai].size(); j++)
				{
					if (neighbor[nai][j] == a)
					{
						naia = j;
						break;
					}
				}
				lh3d -= z[a][b][a][i] - z[a][b][nai][naia];
			}
			sub_model.addConstr(lh3d, GRB_LESS_EQUAL, 1);
		}

	}
	sub_model.set(GRB_IntParam_Threads, 1);
	clock_t build_end = clock();

	passing_time = static_cast<double>(build_end - start) / CLOCKS_PER_SEC;
	optimization_time_limit = total_time_limit - passing_time;
	if (optimization_time_limit <= 0)
	{//reach time limit
		return true;
	}
	sub_model.set(GRB_DoubleParam_TimeLimit, optimization_time_limit);
	sub_model.update();

	sub_model.optimize();
	ofstream fout;
	fout.open(fn, std::ios_base::app);
	fout << "Lazy_status:" << sub_model.get(GRB_IntAttr_Status) << endl;
	if (sub_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		vector<char> yisD;
		yisD.resize(instance.vertices.size(), 'F');
		int headD;//store a vertex in D as the root of spanning tree
		int countD = 0;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			double temp_solution_y = y[i].get(GRB_DoubleAttr_X);
			if (temp_solution_y > 0.5)
			{
				yisD[i] = 'T';
				headD = i;
				countD++;
			}
		}
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (yisD[i] == 'F')
			{//not in D, add arbitrary edge to connect to D
				for (int j = 0; j < neighbor[i].size(); j++)
				{
					int jj = neighbor[i][j];
					if (yisD[jj] == 'T')
					{
						eta_EST.push_back(instance.two_d[i][jj]);
						break;
					}
				}
			}
		}
		//find a spanning tree in D
		int tailD = get_ST(instance, neighbor, headD, yisD, eta_EST);
		if ((tailD == -1) || instance.vertices.size() - 2 > countD)
		{
			(*eta_ST) = 'T';
		}
		else
		{
			bool is_arc = check_full_cycle(instance, remain_edges, yisD, countD, headD, tailD, eta_EST);
			if (!is_arc)
			{
				*eta_ST = 'T';
			}
			else
			{
				*eta_ST = check_arc_cycle(instance, remain_edges, yisD, countD, headD, tailD, eta_EST);
			}
		}
		delete[] y;
		delete[] z;
		delete sub_env;
		return false;
	}
	else if (sub_model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
	{
		delete[] y;
		delete[] z;
		delete sub_env;
		return true;
	}
	else
	{
		delete[] y;
		delete[] z;
		delete sub_env;
		return false;
	}
}

int get_ST(CDSI_instance& instance, vector<vector<int>>& neighbor, int D, vector<char>& D2, vector<int>& STE)
{//get a spanning tree on D with breadth first search strategy
 //D is an abitrary vertex in D used as the root of spanning tree
 //if D is a path, return tail of D, else return -1
	vector<char> accessD;
	accessD.resize(instance.vertices.size(), 'F');
	queue<int> searchQ;//the queue for searching
	accessD[D] = 'T';
	searchQ.push(D);
	bool ispath = true;
	int tail = D;//record the last accessed vertex
	while (!searchQ.empty())
	{
		int tV = searchQ.front();
		searchQ.pop();
		int addcount = 0;
		for (int i = 0; i < neighbor[tV].size(); i++)
		{
			int nti = neighbor[tV][i];
			if (D2[nti] == 'T' && accessD[nti] == 'F')
			{
				searchQ.push(nti);
				tail = nti;
				accessD[nti] = 'T';
				STE.push_back(instance.two_d[tV][nti]);
				addcount++;
				if (addcount > 1)
				{
					ispath = false;
				}
			}
		}
	}
	if (ispath)
	{
		return tail;
	}
	else
	{
		return -1;
	}
}

bool check_full_cycle(CDSI_instance& instance, vector<int>& remain_edges, vector<char>& D2, int D_size, int headD, int tailD, vector<int> EST)
{//return true if the est can make full cycle,D is assumed as a path
	if (D_size == instance.vertices.size())
	{
		if (remain_edges[instance.two_d[headD][tailD]] == 1)
		{
			return true;
		}
		else
			return false;
	}
	else if (D_size == instance.vertices.size() - 1)
	{
		int VnotinD;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (D2[i] == 'F')
			{
				VnotinD = i;
				break;
			}
		}
		int Vt = instance.two_d[VnotinD][tailD];
		int Vh = instance.two_d[VnotinD][headD];
		if (Vt == -1 || Vh == -1)
		{
			return false;
		}
		if (remain_edges[Vt] == 0 || remain_edges[Vh] == 0)
		{
			return false;
		}
		for (int i = 0; i < EST.size(); i++)
		{
			if (EST[i] == Vt || EST[i] == Vh)
			{
				return true;
			}
		}
		return false;
	}
	else
	{//|D|=|V|-2
		int V1 = -1, V2 = -1;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (D2[i] == 'F')
			{
				if (V1 == -1)
				{
					V1 = i;
				}
				else
				{
					V2 = i;
					break;
				}
			}
		}
		if (instance.two_d[V1][V2] == -1 || remain_edges[instance.two_d[V1][V2]] == 0)
		{
			return false;
		}
		int V1t = instance.two_d[V1][tailD];
		int V1h = instance.two_d[V1][headD];
		int V2t = instance.two_d[V2][tailD];
		int V2h = instance.two_d[V2][headD];
		if ((V1t == -1 || V2h == -1) && (V2t == 1 || V1h == -1))
		{
			return false;
		}
		if ((remain_edges[V1t] == 0 || remain_edges[V2h] == 0) && (remain_edges[V2t] == 0 || remain_edges[V1h] == 0))
		{
			return false;
		}
		for (int i = 0; i < EST.size(); i++)
		{
			if (EST[i] == V1t)
			{
				for (int j = i + 1; j < EST.size(); j++)
				{
					if (EST[j] == V2h)
					{
						return true;
					}
				}
				return false;
			}
			else if (EST[i] == V1h)
			{
				for (int j = i + 1; j < EST.size(); j++)
				{
					if (EST[j] == V2t)
					{
						return true;
					}
				}
				return false;
			}
			else if (EST[i] == V2h)
			{
				for (int j = i + 1; j < EST.size(); j++)
				{
					if (EST[j] == V1t)
					{
						return true;
					}
				}
				return false;
			}
			else if (EST[i] == V2t)
			{
				for (int j = i + 1; j < EST.size(); j++)
				{
					if (EST[j] == V1h)
					{
						return true;
					}
				}
				return false;
			}
		}
		return false;
	}
}

char check_arc_cycle(CDSI_instance& instance, vector<int>& remain_edges, vector<char>& D2, int D_size, int headD, int tailD, vector<int> EST)
{//return 'A' if the est is an arc, 'C' if the est is a cycle
 //we add the edge to close the cycle if it is a cycle
	if (D_size == instance.vertices.size())
	{
		for (int i = 0; i < EST.size(); i++)
		{
			int h = instance.one_d[EST[i]].head;
			int t = instance.one_d[EST[i]].tail;
			if (instance.total_weight - instance.vertices[h].vertex_weight - instance.vertices[t].vertex_weight >= instance.desired_clique_weight)
			{
				return 'A';
			}
		}
		if (instance.total_weight - instance.vertices[headD].vertex_weight - instance.vertices[tailD].vertex_weight >= instance.desired_clique_weight)
		{
			return 'A';
		}
		else
		{
			EST.push_back(instance.two_d[headD][tailD]);
			return 'C';
		}
	}
	else if (D_size == instance.vertices.size() - 1)
	{
		int VnotinD;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (D2[i] == 'F')
			{
				VnotinD = i;
				break;
			}
		}
		for (int i = 0; i < EST.size(); i++)
		{
			int h = instance.one_d[EST[i]].head;
			int t = instance.one_d[EST[i]].tail;
			if (instance.total_weight - instance.vertices[h].vertex_weight - instance.vertices[t].vertex_weight >= instance.desired_clique_weight)
			{
				return 'A';
			}
			if (h == VnotinD)
			{
				if (t = headD)
				{
					if (instance.total_weight - instance.vertices[tailD].vertex_weight - instance.vertices[VnotinD].vertex_weight >= instance.desired_clique_weight)
					{
						return 'A';
					}
					else
					{
						EST.push_back(instance.two_d[tailD][VnotinD]);
						return 'C';
					}
				}
				else
				{
					if (instance.total_weight - instance.vertices[headD].vertex_weight - instance.vertices[VnotinD].vertex_weight >= instance.desired_clique_weight)
					{
						return 'A';
					}
					else
					{
						EST.push_back(instance.two_d[headD][VnotinD]);
						return 'C';
					}
				}
			}
			else if (t == VnotinD)
			{
				if (h = headD)
				{
					if (instance.total_weight - instance.vertices[tailD].vertex_weight - instance.vertices[VnotinD].vertex_weight >= instance.desired_clique_weight)
					{
						return 'A';
					}
					else
					{
						EST.push_back(instance.two_d[tailD][VnotinD]);
						return 'C';
					}
				}
				else
				{
					if (instance.total_weight - instance.vertices[headD].vertex_weight - instance.vertices[VnotinD].vertex_weight >= instance.desired_clique_weight)
					{
						return 'A';
					}
					else
					{
						EST.push_back(instance.two_d[headD][VnotinD]);
						return 'C';
					}
				}
			}
		}
	}
	else
	{//|D|=|V|-2
		int V1 = -1, V2 = -1;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (D2[i] == 'F')
			{
				if (V1 == -1)
				{
					V1 = i;
				}
				else
				{
					V2 = i;
					break;
				}
			}
		}
		if (instance.total_weight - instance.vertices[V1].vertex_weight - instance.vertices[V2].vertex_weight >= instance.desired_clique_weight)
		{
			return 'A';
		}
		int V1t = instance.two_d[V1][tailD];
		int V1h = instance.two_d[V1][headD];
		int V2t = instance.two_d[V2][tailD];
		int V2h = instance.two_d[V2][headD];
		for (int i = 0; i < EST.size(); i++)
		{
			int h = instance.one_d[EST[i]].head;
			int t = instance.one_d[EST[i]].tail;
			if (instance.total_weight - instance.vertices[h].vertex_weight - instance.vertices[t].vertex_weight >= instance.desired_clique_weight)
			{
				return 'A';
			}
			if (EST[i] == V1t)
			{
				if (instance.total_weight - instance.vertices[V2].vertex_weight - instance.vertices[headD].vertex_weight >= instance.desired_clique_weight)
				{
					return 'A';
				}
				else
				{
					EST.push_back(instance.two_d[headD][V2]);
					return 'C';
				}
			}
			else if (EST[i] == V1h)
			{
				if (instance.total_weight - instance.vertices[V2].vertex_weight - instance.vertices[tailD].vertex_weight >= instance.desired_clique_weight)
				{
					return 'A';
				}
				else
				{
					EST.push_back(instance.two_d[tailD][V2]);
					return 'C';
				}
			}
			else if (EST[i] == V2h)
			{
				if (instance.total_weight - instance.vertices[V1].vertex_weight - instance.vertices[tailD].vertex_weight >= instance.desired_clique_weight)
				{
					return 'A';
				}
				else
				{
					EST.push_back(instance.two_d[tailD][V1]);
					return 'C';
				}
			}
			else if (EST[i] == V2t)
			{
				if (instance.total_weight - instance.vertices[V1].vertex_weight - instance.vertices[headD].vertex_weight >= instance.desired_clique_weight)
				{
					return 'A';
				}
				else
				{
					EST.push_back(instance.two_d[headD][V1]);
					return 'C';
				}
			}
		}
	}
}

int readdata(CDSI_instance& instance, string filename)
{
	ifstream fin;
	fin.open(filename);
	string str3;
	int number_of_points, number_of_edges;
	fin >> number_of_points;
	fin >> number_of_edges;
	fin >> instance.desired_clique_weight;
	instance.two_d.resize(number_of_points);
	instance.vertices.resize(number_of_points);
	instance.neighbor.resize(number_of_points);
	for (int i = 0; i < number_of_points; i++)
	{
		instance.two_d[i].resize(number_of_points, -1);
	}
	getline(fin, str3);//finish current line
	for (int i = 0; i < number_of_points; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		ss0 >> instance.vertices[i].vertex_weight;
		instance.total_weight += instance.vertices[i].vertex_weight;
	}
	int ondn = 0;
	for (int i = 0; i < number_of_edges; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		int h, t;
		ss0 >> h;
		ss0 >> t;
		h--;
		t--;
		if (instance.two_d[h][t] == -1)
		{
			double temp_cost;
			ss0 >> temp_cost;
			instance.one_d.push_back(edge(h, t, temp_cost));
			instance.two_d[h][t] = ondn++;
			instance.two_d[t][h] = instance.two_d[h][t];
			instance.neighbor[h].push_back(t);
			instance.neighbor[t].push_back(h);
		}
	}
	instance.total_edge = instance.vertices.size() * (instance.vertices.size() - 1) / 2;
	return instance.one_d.size();
}
