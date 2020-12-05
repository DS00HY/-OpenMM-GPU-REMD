#include "mpiplus.h"
#include <stdarg.h>
//#define DEBUG
//#undef DEBUG
static int MyDebugPrintf(const char* format, ...)
{
	va_list argPtr;
	int     count;
	FILE*   fp;
	va_start(argPtr, format);                  /*  获取可变参数列表  */
	fflush(stdout);                            /*  强制刷新输出缓冲区  */
	fp = fopen("debug.log","a");
	count = vfprintf(fp, format, argPtr);  /*  将信息输出到标准输出流设备stdout  ，错误流为stderr*/
	fclose(fp);
	va_end(argPtr);                            /*  可变参数列表结束  */
}

#ifdef DEBUG   /*  如果定义了插桩信息宏，就将调试信息指向调试函数  */
#define DebugPrintf  MyDebugPrintf

#else           /*  如果未定义插桩信息宏，那么就将调试信息指向空NOP  */
#define DebugPrintf

#endif

void init_states(double*& states, int s_size)
{
	states = (double*)malloc(s_size * sizeof(double));
	if (states == NULL) {
		throw "states malloc error!!";
		return;
		// 如何报错
	}
}
void init_partners(std::vector<std::pair<int, int> >& partners, int comm_size, int nowitera, int ex_kind)
{
	/*
	  init the exchange partners, odd exchange or even exchange

	  * partners: save the exchange partners
	  * comm_size: size of partners
	  * isodd: odd exchange or even exchange
	  * ex_kind: exchange kind

	*/
	bool  isodd = nowitera % 2;
	partners.clear();
	//DebugPrintf("C++ pre partner %d :", isodd);
	if (ex_kind == 1) {
		if (isodd) {
			for (int i = 0; i < comm_size; i += 2) {
				partners.push_back(std::make_pair(i, i + 1));
				//DebugPrintf("(%d,%d) ", i, i + 1);
			}
			//DebugPrintf("end \n");
		}
		else {
			for (int i = 0; i < comm_size; i += 2) {
				if (i == 0) {
					partners.push_back(std::make_pair(i, comm_size - 1));
					//DebugPrintf("(%d,%d) ", i, comm_size - 1);
					continue;
				}
				partners.push_back(std::make_pair(i - 1, i));
				//DebugPrintf("(%d,%d) ", i - 1, i);
			}
			//DebugPrintf("end \n");
		}
	}
	else if (ex_kind == 2) {
		throw "Exchange kind not exist!";
		return;
	}
	else {
		throw "Exchange kind not exist!";
		return;
	}
}


long getExchange(MPI_Comm comm, int nowitera, double* replica_state, int state_size, long* paramlist, int param_size, int ex_kind, int md_kind) {

	int  comm_rank;
	long idx;
	MPI_Comm_rank(comm, &comm_rank);

	if (comm_rank == 0) {
		idx = master(comm, nowitera, replica_state, state_size, paramlist, ex_kind, md_kind);

	}
	else {
		idx = replica(comm, replica_state, state_size);
	}
	return idx;
}

long master(MPI_Comm comm, int nowitera, double* replica_state, int state_size, long*& paramlist, int ex_kind, int md_kind)
{
	/*
	  master in python call this method

	  * n_iterations:
	  * state_size： exchangesize
	  * ex_kind: exchange with who
	  * md_kind: md/mc ? md+mc
	*/
	double* states;
	int s_size, comm_size, comm_rank, n_replica;
	long idx;
	std::vector<std::pair<int, int> > partners;
	std::vector<std::pair<int, int> > ex_partner_id;
	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &comm_rank);

	/*init parameter list*/
	DebugPrintf("DEBUG[c++]: master get size=%d, function is %d, %d, %d, %d\n", comm_size, nowitera, state_size, ex_kind, md_kind);

	s_size = state_size * comm_size;
	init_states(states, s_size);
	/*init partners*/
	init_partners(partners, comm_size, nowitera, ex_kind);

	/*gather states;*/
	MPI_Gather(replica_state, state_size, MPI_DOUBLE, states, state_size, MPI_DOUBLE, 0, comm);

	DebugPrintf("DEBUG[c++]:(%d) preparamlist = ", nowitera);
	for (int i = 0; i < comm_size; i++) {
		DebugPrintf("%d ", paramlist[i]);
	}
	DebugPrintf("\nDEBUG[c++]:(%d) states = ", nowitera);
	for (int i = 0; i < s_size; i++) {
		DebugPrintf("%f ", states[i]);		
	}
	DebugPrintf("\n");
	/*判断交换*/
	acceptance_exchange(states, state_size, partners, ex_partner_id, md_kind);
	DebugPrintf("DEBUG[c++]:(%d) ", nowitera);
	for (int i = 0; i < ex_partner_id.size(); i++) {
		DebugPrintf("ex(%d,%d) ", ex_partner_id[i].first, ex_partner_id[i].second);
	}
	DebugPrintf("\nDEBUG[c++]:(%d) exsize is %d now is exchange\n", nowitera, ex_partner_id.size());


	/*交换*/
	exchange(paramlist, ex_partner_id);

	DebugPrintf("DEBUG[c++]:(%d) now newparamlist = ", nowitera);
	for (int i = 0; i < comm_size; i++) {
		DebugPrintf("%d ", paramlist[i]);
	}
	DebugPrintf("\n");

	//分配 paramlist idx
	MPI_Scatter(paramlist, 1, MPI_LONG, &idx, 1, MPI_LONG, 0, comm);

	writeinfo(nowitera, partners);

	free(states);
	partners.clear();
	ex_partner_id.clear();

	return idx;
}

long replica(MPI_Comm comm, double* replica_state, int state_size)
{
	/*
	   replica send states to master （Python call)
	*/
	int* paramlist = NULL;
	long idx;
	double* states = NULL;
	int  comm_rank;
	MPI_Comm_rank(comm, &comm_rank);
	//printf("C++ now task %d gather\n",comm_rank);
	MPI_Gather(replica_state, state_size, MPI_DOUBLE, states, state_size, MPI_DOUBLE, 0, comm);

	//
	MPI_Scatter(paramlist, 1, MPI_LONG, &idx, 1, MPI_LONG, 0, comm);
	//printf("C++ now task %d get param = %ld\n", comm_rank, idx);
	return idx;
}




void exchange(long*& paramlist, std::vector<std::pair<int, int> > ex_partner_id)
{
	/*
	  exchange the paramlist id between exchange partners

	  * paramlist:
	  * ex_partner_id:
	*/
	for (int i = 0; i < ex_partner_id.size(); i++) {
		int a = ex_partner_id[i].first, b = ex_partner_id[i].second;
		int tmp;
		tmp = paramlist[a];
		paramlist[a] = paramlist[b];
		paramlist[b] = tmp;
	}
}


void acceptance_exchange(double* states, int state_size, std::vector<std::pair<int, int> > partners, std::vector<std::pair<int, int> >& ex_partner_id, int md_kind)
{
	/*
	  get replica pairs that need to exchange ，return the exchange pair by "ex_partner_id"

	  * states: all replicas states
	  * state_size:
	  * partners: the id pair of exchange replica partners
	  * ex_partner_id: the relica need to exchange
	  * md_kind:

	*/
	ex_partner_id.clear();
	srand((time(NULL)));
	//DebugPrintf("DEBUG: now in accc  size is %d , md_kind is \n", partners.size(), md_kind);
	for (int i = 0; i < partners.size(); i++) {
		int a = partners[i].first, b = partners[i].second;
		if (md_kind == 1) {
			DebugPrintf("DEBUG[c++]: now judge (%d,%d): \n",a , b);
			bool acc = get_acceptance(states + (a * state_size), states + (b * state_size)); //  // 
			//printf("%d with %d is %d\n", a, b, acc);
			//acc = true;
			if (acc) {
				ex_partner_id.push_back(std::make_pair(a, b));
			}
		}
		else {
			throw "The kind of md is error!\n";
		}
	}

}

bool get_acceptance(double* state_a, double* state_b)
{
	/*
	  calculate the exchange probability ,if ok return true, else returen false.

	  * state_a: state of replica A
	  * state_b:

	  return whether or not to exchange

	*/
	double ua = state_a[0], ub = state_b[0], ta = state_a[1], tb = state_b[1], r = 8.3138462;
	//double delta = ua / (r * tb) + ub / (r * ta) - ua / (r * ta) - ub / (r * tb);
	//delta = -delta;
	//delta *= 1000;
	double delta = ((ua - ub) * (1 / (r * ta) - 1 / (r * tb))) * 1000;
	double randx = rand() / double(RAND_MAX);
	DebugPrintf("            ua = %f,ta= %f,ub= %f, tb= %f\n            randx = %f < exp(%f) = %f ?\n",  ua, ta, ub, tb, randx, delta, exp(delta));
	if (delta >= 0.0 || randx < exp(delta))
		return true;
	else
		return false;
}



void writeinfo(int nowitera, std::vector<std::pair<int, int> >& partners)
{
	/*

	*/
	if (nowitera == 0) {
		std::ofstream outfile("ex_pair.out", std::ios::out);
		if (!outfile.is_open())
		{
			throw "Outfile can't open!";
		}
		outfile << "now_iteration";
		for (int i = 0; i < partners.size(); i++) {
			outfile << "  replica" << i << " with " << "replica" << i + 1;
		}
		outfile << std::endl;
		outfile.close();
	}
	std::ofstream outfile("ex_pair.out", std::ios::app);
	if (!outfile.is_open())
	{
		throw "Outfile can't open!";
	}
	//DebugPrintf("now write\n");
	outfile << nowitera << ":";
	for (int i = 0; i < partners.size(); i++) {
		outfile << "  " << partners[i].first << " " << partners[i].second;
	}
	outfile << std::endl;
	outfile.close();

}
