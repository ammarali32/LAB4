#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include<iostream>
#include<chrono>
#include<fstream>
using namespace std;
#define INFINITY 1000000000


int p;
int* sendcounts;
int* sendcounts2;
int* displs;
int* displs2;
bool fix = 0;
int originalN;
void Gen_matrix(int LocAdjMatrix[], int n, int nPartition, int my_rank, MPI_Comm comm) {
	int* mat = NULL, i, j;
	if (my_rank == 0) {
		mat = (int*)malloc(n * n * sizeof(int));
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if ((fix && i < originalN && j<originalN) || !fix)
					mat[i * n + j] = rand() % 10000 + 1;
				else {
					mat[i * n + j] = INFINITY;
				}
	}
	MPI_Scatterv(mat, sendcounts, displs,MPI_INT, LocAdjMatrix, sendcounts[my_rank], MPI_INT, 0, comm);
	//MPI_Scatter(mat, 1, MPIblock, LocAdjMatrix, n * nPartition, MPI_INT, 0, comm);
	if (my_rank == 0) free(mat);
}



void Init(int LocAdjMatrix[], int LocParent[], int LocCost[], int locVisited[],
	int my_rank, int nPartition) {
	int loc_v;
	if (my_rank == 0)
		locVisited[0] = 1;
	else
		locVisited[0] = 0;
	for (loc_v = 1; loc_v < nPartition; loc_v++)
		locVisited[loc_v] = 0;
	for (loc_v = 0; loc_v < nPartition; loc_v++) {
		LocCost[loc_v] = LocAdjMatrix[my_rank * nPartition + loc_v];
		LocParent[loc_v] = 0;
	}
}

int go_tot(int LocCost[], int locVisited[], int nPartition) {
	int loc_u = -1, loc_v;
	int shortest_dist = INFINITY;
	for (loc_v = 0; loc_v < nPartition; loc_v++) {
		if (!locVisited[loc_v]) {
			if (LocCost[loc_v] < shortest_dist) {
				shortest_dist = LocCost[loc_v];
				loc_u = loc_v;
			}
		}
	}
	return loc_u;
}


void Dijkstra(int LocAdjMatrix[], int LocCost[], int LocParent[], int nPartition, int n,MPI_Comm comm) {

	int i, loc_v, loc_u, glbl_u, new_dist, my_rank, dist_glbl_u;
	int* locVisited;
	int myBest[2];
	int glBest[2];
	MPI_Comm_rank(comm, &my_rank);
	locVisited = (int*)malloc(nPartition * sizeof(int));
	Init(LocAdjMatrix, LocParent, LocCost, locVisited, my_rank, nPartition);
	for (i = 0; i < n - 1; i++) {
		loc_u = go_tot(LocCost, locVisited, nPartition);
		if (loc_u != -1) {
			myBest[0] = LocCost[loc_u];
			myBest[1] = loc_u + my_rank * nPartition;
		}
		else {
			myBest[0] = INFINITY;
			myBest[1] = -1;
		}
		MPI_Allreduce(myBest, glBest, 1, MPI_2INT, MPI_MINLOC, comm);

		dist_glbl_u = glBest[0];
		glbl_u = glBest[1];
		if (glbl_u == -1)
			break;
		if ((glbl_u / nPartition) == my_rank) {
			loc_u = glbl_u % nPartition;
			locVisited[loc_u] = 1;
		}

		for (loc_v = 0; loc_v < nPartition; loc_v++) {
			if (!locVisited[loc_v]) {

				new_dist = dist_glbl_u + LocAdjMatrix[glbl_u * nPartition + loc_v];
				if (new_dist < LocCost[loc_v]) {
					LocCost[loc_v] = new_dist;
					LocParent[loc_v] = glbl_u;
				}

			}
		}

	}
	free(locVisited);
}


 void printShortestPaths(int cost[], int n) {
	 for ( int i= 0 ; i<n;i++){
        if (cost[i] == INFINITY)cout<<"INF ";
        else cout<<cost[i]<<" ";
	 }
	 printf("\n");
 }

int main(int argc, char** argv) {
	ifstream inf("test.csv");

	char c;
	string ss; ss.clear();
	while (inf.get(c))
	{
		ss.push_back(c);
	}
	inf.close();
	ofstream of("test.csv");
	of << ss << endl;
	int* LocAdjMatrix, * LocCost, * LocParent, * cost = NULL, * parent = NULL;
	int my_rank, nPartition, n;
	MPI_Comm comm;

	MPI_Init(NULL, NULL);
	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &my_rank);
	MPI_Comm_size(comm, &p);
	n = atoi(argv[1]);
	if (n < p) {
		originalN = n;
		n = p;
		fix = 1;
	}
	nPartition = n / p;
	sendcounts = new int[p];
	sendcounts2 = new int[p];
	displs = new int[p];
	displs2 = new int[p];
	int rem = n % p;
	int sum = 0,sum2=0;
	for (int i = 0; i < p; i++) {
		sendcounts[i] = nPartition * n;
		if (rem > 0) {
			sendcounts[i] += n;
			rem--;
		}
		displs[i] = sum;
		displs2[i] = sum2;
		sum += sendcounts[i];
		sendcounts2[i] = sendcounts[i] / n;
		sum2 += sendcounts2[i];

	}
	nPartition = sendcounts[my_rank] / n;
	LocAdjMatrix = (int *)malloc(n * nPartition * sizeof(int));
	LocCost = (int*)malloc(nPartition * sizeof(int));
	LocParent = (int*)malloc(nPartition * sizeof(int));

	if (my_rank == 0) {
		cost = (int*)malloc(n * sizeof(int));
		parent = (int*)malloc(n * sizeof(int));
	}
	Gen_matrix(LocAdjMatrix, n, nPartition, my_rank, comm);
	auto start = std::chrono::steady_clock::now();
	Dijkstra(LocAdjMatrix, LocCost, LocParent, nPartition, n, comm);
	auto finish = std::chrono::steady_clock::now();
	MPI_Gatherv(LocCost, nPartition, MPI_INT, cost,sendcounts2,displs2, MPI_INT, 0, comm);
	MPI_Gatherv(LocParent, nPartition, MPI_INT, parent, sendcounts2, displs2, MPI_INT, 0, comm);


	if (my_rank == 0) {


		cout << (std::chrono::duration<double, std::milli>(finish - start).count()) << ',' << n << ',' << p;
		//printShortestPaths(cost, n);
		free(cost);
		free(parent);
	}
	free(LocAdjMatrix);
	free(LocParent);
	free(LocCost);
	MPI_Finalize();
	return 0;
}



