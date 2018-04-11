#pragma once
#include<iostream>
#include<string>
using namespace std;

// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           100      /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    61000   /* run-time in MCS     */
#define K           0.1      /* temperature */
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215



class antiHebbian {
public:
	int vexnum;   //图的顶点个数
	int edge;     //图的边数
	int **arc;   //邻接矩阵
	int ** dis;   //记录各个顶点最短路径的信息
	int ** path;  //记录各个最短路径的信息

	int defector1, cooperator1, defector2, cooperator2;
	double score;
	double b;
	double delta;
	int propor;						//初始化邻接矩阵的比例
	const double Upper = 1.2;
	const double Lower = 0.8;
	const double EPS = 1e-8;

	typedef int       tomb1[SIZE];
	typedef long int  tomb3[SIZE][SIZE];
	typedef double    tomb4[SIZE];
	typedef long double  tomb5[SIZE][SIZE];

	tomb1 player_s;           /* matrix ,containing players strategies */
	tomb1 player_sx;           /* matrix ,containing players strategies */
	tomb1 player_a;           /* matrix ,containing players strategies */
	tomb3 player_n;           /* matrix, containing players neighbours */
	tomb5 player_we;           /* matrix, containing players contribute */
	tomb3 player_adjacency;     //邻接矩阵
	tomb1 player_degree;        //每个节点的度
	tomb4 player_payoff;        //收益矩阵
	tomb4 player_delta;         //收益与场均收益的差值



public:
	//构造函数
	antiHebbian(double b, int propor);
	//析构函数
	~antiHebbian();
	// 判断我们每次输入的的边的信息是否合法
	//顶点从1开始编号
	bool check_edge_value(int start, int end, int weight);
	//创建图
	void createGraph(int);
	//打印邻接矩阵
	void print();
	//求最短路径
	void Floyd();
	//打印最短路径
	void print_path();

	void prodgraph(void);      /* creates host graph                    */
	void initial(void);        /* initial state                         */
	void update(void);
	void tongji(void);
	double payoff(int i);       // 个体i收益计算
	int degree(int i);          // 计算节点 i 的度
	void payoffCalculate(void); // 计算收益，放入矩阵player_payoff中
	void degreeCalculate(void); // 计算节点度，放入player_degree中
	void createLink(void);      //建立新的连接
	double max(double a, double b);
	double min(double a, double b);
	void initialAdjacency(void);//初始化邻接矩阵
};