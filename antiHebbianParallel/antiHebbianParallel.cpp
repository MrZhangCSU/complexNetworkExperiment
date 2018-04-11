// standard include
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ppl.h>
//#include <windows.h>
using namespace concurrency;
using namespace std;

// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           100      /* lattice size                   */		//一共的个体
#define SIZE        1000//(L*L)    /* number of sites                */
#define MC_STEPS    61000   /* run-time in MCS     */
#define K           0.1      /* temperature */
#define MEMORY		2		// Memory of strategy  2,7,30
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215

int defector1, cooperator1, defector2, cooperator2;
double score;
double b;
double delta;
int propor;						//初始化邻接矩阵的比例
const double Upper = 1.2;
const double Lower = 0.8;
const double EPS = 1e-8;
bool floydFlagPass = false;
double ratioOfCoperator = 0;	//平均场中合作者比例
double ratioOfDefector = 0;		//平均场中背叛者比例
double meanFieldCooperatePayoff = 0.0;	//平均场合作者收益
double meanFieldDefectePayoff = 0.0;	//平均场背叛者收益
double averageDegree = 0.0;


typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][SIZE];
typedef double    tomb4[SIZE];
typedef long double  tomb5[SIZE][SIZE];
typedef long double  tomb6[SIZE][MEMORY];

tomb1 player_s;           /* matrix ,containing players strategies */
//tomb1 player_sx;           /* matrix ,containing players strategies */
tomb1 player_a;           /* matrix ,containing players strategies */
tomb3 player_n;           /* matrix, containing players neighbours */
//tomb5 player_we;           /* matrix, containing players contribute */
tomb3 player_adjacency;     //邻接矩阵
tomb1 player_degree;        //每个节点的度
tomb4 player_payoff;        //收益矩阵
tomb4 player_payoff_virtual;//虚拟收益矩阵
tomb4 player_delta;         //收益与场均收益的差值
tomb3 player_dis;			//弗洛伊德算法 D 矩阵
tomb3 player_path;			//弗洛伊德算法 P 矩阵
tomb6 player_memory;		//玩家的记忆矩阵

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


ofstream outfile1;
ofstream outfile2;
ofstream outfile3;
ofstream outfile4;


/*************************** RNG procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti = NN + 1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{
	int i;
	for (i = 0; i<NN; i++) {
		mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
	}
	mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{
	int i; for (i = 0; i<NN; i++) mt[i] = seed_array[i]; mti = NN;
}
double genrand()
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	if (mti >= NN)
	{
		int kk;
		if (mti == NN + 1) sgenrand(4357);
		for (kk = 0; kk<NN - MM; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MM] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk<NN - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MM - NN)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[NN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[NN - 1] = mt[MM - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
	return y;
}

double randf() { return ((double)genrand() * 2.3283064370807974e-10); }
long randi(unsigned long LIM) { return((unsigned long)genrand() % LIM); }

/********************** END of RNG ************************************/

double max(double a, double b)
{
	if (a>b)
		return a;
	else
		return b;
}

double min(double a, double b)
{
	if (a>b)
		return b;
	else
		return a;
}

/*********************初始化空间个体中的策略**************************/
void initial(void)
{
	cooperator1 = 0;
	defector1 = 0;
	//for (int i = 0; i < SIZE; i++)
	//{
	//	player_s[i] = (int)randi(2);
	//	switch (player_s[i])
	//	{
	//	case 0: cooperator1++;  break;
	//	case 1: defector1++;    break;
	//	}
	//}
	parallel_for(int(0), SIZE, [&](int i)
	{
		player_s[i] = (int)randi(2);
		switch (player_s[i])
		{
		case 0: cooperator1++;  break;
		case 1: defector1++;    break;
		}
	});

	parallel_for(int(0), SIZE, [&](int i)
	{
		for(int j=0; j<MEMORY; j++)
		{
			player_memory[i][j] = player_s[i];
		}
	});
}

/*************初始化邻接矩阵**************/
void initialAdjacency(int proportion)
{
	//int i, j;
	//int player1, player2;
	//for (int i = 0; i < SIZE; i++)
	//{
	//	for (int j = 0; j < SIZE; j++)
	//	{
	//		player_adjacency[i][j] = 0;

	//	}
	//}
	parallel_for(int(0), SIZE, [&](int i)
	{
		parallel_for(int(0), SIZE, [&](int j)
		{
			player_adjacency[i][j] = 0;
		});
	});
	//for (int i = 0; i < SIZE; i++)
	//{
	//	for (int j = 0; j < SIZE; j++)
	//	{
	//		if (i != j)
	//		{
	//			if (randi(50000) < proportion)
	//			{
	//				player_adjacency[i][j] = 1;
	//				player_adjacency[j][i] = 1;
	//			}
	//			else
	//			{
	//				player_adjacency[i][j] = 0;
	//				player_adjacency[j][i] = 0;
	//			}
	//		}
	//	}
	//}
	parallel_for(int(0), SIZE, [&](int i)
	{
		parallel_for(int(0), SIZE, [&](int j)
		{
			if (i != j)
			{
				if (randi(100*SIZE) < 8*proportion)							//冯诺依曼结构为5，高斯为9
				//if(randi(SIZE) < 9)
				{
					player_adjacency[i][j] = 1;
					player_adjacency[j][i] = 1;
				}			
				else
				{
					player_adjacency[i][j] = 0;
					player_adjacency[j][i] = 0;
				}			
			}
		});
	});

	//检测连通性
	int countNum = 0;
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			if (player_adjacency[i][j] != 0)
				countNum++;
		}
		if (countNum == 0)
			cout << "在第" << i << "行发现孤立点" << endl;
		countNum = 0;
	}
}

/************************* 个体收益计算 **************************/
double payoff(int i)
{
	int star1, star2;
	int player1, player2;
	double score = 0.0;
	player1 = i;
	star1 = player_s[player1];
	for (int j = 0; j<SIZE; j++)
	{
		if (j != player1)
		{
			player2 = j;
			star2 = player_s[player2];
			switch (2 * star1 + star2)
			{
			case 0:score = score + player_adjacency[player1][player2] * 1.0; break;
			case 1:score = score + player_adjacency[player1][player2] * 0.0; break;
			case 2:score = score + player_adjacency[player1][player2] * b; break;
			case 3:score = score + player_adjacency[player1][player2] * 0.0; break;
			}
		}
	}
	return score;
}

/************** 个体虚拟收益计算*********************/
double payoff_virtual(int i)
{
	int star1, star2;
	int player1, player2;
	double score = 0.0;
	player1 = i;
	star1 = player_s[player1];
	for (int j = 0; j<SIZE; j++)
	{
		if (j != player1)
		{
			player2 = j;
			star2 = player_s[player2];
			switch (2 * star1 + star2)
			{
				case 0:score = score + player_adjacency[player1][player2] * b; break;
				case 1:score = score + player_adjacency[player1][player2] * 0.0; break;
				case 2:score = score + player_adjacency[player1][player2] * 1.0; break;
				case 3:score = score + player_adjacency[player1][player2] * 0.0; break;
			}
		}
	}
	return score;
}
/************** 计算每个节点的收益和虚拟收益 ******************/
void payoffCalculate(void)
{
	//for (int i = 0; i < SIZE; i++)
	//{
	//	player_payoff[i] = payoff(i);

	//}
	parallel_for(int(0), SIZE, [&](int j)
	{
		player_payoff[j] = payoff(j);
	});
	
	parallel_for(int(0), SIZE, [&](int i)
	{
		player_payoff_virtual[i] = payoff_virtual(i);
	});

	//更新记忆矩阵
	static int count = 0;
	parallel_for(int(0), SIZE, [&](int k)
	{
		if(player_payoff[k] < player_payoff_virtual[k])
		{
			if(player_s[k] == 0)
				player_memory[k][count] = 1;
			else
				player_memory[k][count] = 0;
		}
		else{
			player_memory[k][count] = player_s[k];
		}
	});
		// for(int i=0; i<SIZE; i++)
		// {
		// 	if(player_payoff[i] < player_payoff_virtual[i])
		// 	{
		// 		if(player_s[i] == 0)
		// 			player_memory[i][count] = 1;
		// 		else
		// 			player_memory[i][count] = 0;
		// 	}
		// 	else
		// 	{
		// 		player_memory[i][count] = player_s[i];
		// 	}
		// }
		count++;
	if(count >= MEMORY)
		count = 0;
}

/*************** 计算该节点i的度 *********************/
int degree(int i)
{
	int player1;
	int k, j;
	k = 0;
	player1 = i;
	for (j = 0; j<SIZE; j++)
	{
		if (player_adjacency[player1][j] != 0)
		{
			k++;
		}
	}
	return k;
}

void degreeCalculate(void)
{
	//for (int i = 0; i < SIZE; i++)
	//{
	//	player_degree[i] = degree(i);
	//}
	parallel_for(int(0), SIZE, [&](int i)
	{
		player_degree[i] = degree(i);
	});
}

/****************建立新连接*********************/
void createLink(void)
{
	//parallel_for(int(0), SIZE, [&](int num)
	//{
	//	int player1, player2;
	//	double femi = 0.0;
	//	player1 = num;
	//	if (player_degree[player1] < SIZE)
	//	{
	//		if (player_degree[player1] < 50)								//增加节点度上限
	//		{
	//			do
	//			{
	//				do
	//				{
	//					player2 = (int)randi(SIZE);
	//				} while (player1 == player2);
	//			} while (player_adjacency[player1][player2] != 0);
	//			if (player_degree[player2] < 50)
	//			{
	//				femi = abs(player_payoff[player1] - player_payoff[player2]) / (max(player_degree[player1], player_degree[player2])*b + EPS);
	//				if (femi > randf())
	//				{
	//					player_adjacency[player1][player2] = 1;
	//					player_adjacency[player2][player1] = 1;
	//					player_degree[player1] += 1;
	//					player_degree[player2] += 1;
	//				}
	//			}
	//		}
	//	}
	//});
	/*********************************并行计算*****************************************/

	int i;
	int player1, player2;
	double femi = 0.0;
	double sumGain = 0.0;
	double averageGain = 0.0;
	double sumDelta = 0.0;	
	for (int num = 0; num < SIZE; num++)
	{
		player1 = num;
		//if (player_degree[player1] >= 50)
		//	continue;
		do
		{
			do
			{
				player2 = (int)randi(SIZE);
			} while (player1 == player2);
		} while (player_adjacency[player1][player2] != 0);
		//if (player_degree[player2] >= 50)
		//	continue;
		femi = abs(player_payoff[player1] - player_payoff[player2]) / (max(player_degree[player1], player_degree[player2])*b + EPS);
		if (femi > randf())
		{
			player_adjacency[player1][player2] = 1;
			player_adjacency[player2][player1] = 1;
			player_degree[player1] += 1;
			player_degree[player2] += 1;
		}
	}
	//同步更新方式
}

/****************************** 策略更新 **************************/
void update(void)
{
	for (int i = 0; i < SIZE; i++)
	{
		double average = 0;
		for (int j = 0; j < MEMORY; j++)
		{
			average += player_memory[i][j];
		}
		average /= MEMORY;
		if (average > randf()) 
			player_s[i] = 0;		//	选择合作
		else
			player_s[i] = 1;
	}
	//for (int i = 0; i < SIZE; i++)
	//{
	//	int j, k, choice;
	//	int player1, player2;
	//	int degree, degree1, degree2;
	//	int numberOfNeighbour = 0;                            //统计邻居有多少个
	//	int numberOfCount = 0;                                //统计数到了第几个邻居
	//	double numberOfRandom = 0.0;                          //生成的随机数
	//	double score1, score2, dscore;
	//	double femi;

	//	player1 = i;
	//	for (j = 0; j<SIZE; j++)                           //统计有多少个邻居
	//	{
	//		if (player_adjacency[player1][j] != 0)
	//		{
	//			numberOfNeighbour++;
	//		}
	//	}
	//	if (numberOfNeighbour != 0)
	//	{
	//		numberOfRandom = randi(numberOfNeighbour) + 1;  //随机挑选邻居
	//		for (k = 0; k<SIZE; k++)
	//		{
	//			if (player_adjacency[player1][k] != 0)
	//			{
	//				numberOfCount++;
	//				if (numberOfCount >= numberOfRandom)
	//				{
	//					choice = k; break;                   //挑选好之后退出循环
	//				}
	//			}
	//		}
	//		player2 = choice;

	//		degree1 = player_degree[player1];              //找到度大的节点
	//		degree2 = player_degree[player2];
	//		if (degree1 > degree2)
	//			degree = degree1;
	//		else
	//			degree = degree2;
	//		score1 = player_payoff[player1];
	//		score2 = player_payoff[player2];
	//		if (player_s[player1] != player_s[player2])
	//		{
	//			dscore = score2 / degree - score1 / degree;				//在这一步中已经确定了，为负时候肯定不变，为正再看概率
	//			femi = dscore / b;
	//			if (femi>randf()) player_s[player1] = player_s[player2];
	//		}
	//	}
	//}


	//parallel_for(int(0), SIZE, [&](int i)
	//{
	//	int j, k, choice;
	//	int player1, player2;
	//	int degree, degree1, degree2;
	//	int numberOfNeighbour = 0;                            //统计邻居有多少个
	//	int numberOfCount = 0;                                //统计数到了第几个邻居
	//	double numberOfRandom = 0.0;                          //生成的随机数
	//	double score1, score2, dscore;
	//	double femi;

	//	player1 = i;
	//	for (j = 0; j<SIZE; j++)                           //统计有多少个邻居
	//	{
	//		if (player_adjacency[player1][j] != 0)
	//		{
	//			numberOfNeighbour++;
	//		}
	//	}
	//	if (numberOfNeighbour != 0)
	//	{
	//		numberOfRandom = randi(numberOfNeighbour) + 1;  //随机挑选邻居
	//		for (k = 0; k<SIZE; k++)
	//		{
	//			if (player_adjacency[player1][k] != 0)
	//			{
	//				numberOfCount++;
	//				if (numberOfCount >= numberOfRandom)
	//				{
	//					choice = k; break;                   //挑选好之后推出循环
	//				}
	//			}
	//		}
	//		player2 = choice;

	//		degree1 = player_degree[player1];              //找到度大的节点
	//		degree2 = player_degree[player2];
	//		if (degree1 > degree2)
	//			degree = degree1;
	//		else
	//			degree = degree2;
	//		score1 = player_payoff[player1];
	//		score2 = player_payoff[player2];
	//		if (player_s[player1] != player_s[player2])
	//		{
	//			dscore = score2 / degree - score1 / degree;				//在这一步中已经确定了，为负时候肯定不变，为正再看概率
	//			femi = dscore / b;
	//			if (femi>randf()) player_s[player1] = player_s[player2];
	//		}
	//	}
	//});
}

/************************* 统计群体中各种策略的数目 **************************/
void tongji(void)
{
	int i;
	cooperator1 = 0;
	defector1 = 0;
	ratioOfCoperator = 0.0;
	ratioOfDefector = 0.0;
	meanFieldCooperatePayoff = 0.0;
	meanFieldDefectePayoff = 0.0;
	for (i = 0; i<SIZE; i++)
	{
		switch (player_s[i])
		{
		case 0: cooperator1++; break;
		case 1: defector1++; break;
		}
	}//统计合作者与背叛者人数
	ratioOfCoperator = cooperator1 / SIZE;
	ratioOfDefector = defector1 / SIZE;

	averageDegree = 0;
	for (int num = 0; num < SIZE; num++)
	{
		averageDegree += player_degree[num];
	}
	averageDegree = averageDegree / SIZE;		//计算平均度

	meanFieldCooperatePayoff = (ratioOfCoperator * 1 + ratioOfDefector * 0) * averageDegree;
	meanFieldCooperatePayoff = (ratioOfCoperator * b + ratioOfDefector * 0) * averageDegree;

	


}

/************************* 弗洛伊德算法进行网络连通性检测 **************************/
void Floyd(void) 
{
	int row = 0;
	int col = 0;
	for (row = 0; row < SIZE; row++) {
		for (col = 0; col < SIZE; col++) {
			//把矩阵D初始化为邻接矩阵的值
			player_dis[row][col] = player_adjacency[row][col];
			//矩阵中没有连接的点变为距离无穷大
			if (player_dis[row][col] == 0)
				player_dis[row][col] = INT_MAX;
			//矩阵P的初值则为各个边的终点顶点的下标
			player_path[row][col] = col;
		}
	}

	//三重循环，用于计算每个点对的最短路径
	int temp = 0;
	int select = 0;
	for (temp = 0; temp < SIZE; temp++) 
	{
		for (row = 0; row < SIZE; row++) 
		{
			for (col = 0; col < SIZE; col++) 
			{
				//为了防止溢出，所以需要引入一个select值
				select = (player_dis[row][temp] == INT_MAX || player_dis[temp][col] == INT_MAX) ? INT_MAX : (player_dis[row][temp] + player_dis[temp][col]);
				if (player_dis[row][col] > select) {
					//更新我们的D矩阵
					player_dis[row][col] = select;
					//更新我们的P矩阵
					player_path[row][col] = player_path[row][temp];
				}
			}
		}
		//cout<< temp << endl;
	}
}

void checkPath(void)
{
	bool breakFlag = false;				// true: 不通过   false:通过
	int row = 0;
	int col = 0;
	int temp = 0;
	for (row = 0; row < SIZE; row++)
	{
		for (col = 0; col < SIZE; col++)
		{
			if (player_dis[row][col] > SIZE)
			{
				cout << "第" << row << "与第" << col << "个节点断连," << "增加两点之间的连接" << endl;
				player_adjacency[row][col] = 1;
				player_adjacency[col][row] = 1;
				breakFlag = true;
				break;
			}
		}
		if (breakFlag)
			break;		
	}
	if (!breakFlag)
		floydFlagPass = true;
}

///******************************  主函数   ***************************/
//int main()
//{
//	int steps, propor;
//	int countNumber = 0;
//	int i, j;
//	double x, XX, aa, x1, XX1, aa1, averageDegree;
//	outfile1.open("frequency.txt");
//	outfile2.open("average.txt");
//	outfile3.open("degree.csv");
//	outfile4.open("adjacency.csv");
//	if (!outfile1)
//	{
//		cout << "can not open";
//		abort();
//	}
//	if (!outfile2)
//	{
//		cout << "can not open";
//		abort();
//	}
//	if (!outfile3)
//	{
//		cout << "can not open degree.csv";
//		abort();
//	}
//
//	// initialize the random number generation
//	sgenrand(RANDOMIZE);
//
//	cout << "Start the experiment!" << endl;
//	// begins the mc steps
//	int ii;
//	//for (b = 1.0; b<1.5; b = b + 0.1)
//	b = 1.25;
//	{
//		initial();                                      //初始化空间个体中的策略
//		//for (propor = 10; propor<60; propor = propor + 10)    //分别从0-50%的初始状态进行迭代
//		propor = 4;
//		{
//			cout << "The temptation to defect is :" << b << " and the proportion of beginning is:" << propor << endl;
//			aa = 0;
//			aa1 = 0;
//			initialAdjacency(propor);
//			payoffCalculate();
//			degreeCalculate();
//			for (steps = 0; steps<MC_STEPS; steps++)
//			{
//				update();                               //更新策略
//				createLink();                           //添加新连接
//				payoffCalculate();                      //计算节点收益
//				degreeCalculate();                      //计算节点的度(可以作废了)
//				tongji();
//				averageDegree = 0;
//				for (int num = 0; num < SIZE; num++)
//				{
//					averageDegree += player_degree[num];
//				}
//				averageDegree = averageDegree / SIZE;
//				x = (double)cooperator1 / SIZE;
//				x1 = (double)defector1 / SIZE;
//				outfile1 << steps << '\t' << x << '\t' << x1 << '\t' << averageDegree << endl;//C,D
//				cout << steps << '\t' << x << '\t' << x1 << '\t' << averageDegree << endl;
//				if ((x == 1) || (x1 == 1))
//				{
//					XX = x;
//					XX1 = x1;
//					break;
//				}
//
//				for (i = 0; i<SIZE; i++)                   //矩阵满时候停止本轮循环
//				{
//					for (j = 0; j<SIZE; j++)
//					{
//						if (player_adjacency[i][j] != 0)
//							countNumber++;
//					}
//				}
//				if (countNumber == SIZE*SIZE)
//					break;
//
//				if (steps > MC_STEPS - 1001)
//				{
//					aa += x;
//					aa1 += x1;
//					XX = (double)(aa) / 1000;
//					XX1 = (double)(aa1) / 1000;
//				}
//			}
//			outfile2 << propor << '\t' << b << '\t' << XX << '\t' << XX1 << endl;//C，D
//			cout << propor << '\t' << b << '\t' << XX << '\t' << XX1 << endl;
//
//			outfile3 << propor << ',' << b << endl;
//			for (int num = 0; num<SIZE; num++)
//			{
//				outfile3 << player_degree[num] << ',';
//				// cout << player_degree[num]<<'\t';
//			}
//			outfile3 << endl;
//			// cout<<endl;
//
//			for (int num1 = 0; num1 < SIZE; num1++)
//			{
//				for (int num2 = 0; num2 < SIZE; num2++)
//				{
//					outfile4 << player_adjacency[num1][num2] << ',';
//				}
//				outfile4 << endl;
//			}
//			outfile4 << endl;
//		}
//	}
//	outfile1.close();
//	outfile2.close();
//	outfile3.close();
//	outfile4.close();
//	printf("This is the end\n");
//	return 0;
//}
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////测试
int main(void)
{
	int steps, propor;
	int countNumber = 0;
	int i, j;
	double x, XX, aa, x1, XX1, aa1;
	outfile1.open("frequency.txt");
	outfile2.open("average.txt");
	outfile3.open("degree.csv");
	outfile4.open("adjacency.csv");
	if (!outfile1)
	{
		cout << "can not open";
		abort();
	}
	if (!outfile2)
	{
		cout << "can not open";
		abort();
	}
	if (!outfile3)
	{
		cout << "can not open degree.csv";
		abort();
	}
	if (!outfile4)
	{
		cout << "can not open adjacency.csv";
		abort();
	}
	// initialize the random number generation
	sgenrand(RANDOMIZE);
	
	cout << "Start the experiment!" << endl;

	// begins the mc steps
	int ii;
	//for (b = 1.0; b<1.5; b = b + 0.1)
	b = 1.2;
	{
		initial();                                      //初始化空间个体中的策略
														//for (propor = 10; propor<60; propor = propor + 10)    //分别从0-50%的初始状态进行迭代
		propor = 100;
		{
			cout << "The temptation to defect is :" << b << " and the proportion of beginning is:" << propor << endl;
			aa = 0;
			aa1 = 0;
			initialAdjacency(propor);
			
			while (!floydFlagPass)						//弗洛伊德算法实现网络联通
			{
				cout << "Floyd start!" << endl;
				Floyd();				//弗洛伊德算法测试连通性
				cout << "Check path!" << endl;
				checkPath();
			}

			payoffCalculate();
			degreeCalculate();
			tongji();
			for (steps = 0; steps<MC_STEPS; steps++)
			{
				update();                               //更新策略
				createLink();                           //添加新连接
				payoffCalculate();                      //计算节点收益
				//degreeCalculate();                      //计算节点的度(可以作废了)
				tongji();
				//averageDegree = 0;
				//for (int num = 0; num < SIZE; num++)
				//{
				//	averageDegree += player_degree[num];
				//}
				//averageDegree = averageDegree / SIZE;
				x = (double)cooperator1 / SIZE;
				x1 = (double)defector1 / SIZE;
				outfile1 << steps << '\t' << x << '\t' << x1 << '\t' << averageDegree << endl;//C,D
				//cout << steps << '\t' << x << '\t' << x1 << '\t' << averageDegree << endl;
				if ((x == 1) || (x1 == 1))
				{
					XX = x;
					XX1 = x1;
					break;
				}

				for (i = 0; i<SIZE; i++)                   //矩阵满时候停止本轮循环(没什么用了)
				{
					for (j = 0; j<SIZE; j++)
					{
						if (player_adjacency[i][j] != 0)
							countNumber++;
					}
				}
				if (countNumber == SIZE*SIZE)
					break;

				if (steps > MC_STEPS - 1001)
				{
					aa += x;
					aa1 += x1;
					XX = (double)(aa) / 1000;
					XX1 = (double)(aa1) / 1000;
				}
			}
			outfile2 << propor << '\t' << b << '\t' << XX << '\t' << XX1 << endl;//C，D
			cout << propor << '\t' << b << '\t' << XX << '\t' << XX1 << endl;

			outfile3 << propor << ',' << b << endl;
			for (int num = 0; num<SIZE; num++)
			{
				outfile3 << player_degree[num] << ',';
				// cout << player_degree[num]<<'\t';
			}
			outfile3 << endl;
			// cout<<endl;

			//for (int num1 = 0; num1 < SIZE; num1++)								//输出邻接矩阵
			//{
			//	for (int num2 = 0; num2 < SIZE; num2++)
			//	{
			//		outfile4 << player_adjacency[num1][num2] << ',';
			//	}
			//	outfile4 << endl;
			//}
			//outfile4 << endl;
		}
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	outfile4.close();
	printf("This is the end\n");
	return 0;












	//int propor = 20000;
	//initialAdjacency(propor);
	////degreeCalculate();
	//payoffCalculate();
	////for (int i = 0; i < SIZE; i++)
	////{
	////	cout << player_s[i] << '\t';
	////}
	////cout << endl;										//测试策略
	//for (int i = 0; i < SIZE; i++)
	//{
	//	cout << player_payoff[i] << '\t';
	//}
	//cout << endl<< endl;										//测试收益

	//for (int i = 0; i < SIZE; i++)
	//{
	//	for (int j = 0; j < SIZE; j++)
	//	{
	//		cout << player_adjacency[i][j] << '\t';
	//	}
	//	cout << endl;
	//}
	//cout << endl;										//测试邻接矩阵

	//for (int i = 0; i < SIZE; i++)
	//{
	//	cout << player_degree[i] << '\t';
	//}
	//cout << endl;										//测试节点度计算
}