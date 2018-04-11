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
	int vexnum;   //ͼ�Ķ������
	int edge;     //ͼ�ı���
	int **arc;   //�ڽӾ���
	int ** dis;   //��¼�����������·������Ϣ
	int ** path;  //��¼�������·������Ϣ

	int defector1, cooperator1, defector2, cooperator2;
	double score;
	double b;
	double delta;
	int propor;						//��ʼ���ڽӾ���ı���
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
	tomb3 player_adjacency;     //�ڽӾ���
	tomb1 player_degree;        //ÿ���ڵ�Ķ�
	tomb4 player_payoff;        //�������
	tomb4 player_delta;         //�����볡������Ĳ�ֵ



public:
	//���캯��
	antiHebbian(double b, int propor);
	//��������
	~antiHebbian();
	// �ж�����ÿ������ĵıߵ���Ϣ�Ƿ�Ϸ�
	//�����1��ʼ���
	bool check_edge_value(int start, int end, int weight);
	//����ͼ
	void createGraph(int);
	//��ӡ�ڽӾ���
	void print();
	//�����·��
	void Floyd();
	//��ӡ���·��
	void print_path();

	void prodgraph(void);      /* creates host graph                    */
	void initial(void);        /* initial state                         */
	void update(void);
	void tongji(void);
	double payoff(int i);       // ����i�������
	int degree(int i);          // ����ڵ� i �Ķ�
	void payoffCalculate(void); // �������棬�������player_payoff��
	void degreeCalculate(void); // ����ڵ�ȣ�����player_degree��
	void createLink(void);      //�����µ�����
	double max(double a, double b);
	double min(double a, double b);
	void initialAdjacency(void);//��ʼ���ڽӾ���
};