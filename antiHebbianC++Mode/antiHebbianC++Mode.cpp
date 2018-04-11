// standard include
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ppl.h>
#include "RNGProcess.h"
#include "antiHebbianC++Mode.h"
//#include <windows.h>
using namespace concurrency;
using namespace std;

double antiHebbian::max(double a, double b)
{
	if (a>b)
		return a;
	else
		return b;
}

double antiHebbian::min(double a, double b)
{
	if (a>b)
		return b;
	else
		return a;
}

void antiHebbian::initialAdjacency(void)
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
				if (randi(SIZE) < propor)
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

void antiHebbian::tongji(void)
{
	int i;
	cooperator1 = 0;
	defector1 = 0;
	for (i = 0; i<SIZE; i++)
	{
		switch (player_s[i])
		{
		case 0: cooperator1++; break;
		case 1: defector1++; break;
		}
	}
}

void antiHebbian::initial(void)
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
}

antiHebbian::antiHebbian(double b_set, int propor_set)
{
	b = b_set;
	propor = propor_set;
	//this->b = b;
	//this->propor = propor;
}
