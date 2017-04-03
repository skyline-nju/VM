#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include "comn.h"


struct Node
{
	Node();
	double rr(Node *node);
	double rr(Node *node, double a, double b);
	void addV(Node *node);
	void align(Node *node);
	void align(Node *node, double a, double b);
	void move(double noise);

	static Node *ini_rand(Ran *myran);
	static Node *ini_move_left(Ran *myran);
	static Node *ini_copy_snap(double _eta, double _eps, double _rho0, double _Lx, double _Ly, unsigned long long _seed, int _t);


	double x, y, vx, vx0, vy, vy0;
	int cell_idx;
	Node* next;

	static double Lx;
	static double Ly;
	static double v0;
	static double rho_0;
	static int N;
};

inline Node::Node()
{
	x = y = vx = vx0 = vy = vy0;
	cell_idx = 0;
	next = nullptr;
}

inline double Node::rr(Node *node)
{
	double dx = node->x - x;
	double dy = node->y - y;
	return dx*dx + dy*dy;
}

inline double Node::rr(Node *node, double a, double b)
{
	double dx = node->x - x + a;
	double dy = node->y - y + b;
	return dx*dx + dy*dy;
}

inline void Node::addV(Node *node)
{
	vx += node->vx0;
	vy += node->vy0;
	node->vx += vx0;
	node->vy += vy0;
}

#endif

