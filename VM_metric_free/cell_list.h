#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <cstdio>
#include <vector>
//class Node {
//public:
//	Node() { next = NULL; }
//	double x;
//	double y;
//	int cell_idx;
//	Node * next;
//
//	static int Lx;
//	static int Ly;
//	static int N;
//};

template<class Node>
class Cell {
public:
	Cell() { head = NULL; }
	void cell_cell();
	void cell_cell(Cell<Node>*);
	void cell_cell(Cell<Node>*, double, double);
	static void all_pairs(Cell<Node>*cell, double Lx, double Ly);
  static void link_nodes(Cell<Node> *cell, std::vector<Node> &node);
  static void refresh(Cell<Node> *cell, std::vector<Node> &Node);
	//static Cell<Node> *ini(Node *bird, double Lx, double Ly);

	Node* head;

	static int mx;
	static int my;
	static int mm;
	static double l0;
};

template <class Node>
int Cell<Node>::mx;

template <class Node>
int Cell<Node>::my;

template <class Node>
int Cell<Node>::mm;

template <class Node>
double Cell<Node>::l0;

template <class Node>
void Cell<Node>::cell_cell() {
	Node *node1 = head;
	Node *node2;
	while (node1->next) {
		node2 = node1->next;
		do {
			node1->interact(node2);
			node2 = node2->next;
		} while (node2);
		node1 = node1->next;
	}
}

template <class Node>
void Cell<Node>::cell_cell(Cell<Node> *cell) {
	if (cell->head) {
		Node* node1 = head;
		Node* node2;
		do {
			node2 = cell->head;
			do {
				node1->interact(node2);
				node2 = node2->next;
			} while (node2);
			node1 = node1->next;
		} while (node1);
	}
}

template <class Node>
void Cell<Node>::cell_cell(Cell<Node> *cell, double a, double b) {
	if (cell->head) {
		Node* node1 = head;
		Node* node2;
		do {
			node2 = cell->head;
			do {
				node1->interact(node2, a, b);
				node2 = node2->next;
			} while (node2);
			node1 = node1->next;
		} while (node1);
	}
}

template <class Node>
void Cell<Node>::all_pairs(Cell<Node> * cell, double Lx, double Ly) {
	int i, j;
	Cell* p = cell;
	for (j = 0; j <= my - 2; j++) {
		if (p->head) {
			p->cell_cell();
			p->cell_cell(p + 1);
			p->cell_cell(p + mx + mx - 1, -Lx, 0);
			p->cell_cell(p + mx);
			p->cell_cell(p + mx + 1);
		}
		p++;
		for (i = 1; i <= mx - 2; i++) {
			if (p->head) {
				p->cell_cell();
				p->cell_cell(p + 1);
				p->cell_cell(p + mx - 1);
				p->cell_cell(p + mx);
				p->cell_cell(p + mx + 1);
			}
			p++;
		}
		if (p->head) {
			p->cell_cell();
			p->cell_cell(p - mx + 1, Lx, 0);
			p->cell_cell(p + 1, Lx, 0);
			p->cell_cell(p + mx);
			p->cell_cell(p + mx - 1);
		}
		p++;
	}
	if (p->head) {
		p->cell_cell();
		p->cell_cell(p + 1);
		p->cell_cell(cell, 0, Ly);
		p->cell_cell(cell + 1, 0, Ly);
		p->cell_cell(cell + mx - 1, -Lx, Ly);
	}
	p++;
	for (i = 1; i <= mx - 2; i++) {
		if (p->head) {
			p->cell_cell();
			p->cell_cell(p + 1);
			p->cell_cell(cell + i - 1, 0, Ly);
			p->cell_cell(cell + i, 0, Ly);
			p->cell_cell(cell + i + 1, 0, Ly);
		}
		p++;
	}
	if (p->head) {
		p->cell_cell();
		p->cell_cell(p - mx + 1, Lx, 0);
		p->cell_cell(cell + mx - 2, 0, Ly);
		p->cell_cell(cell + mx - 1, 0, Ly);
		p->cell_cell(cell, Lx, Ly);
	}
}

template <class Node>
void Cell<Node>::link_nodes(Cell<Node> *cell, std::vector<Node> &node) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    int col = int((*it).x);
    int row = int((*it).y);
    int idx_cell = (*it).cell_idx = col + mx * row;
    (*it).next = cell[idx_cell].head;
    cell[idx_cell].head = &(*it);
  }
}

template <class Node>
void Cell<Node>::refresh(Cell<Node> *cell, std::vector<Node> &node) {
  for (int i = 0; i < mm; i++) {
    cell[i].head = nullptr;
  }
  Cell::link_nodes(cell, node);
}
#endif