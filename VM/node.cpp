#include "node.h"

using namespace std;
double Node::Lx;
double Node::Ly;
double Node::v0 = 0.5;
double Node::rho_0;
int Node::N;

void Node::align(Node *node) {
  if (rr(node) < 1)
    addV(node);
}

void Node::align(Node *node, double a, double b) {
  if (rr(node, a, b) < 1)
    addV(node);
}

void Node::move(double noise) {
  double tmp = sqrt(vx*vx + vy*vy);
  double c1 = vx / tmp;
  double s1 = vy / tmp;
  double c2 = cos(noise);
  double s2 = sin(noise);
  vx = vx0 = c1 * c2 - s1 * s2;
  vy = vy0 = c1 * s2 + c2 * s1;

  x += v0*vx;
  if (x >= Lx)
    x -= Lx;
  else if (x < 0)
    x += Lx;
  y += v0*vy;
  if (y >= Ly)
    y -= Ly;
  else if (y < 0)
    y += Ly;

}

Node * Node::ini_rand(Ran * myran) {
  N = int(rho_0 * Lx * Ly);
  Node *bird = new Node[N];
  for (int i = 0; i < N; i++) {
    bird[i].x = myran->doub() * Lx;
    bird[i].y = myran->doub() * Ly;
    double theta = myran->doub() * 2 * PI;
    bird[i].vx = bird[i].vx0 = cos(theta);
    bird[i].vy = bird[i].vy0 = sin(theta);
    bird[i].next = nullptr;
  }
  return bird;
}

Node * Node::ini_move_left(Ran *myran) {
  N = int(rho_0 * Lx * Ly);
  Node *bird = new Node[N];
  for (int i = 0; i < N; i++) {
    bird[i].x = myran->doub() * Lx;
    bird[i].y = myran->doub() * Ly;
    bird[i].vx = bird[i].vx0 = 1;
    bird[i].vy = bird[i].vy0 = 0;
    bird[i].next = nullptr;
  }
  return bird;
}
Node * Node::ini_from_snap(double Lx0, double Ly0,
                           const vector<float>& x0,
                           const vector<float>& y0,
                           const vector<float>& theta0) {
  int nrows;
  int ncols;
  // check the compatibility of size of input snapshot
  if (int(Lx) % int(Lx0) == 0 && int(Ly) % int(Ly0) == 0) {
    ncols = int(Lx) / int(Lx0);
    nrows = int(Ly) / int(Ly0);
  } else {
    cout << "Error, the size of input snapshot is not right" << endl;
    exit(1);
  }
  int N0 = x0.size();
  N = N0 * nrows * ncols;
  Node *bird = new Node[N];
  for (int row = 0; row < nrows; row++) {
    double dy = row * Ly0;
    for (int col = 0; col < ncols; col++) {
      double dx = col * Lx0;
      for (int i = 0; i < N0; i++) {
        int j = i + (col + row * ncols) * N0;
        bird[j].x = x0[i] + dx;
        bird[j].y = y0[i] + dy;
        bird[j].vx = bird[j].vx0 = cos(theta0[i]);
        bird[j].vy = bird[j].vy0 = sin(theta0[i]);
        bird[j].next = nullptr;
      }
    }
  }
  return bird;
}
