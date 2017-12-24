#ifndef BIRD_H
#define BIRD_H
#include <fstream>
#include <vector>
#include <CGAL/Vector_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "rand.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Vector_2<K> Vec_2;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, GT, Vb> VbInfo;
typedef CGAL::Periodic_2_triangulation_face_base_2<GT> Fb;
typedef CGAL::Triangulation_data_structure_2<VbInfo, Fb> Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds> PDT;
typedef PDT::Point Point;
typedef std::pair<Point, unsigned int> Pair_P_I;

const double PI = 3.14159265358979;

// class for bird without positional information
struct Bird_wo_pos {
  Bird_wo_pos() {}
  Bird_wo_pos(double u, double v):
    vx(u), vy(v), vx_next(u), vy_next(v), n_neighbor(1) {}
  void move(Point &p) const;
	void vectorial_noise(double eta, Ran &myran);
	void scalar_noise(double eta, Ran &myran);
  static void ini_rand(std::vector<Bird_wo_pos> &birds,
											 std::vector<Pair_P_I> &pos_idx_pair, int nBird,
                       double LX, double LY,  double V0, Ran &myran);

  double vx;
  double vy;
  double vx_next;
  double vy_next;
  int n_neighbor;

  static double v0;
  static double Lx;
  static double Ly;
};

void align(PDT &DT, std::vector<Pair_P_I> &pos_idx_pair, std::vector<Bird_wo_pos> &birds);

void move_with_vectorial_noise(double eta, Ran &myran,
															 std::vector<Pair_P_I> &pos_idx_pair,
															 std::vector<Bird_wo_pos> &bird);

#endif