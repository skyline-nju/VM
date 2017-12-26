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

// const double PI = 3.14159265358979;

// class for bird without positional information
struct Bird_wo_pos {
  Bird_wo_pos() {}
  Bird_wo_pos(double u, double v):
    vx(u), vy(v), vx_next(u), vy_next(v), n_neighbor(1) {}
  void move(Point &p) const;
	void vectorial_noise(double eta, Ran &myran);
	void scalar_noise(double eta, Ran &myran);
	template <class T1, class T2>
  static void ini_rand(std::vector<T1> &birds, std::vector<T2> &pos_idx_pair,
											 int nBird, double LX, double LY,  double V0, Ran &myran);

  double vx;
  double vy;
  double vx_next;
  double vy_next;
  int n_neighbor;

  static double v0;
  static double Lx;
  static double Ly;
};

template <class T1, class T2>
void Bird_wo_pos::ini_rand(std::vector<T1> &birds, std::vector<T2> &pos_idx_pair,
													 int nBird, double LX, double LY, double V0, Ran &myran) {
	birds.reserve(nBird);
	for (unsigned int i = 0; i < nBird; i++) {
		double x = myran.doub() * Lx;
		double y = myran.doub() * Ly;
		pos_idx_pair.push_back(std::make_pair(Point(x, y), i));
		double vx, vy;
		myran.circle_point_picking(vx, vy);
		birds.emplace_back(vx, vy);
	}
	Lx = LX;
	Ly = LY;
	v0 = V0;
}
// set n_neighbor for each bird
template <class T>
void set_n_neighbor(std::vector<T> &bird, unsigned int n) {
  for (auto it = bird.begin(); it != bird.end(); ++it) {
    (*it).n_neighbor = n;
  }
}

// aligning with the nearest neighbors
template <class T1, class T2>
void align(PDT &DT, std::vector<T1> &pos_idx_pairs, std::vector<T2> &birds) {
  DT.insert(pos_idx_pairs.begin(), pos_idx_pairs.end());
  	for (auto fit = DT.periodic_triangles_begin(PDT::UNIQUE);
         fit != DT.periodic_triangles_end(PDT::UNIQUE); ++fit) {
		unsigned int idx0 = fit.get_face()->vertex(0)->info();
		unsigned int idx1 = fit.get_face()->vertex(1)->info();
		unsigned int idx2 = fit.get_face()->vertex(2)->info();
		if (idx0 < idx1) {
			birds[idx0].vx_next += birds[idx1].vx;
			birds[idx0].vy_next += birds[idx1].vy;
			birds[idx1].vx_next += birds[idx0].vx;
			birds[idx1].vy_next += birds[idx0].vy;
			birds[idx0].n_neighbor++;
			birds[idx1].n_neighbor++;
		}
		if (idx1 < idx2) {
			birds[idx1].vx_next += birds[idx2].vx;
			birds[idx1].vy_next += birds[idx2].vy;
			birds[idx2].vx_next += birds[idx1].vx;
			birds[idx2].vy_next += birds[idx1].vy;
			birds[idx1].n_neighbor++;
			birds[idx2].n_neighbor++;
		}
		if (idx2 < idx0) {
			birds[idx2].vx_next += birds[idx0].vx;
			birds[idx2].vy_next += birds[idx0].vy;
			birds[idx0].vx_next += birds[idx2].vx;
			birds[idx0].vy_next += birds[idx2].vy;
			birds[idx2].n_neighbor++;
			birds[idx0].n_neighbor++;
		}
	}
	DT.clear();
}

template <class T1, class T2>
void move_with_vectorial_noise(double eta, Ran &myran,
															 std::vector<T1> &pos_idx_pair,
															 std::vector<T2> &birds) {
	unsigned int n = birds.size();
	for (unsigned int i = 0; i < n; i++) {
		// update velocity
		birds[i].vectorial_noise(eta, myran);
		// update coordination
		birds[i].move(pos_idx_pair[i].first);
	}
}

template <class T>
double cal_order_para(const std::vector<T> &bird) {
  double svx = 0;
  double svy = 0;
  for (auto it = bird.cbegin(); it != bird.cend(); ++it) {
    svx += (*it).vx;
    svy += (*it).vy;
  }
  return std::sqrt(svx * svx + svy * svy) / bird.size();
}

template <class T>
double cal_mean_neighbor(const std::vector<T> &bird) {
  double n = 0;
  for (auto it = bird.cbegin(); it != bird.cend(); ++it) {
    n += (*it).n_neighbor;
  }
  return n / bird.size();
}


#endif