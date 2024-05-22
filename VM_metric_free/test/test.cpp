#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> PDT;
typedef CGAL::Triangulation_data_structure_2<
  CGAL::Periodic_2_triangulation_vertex_base_2<GT>,
  CGAL::Periodic_2_triangulation_face_base_2<GT> >  Tds;
typedef Tds::Vertex_handle Vertex_handle;
typedef PDT::Iso_rectangle Iso_rectangle;
typedef PDT::Point          Point;
int main(int argc, char* argv[]) {
  CGAL::Timer t;
  typedef CGAL::Creator_uniform_2<double, Point> Creator;
  CGAL::Random random(7);
  CGAL::Random_points_in_square_2<Point, Creator> in_square(.5, random);
  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  int n = int(Lx * Ly);

  Iso_rectangle domain(0, 0, Lx, Ly);
  std::vector<Point> pts;
  PDT PT1(domain);
  PDT PT2(domain);
  PDT PT3(domain);
  // PDT PT1_copy(domain);

  // Generating n random points
  for (int i = 0 ; i < n ; i++)
  {
    Point p = *in_square;
    in_square++;
    pts.push_back(Point((p.x() + .5) * Lx, (p.y() + .5) * Ly));
  }
  // Standard insertion
  std::vector<Vertex_handle> vh_arr;
  vh_arr.reserve(n);
  t.start();
  for (int i = 0 ; i < n ; i++)
    vh_arr.push_back(PT1.insert(pts[i]));
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (Standard insertion)" << std::endl;
  t.reset();
  // Iterator range insertion using spatial sorting but no dummy points
  t.start();
  PT2.insert(pts.begin(), pts.end()); // third parameter defaults to false
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (with spatial sorting)" << std::endl;
  t.reset();
  // Iterator range insertion using spatial sorting and dummy point heuristic
  t.start();
  PT3.insert(pts.begin(), pts.end(), true);
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (Dummy point heuristic)" << std::endl;



  std::vector<Point> dis;
  for (int i = 0 ; i < n ; i++) {
    Point p = *in_square;
    in_square++;
    dis.push_back(p);
  }

  std::cout << vh_arr[3]->point().x() << std::endl;

  std::vector<Point> pts_new;
  // pts_new.reserve(n);
  for (int i = 0; i < n; i++) {
    auto p_old = vh_arr[i]->point();
    double x = p_old.x() + dis[i].x();
    double y = p_old.y() + dis[i].y();
    if (x < 0) {
      x += Lx;
    } else if (x >= Lx) {
      x -= Lx;
    }
    if (y < 0) {
      y += Ly;
    } else if (y >= Ly) {
      y -= Ly;
    }
    pts_new.push_back(Point(x, y));
  }
  std::cout << pts_new[0].x() << std::endl;

  t.reset();
  t.start();
  for (int i = 0; i < n; i++) {
    // vh_arr[i] = PT1.move_if_no_collision(vh_arr[i], pts_new[i]);
    auto vh_old = vh_arr[i];
    PT1.insert(pts_new[i], vh_old->face());
    PT1.remove(vh_old);
  }
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (move_if_no_collision)" << std::endl;
  std::cout << vh_arr[0]->point().x() << std::endl;

  return 0;
}