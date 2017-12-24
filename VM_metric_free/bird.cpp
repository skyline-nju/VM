#include "bird.h"

double Bird_wo_pos::v0;
double Bird_wo_pos::Lx;
double Bird_wo_pos::Ly;

void Bird_wo_pos::move(Point &p) const {
	double x = p[0] + v0 * vx;
	double y = p[1] + v0 * vy;
	if (x >= Lx) {
		x -= Lx;
	} else if (x < 0) {
		x += Lx;
	}
	if (y >= Ly) {
		y -= Ly;
	} else if (y < 0) {
		y += Ly;
	}
	p += Vec_2(x-p[0], y-p[1]);  
}

void Bird_wo_pos::vectorial_noise(double eta, Ran & myran) {
	double noise_x, noise_y;
	myran.circle_point_picking(noise_x, noise_y);
	double tmp = eta * n_neighbor;
	vx_next += noise_x * tmp;
	vy_next += noise_y * tmp;
	tmp = 1 / sqrt(vx_next * vx_next + vy_next * vy_next);
	vx_next *= tmp;
	vy_next *= tmp;
	vx = vx_next;
	vy = vy_next;
	n_neighbor = 1;
}

void Bird_wo_pos::scalar_noise(double eta, Ran & myran) {
	double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
	double c1 = vx / tmp;
	double s1 = vy / tmp;
	double noise = myran.doub() * eta * 2 * PI;
	double c2 = cos(noise);
	double s2 = sin(noise);
	vx = vx_next = c1 * c2 - s1 * s2;
	vy = vy_next = c1 * s2 + c2 * s1;
}

void Bird_wo_pos::ini_rand(std::vector<Bird_wo_pos> &birds,
													 std::vector<Pair_P_I> &pos_idx_pair, int nBird,
                           double LX, double LY, double V0, Ran &myran) {
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

// aligning with the nearest neighbors
void align(PDT &DT, std::vector<Pair_P_I> &pos_idx_pair, std::vector<Bird_wo_pos> &birds) {
	DT.insert(pos_idx_pair.begin(), pos_idx_pair.end());
	for (auto fit = DT.periodic_triangles_begin(PDT::UNIQUE); fit != DT.periodic_triangles_end(PDT::UNIQUE); ++fit) {
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

void move_with_vectorial_noise(double eta, Ran &myran,
	                             std::vector<Pair_P_I> &pos_idx_pair,
	                             std::vector<Bird_wo_pos> &bird) {
	unsigned int n = bird.size();
	for (unsigned int i = 0; i < n; i++) {
		// update velocity
		bird[i].vectorial_noise(eta, myran);
		// update coordination
		bird[i].move(pos_idx_pair[i].first);
	}
}

