#ifndef GRID_H
#define GRID_H
#include <vector>

template <typename TPar>
class Grid{
public:
  Grid(double Lx, double Ly, double l0);

  void cell_cell(int cell_idx);
  void cell_cell(int cell_idx1, int cell_idx2);
  void cell_cell(int cell_idx1, int cell_idx2, double, double);
  void all_pairs();

  void link_nodes(std::vector<TPar> &p_arr);

  void refresh(std::vector<TPar>& p_arr);

  int counter_particles() const;

  std::vector<TPar*> heads;
protected:
  double Lx_;
  double Ly_;
  double l0_;
  int mm_;
  int mx_;
  int my_;
};

template<typename TPar>
Grid<TPar>::Grid(double Lx, double Ly, double l0) {
  Lx_ = Lx;
  Ly_ = Ly;
  l0_ = l0;
  mx_ = int(Lx / l0);
  my_ = int(Ly / l0);
  mm_ = mx_ * my_;

  heads.reserve(mm_);
  for (int i = 0; i < mm_; i++) {
    heads.push_back(nullptr);
  }
}

template<typename TPar>
void Grid<TPar>::cell_cell(int cell_idx) {
  TPar* node1 = heads[cell_idx];
  TPar* node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      //node1->align(node2);
      node1->asym_align(node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template<typename TPar>
void Grid<TPar>::cell_cell(int cell_idx1, int cell_idx2) {
  if (heads[cell_idx2]) {
    TPar* node1 = heads[cell_idx1];
    TPar* node2;
    do {
      node2 = heads[cell_idx2];
      do {
        //node1->align(node2);
        node1->asym_align(node2);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

template<typename TPar>
void Grid<TPar>::cell_cell(int cell_idx1, int cell_idx2, double a, double b) {
  if (heads[cell_idx2]) {
    TPar* node1 = heads[cell_idx1];
    TPar* node2;
    do {
      node2 = heads[cell_idx2];
      do {
        //node1->align(node2, a, b);
        node1->asym_align(node2, a, b);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

template<typename TPar>
void Grid<TPar>::all_pairs() {
  int i, j;
  for (j = 0; j <= my_ - 2; j++) {
    int idx0 = j * mx_;
    if (heads[idx0]) {
      cell_cell(idx0);
      cell_cell(idx0, idx0 + 1);
      cell_cell(idx0, idx0 + mx_ + mx_ - 1, -Lx_, 0);
      cell_cell(idx0, idx0 + mx_);
      cell_cell(idx0, idx0 + mx_ + 1);
    }
    for (i = 1; i <= mx_ - 2; i++) {
      int idx1 = idx0 + i;
      if (heads[idx1]) {
        cell_cell(idx1);
        cell_cell(idx1, idx1 + 1);
        cell_cell(idx1, idx1 + mx_ - 1);
        cell_cell(idx1, idx1 + mx_);
        cell_cell(idx1, idx1 + mx_ + 1);
      }
    }
    int idx2 = idx0 + mx_ - 1;
    if (heads[idx2]) {
      cell_cell(idx2);
      cell_cell(idx2, idx2 - mx_ + 1, Lx_, 0);
      cell_cell(idx2, idx2 + 1, Lx_, 0);
      cell_cell(idx2, idx2 + mx_);
      cell_cell(idx2, idx2 + mx_ - 1);
    }
  }
  int idx0 = (my_ - 1) * mx_;
  if (heads[idx0]) {
    cell_cell(idx0);
    cell_cell(idx0, idx0 + 1);
    cell_cell(idx0, 0, 0, Ly_);
    cell_cell(idx0, 1, 0, Ly_);
    cell_cell(idx0, mx_ - 1, -Lx_, Ly_);
  }
  for (i = 1; i <= mx_ - 2; i++) {
    int idx1 = idx0 + i;
    if (heads[idx1]) {
      cell_cell(idx1);
      cell_cell(idx1, idx1 + 1);
      cell_cell(idx1, i - 1, 0, Ly_);
      cell_cell(idx1, i, 0, Ly_);
      cell_cell(idx1, i + 1, 0, Ly_);
    }
  }
  int idx2 = idx0 + mx_ - 1;
  if (heads[idx2]) {
    cell_cell(idx2);
    cell_cell(idx2, idx2 - mx_ + 1, Lx_, 0);
    cell_cell(idx2, mx_ - 2, 0, Ly_);
    cell_cell(idx2, mx_ - 1, 0, Ly_);
    cell_cell(idx2, 0, Lx_, Ly_);
  }
}

template<typename TPar>
void Grid<TPar>::link_nodes(std::vector<TPar>& p_arr) {
  for (auto& p : p_arr) {
    int col = int(p.x / l0_);
    if (col >= mx_) {
      col -= mx_;
    } else if (col < 0) {
      col += mx_;
    }

    int row = int(p.y / l0_);
    if (row >= my_) {
      row -= my_;
    } else if (row < 0) {
      row += my_;
    }

    int idx = col + row * mx_;
    p.next = heads[idx];
    heads[idx] = &p;
  }
}

template<typename TPar>
void Grid<TPar>::refresh(std::vector<TPar>& p_arr) {
  for (auto& head : heads) {
    head = nullptr;
  }
  link_nodes(p_arr);
}


template<typename TPar>
int Grid<TPar>::counter_particles() const {
  int n = 0;

  for (int i = 0; i < mm_; i++) {
    TPar* p = heads[i];
    while (p) {
      n++;
      p = p->next;
    }
  }
  return n;
}
#endif

