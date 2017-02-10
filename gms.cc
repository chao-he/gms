#include <vector>
#include <string>
#include <stdint.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <fstream>
#include <thread>
#include <math.h>
#include "geo.h"
#include "util.h"

static const int kMaxDistance = 1<<15;

const static int kCountBits = 12;
const static int kMaxCount = (1<<kCountBits) - 1;

#define POS(id) ((id)>>kCountBits)
#define CNT(id) ((id) % (kMaxCount + 1))
#define COMBINE(id, cnt) (id<<kCountBits) + std::min(kMaxCount, (int)cnt)

struct Point {
  uint64_t id;
  int32_t x;
  int32_t y;
  int c;
  double w;

  Point(): id(0), x(0), y(0), c(1) {}

  Point(int x, int y): id(0), x(x), y(y), c(1) { }

  Point(int x, int y, int c):
    id(0), x(x), y(y), c(c) {}

  Point(int x, int y, int c, int64_t id):
    id(id), x(x), y(y), c(c) {}

  int DistanceTo(const Point &other) const {
    int dx = abs(x - other.x);
    int dy = abs(y - other.y);
    if (dx > kMaxDistance)
      return kMaxDistance;
    if (dy > kMaxDistance)
      return kMaxDistance;
    return dx*dx + dy*dy;
  }

  bool operator==(const Point &other) {
    return x == other.x && y == other.y;
  }

  bool operator<(const Point &other) {
    if (x == other.x) return y < other.y;
    else return x < other.x;
  }
};


inline uint64_t XY2D(const Point &pt) {
  S2CellId id = S2CellId::FromFaceIJ(0, pt.x, pt.y);
  return COMBINE(id.pos(), pt.c);
}

inline void D2XY(uint64_t pos, Point &pt) {
  S2CellId id = S2CellId::FromFacePosLevel(0, POS(pos), S2CellId::kMaxLevel);
  id.ToFaceIJOrientation(&pt.x, &pt.y, NULL);
  pt.c = CNT(pos);
}

class ApproximalNeighbor {
  public:
    void Insert(const Point &pt) {
      uint64_t d = XY2D(pt);
      points1d.push_back(d);
    }

    void Insert(const std::vector<Point> &points) {

      points1d.clear();
      points1d.reserve(points.size());

      for(const auto &pt: points) {
        Insert(pt);
      }

      // 合并相同的网格
      sort(points1d.begin(), points1d.end());

      std::vector<uint64_t> points1d_dup;

      for(auto iter = points1d.begin(); iter != points1d.end(); ) {
        uint64_t d = POS(*iter);
        uint64_t c = CNT(*iter);

        auto next = iter + 1;
        while (next != points1d.end() && d == POS(*next)) {
          c += CNT(*next);
          ++ next;
        }
        iter = next;
        points1d_dup.push_back(COMBINE(d, c));
      }
      std::swap(points1d, points1d_dup);
    }


    void Find(const Point &where, int d, std::vector<Point> &points) const {

      points.clear();

      uint64_t target = XY2D(where);
      auto iter = std::lower_bound(points1d.begin(), points1d.end(), target);

      size_t start = 0, end = points1d.size();
      if (iter != points1d.end()) {
        start = std::max<int>(0, iter - points1d.begin() - 100);
        end = std::min(start + 200, points1d.size());

        Point pt;
        int d2 = d*d;
        for(size_t i = start; i < end; ++ i) {
          D2XY(points1d[i], pt);
          if (where.DistanceTo(pt) < d2) {
            points.push_back(pt);
          }
        }
      }
    }
  private:
    std::vector<uint64_t> points1d;
};

int bandwidth{200};
int innerloop{10};
int outerloop{10};
int maxSamples{500};


template <typename Iter, typename KNN>
void Fit1(Iter iter, Iter end, const KNN& nn) {
  double R = 5 * bandwidth;
  double H = 2 * bandwidth * bandwidth;

  while (iter != end) {
    auto next = iter + 1;
    Point &pt = *iter;
    while (next != end && *next == *iter) {
      ++ next;
    }
    if (pt.c < maxSamples) {
      std::vector<Point> nears;
      for (int j = 0; j < innerloop; ++ j) {
        nears.clear();
        nn.Find(pt, R, nears);
        if (nears.empty()) {
          break;
        }
        double n = 0, x = 0, y = 0;
        for(auto &qt: nears) {
          qt.w = 0.5 * qt.c * Exp(-qt.DistanceTo(pt)/H);
          n += qt.w;
        }

        if (n < 1e-8) {
          break;
        }

        for(auto &qt: nears) {
          double w = qt.w / n;
          x += w * qt.x;
          y += w * qt.y;
        }
        int d = fabs(pt.x - x) + fabs(pt.x - x);
        pt.x = x;
        pt.y = y;
        if (d < 3) {
          break;
        }
      }
    }

    for(; iter != next; ++ iter) {
      iter->x = pt.x;
      iter->y = pt.y;
    }
  }
}

void Pack(std::vector<Point>& points) {
  std::vector<Point> dup;
  std::sort(points.begin(), points.end());
  for(auto iter = points.begin(); iter != points.end(); ) {
    auto next = iter + 1;
    while (next != points.end() && *iter == *next) {
      iter->c += next->c;
      ++ next;
    }
    dup.emplace_back(*iter);
    iter = next;
  }
  std::cerr << "Pack " << dup.size() << " -> " << points.size() << " points\n";
  std::swap(points, dup);
}

void Fit(std::vector<Point>& points) {
  ApproximalNeighbor nn;
  std::sort(points.begin(), points.end()); 
  for (int i = 0; i < outerloop; ++ i) {
    nn.Insert(points);
    if (points.size() < 100000) {
      Fit1(points.begin(), points.end(), nn);
    } else {
      const int n = 16;
      std::vector<std::thread> workers;
      int batch = points.size() / n;
      for (int i = 0; i < n; ++ i) {
        workers.emplace_back(std::thread([&points, &nn, i, n, batch] { 
              Fit1(points.begin() + batch * i,
                  i == n - 1 ? points.end() : points.begin() + (i + 1) * batch,
                  nn); }));
      }
      for (int i = 0; i < n; ++ i) {
        workers[i].join();
      }
    }
    Pack(points);
  }
}

void ReadPoints(const std::string& file, std::vector<Point>& points) {
  std::string line;
  std::ifstream ifs(file);
  double c,lat,lng,x,y;
  while(ifs>>lat>>lng>>c) {
    geo::LatLng2UTM(lat, lng, x, y);
    points.push_back(Point(x, y, c));
  }
  ifs.close();
}

int main(int argc, char** argv) {
  int c = -1;
  while (-1 != (c = getopt(argc, argv, "b:l:"))) {
    switch(c) {
      case 'b':
        bandwidth = atoi(optarg);
        break;
      case 'l':
        outerloop = atoi(optarg);
        break;
    }
  }
  argc -= optind;
  argv += optind;

  // fprintf(stderr, "%s %d %d\n", argv[0], gms.outerloop, gms.bandwidth);

  std::vector<Point> points;
  ReadPoints(argv[0], points);
  Fit(points);
  for (auto& pt: points) {
    double lat, lng;
    geo::UTM2LatLng(pt.x, pt.y, lat, lng);
    printf("%.5f\t%.5f\t%d\n", lat, lng, pt.c);
  }
  return 0;
}
