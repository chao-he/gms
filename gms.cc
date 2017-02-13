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
static const int kCountBits = 12;
static const int kMaxCount = (1<<kCountBits) - 1;

static int bandwidth{200};
static int innerloop{10};
static int outerloop{10};
static int max_samples{1000};
static int num_thread{16};

#define POS(id) ((id)>>kCountBits)
#define CNT(id) ((id) % (kMaxCount + 1))
#define COMBINE(id, cnt) (id<<kCountBits) + std::min(kMaxCount, (int)cnt)

struct Point {
  int32_t x;
  int32_t y;
  int32_t c;
  double w;

  Point(): x(0), y(0), c(1) {}
  Point(int x, int y): x(x), y(y), c(1) { }
  Point(int x, int y, int c): x(x), y(y), c(c) {}

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
    void Insert(const std::vector<Point> &points) {
      points1d.clear();
      points1d.reserve(points.size());
      for(const auto &pt: points) {
        uint64_t d = XY2D(pt);
        points1d.emplace_back(d);
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
        points1d_dup.emplace_back(COMBINE(d, c));
      }
      std::swap(points1d, points1d_dup);
    }


    void Find(const Point &where, int d, int n, std::vector<Point>& near) const {
      near.clear();
      uint64_t target = XY2D(where);
      auto iter = std::lower_bound(points1d.begin(), points1d.end(), target);
      size_t start = 0, end = points1d.size();
      if (iter != points1d.end()) {
        start = std::max<int>(0, iter - points1d.begin() - n);
        end = std::min(start + 2 * n, points1d.size());

        Point pt;
        int d2 = d*d;
        for(size_t i = start; i < end; ++ i) {
          D2XY(points1d[i], pt);
          if (where.DistanceTo(pt) < d2) {
            near.emplace_back(pt);
          }
        }
      }
    }
  private:
    std::vector<uint64_t> points1d;
};


template <typename Iter, typename KNN>
void Fit1(Iter iter, Iter end, const KNN& nn) {
  double R = 5 * bandwidth;
  double H = 2 * bandwidth * bandwidth;
  int N = 100;

  while (iter != end) {
    auto next = iter + 1;
    Point &pt = *iter;
    while (next != end && *next == *iter) {
      ++ next;
    }
    if (pt.c < max_samples) {
      std::vector<Point> nears;
      for (int j = 0; j < innerloop; ++ j) {
        nn.Find(pt, R, N, nears);
        if (nears.empty()) {
          break;
        }
        double n = 0, x = 0, y = 0;
        for(auto &qt: nears) {
          qt.w = 0.5 * qt.c * Exp(-qt.DistanceTo(pt) / H);
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
        pt.x = x;
        pt.y = y;

        if ((pt.x - x) * (pt.x - x) + (pt.y - y) * (pt.y - y) < 5) {
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
      std::vector<std::thread> workers;
      int n = num_thread;
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

void DumpPoints(const std::vector<Point>& points, const std::string& ofname) {
  char of[200];
  snprintf(of, sizeof(of), "%s-b_%d-l_%d-i_%d-m_%d", ofname.c_str(),
      bandwidth, outerloop, innerloop, max_samples);
  FILE* output = fopen(of, "w");
  for (const auto& pt: points) {
    double lat, lng;
    geo::UTM2LatLng(pt.x, pt.y, lat, lng);
    fprintf(output, "%.5f\t%.5f\t%d\n", lat, lng, pt.c);
  }
  fclose(output);
}

int main(int argc, char** argv) {
  int c = -1;
  std::string output;
  while (-1 != (c = getopt(argc, argv, "b:i:m:n:l:"))) {
    switch(c) {
      case 'b':
        bandwidth = atoi(optarg);
        break;
      case 'i':
        innerloop = atoi(optarg);
        break;
      case 'l':
        outerloop = atoi(optarg);
        break;
      case 'm':
        max_samples = atoi(optarg);
        break;
      case 'n':
        num_thread = atoi(optarg);
        break;
      case 'o':
        output = optarg;
        break;
      default:
        break;
    }
  }
  argc -= optind;
  argv += optind;

  if (argc == 0) {
    printf("Usage: gms -b bandwidth -i innerloop -l outerloop -m max_samples -n threads -o output input\n");
    return 0;
  }
  // fprintf(stderr, "%s %d %d\n", argv[0], gms.outerloop, gms.bandwidth);

  std::vector<Point> points;
  ReadPoints(argv[0], points);
  Fit(points);
  DumpPoints(points, output);
  return 0;
}
