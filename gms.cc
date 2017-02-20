#include <unistd.h>
#include <sys/syscall.h>
#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <thread>
#include <iostream>
#include <fstream>
#include "geo.h"
#include "util.h"

static int bandwidth = 30;
static int innerloop = 10;
static int outerloop = 50;
static int max_samples = 200;
static int num_thread = 40;
static int MIN_BATCH_SIZE = 1e6;
static int MIN_PACK_DISTANCE = 2;

struct Point1D {
  uint64_t id;
  uint32_t c;
  Point1D(uint64_t id, uint32_t c): id(id), c(c) {}
};

struct PointLL {
  double lat;
  double lng;
  uint32_t c;
  PointLL(double lat, double lng, uint32_t c):
    lat(lat), lng(lng), c(c) {}
};

struct PointXY {
  double x;
  double y;
  uint32_t c;
  PointXY(double x, double y, uint32_t c):
    x(x), y(y), c(c) {}
};

class Point {
  public:
    explicit Point(const Point1D& pt)
      : id_(pt.id), c_(pt.c), x_(0), y_(0)
    {
      double lat = 0, lng = 0;
      geo::Id2LatLng(id_, &lat, &lng, NULL);
      geo::LatLng2UTM(lat, lng, x_, y_);
    }

    explicit Point(const PointLL& pt)
      : id_(0), c_(pt.c), x_(0), y_(0)
    {
      geo::LatLng2UTM(pt.lat, pt.lng, x_, y_);
      id_ = geo::LatLng2Id(pt.lat, pt.lng, 30);
    }

    explicit Point(const PointXY& pt)
      : id_(0), c_(pt.c), x_(pt.x), y_(pt.y)
    {
      set_xy(pt.x, pt.y); 
    }

    inline uint64_t id() const { return id_; }
    inline uint32_t c() const { return c_; }
    inline double x() const { return x_; }
    inline double y() const { return y_; }
    void set_c(int c) {
      c_ = c;
    }
    void set_xy(double x, double y) {
      double lat = 0, lng = 0;
      geo::UTM2LatLng(x, y, lat, lng);
      id_ = geo::LatLng2Id(lat, lng, 30);
      x_ = x;
      y_ = y;
    }

  private:
    uint64_t id_;
    uint32_t c_;
    double x_;
    double y_;
};

struct PointNearBy {
  Point1D pt;
  double d;
  PointNearBy(const Point1D& pt, double d): pt(pt), d(d) {}
};

bool operator<(const Point1D& a, const Point1D& b) { return a.id < b.id; }
bool operator==(const Point1D& a, const Point1D& b) { return a.id == b.id; }

bool operator<(const Point& a, const Point& b) { return a.id() < b.id() || (a.id() == b.id() && a.c() > b.c()); }
bool operator==(const Point& a, const Point& b) { return a.id() == b.id(); }


class ApproximalNeighbor {
  public:
    template <typename iterator>
    ApproximalNeighbor(iterator beg, iterator end) {
     for(auto iter = beg; iter != end; ++ iter) {
       points1d_.emplace_back(iter->id(), iter->c());
      }
    }

    void Find(const Point &where, int d, int n, std::vector<PointNearBy>& near) const {
      Point1D target(where.id(), where.c());
      auto iter = std::lower_bound(points1d_.begin(), points1d_.end(), target);
      size_t start = 0, end = points1d_.size();
      if (iter != points1d_.end()) {
        start = std::max<int>(0, iter - points1d_.begin() - n);
        end = std::min(start + 2 * n, points1d_.size());

        int d2 = d*d;
        for(size_t i = start; i < end; ++ i) {
          auto d = geo::DistanceOfId(where.id(), points1d_[i].id);
          if (d >= d2) continue;
          near.emplace_back(points1d_[i], d);
        }
      }
    }

    void Project(const Point &where, int d, int n, double H, std::vector<PointNearBy>& near) const {
      Find(where, d, n, near);
      for (auto& p: near) {
        p.d = 0.5 * p.pt.c * Exp(-p.d / H);
      }
      double z = 1e-15;
      for (auto& p: near) {
        z += p.d;
      }
      for (auto& p: near) {
        p.d /= z;
      }
    }

  private:
    std::vector<Point1D> points1d_;
};


void Fit1(std::vector<Point>::iterator beg, std::vector<Point>::iterator end) {
  double R = 5 * bandwidth;
  double H = 2 * bandwidth * bandwidth;
  int N = 100;

  ApproximalNeighbor nn(beg, end);
  for (auto iter = beg; iter != end; ++ iter) {
    std::vector<PointNearBy> nears;
    Point &pt = *iter;
    for (int j = 0; j < innerloop; ++ j) {
      nears.clear();
      nn.Project(pt, R, N, H, nears);
      if (nears.empty()) {
        break;
      }
      double x = 0, y = 0;
      for (auto& p: nears) {
        Point q(p.pt);
        x += p.d * q.x();
        y += p.d * q.y();
      }
      pt.set_xy(x, y);
    }
  }
}

void Pack(std::vector<Point>& points) {
  std::vector<Point> dup;
  std::sort(points.begin(), points.end());
  for(auto iter = points.begin(), next = iter; iter != points.end(); iter = next) {
    next = iter + 1;
    double c = iter->c();
    while (next != points.end() && geo::DistanceOfId(next->id(), iter->id()) < MIN_PACK_DISTANCE) {
      c += next->c();
      ++ next;
    }
    iter->set_c(c);
    dup.emplace_back(*iter);
  }
  std::cerr << "Pack " << points.size() << " -> " << dup.size() << " points\n";
  std::swap(points, dup);
}

void Fit(std::vector<Point>& points) {
  std::sort(points.begin(), points.end()); 
  Pack(points);

  for (int i = 0; i < outerloop; ++ i) {
    std::vector<std::thread> workers;
    int n = std::min<int>(num_thread, ceil(points.size() / (double)MIN_BATCH_SIZE));
    int b = points.size() / n;
    for (int i = 0; i < n; ++ i) {
      workers.emplace_back(std::thread([&points, i, n, b] { 
            Fit1(points.begin() + b * i,
                i == n - 1 ? points.end() : points.begin() + (i + 1) * b); }));
    }
    for (int i = 0; i < n; ++ i) {
      workers[i].join();
    }
    Pack(points);
  }
}

void ReadPoints(const std::string& file, std::vector<Point>& points) {
  double lat = 0, lng = 0;
  std::ifstream ifs(file);
  while(ifs>>lat>>lng) {
    points.emplace_back(PointLL(lat, lng, 1));
  }
  ifs.close();
}

static int output_format = 0;

void DumpPoints(const std::vector<Point>& points, const std::string& ofname) {
  char of[200];
  snprintf(of, sizeof(of), "%s-b_%d-l_%d-i_%d-m_%d", ofname.c_str(),
      bandwidth, outerloop, innerloop, max_samples);
  FILE* output = fopen(of, "w");
  for (auto& pt: points) {
    switch(output_format) {
      case 0:
        fprintf(output, "%llu\t%u\n", pt.id(), pt.c());
        break;
      case 1:
        {
          double lat = 0, lng = 0;
          geo::Id2LatLng(pt.id(), &lat, &lng, NULL);
          fprintf(output, "%.5f\t%.5f\t%u\n", lat, lng, pt.c());
        }
        break;
    }
  }
  fclose(output);
}

int main(int argc, char** argv) {
  int c = -1;
  std::string output;
  while (-1 != (c = getopt(argc, argv, "b:i:m:n:l:o:d"))) {
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
      case 'd':
        output_format = 1;
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
