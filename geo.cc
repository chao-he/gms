#include <algorithm>
#include <cmath>
#include "geometry/s2/s2cellid.h"
#include "geometry/s2/s2latlng.h"
#include "geo.h"

namespace geo {

  static const double kEarthRadius = 6378137;
  static const double kEarthCircumferenceMeters = 2 * M_PI * kEarthRadius; //40075017;

  inline double EarthMetersToRadians(double meters) {
    return meters / kEarthRadius;
  }

  inline double RadiansToEarthMeters(double radians) {
    return radians * kEarthRadius;
  }

  double Distance(double latitude1, double lngitude1, double latitude2, double lngitude2) {
    double dx = sin((latitude1 - latitude2) * M_PI / 360);
    double dy = sin((lngitude1 - lngitude2) * M_PI / 360);
    return 2 * kEarthRadius * asin(sqrt(dx*dx + cos(latitude1*M_PI/180)*cos(latitude2*M_PI/180)*dy*dy));
  }

  int DistanceOfId(uint64_t c1, uint64_t c2) {
    return RadiansToEarthMeters(
        S2CellId(c1).ToLatLng().GetDistance(
          S2CellId(c2).ToLatLng()).radians());
  }

  u64 LatLng2Id(double lat, double lng, int level) {
    return S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng)).parent(level).id();
  }

  u64 LatLng2Id(i32 lat, i32 lng, int level) {
    return S2CellId::FromLatLng(S2LatLng::FromE6(lat, lng)).parent(level).id();
  }

  void Id2LatLng(u64 id, double *lat, double *lng, int *level) {
    S2CellId sid(id);
    S2LatLng ll = sid.ToLatLng();
    if (lat) *lat = ll.lat().degrees();
    if (lng) *lng = ll.lng().degrees();
    if (level) *level = sid.level();
  }

  static const double MC2LL[][10] = {
    {
      1.410526172116255e-8,
      0.00000898305509648872,
      -1.9939833816331,
      200.9824383106796,
      -187.2403703815547,
      91.6087516669843,
      -23.38765649603339,
      2.57121317296198,
      -0.03801003308653,
      17337981.2
    },
    {
      -7.435856389565537e-9,
      0.000008983055097726239,
      -0.78625201886289,
      96.32687599759846,
      -1.85204757529826,
      -59.36935905485877,
      47.40033549296737,
      -16.50741931063887,
      2.28786674699375,
      10260144.86
    },
    {
      -3.030883460898826e-8,
      0.00000898305509983578,
      0.30071316287616,
      59.74293618442277,
      7.357984074871,
      -25.38371002664745,
      13.45380521110908,
      -3.29883767235584,
      0.32710905363475,
      6856817.37
    },
    {
      -1.981981304930552e-8,
      0.000008983055099779535,
      0.03278182852591,
      40.31678527705744,
      0.65659298677277,
      -4.44255534477492,
      0.85341911805263,
      0.12923347998204,
      -0.04625736007561,
      4482777.06
    },
    {
      3.09191371068437e-9,
      0.000008983055096812155,
      0.00006995724062,
      23.10934304144901,
      -0.00023663490511,
      -0.6321817810242,
      -0.00663494467273,
      0.03430082397953,
      -0.00466043876332,
      2555164.4
    },
    {
      2.890871144776878e-9,
      0.000008983055095805407,
      -3.068298e-8,
      7.47137025468032,
      -0.00000353937994,
      -0.02145144861037,
      -0.00001234426596,
      0.00010322952773,
      -0.00000323890364,
      826088.5
    }
  };

  static const double LL2MC[][10] = {
    {
      -0.0015702102444,
      111320.7020616939,
      1704480524535203,
      -10338987376042340,
      26112667856603880,
      -35149669176653700,
      26595700718403920,
      -10725012454188240,
      1800819912950474,
      82.5
    },
    {
      0.0008277824516172526,
      111320.7020463578,
      647795574.6671607,
      -4082003173.641316,
      10774905663.51142,
      -15171875531.51559,
      12053065338.62167,
      -5124939663.577472,
      913311935.9512032,
      67.5
    },
    {
      0.00337398766765,
      111320.7020202162,
      4481351.045890365,
      -23393751.19931662,
      79682215.47186455,
      -115964993.2797253,
      97236711.15602145,
      -43661946.33752821,
      8477230.501135234,
      52.5
    },
    {
      0.00220636496208,
      111320.7020209128,
      51751.86112841131,
      3796837.749470245,
      992013.7397791013,
      -1221952.21711287,
      1340652.697009075,
      -620943.6990984312,
      144416.9293806241,
      37.5
    },
    {
      -0.0003441963504368392,
      111320.7020576856,
      278.2353980772752,
      2485758.690035394,
      6070.750963243378,
      54821.18345352118,
      9540.606633304236,
      -2710.55326746645,
      1405.483844121726,
      22.5
    },
    {
      -0.0003218135878613132,
      111320.7020701615,
      0.00369383431289,
      823725.6402795718,
      0.46104986909093,
      2351.343141331292,
      1.58060784298199,
      8.77738589078284,
      0.37238884252424,
      7.45
    }  
  };

  // BMap.Projection.convertLL2MC(new BMap.Point(113.768261,23.036282))
  // => 12664762.69, 2619536
  //
  // const double LLBAND[] = {75, 60, 45, 30, 15, 0};

  #define CHECK_LATLNG(lat, lng) \
  if (lat >= 90 || lat <= -90 || lng >= 180 || lng <= -180) return false

  bool LatLng2UTM(double lat, double lng, double &x, double &y) {
    CHECK_LATLNG(lat, lng);
    int band  = int(std::floor(fabs(lat)));
    if (band >= 75) band = 0;
    else band = 5 - band/15;

    const double *C = LL2MC[band];
    const double alpha = fabs(lat) / C[9];

    // std::cout<<"band = "<<band<<", "<<C[0]<<"+"<<C[1]<<"*"<<fabs(lng)<<"="<<x<<std::endl;
    x = C[0] + C[1] * fabs(lng);
    y = C[2] + C[3] * alpha
      + C[4] * std::pow(alpha, 2)
      + C[5] * std::pow(alpha, 3)
      + C[6] * std::pow(alpha, 4)
      + C[7] * std::pow(alpha, 5)
      + C[8] * std::pow(alpha, 6);

    if (lng < 0) x *= -1;
    if (lat < 0) y *= -1;

    return true;
  }

  static const double MCBAND[] = {
    12890594.86, 8362377.87, 5591021, 3481989.83, 1678043.12, 0 };

  bool UTM2LatLng(double x, double y, double &lat, double &lng) {
    int band = 0;
    for(int i = 0; i < 6; ++ i) {
      if (fabs(y) >= MCBAND[i]) { band = i; break; }
    }
    const double *C = MC2LL[band];
    const double alpha = fabs(y) / C[9];

    // std::cout<<C[0]<<"+"<<C[1]<<"*"<<fabs(x)<<std::endl;
    //
    lng = C[0] + C[1] * fabs(x);
    lat = C[2] + C[3] * alpha
      + C[4] * std::pow(alpha, 2)
      + C[5] * std::pow(alpha, 3)
      + C[6] * std::pow(alpha, 4)
      + C[7] * std::pow(alpha, 5)
      + C[8] * std::pow(alpha, 6);
    if (x < 0) lng *= -1;
    if (y < 0) lat *= -1;
    return true;
  }

};

