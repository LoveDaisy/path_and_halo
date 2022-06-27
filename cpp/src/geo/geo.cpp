#include "geo/geo.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include "core/math.hpp"
#include "core/types.hpp"
#include "util/log.hpp"

namespace halo_pm {

void Ll2Xyz(const Vec2f* ll, Vec3f* xyz,                    // input & output, ll in degree
            size_t num,                                     // data number
            size_t ll_step_bytes, size_t xyz_step_bytes) {  // step of input & output
  assert(ll_step_bytes % sizeof(float) == 0);
  assert(xyz_step_bytes % sizeof(float) == 0);

  const auto* p = reinterpret_cast<const float*>(ll);
  auto* q = reinterpret_cast<float*>(xyz);
  for (size_t i = 0; i < num; i++) {
    q[0] = std::cos(p[1] * kDegree2Rad) * std::cos(p[0] * kDegree2Rad);
    q[1] = std::cos(p[1] * kDegree2Rad) * std::sin(p[0] * kDegree2Rad);
    q[2] = std::sin(p[1] * kDegree2Rad);

    p += ll_step_bytes == 0 ? 2 : ll_step_bytes / sizeof(float);
    q += xyz_step_bytes == 0 ? 3 : xyz_step_bytes / sizeof(float);
  }
}

Vec3f Ll2Xyz(const Vec2f& ll) {
  Vec3f xyz;
  Ll2Xyz(&ll, &xyz);
  return xyz;
}


void Xyz2Ll(const Vec3f* xyz, Vec2f* ll,                    // input & output, ll in degree, xyz normalized
            size_t num,                                     // data number
            size_t xyz_step_bytes, size_t ll_step_bytes) {  // step of input & output
  assert(xyz_step_bytes % sizeof(float) == 0);
  assert(ll_step_bytes % sizeof(float) == 0);

  const auto* p_ptr = reinterpret_cast<const char*>(xyz);
  auto* q_ptr = reinterpret_cast<char*>(ll);
  for (size_t i = 0; i < num; i++) {
    const auto& curr_xyz = *reinterpret_cast<const Vec3f*>(p_ptr);
    auto& curr_ll = *reinterpret_cast<Vec2f*>(q_ptr);

    curr_ll.y() = std::asin(curr_xyz.z()) * kRad2Degree;
    curr_ll.x() = std::atan2(curr_xyz.y(), curr_xyz.x()) * kRad2Degree;

    p_ptr += xyz_step_bytes == 0 ? 3 * sizeof(float) : xyz_step_bytes;
    q_ptr += ll_step_bytes == 0 ? 2 * sizeof(float) : ll_step_bytes;
  }
}

Vec2f Xyz2Ll(const Vec3f& xyz) {
  Vec2f ll;
  Xyz2Ll(&xyz, &ll);
  return ll;
}


void Llr2Mat(const Vec3f* llr, Mat3x3f* mat,                  // input & output, llr in degree
             size_t num,                                      // data number
             size_t llr_step_bytes, size_t mat_step_bytes) {  // step of input & output
  assert(llr_step_bytes % sizeof(float) == 0);
  assert(mat_step_bytes % sizeof(float) == 0);
  if (llr_step_bytes == 0) {
    llr_step_bytes = 3 * sizeof(float);
  }
  if (mat_step_bytes == 0) {
    mat_step_bytes = 9 * sizeof(float);
  }

  for (size_t i = 0; i < num; i++) {
    const auto& p = *reinterpret_cast<const Vec3f*>(reinterpret_cast<const uint8_t*>(llr) + i * llr_step_bytes);
    auto c1 = std::cos(p(0) * kDegree2Rad);
    auto s1 = std::sin(p(0) * kDegree2Rad);
    auto c2 = std::cos(p(1) * kDegree2Rad);
    auto s2 = std::sin(p(1) * kDegree2Rad);
    auto c3 = std::cos(p(2) * kDegree2Rad);
    auto s3 = std::sin(p(2) * kDegree2Rad);

    auto& q = *reinterpret_cast<Mat3x3f*>(reinterpret_cast<uint8_t*>(mat) + i * mat_step_bytes);
    q(0, 0) = -s1 * c3 - c1 * s2 * s3;
    q(0, 1) = s1 * s3 - c1 * s2 * c3;
    q(0, 2) = c1 * c2;
    q(1, 0) = c1 * c3 - s1 * s2 * s3;
    q(1, 1) = -c1 * s3 - s1 * s2 * c3;
    q(1, 2) = s1 * c2;
    q(2, 0) = c2 * s3;
    q(2, 1) = c2 * c3;
    q(2, 2) = s2;
  }
}


Quatf Llr2Quat(const Vec3f& llr) {
  auto lon = kPi / 2.0f + llr.x() * kDegree2Rad;
  auto lat = kPi / 2.0f - llr.y() * kDegree2Rad;
  auto roll = llr.z() * kDegree2Rad;
  Quatf q_lon{ std::cos(lon / 2.0f), 0.0f, 0.0f, std::sin(lon / 2.0f) };
  Quatf q_lat{ std::cos(lat / 2.0f), std::sin(lat / 2.0f), 0.0f, 0.0f };
  Quatf q_roll{ std::cos(roll / 2.0f), 0.0f, 0.0f, std::sin(roll / 2.0f) };
  return q_lon * q_lat * q_roll;
}


Vec3f Quat2Llr(const Quatf& quat) {
  Curve3f basis{ { 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f } };
  for (auto& b : basis) {
    b = quat * b;
  }
  // rotated x: [..., ..., cos(lat)sin(roll)]
  // rotated y: [..., ..., cos(lat)cos(roll)]
  // rotated z: [cos(lat)cos(lon), cos(lat)sin(lon), sin(lat)]

  float lat = std::asin(basis[2].z()) * kRad2Degree;
  float lon = std::atan2(basis[2].y(), basis[2].x()) * kRad2Degree;
  float roll = std::atan2(basis[0].z(), basis[1].z()) * kRad2Degree;
  if (roll < 0) {
    roll += 360.0f;
  }

  return Vec3f{ lon, lat, roll };
}


Vec3f Quat2Llr(const Vec4f& qv) {
  Quatf q{ qv(0), qv(1), qv(2), qv(3) };
  return Quat2Llr(q);
}


AxisPdf::AxisPdf() : AxisPdf(0, -1, 0, -1) {}

AxisPdf::AxisPdf(float zenith_mean, float zenith_std) : AxisPdf(zenith_mean, zenith_std, 0, -1) {}

AxisPdf::AxisPdf(float zenith_mean, float zenith_std, float roll_mean, float roll_std)
    : zenith_mean_(zenith_mean), zenith_std_(zenith_std), roll_mean_(roll_mean), roll_std_(roll_std), zenith_int_c_(-1),
      roll_int_c_(-1) {
  zenith_int_c_ = 0.0f;
  auto dx = 180.0f / 50000.0f;
  for (float x = -90.0f; x <= 90.0f; x += dx) {
    auto x1 = (90.0f - x - zenith_mean_);
    auto p = std::exp(-x1 * x1 / 2.0f / (zenith_std_ * zenith_std_)) / zenith_std_;
    zenith_int_c_ += (p * std::cos(x * kDegree2Rad));
  }
  zenith_int_c_ *= dx * kPi / 180.0f * 2.0f * kPi;

  roll_int_c_ = 0.0f;
  dx = 360.0f / 50000.0f;
  for (float x = 0.0f; x <= 360.0f; x += dx) {
    auto x1 = (x - roll_mean_);
    auto p = std::exp(-x1 * x1 / 2.0f / (roll_std_ * roll_std_)) / roll_std_;
    roll_int_c_ += p;
  }
  roll_int_c_ *= dx * kPi / 180.0f;
}

float AxisPdf::operator()(const Vec3f& llr) const {
  float a = 1.0f / (4 * kPi);
  float r = 1.0f / (2 * kPi);

  if (zenith_std_ > 0) {
    auto lat = llr(1);
    if (lat > 90.0f) {
      lat = 180.0f - lat;
    }
    if (lat < -90.0f) {
      lat = -180.0f - lat;
    }
    auto z = 90.0f - lat - zenith_mean_;
    a = std::exp(-z * z / 2.0f / (zenith_std_ * zenith_std_)) / zenith_std_ / zenith_int_c_;
  }

  if (roll_std_ > 0) {
    auto roll = llr(2);
    roll -= std::floor(roll / 360.0f) * 360.0f;
    auto z = roll - roll_mean_;
    a = std::exp(-z * z / 2.0f / (roll_std_ * roll_std_)) / roll_std_ / roll_int_c_;
  }

  return a * r;
}


void RotateByQuat(const Quatf& quat, const Vec3f* xyz0, Vec3f* xyz1,  // input & output
                  size_t num,                                         // data number
                  size_t xyz0_step_bytes, size_t xyz1_step_bytes) {   // steps
  // Quaternion rotation. See wiki link for details: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
  // q =    [cos(theta/2),  [ux, uy, uz] * sin(theta/2)]
  // p =    [0,             [px, py, pz]]
  // q^-1 = [cos(theta/2), -[ux, uy, uz] * sin(theta/2)]
  // p' =   q * p * q^-1

  assert(xyz0_step_bytes % sizeof(float) == 0);
  assert(xyz1_step_bytes % sizeof(float) == 0);
  if (xyz0_step_bytes == 0) {
    xyz0_step_bytes = 3 * sizeof(float);
  }
  if (xyz1_step_bytes == 0) {
    xyz1_step_bytes = 3 * sizeof(float);
  }

  Quatf tmp{};

  for (size_t i = 0; i < num; i++) {
    const auto& p0 = *reinterpret_cast<const Vec3f*>(reinterpret_cast<const uint8_t*>(xyz0) + i * xyz0_step_bytes);
    tmp.x() = quat.w() * p0(0) + quat.y() * p0(2) - quat.z() * p0(1);
    tmp.y() = quat.w() * p0(1) - quat.x() * p0(2) + quat.z() * p0(0);
    tmp.z() = quat.w() * p0(2) + quat.x() * p0(1) - quat.y() * p0(0);
    tmp.w() = -quat.x() * p0(0) - quat.y() * p0(1) - quat.z() * p0(2);

    LOG_DEBUG("tmp=[% .4f,% .4f,% .4f,% .4f]", tmp.x(), tmp.y(), tmp.z(), tmp.w());

    auto& p1 = *reinterpret_cast<Vec3f*>(reinterpret_cast<uint8_t*>(xyz1) + i * xyz1_step_bytes);
    p1(0) = -tmp.w() * quat.x() + tmp.x() * quat.w() - tmp.y() * quat.z() + tmp.z() * quat.y();
    p1(1) = -tmp.w() * quat.y() + tmp.x() * quat.z() + tmp.y() * quat.w() - tmp.z() * quat.x();
    p1(2) = -tmp.w() * quat.z() - tmp.x() * quat.y() + tmp.y() * quat.x() + tmp.z() * quat.w();
    LOG_DEBUG("xyz1=%s", ObjLogFormatter<Vec3f>{ p1 }.Format());
  }
}


Curve2f Polygon2DIntersection(const Curve2f& p1, const Curve2f& p2) {
  // Calculate intersection of two *CONVEX* 2D polygon.
  // Iterate every edge of polygon 2, and cut off outer part of polygon 1.

  // First find out which side is inner side. Those have same sign to z0 considered as inner.
  auto d1 = p2[1] - p2[0];
  auto d2 = p2[2] - p2[1];
  float z0 = d1.x() * d2.y() - d1.y() * d2.x();

  // Start cutting
  Curve2f p_int;
  Curve2f new_p1 = p1;
  for (size_t i1 = 0; i1 < p2.size(); i1++) {
    size_t i2 = (i1 + 1) % p2.size();
    auto edge2 = p2[i2] - p2[i1];

    // Find out which point of vtx1 lies inner side,
    // and deal with edge crossing.
    p_int.clear();
    for (size_t j1 = 0; j1 < new_p1.size(); j1++) {
      size_t j2 = (j1 + 1) % new_p1.size();

      auto e1 = new_p1[j1] - p2[i2];
      auto e2 = new_p1[j2] - p2[i2];
      float z1 = edge2.x() * e1.y() - edge2.y() * e1.x();
      float z2 = edge2.x() * e2.y() - edge2.y() * e2.x();

      if (z1 * z0 > -kDefaultFloatEps && z2 * z0 > -kDefaultFloatEps) {
        // Two inner (at edge) points. Record the first one
        p_int.emplace_back(new_p1[j1]);
      } else if (z1 * z0 > -kDefaultFloatEps && z2 * z0 < kDefaultFloatEps) {
        // Inner (at edge) --> outer. Record the first point and the intersection point.
        p_int.emplace_back(new_p1[j1]);
        Mat<float, 2, 2> m;
        m << p2[i2] - p2[i1], new_p1[j2] - new_p1[j1];
        Vec2f s = m.partialPivLu().solve(p2[i1] - new_p1[j1]);
        p_int.emplace_back(s(1) * (new_p1[j2] - new_p1[j1]) + new_p1[j1]);
      } else if (z1 * z0 < kDefaultFloatEps && z2 * z0 > -kDefaultFloatEps) {
        // Outer --> inner (at edge). Record the intersection point.
        Mat<float, 2, 2> m;
        m << p2[i2] - p2[i1], new_p1[j2] - new_p1[j1];
        Vec2f s = m.partialPivLu().solve(p2[i1] - new_p1[j1]);
        p_int.emplace_back(s(1) * (new_p1[j2] - new_p1[j1]) + new_p1[j1]);
      }
    }
    new_p1.swap(p_int);
  }

  return new_p1;
}


std::tuple<Curve3f, float>  // (projected intersect polygon, area factor)
ProjectPolygonIntersection(const Curve3f& p1, const Curve3f& p2, const Vec3f& dir) {
  // The vertices of project-to polygon p2 must be coplannar.
  assert(p1.size() > 2);
  assert(p2.size() > 2);

  // Say original point is p0, projection ray is r, target polygon is qi, i = 1, 2, ...
  // and projection of p0 is pq, then we will see:
  // (pq - q1) - (p0 - q1) = t * r  ......(1)
  // In order to compute area we need to express pq in 2D form:
  // pq - q1 = [q2 - q1, q3 - q1] * uv = Bq * uv, where Bq is 3*2 and uv is 2*1. So (1) becomes
  // Bq * uv - p0 + q1 = t * r, or write in matrix:
  // [Bq, -r] * [uv; t] = p0 - q1  ......(2)

  Mat<float, 3, 3> mat1;
  mat1 << p2[1] - p2[0], p2[2] - p2[0], -dir;  // [q2 - q1, q3 - q1, -r]
  Curve2f uv1;
  for (const auto& p : p1) {
    Vec3f uvt = mat1.fullPivHouseholderQr().solve(p - p2[0]);
    uv1.emplace_back(Vec2f{ uvt(0), uvt(1) });
  }

  Mat<float, 3, 2> mat2;
  mat2 << p2[1] - p2[0], p2[2] - p2[0];
  Curve2f uv2;
  for (const auto& p : p2) {
    Vec2f uv = mat2.fullPivHouseholderQr().solve(p - p2[0]);
    uv2.emplace_back(uv);
  }

  auto uv_int = Polygon2DIntersection(uv1, uv2);
  float area_from = 0.0f;
  for (size_t i = 0; i < uv1.size(); i++) {
    size_t j = (i + 1) % uv1.size();
    area_from += (uv1[i].x() * uv1[j].y() - uv1[j].x() * uv1[i].y());
  }
  float area_int = 0.0f;
  for (size_t i = 0; i < uv_int.size(); i++) {
    size_t j = (i + 1) % uv_int.size();
    area_int += (uv_int[i].x() * uv_int[j].y() - uv_int[j].x() * uv_int[i].y());
  }
  Curve3f p_int;
  for (const auto& uv : uv_int) {
    Vec3f curr_p = (p2[1] - p2[0]) * uv.x() + (p2[2] - p2[0]) * uv.y();
    p_int.emplace_back(curr_p);
  }

  return { p_int, area_int / area_from };
}

}  // namespace halo_pm