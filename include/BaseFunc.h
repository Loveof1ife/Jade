#pragma once

using namespace acamcad;
using namespace polymesh;

namespace Jade{

    // have the function implementations of Cotangent, AngleBetween, and CalCircumCenter in the header file BaseTriFunc.h,
    // which is leading to multiple definitions of these functions.
    // When a header file with function definitions is included in multiple translation units (like in Curvature.cpp and main.cpp), the linker encounters multiple definitions of the same function.
    // Use inline: If the functions in BaseTriFunc.h are small and you want them to be defined in the header file
    // mark them as inline to prevent multiple definitions during linking.

    inline double Cotangent(const MVector3 &v1, const MVector3 &v2) {
        double cosTheta = v1.dot(v2) / (v1.norm() * v2.norm());
        double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        return cosTheta / sinTheta;
    }

    inline double AngleBetween(const MVector3 &a, const MVector3 &b) {
        double dotProduct = a.dot(b);
        double normsProduct = std::sqrt(a.norm() * b.norm());
        double cosTheta = dotProduct / normsProduct;
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // 防止超出 acos 的输入范围
        return std::acos(cosTheta); // 返回弧度值
    }

    inline MVector3 CalCircumCenter(const MVector3 &p1, const MVector3 &p2, const MVector3 &p3) {
        MVector3 edge1 = p2 - p1;
        MVector3 edge2 = p3 - p1;
        MVector3 v_cross = cross(edge1, edge2);
        // Norm squared of the cross product (magnitude squared of the area of the parallelogram)
        double denominator = 2.0 * norm(v_cross);

        // If the denominator is too small, the points are collinear (handle degenerate case)
        if (denominator < 1e-8) {
            return (p1 + p2 + p3) / 3.0; // Return centroid for degenerate case
        }
        // Compute the squared lengths of the sides
        double d1 = edge1.norm(); // squared length of p2 - p1
        double d2 = edge2.norm(); // squared length of p3 - p1

        // Compute the circumvent using determinant form
        MVector3 circumcenter = ((edge2 * d1) - (edge1 * d2)).cross(v_cross) / denominator;

        // Shift to the correct location relative to p1
        return circumcenter + p1;
    }

    inline  double Gaussian(const MVector3& deltaVector, const double& sigma)
    {
        double delta = (deltaVector).norm();
        return exp((-delta * delta) / (2 *  sigma * sigma) );
    }

    inline double Gaussian(const MVector3& ni, const MVector3& nj, const double& sigma)
    {
        double delta_n = (ni - nj).norm();
        return exp(sqrt(delta_n) / 2 / sqrt(sigma));
    }
    inline double Gaussian(const MPoint3& ci, const MPoint3& cj,  const double& sigma){
        double delta_c = (ci - cj).norm();
        return exp(sqrt(delta_c) / 2 / sqrt(sigma));
    }

}