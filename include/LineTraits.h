#ifndef LINETRAITS_H_
#define LINETRAITS_H_
#include "Utility.h"

template<typename LineT> struct line_traits;
template<> struct line_traits<Segment2D> {
    // For isPointOnLine2D
    static bool validTScaled(double tScaled, double lenSq) noexcept {
        return (tScaled >= -epsilon) && (tScaled <= lenSq + epsilon);
    }

    // For isLinesIntersection2D
    static bool validIntersect(double tScaled, double uScaled, double denom) noexcept {
        if (denom > 0.0) {                           
            return (tScaled >= -epsilon) && (tScaled <= denom + epsilon) &&
                (uScaled >= -epsilon) && (uScaled <= denom + epsilon);
        }
        else {                                     
            return (tScaled <= epsilon) && (tScaled >= denom - epsilon) &&
                (uScaled <= epsilon) && (uScaled >= denom - epsilon);
        }
    }

};

template<> struct line_traits<Ray2D> {
    // For isPointOnLine2D
    static bool validTScaled(double tScaled, double) noexcept {
        return tScaled >= -epsilon;
    }

    // For isLinesIntersection2D
    static bool validIntersect(double tScaled, double uScaled, double denom) noexcept {
        if (denom > 0.0)
            return (tScaled >= -epsilon) && (uScaled >= -epsilon);
        else
            return (tScaled <= epsilon) && (uScaled <= epsilon);
    }
};

template<> struct line_traits<Line2D> {
    static bool validTScaled(double, double) noexcept { return true; }
    static bool validIntersect(double, double, double) noexcept { return true; }
};

template<> struct line_traits<Segment3D> {
    // For isPointOnLine3D
    static bool validTScaled(double tScaled, double lenSq) noexcept {
        return (tScaled >= -epsilon) && (tScaled <= lenSq + epsilon);
    }

    // For computeLinesDistance
    static void clampParams(double& t, double& u) noexcept {
        t = std::clamp(t, 0.0, 1.0);
        u = std::clamp(u, 0.0, 1.0);
    }

    // For isLinePolygonIntersection3D ºÍ isLinePolyhedronIntersection
    static bool validT(double t) noexcept {
        return t >= -epsilon && t <= epsilon_plus_1;
    }
    using Line2DType = Segment2D;
};

template<> struct line_traits<Ray3D> {
    static bool validTScaled(double tScaled, double) noexcept {
        return tScaled >= -epsilon;
    }

    static void clampParams(double& t, double& u) noexcept {
        t = std::max(t, 0.0);
        u = std::max(u, 0.0);
    }

    static bool validT(double t) noexcept {
        return t >= -epsilon;
    }
    using Line2DType = Ray2D;

};

template<> struct line_traits<Line3D> {
    static bool validTScaled(double, double) noexcept { return true; }
    static void clampParams(double&, double&) noexcept {}
    static bool validT(double) noexcept { return true; }
    using Line2DType = Line2D;
};

#endif