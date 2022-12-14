/***
 * @Description: Geometry Define
 * @Author: jwimd chenjiewei@zju.edu.cn
 * @Date: 2022-11-12 19:16:44
 * @LastEditors: jwimd chenjiewei@zju.edu.cn
 * @LastEditTime: 2022-12-05 23:22:38
 * @FilePath: /GDB-HW6/src/Geometry.cpp
 */
#include "Geometry.h"
#include <cmath>
#include <memory>
#include <GL/freeglut.h>

#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>

#define NOT_IMPLEMENT -1.0

namespace hw6
{

    /*
     * Envelope functions
     */
    bool Envelope::contain(double x, double y) const
    {
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }

    /***
     * @Description: Calculate evenlope contain
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &envelope
     * @return {*}
     */
    bool Envelope::contain(const Envelope &envelope) const
    {
        return (minX <= envelope.getMinX() && minY <= envelope.getMinY() && maxX >= envelope.getMaxX() && maxY >= envelope.getMaxY());
    }

    /***
     * @Description: Calculate evenlope intersect
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &envelope
     * @return {*}
     */
    bool Envelope::intersect(const Envelope &envelope) const
    {
        if (this->contain(envelope.getMinX(), envelope.getMinY()))
            return true;
        if (this->contain(envelope.getMinX(), envelope.getMaxY()))
            return true;
        if (this->contain(envelope.getMaxX(), envelope.getMinY()))
            return true;
        if (this->contain(envelope.getMaxX(), envelope.getMaxY()))
            return true;

        if (envelope.contain(this->getMinX(), this->getMinY()))
            return true;
        if (envelope.contain(this->getMinX(), this->getMaxY()))
            return true;
        if (envelope.contain(this->getMaxX(), this->getMinY()))
            return true;
        if (envelope.contain(this->getMaxX(), this->getMaxY()))
            return true;

        if (this->getMaxX() >= envelope.getMaxX() && this->getMaxY() <= envelope.getMaxY() && this->getMinX() <= envelope.getMinX() && this->getMinY() >= envelope.getMinY())
            return true;
        if (this->getMaxX() <= envelope.getMaxX() && this->getMaxY() >= envelope.getMaxY() && this->getMinX() >= envelope.getMinX() && this->getMinY() <= envelope.getMinY())
            return true;

        return false;
    }

    /***
     * @Description: Calculate evenlope UNION
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &envelope
     * @return {*}
     */
    Envelope Envelope::unionEnvelope(const Envelope &envelope) const
    {
        return Envelope(std::min(this->minX, envelope.getMinX()), std::max(this->maxX, envelope.getMaxX()), std::min(this->minY, envelope.getMinY()), std::max(this->maxY, envelope.getMaxY()));
    }

    /***
     * @Description: Return min Distance Between point(x,y)
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {double} x
     * @param {double} y
     * @return {*}
     */
    double Envelope::minDistance(double x, double y) const
    {
        if (contain(x, y))
            return 0;

        if (x >= minX && x <= maxX)
            return std::min(fabs(y - maxY), fabs(y - minY));
        else if (y >= minY && y <= maxY)
            return std::min(fabs(x - maxX), fabs(x - minX));

        std::unique_ptr<Point> p1(new Point(minX, minY));
        std::unique_ptr<Point> p2(new Point(minX, maxY));
        std::unique_ptr<Point> p3(new Point(maxX, minY));
        std::unique_ptr<Point> p4(new Point(maxX, maxY));

        return std::min(std::min(Point(x, y).distance(p1.get()), Point(x, y).distance(p2.get())), std::min(Point(x, y).distance(p3.get()), Point(x, y).distance(p4.get())));
    }

    /***
     * @Description: Return max Distance Between point(x,y)
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {double} x
     * @param {double} y
     * @return {*}
     */
    double Envelope::maxDistance(double x, double y) const
    {
        std::unique_ptr<Point> p1(new Point(minX, minY));
        std::unique_ptr<Point> p2(new Point(minX, maxY));
        std::unique_ptr<Point> p3(new Point(maxX, minY));
        std::unique_ptr<Point> p4(new Point(maxX, maxY));

        return std::max(std::max(Point(x, y).distance(p1.get()), Point(x, y).distance(p2.get())), std::max(Point(x, y).distance(p3.get()), Point(x, y).distance(p4.get())));
    }

    void Envelope::draw() const
    {
        glBegin(GL_LINE_STRIP);

        glVertex2d(minX, minY);
        glVertex2d(minX, maxY);
        glVertex2d(maxX, maxY);
        glVertex2d(maxX, minY);
        glVertex2d(minX, minY);

        glEnd();
    }

    /*
     * Points functions
     */
    double Point::distance(const Point *point) const
    {
        return sqrt((x - point->x) * (x - point->x) +
                    (y - point->y) * (y - point->y));
    }

    /***
     * @Description: Distance between line and point
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {LineString} *line
     * @return {*}
     */
    double Point::distance(const LineString *line) const
    {
        double mindist = line->getPointN(0).distance(this);
        for (size_t i = 0; i < line->numPoints() - 1; ++i)
        {
            double dist = 0;
            double x1 = line->getPointN(i).getX();
            double y1 = line->getPointN(i).getY();
            double x2 = line->getPointN(i + 1).getX();
            double y2 = line->getPointN(i + 1).getY();
            // Task calculate the distance between Point P(x, y) and Line [P1(x1,
            // y1), P2(x2, y2)] (less than 10 lines)
            glm::vec2 vec_p1_p((x - x1), (y - y1));
            glm::vec2 vec_p2_p((x - x2), (y - y2));
            glm::vec2 vec_p1_p2((x2 - x1), (y2 - y1));

            if (glm::length(vec_p1_p) == 0 || glm ::angle(vec_p1_p, vec_p1_p2) > M_PI / 2)
                dist = glm::length(vec_p1_p);
            else if (glm::length(vec_p2_p) == 0 || glm::angle(-vec_p1_p2, vec_p2_p) > M_PI / 2)
                dist = glm::length(vec_p2_p);
            else
                dist = fabs(x * y1 + y * x2 + x1 * y2 - x * y2 - y1 * x2 - x1 * y) / glm::length(vec_p1_p2);

            if (dist < mindist)
                mindist = dist;
        }
        return mindist;
    }

    /***
     * @Description: Distance between point and polygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Polygon} *polygon
     * @return {*}
     */
    double Point::distance(const Polygon *polygon) const
    {
        LineString line = polygon->getExteriorRing();
        size_t n = line.numPoints();

        std::vector<LineString> innerRings = polygon->getInnerRings();

        for (std::vector<LineString>::iterator ring = innerRings.begin(); ring < innerRings.end(); ++ring)
        {
            for (size_t i = 0; i < ring->numPoints() - 1; ++i)
            {
                double x1 = line.getPointN(i).getX();
                double y1 = line.getPointN(i).getY();
                double x2 = line.getPointN(i + 1).getX();
                double y2 = line.getPointN(i + 1).getY();

                if ((x1 <= x || x2 <= x) &&
                    ((y1 >= y && y2 <= y && (x1 - x2) / (y1 - y2) * (y - y2) + x2 <= x) || (y1 <= y && y2 >= y && (x1 - x2) / (y1 - y2) * (y - y2) + x2 <= x)))
                    return this->distance(ring.base());
            }
        }

        double mindist = 0.0;

        if (!polygon->contains(this))
            mindist = this->distance(&line);
        return mindist;
    }

    /***
     * @Description: Distance between point and mutipoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPolygon
     * @return {*}
     */
    double Point::distance(const MutiPoint *mutiPoint) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPoint->numPoints(); ++i)
        {
            if (mindst < 0 || mutiPoint->getPointN(i).distance(this) < mindst)
                mindst = mutiPoint->getPointN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Calculate distance between MutiLineString and point
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiLineString} *mutiLineString
     * @return {*}
     */
    double Point::distance(const MutiLineString *mutiLineString) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiLineString->numLineStrings(); ++i)
        {
            if (mindst < 0 || mutiLineString->getLineStringN(i).distance(this) < mindst)
                mindst = mutiLineString->getLineStringN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Calculate distance between MutiPolygon and point
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPolygon} *mutiPolygon
     * @return {*}
     */
    double Point::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPolygon->numPolygons(); ++i)
        {
            if (mindst < 0 || mutiPolygon->getPolygonN(i).distance(this) < mindst)
                mindst = mutiPolygon->getPolygonN(i).distance(this);
        }

        return mindst;
    }

    bool Point::intersects(const Envelope &rect) const
    {
        return (x >= rect.getMinX()) && (x <= rect.getMaxX()) &&
               (y >= rect.getMinY()) && (y <= rect.getMaxY());
    }

    void Point::draw() const
    {
        glBegin(GL_POINTS);
        glVertex2d(x, y);
        glEnd();
    }

    /*
     * LineString functions
     */
    void LineString::constructEnvelope()
    {
        double minX, minY, maxX, maxY;
        maxX = minX = points[0].getX();
        maxY = minY = points[0].getY();
        for (size_t i = 1; i < points.size(); ++i)
        {
            maxX = std::max(maxX, points[i].getX());
            maxY = std::max(maxY, points[i].getY());
            minX = std::min(minX, points[i].getX());
            minY = std::min(minY, points[i].getY());
        }
        envelope = Envelope(minX, maxX, minY, maxY);
    }

    /***
     * @Description: To Calculate the distance between two LineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {LineString} *line
     * @return {*}
     */
    double LineString::distance(const LineString *line) const
    {
        double mindst = -1.0;

        if (this->intersects(line->getEnvelope())) // if envelope intersert
        {
            // calculate intersects
            if (this->intersects(line))
                return 0.0;
        }

        for (size_t i = 1; i < this->numPoints(); ++i)
        {
            if (mindst < 0 || this->points[i].distance(line) < mindst)
                mindst = this->points[i].distance(line);
        }

        for (size_t i = 1; i < line->numPoints(); ++i)
        {
            if (mindst < 0 || line->getPointN(i).distance(this) < mindst)
                mindst = line->getPointN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Calculate Distance between LineString and Polygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Polygon} *polygon
     * @return {*}
     */
    double LineString::distance(const Polygon *polygon) const
    {
        double mindst = -1.0;

        for (size_t i = 1; i < this->numPoints(); ++i)
            if (mindst < 0 || this->points[i].distance(polygon) < mindst)
                mindst = this->points[i].distance(polygon);

        return mindst;
    }

    /***
     * @Description: Distance between lineString and mutipoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPolygon
     * @return {*}
     */
    double LineString::distance(const MutiPoint *mutiPoint) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPoint->numPoints(); ++i)
        {
            if (mindst < 0 || mutiPoint->getPointN(i).distance(this) < mindst)
                mindst = mutiPoint->getPointN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Distance between lineString and mutilineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPolygon
     * @return {*}
     */
    double LineString::distance(const MutiLineString *mutiLineString) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiLineString->numLineStrings(); ++i)
        {
            if (mindst < 0 || mutiLineString->getLineStringN(i).distance(this) < mindst)
                mindst = mutiLineString->getLineStringN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description:  Distance between lineString and mutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPolygon} *mutiPolygon
     * @return {*}
     */
    double LineString::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPolygon->numPolygons(); ++i)
        {
            if (mindst < 0 || mutiPolygon->getPolygonN(i).distance(this) < mindst)
                mindst = mutiPolygon->getPolygonN(i).distance(this);
        }

        return mindst;
    }

    typedef int OutCode;

    const int INSIDE = 0; // 0000
    const int LEFT = 1;   // 0001
    const int RIGHT = 2;  // 0010
    const int BOTTOM = 4; // 0100
    const int TOP = 8;    // 1000

    // Compute the bit code for a point (x, y) using the clip rectangle
    // bounded diagonally by (xmin, ymin), and (xmax, ymax)
    // ASSUME THAT xmax, xmin, ymax and ymin are global constants.
    OutCode ComputeOutCode(double x, double y, double xmin, double xmax,
                           double ymin, double ymax)
    {
        OutCode code;

        code = INSIDE; // initialised as being inside of [[clip window]]

        if (x < xmin) // to the left of clip window
            code |= LEFT;
        else if (x > xmax) // to the right of clip window
            code |= RIGHT;
        if (y < ymin) // below the clip window
            code |= BOTTOM;
        else if (y > ymax) // above the clip window
            code |= TOP;

        return code;
    }

    // Cohen�CSutherland clipping algorithm clips a line from
    // P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
    // diagonal from (xmin, ymin) to (xmax, ymax).
    bool intersectTest(double x0, double y0, double x1, double y1, double xmin,
                       double xmax, double ymin, double ymax)
    {
        // compute outcodes for P0, P1, and whatever point lies outside the clip
        // rectangle
        OutCode outcode0 = ComputeOutCode(x0, y0, xmin, xmax, ymin, ymax);
        OutCode outcode1 = ComputeOutCode(x1, y1, xmin, xmax, ymin, ymax);
        bool accept = false;

        while (true)
        {
            if (!(outcode0 | outcode1))
            {
                // bitwise OR is 0: both points inside window; trivially accept and
                // exit loop
                accept = true;
                break;
            }
            else if (outcode0 & outcode1)
            {
                // bitwise AND is not 0: both points share an outside zone (LEFT,
                // RIGHT, TOP, or BOTTOM), so both must be outside window; exit loop
                // (accept is false)
                break;
            }
            else
            {
                // failed both tests, so calculate the line segment to clip
                // from an outside point to an intersection with clip edge
                double x, y;

                // At least one endpoint is outside the clip rectangle; pick it.
                OutCode outcodeOut = outcode0 ? outcode0 : outcode1;

                // Now find the intersection point;
                // use formulas:
                //   slope = (y1 - y0) / (x1 - x0)
                //   x = x0 + (1 / slope) * (ym - y0), where ym is ymin or ymax
                //   y = y0 + slope * (xm - x0), where xm is xmin or xmax
                // No need to worry about divide-by-zero because, in each case, the
                // outcode bit being tested guarantees the denominator is non-zero
                if (outcodeOut & TOP)
                { // point is above the clip window
                    x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
                    y = ymax;
                }
                else if (outcodeOut & BOTTOM)
                { // point is below the clip window
                    x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
                    y = ymin;
                }
                else if (outcodeOut &
                         RIGHT)
                { // point is to the right of clip window
                    y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
                    x = xmax;
                }
                else if (outcodeOut &
                         LEFT)
                { // point is to the left of clip window
                    y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
                    x = xmin;
                }

                // Now we move outside point to intersection point to clip
                // and get ready for next pass.
                if (outcodeOut == outcode0)
                {
                    x0 = x;
                    y0 = y;
                    outcode0 = ComputeOutCode(x0, y0, xmin, xmax, ymin, ymax);
                }
                else
                {
                    x1 = x;
                    y1 = y;
                    outcode1 = ComputeOutCode(x1, y1, xmin, xmax, ymin, ymax);
                }
            }
        }
        return accept;
    }

    bool LineString::intersects(const Envelope &rect) const
    {
        double xmin = rect.getMinX();
        double xmax = rect.getMaxX();
        double ymin = rect.getMinY();
        double ymax = rect.getMaxY();

        for (size_t i = 1; i < points.size(); ++i)
            if (intersectTest(points[i - 1].getX(), points[i - 1].getY(),
                              points[i].getX(), points[i].getY(), xmin, xmax, ymin,
                              ymax))
                return true;
        return false;
    }

    /***
     * @Description: Intersection between LineString(Update vertice method)
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {LineString} *line
     * @return {bool for intersects}
     */
    bool LineString::intersects(const LineString *line) const
    {
        for (size_t i = 0; i < line->numPoints() - 1; ++i)
        {
            double Ax, Ay, Bx, By;
            Ax = line->getPointN(i).getX();
            Ay = line->getPointN(i).getY();
            Bx = line->getPointN(i + 1).getX();
            By = line->getPointN(i + 1).getY();

            for (size_t i = 0; i < this->numPoints() - 1; ++i)
            {
                double Cx, Cy, Dx, Dy;
                Cx = this->getPointN(i).getX();
                Cy = this->getPointN(i).getY();
                Dx = this->getPointN(i + 1).getX();
                Dy = this->getPointN(i + 1).getY();

                if (std::max(Cx, Dx) < std::min(Ax, Bx) || std::max(Cy, Dy) < std::min(Ay, By) ||
                    std::max(Ax, Bx) < std::min(Cx, Dx) || std::max(Ay, By) < std::min(Cy, Dy))
                    continue;

                double z1 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Cx - Ax, Cy - Ay, 0.0))).z;
                double z2 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Dx - Ax, Dy - Ay, 0.0))).z;
                double z3 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Ax - Cx, Ay - Cy, 0.0))).z;
                double z4 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Bx - Ax, By - Ay, 0.0))).z;

                if (z1 * z2 <= 0 && z3 * z4 <= 0)
                    return true;
            }
        }

        return false;
    }

    void LineString::draw() const
    {
        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < points.size(); ++i)
            glVertex2d(points[i].getX(), points[i].getY());
        glEnd();
    }

    void LineString::print() const
    {
        std::cout << "LineString(";
        for (size_t i = 0; i < points.size(); ++i)
        {
            if (i != 0)
                std::cout << ", ";
            std::cout << points[i].getX() << " " << points[i].getY();
        }
        std::cout << ")";
    }

    /*
     * Polygon
     */
    /***
     * @Description: Update innerRing
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Polygon} *polygon
     * @return {*}
     */
    double Polygon::distance(const Polygon *polygon) const
    {

        double mindst = exteriorRing.distance(polygon);
        mindst = std::min(mindst, polygon->getExteriorRing().distance(this));

        for (size_t i = 0; i < innerRings.size(); i++)
            mindst = std::min(mindst, innerRings[i].distance(polygon));

        for (size_t i = 0; i < innerRings.size(); i++)
            mindst = std::min(mindst, polygon->getInnerRings()[i].distance(this));

        return mindst;
    }

    /***
     * @Description: Distance between polygon and mutipoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPolygon
     * @return {*}
     */
    double Polygon::distance(const MutiPoint *mutiPoint) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPoint->numPoints(); ++i)
        {
            if (mindst < 0 || mutiPoint->getPointN(i).distance(this) < mindst)
                mindst = mutiPoint->getPointN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Distance between polygon and mutilinestring
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPolygon
     * @return {*}
     */
    double Polygon::distance(const MutiLineString *mutiLineString) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiLineString->numLineStrings(); ++i)
        {
            if (mindst < 0 || mutiLineString->getLineStringN(i).distance(this) < mindst)
                mindst = mutiLineString->getLineStringN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Distance between polygon and mutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPolygon} *mutiPolygon
     * @return {*}
     */
    double Polygon::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPolygon->numPolygons(); ++i)
        {
            if (mindst < 0 || mutiPolygon->getPolygonN(i).distance(this) < mindst)
                mindst = mutiPolygon->getPolygonN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Intersection between Evenlope and Polygon(Version without considering innerRing)
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &rect
     * @return {*}
     */
    bool Polygon::intersects(const Envelope &rect) const
    {
        // point inside
        for (size_t i = 0; i < getExteriorRing().numPoints(); ++i)
            if (rect.contain(getExteriorRing().getPointN(i).getX(), getExteriorRing().getPointN(i).getY()))
                return true;

        std::vector<Point> points;
        points.push_back(Point(rect.getMinX(), rect.getMinY()));
        points.push_back(Point(rect.getMaxX(), rect.getMinY()));
        points.push_back(Point(rect.getMaxX(), rect.getMaxY()));
        points.push_back(Point(rect.getMinX(), rect.getMaxY()));

        for (size_t i = 0; i < points.size(); ++i)
        {
            if (contains(&points[i]))
                return true;
        }

        // line intersect

        LineString boundary(points);
        if (getExteriorRing().intersects(&boundary))
            return true;

        return false;
    }

    /***
     * @Description: Update vertice method to calculate
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Point} *point
     * @return {*}
     */
    bool Polygon::contains(const Point *point) const
    {
        // Task whether Point P(x, y) is within Polygon
        LineString line = getExteriorRing();
        size_t n = line.numPoints();

        double x = point->getX();
        double y = point->getY();

        bool inPolygon = false;

        for (size_t i = 0; i < line.numPoints() - 1; ++i)
        {
            double Cx = line.getPointN(i).getX();
            double Cy = line.getPointN(i).getY();
            double Dx = line.getPointN(i + 1).getX();
            double Dy = line.getPointN(i + 1).getY();

            double Ax = x;
            double Ay = y;
            double Bx = std::max(Cx, Dx);
            double By = y;

            if (Cx >= Ax && Cy == Ay)
                inPolygon = !inPolygon;

            if (Bx < x)
                continue;

            if (std::max(Cx, Dx) < std::min(Ax, Bx) || std::max(Cy, Dy) < std::min(Ay, By) ||
                std::max(Ax, Bx) < std::min(Cx, Dx) || std::max(Ay, By) < std::min(Cy, Dy))
                continue;

            double z1 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Cx - Ax, Cy - Ay, 0.0))).z;
            double z2 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Dx - Ax, Dy - Ay, 0.0))).z;
            double z3 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Ax - Cx, Ay - Cy, 0.0))).z;
            double z4 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Bx - Ax, By - Ay, 0.0))).z;

            if (z1 * z2 < 0 && z3 * z4 < 0)
                inPolygon = !inPolygon;
        }

        std::vector<LineString> innerRings = getInnerRings();

        for (std::vector<LineString>::iterator ring = innerRings.begin(); ring < innerRings.end(); ++ring)
        {
            for (size_t i = 0; i < ring->numPoints() - 1; ++i)
            {
                double Cx = line.getPointN(i).getX();
                double Cy = line.getPointN(i).getY();
                double Dx = line.getPointN(i + 1).getX();
                double Dy = line.getPointN(i + 1).getY();

                double Ax = x;
                double Ay = y;
                double Bx = std::max(Cx, Dx);
                double By = y;

                if (Cx >= Ax && Cy == Ay)
                    inPolygon = !inPolygon;

                if (Bx < x)
                    continue;

                if (std::max(Cx, Dx) < std::min(Ax, Bx) || std::max(Cy, Dy) < std::min(Ay, By) ||
                    std::max(Ax, Bx) < std::min(Cx, Dx) || std::max(Ay, By) < std::min(Cy, Dy))
                    continue;

                double z1 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Cx - Ax, Cy - Ay, 0.0))).z;
                double z2 = glm::normalize(glm::cross(glm::vec3(Bx - Ax, By - Ay, 0.0), glm::vec3(Dx - Ax, Dy - Ay, 0.0))).z;
                double z3 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Ax - Cx, Ay - Cy, 0.0))).z;
                double z4 = glm::normalize(glm::cross(glm::vec3(Dx - Cx, Dy - Cy, 0.0), glm::vec3(Bx - Ax, By - Ay, 0.0))).z;

                if (z1 + z2 == 0 && z1 != z2 && z3 + z4 == 0 && z3 != z4)
                    inPolygon = !inPolygon;
            }
        }

        return inPolygon;
    }

    /***
     * @Description: Update innerRings
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void Polygon::draw() const
    {
        exteriorRing.draw();
        for (size_t i = 0; i < this->innerRings.size(); i++)
        {
            innerRings[i].draw();
        }
    }

    /***
     * @Description: Update innerRings
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void Polygon::print() const
    {
        std::cout << "Polygon(";
        std::cout << "(";
        for (size_t i = 0; i < exteriorRing.numPoints(); ++i)
        {
            if (i != 0)
                std::cout << ", ";
            Point p = exteriorRing.getPointN(i);
            std::cout << p.getX() << " " << p.getY();
        }
        std::cout << ")";

        for (size_t i = 0; i < innerRings.size(); ++i)
        {
            if (i != innerRings.size() - 1)
                std::cout << ", ";

            std::cout << "(";
            for (size_t j = 0; j < innerRings[i].numPoints(); ++j)
            {
                if (j != 0)
                    std::cout << ", ";
                Point p = innerRings[i].getPointN(j);
                std::cout << p.getX() << " " << p.getY();
            }

            std::cout << ")";
        }
        std::cout << ")";
    }

    /***
     * @Description: Construct Evenlope For MutiPoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPoint::constructEnvelope()
    {
        double minX, minY, maxX, maxY;
        maxX = minX = points[0].getX();
        maxY = minY = points[0].getY();
        for (size_t i = 1; i < points.size(); ++i)
        {
            maxX = std::max(maxX, points[i].getX());
            maxY = std::max(maxY, points[i].getY());
            minX = std::min(minX, points[i].getX());
            minY = std::min(minY, points[i].getY());
        }
        envelope = Envelope(minX, maxX, minY, maxY);
    }

    /***
     * @Description: Calculate Distance between MutiPoints
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPoint
     * @return {*}
     */
    double MutiPoint::distance(const MutiPoint *mutiPoint) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < numPoints(); ++i)
            if (mindst < 0 || this->getPointN(i).distance(mutiPoint) < mindst)
                mindst = this->getPointN(i).distance(mutiPoint);

        return mindst;
    }

    /***
     * @Description: Calculate Distance between MutiPoint and MutiLineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPoint} *mutiPoint
     * @return {*}
     */
    double MutiPoint::distance(const MutiLineString *mutiLineString) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiLineString->numLineStrings(); ++i)
        {
            if (mindst < 0 || mutiLineString->getLineStringN(i).distance(this) < mindst)
                mindst = mutiLineString->getLineStringN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Calculate Distance between MutiPoint and MutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPolygon} *mutiPolygon
     * @return {*}
     */
    double MutiPoint::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPolygon->numPolygons(); ++i)
        {
            if (mindst < 0 || mutiPolygon->getPolygonN(i).distance(this) < mindst)
                mindst = mutiPolygon->getPolygonN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Intersect for MutiPoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &rect
     * @return {*}
     */
    bool MutiPoint::intersects(const Envelope &rect) const
    {
        for (size_t i = 0; i < numPoints(); ++i)
            if (getPointN(i).intersects(rect))
                return true;
        return false;
    }

    /***
     * @Description: Drawing MutiPoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPoint::draw() const
    {
        for (size_t i = 0; i < numPoints(); ++i)
            getPointN(i).draw();
    }

    /***
     * @Description: Print MutiPoint
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPoint::print() const
    {
        std::cout << "MutiPoint(";
        for (size_t i = 0; i < points.size(); ++i)
        {
            if (i != 0)
                std::cout << ", ";
            std::cout << points[i].getX() << " " << points[i].getY();
        }
        std::cout << ")";
    }

    /***
     * @Description: Construct Evenlope For MutiLineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiLineString::constructEnvelope()
    {

        double minX, minY, maxX, maxY;
        maxX = minX = getLineStringN(0).getPointN(0).getX();
        maxY = minY = getLineStringN(0).getPointN(0).getY();
        for (size_t index = 0; index < numLineStrings(); ++index)
        {
            for (size_t i = 1; i < getLineStringN(index).numPoints(); ++i)
            {
                maxX = std::max(maxX, getLineStringN(index).getPointN(i).getX());
                maxY = std::max(maxY, getLineStringN(index).getPointN(i).getY());
                minX = std::min(minX, getLineStringN(index).getPointN(i).getX());
                minY = std::min(minY, getLineStringN(index).getPointN(i).getY());
            }
        }
        envelope = Envelope(minX, maxX, minY, maxY);
    }

    /***
     * @Description: Calculate Distance between MutiLineStrings
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiLineString} *mutiLineString
     * @return {*}
     */
    double MutiLineString::distance(const MutiLineString *mutiLineString) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < numLineStrings(); ++i)
            if (mindst < 0 || this->getLineStringN(i).distance(mutiLineString) < mindst)
                mindst = this->getLineStringN(i).distance(mutiLineString);

        return mindst;
    }

    /***
     * @Description: Calculate Distance between MutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {MutiPolygon} *mutiPolygon
     * @return {*}
     */
    double MutiLineString::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < mutiPolygon->numPolygons(); ++i)
        {
            if (mindst < 0 || mutiPolygon->getPolygonN(i).distance(this) < mindst)
                mindst = mutiPolygon->getPolygonN(i).distance(this);
        }

        return mindst;
    }

    /***
     * @Description: Intersect for MutiLineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &rect
     * @return {*}
     */
    bool MutiLineString::intersects(const Envelope &rect) const
    {
        for (size_t i = 0; i < numLineStrings(); ++i)
            if (getLineStringN(i).intersects(rect))
                return true;
        return false;
    }

    /***
     * @Description: Drawing mutiLineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiLineString::draw() const
    {
        for (size_t i = 0; i < numLineStrings(); i++)
        {
            getLineStringN(i).draw();
        }
    }

    /***
     * @Description: Print MutiLineString
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiLineString::print() const
    {
        std::cout << "MutiLineString(";
        for (size_t index = 0; index < numLineStrings(); ++index)
        {
            if (index != 0)
                std::cout << ',';
            std::cout << "(";
            for (size_t i = 0; i < getLineStringN(index).numPoints(); ++i)
            {
                if (i != 0)
                    std::cout << ", ";
                std::cout << getLineStringN(index).getPointN(i).getX() << " " << getLineStringN(index).getPointN(i).getY();
            }
            std::cout << ")";
        }
        std::cout << ")";
    }

    /***
     * @Description: Construct Evenlope For MutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPolygon::constructEnvelope()
    {
        double minX, minY, maxX, maxY;
        maxX = minX = getPolygonN(0).getExteriorRing().getPointN(0).getX();
        maxY = minY = getPolygonN(0).getExteriorRing().getPointN(0).getY();
        for (size_t index = 0; index < numPolygons(); ++index)
        {
            for (size_t i = 1; i < getPolygonN(index).getExteriorRing().numPoints(); ++i)
            {
                maxX = std::max(maxX, getPolygonN(index).getExteriorRing().getPointN(i).getX());
                maxY = std::max(maxY, getPolygonN(index).getExteriorRing().getPointN(i).getY());
                minX = std::min(minX, getPolygonN(index).getExteriorRing().getPointN(i).getX());
                minY = std::min(minY, getPolygonN(index).getExteriorRing().getPointN(i).getY());
            }
        }
        envelope = Envelope(minX, maxX, minY, maxY);
    }

    double MutiPolygon::distance(const MutiPolygon *mutiPolygon) const
    {
        double mindst = -1.0;

        for (size_t i = 0; i < numPolygons(); ++i)
            if (mindst < 0 || this->getPolygonN(i).distance(mutiPolygon) < mindst)
                mindst = this->getPolygonN(i).distance(mutiPolygon);

        return mindst;
    }

    /***
     * @Description: Intersect for MutiPloygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @param {Envelope} &rect
     * @return {*}
     */
    bool MutiPolygon::intersects(const Envelope &rect) const
    {
        for (size_t i = 0; i < numPolygons(); ++i)
            if (getPolygonN(i).intersects(rect))
                return true;
        return false;
    }

    /***
     * @Description: Drawing mutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPolygon::draw() const
    {
        for (size_t i = 0; i < numPolygons(); i++)
        {
            getPolygonN(i).draw();
        }
    }

    /***
     * @Description: Print MutiPolygon
     * @Author: jwimd chenjiewei@zju.edu.cn
     * @msg: None
     * @return {*}
     */
    void MutiPolygon::print() const
    {
        std::cout << "MutiPolygon(";
        for (size_t index = 0; index < numPolygons(); ++index)
        {
            if (index != 0)
                std::cout << ',';
            std::cout << "(";
            std::cout << "(";
            for (size_t i = 0; i < getPolygonN(index).getExteriorRing().numPoints(); ++i)
            {
                if (i != 0)
                    std::cout << ", ";
                Point p = getPolygonN(index).getExteriorRing().getPointN(i);
                std::cout << p.getX() << " " << p.getY();
            }
            std::cout << ")";

            for (size_t i = 0; i < getPolygonN(index).getInnerRings().size(); ++i)
            {
                if (i != getPolygonN(index).getInnerRings().size() - 1)
                    std::cout << ", ";

                std::cout << "(";
                for (size_t j = 0; j < getPolygonN(index).getInnerRings()[i].numPoints(); ++j)
                {
                    if (j != 0)
                        std::cout << ", ";
                    Point p = getPolygonN(index).getInnerRings()[i].getPointN(j);
                    std::cout << p.getX() << " " << p.getY();
                }

                std::cout << ")";
            }
            std::cout << ")";
        }
        std::cout << ")";
    }

} // namespace hw6
