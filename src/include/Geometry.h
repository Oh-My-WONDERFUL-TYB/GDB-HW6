/***
 * @Description:Update Some Require For Q1 and Q2
 * @Author: jwimd chenjiewei@zju.edu.cn
 * @Date: 2022-11-12 19:16:44
 * @LastEditors: jwimd chenjiewei@zju.edu.cn
 * @LastEditTime: 2022-12-05 23:22:00
 * @FilePath: /GDB-HW6/src/include/Geometry.h
 */
#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

#include <algorithm>
#include <iostream>
#include <vector>

namespace hw6
{

  class Point;
  class LineString;
  class Polygon;
  class MutiPoint;
  class MutiLineString;
  class MutiPolygon;

  class Envelope
  {
  private:
    double minX;
    double minY;
    double maxX;
    double maxY;

  public:
    Envelope() : minX(0), minY(0), maxX(0), maxY(0) {}
    Envelope(double minX, double maxX, double minY, double maxY)
        : minX(minX), minY(minY), maxX(maxX), maxY(maxY) {}

    double getMinX() const { return minX; }
    double getMinY() const { return minY; }
    double getMaxX() const { return maxX; }
    double getMaxY() const { return maxY; }

    double getWidth() const { return maxX - minX; }
    double getHeight() const { return maxY - minY; }
    double getArea() const { return getWidth() * getHeight(); }

    bool contain(double x, double y) const;

    void draw() const;

    void print() const
    {
      std::cout << "Envelope( " << minX << " " << maxX << " " << minY << " "
                << maxY << ") ";
    }

    bool operator==(const Envelope &t1) const
    {
      return (minX == t1.minX && minY == t1.minY && maxX == t1.maxX &&
              maxY == t1.maxY);
    }
    bool operator!=(const Envelope &t1) const { return !(*this == t1); }

    bool contain(const Envelope &envelope) const;
    bool intersect(const Envelope &envelope) const;
    Envelope unionEnvelope(const Envelope &envelope) const;
  };

  /*
   * Geometry hierarchy
   */
  class Geometry
  {
  protected:
    Envelope envelope;

  public:
    Geometry() {}
    virtual ~Geometry() {}

    const Envelope &getEnvelope() const { return envelope; }

    virtual void constructEnvelope() = 0;
    virtual double distance(const Geometry *geom) const
    {
      return geom->distance(this);
    } // Euclidean distance
    virtual double distance(const Point *point) const = 0;
    virtual double distance(const LineString *line) const = 0;
    virtual double distance(const Polygon *polygon) const = 0;
    virtual double distance(const MutiPoint *mutiPoint) const = 0;
    virtual double distance(const MutiLineString *mutiLineString) const = 0;
    virtual double distance(const MutiPolygon *mutiPolygon) const = 0;

    virtual bool intersects(const Envelope &rect) const = 0;
    virtual void draw() const = 0;
    virtual void print() const = 0;
  };

  class Point : public Geometry
  {
  private:
    double x;
    double y;

  public:
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) { constructEnvelope(); }
    virtual ~Point() {}

    double getX() const { return x; }
    double getY() const { return y; }

    virtual void constructEnvelope() { envelope = Envelope(x, x, y, y); }

    // Euclidean distance
    virtual double distance(const Point *point) const;
    virtual double distance(const LineString *line) const;
    virtual double distance(const Polygon *polygon) const;
    virtual double distance(const MutiPoint *mutiPoint) const;
    virtual double distance(const MutiLineString *mutiLineString) const;
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    virtual void draw() const;

    virtual void print() const
    {
      std::cout << "Point(" << x << " " << y << ")";
    }
  };

  class LineString : public Geometry
  {
  private:
    std::vector<Point> points;

  public:
    LineString() {}
    LineString(std::vector<Point> &pts) : points(pts) { constructEnvelope(); }
    virtual ~LineString() {}

    size_t numPoints() const { return points.size(); }
    Point getStartPoint() const { return points.front(); }
    Point getEndPoint() const { return points.back(); }
    Point getPointN(size_t n) const { return points[n]; }

    virtual void constructEnvelope();

    // Euclidean distance
    virtual double distance(const Point *point) const
    {
      return point->distance(this);
    }
    virtual double distance(const LineString *line) const;
    virtual double distance(const Polygon *polygon) const;
    virtual double distance(const MutiPoint *mutiPoint) const;
    virtual double distance(const MutiLineString *mutiLineString) const;
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    bool intersects(const LineString *line) const;

    virtual void draw() const;

    virtual void print() const;
  };

  class Polygon : public Geometry
  {
  private:
    LineString exteriorRing;
    std::vector<LineString> innerRings;

  public:
    Polygon() {}
    Polygon(LineString &ering) : exteriorRing(ering) { constructEnvelope(); }

    //构建内环
    Polygon(LineString &ering, std::vector<LineString> &irings) : exteriorRing(ering), innerRings(irings) { constructEnvelope(); }

    virtual ~Polygon() {}

    LineString getExteriorRing() const { return exteriorRing; }

    std::vector<LineString> getInnerRings() const { return innerRings; }

    virtual void constructEnvelope()
    {
      envelope = exteriorRing.getEnvelope();
    }

    // Euclidean distance
    virtual double distance(const Point *point) const
    {
      return point->distance(this);
    }
    virtual double distance(const LineString *line) const
    {
      return line->distance(this);
    }
    virtual double distance(const Polygon *polygon) const;
    virtual double distance(const MutiPoint *mutiPoint) const;
    virtual double distance(const MutiLineString *mutiLineString) const;
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    virtual void draw() const;

    virtual void print() const;
  };

  class MutiPoint : public Geometry
  {
  private:
    std::vector<Point> points;

  public:
    MutiPoint() {}
    MutiPoint(std::vector<Point> _points) : points(_points) { constructEnvelope(); }
    virtual ~MutiPoint() {}

    size_t numPoints() const { return points.size(); }
    Point getStartPoint() const { return points.front(); }
    Point getEndPoint() const { return points.back(); }
    Point getPointN(size_t n) const { return points[n]; }

    virtual void constructEnvelope();

    // Euclidean distance
    virtual double distance(const Point *point) const { return point->distance(this); }
    virtual double distance(const LineString *line) const { return line->distance(this); }
    virtual double distance(const Polygon *polygon) const { return polygon->distance(this); }
    virtual double distance(const MutiPoint *mutiPoint) const;
    virtual double distance(const MutiLineString *mutiLineString) const;
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    virtual std::vector<Point> spilt() const { return points; }

    virtual void draw() const;

    virtual void print() const;
  };

  class MutiLineString : public Geometry
  {
  private:
    std::vector<LineString> lineStrings;

  public:
    MutiLineString() {}
    MutiLineString(std::vector<LineString> _lineStrings) : lineStrings(_lineStrings) { constructEnvelope(); }
    virtual ~MutiLineString() {}

    size_t numLineStrings() const { return lineStrings.size(); }
    LineString getStartLineString() const { return lineStrings.front(); }
    LineString getEndLineString() const { return lineStrings.back(); }
    LineString getLineStringN(size_t n) const { return lineStrings[n]; }

    virtual void constructEnvelope();

    // Euclidean distance
    virtual double distance(const Point *point) const { return point->distance(this); }
    virtual double distance(const LineString *line) const { return line->distance(this); }
    virtual double distance(const Polygon *polygon) const { return polygon->distance(this); }
    virtual double distance(const MutiPoint *mutiPoint) const { return mutiPoint->distance(this); }
    virtual double distance(const MutiLineString *mutiLineString) const;
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    virtual std::vector<LineString> spilt() const { return lineStrings; }

    virtual void draw() const;

    virtual void print() const;
  };

  class MutiPolygon : public Geometry
  {
  private:
    std::vector<Polygon> polygons;

  public:
    MutiPolygon() {}
    MutiPolygon(std::vector<Polygon> _polygons) : polygons(_polygons) { constructEnvelope(); }
    virtual ~MutiPolygon() {}

    size_t numPolygons() const { return polygons.size(); }
    Polygon getStartPolygon() const { return polygons.front(); }
    Polygon getEndPolygon() const { return polygons.back(); }
    Polygon getPolygonN(size_t n) const { return polygons[n]; }

    virtual void constructEnvelope();

    // Euclidean distance
    virtual double distance(const Point *point) const { return point->distance(this); }
    virtual double distance(const LineString *line) const { return line->distance(this); }
    virtual double distance(const Polygon *polygon) const { return polygon->distance(this); }
    virtual double distance(const MutiPoint *mutiPoint) const { return mutiPoint->distance(this); }
    virtual double distance(const MutiLineString *mutiLineString) const { return mutiLineString->distance(this); }
    virtual double distance(const MutiPolygon *mutiPolygon) const;

    // intersection test with the envelope for range query
    virtual bool intersects(const Envelope &rect) const;

    virtual std::vector<Polygon> spilt() const { return polygons; }

    virtual void draw() const;

    virtual void print() const;
  };

} // namespace hw6

#endif