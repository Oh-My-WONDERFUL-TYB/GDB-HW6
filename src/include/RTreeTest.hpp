#ifndef RTREE_TEST_H_INCLUDED
#define RTREE_TEST_H_INCLUDED

#include "Common.h"
#include "RTreeSrc.hpp"

using namespace hw6;

extern int mode;
extern std::vector<Geometry *> readGeom(const char *filename);
extern std::vector<std::string> readName(const char *filename);
extern void transformValue(double &res, const char *format);
extern void wrongMessage(Envelope e1, Envelope e2, bool cal);
extern void wrongMessage(const Point &pt1, const Point &pt2, double dis,
                         double res);
extern void wrongMessage(Envelope e1, Envelope e2, Envelope cal, Envelope res);

namespace hw6
{

    template <uint8_t M>
    void RTree<M>::test(int t)
    {
        using namespace std;

        std::cout << "*********************Start*********************" << std::endl;
        if (t == TEST1)
        {
            std::cout << "????1: Envelope Contain, Intersect, and Union" << endl;

            int failedCase = 0;
            Envelope e1(-1, 1, -1, 1);
            vector<Envelope> tests;
            tests.push_back(Envelope(-0.5, 0.5, -0.5, 0.5));
            tests.push_back(Envelope(-0.5, 0.5, 0.5, 1.5));
            tests.push_back(Envelope(0.5, 1.5, -0.5, 0.5));
            tests.push_back(Envelope(-1.5, -0.5, -1.5, -0.5));
            tests.push_back(Envelope(-2, -1, -0.5, 0.5));
            tests.push_back(Envelope(1, 1.5, 1, 1.5));
            tests.push_back(Envelope(-2, -1.5, -0.5, 0.5));
            tests.push_back(Envelope(-0.5, 0.5, 1.5, 2));
            tests.push_back(Envelope(-2, -1.5, 0.5, 1.5));
            tests.push_back(Envelope(0.5, 1.5, 1.5, 2));

            for (size_t i = 0; i < tests.size(); ++i)
            {
                if (e1.contain(tests[i]) != (i == 0))
                {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], (i != 0));
                }
                if (tests[i].contain(e1) == true)
                {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, true);
                }
            }
            cout << "Envelope Contain: " << tests.size() * 2 - failedCase << " / "
                 << tests.size() * 2 << " tests are passed" << endl;

            failedCase = 0;
            for (size_t i = 0; i < tests.size(); ++i)
            {
                if (e1.intersect(tests[i]) != (i < 6))
                {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], (i < 6));
                }
                if (tests[i].intersect(e1) != (i < 6))
                {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, (i < 6));
                }
            }
            cout << "Envelope Intersect: " << tests.size() * 2 - failedCase << " / "
                 << tests.size() * 2 << " tests are passed" << endl;

            failedCase = 0;
            vector<Envelope> results;
            results.push_back(Envelope(-1, 1, -1, 1));
            results.push_back(Envelope(-1, 1, -1, 1.5));
            results.push_back(Envelope(-1, 1.5, -1, 1));
            results.push_back(Envelope(-1.5, 1, -1.5, 1));
            results.push_back(Envelope(-2, 1, -1, 1));
            results.push_back(Envelope(-1, 1.5, -1, 1.5));
            results.push_back(Envelope(-2, 1, -1, 1));
            results.push_back(Envelope(-1, 1, -1, 2));
            results.push_back(Envelope(-2, 1, -1, 1.5));
            results.push_back(Envelope(-1, 1.5, -1, 2));
            for (size_t i = 0; i < tests.size(); ++i)
            {
                if (e1.unionEnvelope(tests[i]) != results[i])
                {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], e1.unionEnvelope(tests[i]),
                                 results[i]);
                }
                if (tests[i].unionEnvelope(e1) != results[i])
                {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, e1.unionEnvelope(tests[i]),
                                 results[i]);
                }
            }
            cout << "Envelope Union: " << tests.size() * 2 - failedCase << " / "
                 << tests.size() * 2 << " tests are passed" << endl;
        }
        else if (t == TEST2)
        {
            cout << "????2: Distance between Point and LineString" << endl;

            vector<Point> points;
            points.push_back(Point(0, 0));
            points.push_back(Point(10, 10));
            LineString line(points);

            points.push_back(Point(-10, -10));
            points.push_back(Point(20, 20));
            points.push_back(Point(5, 5));
            points.push_back(Point(10, 0));
            points.push_back(Point(10, -10));
            points.push_back(Point(0, 10));
            points.push_back(Point(0, 20));
            points.push_back(Point(20, 0));

            double dists[] = {0, 0, 14.1421, 14.1421, 0,
                              7.07107, 14.1421, 7.07107, 14.1421, 14.1421};

            int failedCase = 0;
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dist = points[i].distance(&line);
                if (fabs(dist - dists[i]) > 0.0001)
                {
                    failedCase += 1;
                    cout << "Your answer is " << dist << " for test between ";
                    line.print();
                    cout << " and ";
                    points[i].print();
                    cout << ", but the answer is " << dists[i] << endl;
                }
            }
            cout << "Distance between Point and LineString: "
                 << points.size() - failedCase << " / " << points.size()
                 << " tests are passed" << endl;
        }
        else if (t == TEST3)
        {
            cout << "????3: Distance between Point and Polygon" << endl;

            vector<Point> points;
            points.push_back(Point(5, 0));
            points.push_back(Point(3, 6));
            points.push_back(Point(2, 4));
            points.push_back(Point(-2, 4));
            points.push_back(Point(-3, 5));
            points.push_back(Point(-5, 0));
            points.push_back(Point(0, -3));
            points.push_back(Point(5, 0));
            LineString line(points);
            Polygon poly(line);

            points.clear();
            points.push_back(Point(5, 4));
            points.push_back(Point(3, 4));
            points.push_back(Point(0, 4));
            points.push_back(Point(-3, 4));
            points.push_back(Point(-5, 4));
            points.push_back(Point(5, 5));
            points.push_back(Point(3, 5));
            points.push_back(Point(0, 5));
            points.push_back(Point(-3, 5));
            points.push_back(Point(0, 0));

            double dists[] = {1.26491, 0, 0, 0, 1.48556, 1.58114, 0, 1, 0, 0};

            int failedCase = 0;
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dist = points[i].distance(&poly);
                if (fabs(dist - dists[i]) > 0.00001)
                {
                    failedCase += 1;
                    cout << "Your answer is " << dist << " for test between ";
                    poly.print();
                    cout << " and ";
                    points[i].print();
                    cout << ", but the answer is " << dists[i] << endl;
                }
            }
            cout << "Distance between Point and Polygon: "
                 << points.size() - failedCase << " / " << points.size()
                 << " tests are passed" << endl;
        }
        else if (t == TEST4)
        {
            cout << "????4: RTree Construction" << endl;
            int ncase, cct;
            ncase = cct = 2;

            RTree<8> rtree;
            vector<Geometry *> geom = readGeom("../data/station");
            vector<Feature> features;

            for (size_t i = 0; i < geom.size(); ++i)
                features.push_back(Feature("", geom[i]));

            rtree.constructTree(features);

            int height, interiorNum, leafNum;
            rtree.countHeight(height);
            rtree.countNode(interiorNum, leafNum);

            if (!(height == 4 && interiorNum == 18 && leafNum == 86))
            {
                cout << "Case 1: "
                     << "Your answer is height: " << height
                     << ", interiorNum: " << interiorNum << ", leafNum: " << leafNum
                     << ". One possible answer is height: 4, interiorNum: "
                        "18, "
                        "leafNum: 86\n";
                --cct;
            }

            features.clear();
            for (size_t i = 0; i < geom.size(); ++i)
                delete geom[i];
            geom.clear();

            vector<Geometry *> geom2 = readGeom("../data/highway");
            vector<Feature> features2;
            RTree<8> rtree2;

            for (size_t i = 0; i < geom2.size(); ++i)
                features2.push_back(Feature("", geom2[i]));

            rtree2.constructTree(features2);

            int height2, interiorNum2, leafNum2;
            rtree2.countHeight(height2);
            rtree2.countNode(interiorNum2, leafNum2);

            if (!(height2 == 6 && interiorNum2 == 484 && leafNum2 == 2305))
            {
                cout << "Case 2: "
                     << "Your answer is height: " << height2
                     << ", interiorNum: " << interiorNum2
                     << ", leafNum: " << leafNum2
                     << ". One possible answer is height: 6, interiorNum: "
                        "484, leafNum: 2305\n";
                --cct;
            }

            features2.clear();
            for (size_t i = 0; i < geom2.size(); ++i)
                delete geom2[i];
            geom2.clear();

            // cout << "RTree Construction: " << cct << " / " << ncase
            //      << " tests are passed" << endl;
        }
        else if (t == TEST5)
        {
            cout << "????5: Distance Between Point And Polygon With Inner Rings" << endl;

            vector<Point> points;
            points.push_back(Point(5, 0));
            points.push_back(Point(3, 6));
            points.push_back(Point(2, 4));
            points.push_back(Point(-2, 4));
            points.push_back(Point(-3, 5));
            points.push_back(Point(-5, 0));
            points.push_back(Point(0, -3));
            points.push_back(Point(5, 0));
            LineString ex(points);
            vector<LineString> in;
            points.clear();
            points.push_back(Point(3, 1));
            points.push_back(Point(3, 3));
            points.push_back(Point(1, 3));
            points.push_back(Point(1, 1));
            points.push_back(Point(3, 1));
            in.push_back(LineString(points));
            points.clear();
            points.push_back(Point(-1, -1));
            points.push_back(Point(-1, 2));
            points.push_back(Point(-2, 2));
            points.push_back(Point(-2, -1));
            points.push_back(Point(-1, -1));
            in.push_back(LineString(points));
            Polygon poly(ex, in);

            points.clear();
            points.push_back(Point(5, 4));
            points.push_back(Point(3, 4));
            points.push_back(Point(0, 4));
            points.push_back(Point(-3, 4));
            points.push_back(Point(-5, 4));
            points.push_back(Point(5, 5));
            points.push_back(Point(3, 5));
            points.push_back(Point(0, 5));
            points.push_back(Point(-3, 5));
            points.push_back(Point(0, 0));
            points.push_back(Point(3, 3));
            points.push_back(Point(2, 2));
            points.push_back(Point(-1.5, 1));
            points.push_back(Point(5, 2));
            points.push_back(Point(-5, 2));

            double dists[] = {1.26491, 0, 0, 0, 1.48556, 1.58114, 0, 1, 0, 0, 0, 1, 0.5, 0.63246, 0.74278};

            int failedCase = 0;
            for (size_t i = 0; i < points.size(); ++i)
            {
                double dist = points[i].distance(&poly);
                if (fabs(dist - dists[i]) > 0.00001)
                {
                    failedCase += 1;
                    cout << "Your answer is " << dist << " for test between ";
                    poly.print();
                    cout << " and ";
                    points[i].print();
                    cout << ", but the answer is " << dists[i] << endl;
                }
            }
            cout << "Distance Between Point And Polygon With Inner Rings: "
                 << points.size() - failedCase << " / " << points.size()
                 << " tests are passed" << endl;
        }
        else if (t == TEST8)
        {
            cout << "????8: RTreeAnalysis" << endl;
            analyse();
        }

        cout << "**********************End**********************" << endl;
    }

    template <uint8_t I, uint8_t Last, uint8_t Step>
    void forConstCapAnalyseRTree(const std::vector<Feature> &features)
    {
        if constexpr (I <= Last)
        {
            RTree<I> rtree;
            rtree.constructTree(features);

            clock_t start_time = clock();
            rtree.constructTree(features);
            clock_t end_time = clock();

            int height = 0, interiorNum = 0, leafNum = 0;
            rtree.countHeight(height);
            rtree.countNode(interiorNum, leafNum);

            std::cout << "M " << (size_t)I << "\n";
            std::cout << "Height: " << height
                      << " \tInterior node number: " << interiorNum
                      << " \tLeaf node number: " << leafNum << "\n";
            std::cout << "Construction time: " << (end_time - start_time) / 1000.0 / 1000.0 << "s"
                      << std::endl;

            double x, y;

            start_time = clock();
            for (int i = 0; i < 100000; ++i)
            {
                x = -((rand() % 225) / 10000.0 + 73.9812);
                y = (rand() % 239) / 10000.0 + 40.7247;
                std::vector<Feature> candidateFeatures;
                rtree.NNQuery(x, y, candidateFeatures);

                std::unique_ptr<Point> p(new Point(x, y));
                double distance = candidateFeatures[0].getGeom()->distance(p.get());
                Feature nearestFeature = candidateFeatures[0];
                for (size_t i = 1; i < candidateFeatures.size(); ++i)
                {
                    if (candidateFeatures[i].getGeom()->distance(p.get()) < distance)
                    {
                        distance = candidateFeatures[i].getGeom()->distance(p.get());
                        nearestFeature = candidateFeatures[i];
                    }
                }
            }
            end_time = clock();
            std::cout << "NNQuery time: " << (end_time - start_time) / 1000.0 / 1000.0 << "s"
                      << std::endl
                      << std::endl;

            forConstCapAnalyseRTree<I + Step, Last, Step>(features);
        }
    }

    template <uint8_t M>
    void RTree<M>::analyse()
    {
        using namespace std;

        vector<Feature> features;
        vector<Geometry *> geom = readGeom("../data/taxi");
        vector<string> name = readName("../data/taxi");

        features.clear();
        features.reserve(geom.size());
        for (size_t i = 0; i < geom.size(); ++i)
            features.push_back(Feature(name[i], geom[i]));

        cout << "taxi number: " << geom.size() << endl;

        srand(time(nullptr));

        forConstCapAnalyseRTree<70, 200, 10>(features);
    }

} // namespace hw6

#endif // !RTREE_TEST_H_INCLUDED