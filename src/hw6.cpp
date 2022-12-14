/***
 * @Description:main entry
 * @Author: jwimd chenjiewei@zju.edu.cn
 * @Date: 2022-11-11 17:50:56
 * @LastEditors: jwimd chenjiewei@zju.edu.cn
 * @LastEditTime: 2022-11-11 18:38:26
 * @FilePath: /GDB-HW6/src/hw6.cpp
 */

//

#include "Common.h"
#include "Geometry.h"
#include "shapelib/shapefil.h"

#include <GL/freeglut.h>

#include <cstdio>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>
#include <iomanip>

#ifdef USE_RTREE
#include "RTree.h"
using TreeTy = hw6::RTree<8>;
#else
#include "QuadTree.h"
using TreeTy = hw6::QuadTree;
#endif

using namespace std;

int screenWidth = 640;
int screenHeight = 480;

double pointSize = 2.0;

int mode;

vector<hw6::Feature> features;
vector<hw6::Feature> roads;
vector<hw6::Feature> polygons;
bool showRoad = true;
bool showChina = false;

unique_ptr<hw6::Tree> pointTree;
unique_ptr<hw6::Tree> roadTree;
unique_ptr<hw6::Tree> polygonTree;
bool showTree = false;

hw6::Feature nearestFeature;

bool firstPoint = true;
hw6::Point corner[2];
hw6::Envelope selectedRect;
vector<hw6::Feature> selectedFeatures;

/*
 * shapefile
 */
vector<string> readName(const char *filename)
{
    DBFHandle file = DBFOpen(filename, "r");

    vector<string> res;
    int cct = DBFGetRecordCount(file);
    res.reserve(cct);
    for (int i = 0; i < cct; ++i)
    {
        string a = DBFReadStringAttribute(file, i, 0);
        res.push_back(a);
    }

    DBFClose(file);

    return res;
}

vector<hw6::Geometry *> readGeom(const char *filename)
{
    SHPHandle file = SHPOpen(filename, "r");

    int pnEntities, pnShapeType;
    double padfMinBound[4], padfMaxBound[4];
    SHPGetInfo(file, &pnEntities, &pnShapeType, padfMinBound, padfMaxBound);

    vector<hw6::Point> points;
    vector<hw6::Geometry *> geoms;
    geoms.reserve(pnEntities);
    switch (pnShapeType)
    {
    case SHPT_POINT:
        for (int i = 0; i < pnEntities; ++i)
        {
            SHPObject *pt = SHPReadObject(file, i);
            geoms.push_back(new hw6::Point(pt->padfY[0], pt->padfX[0]));
            SHPDestroyObject(pt);
        }
        break;

    case SHPT_ARC:
        for (int i = 0; i < pnEntities; ++i)
        {
            points.clear();
            SHPObject *pt = SHPReadObject(file, i);
            for (int j = 0; j < pt->nVertices; ++j)
            {
                points.push_back(hw6::Point(pt->padfY[j], pt->padfX[j]));
            }
            SHPDestroyObject(pt);
            geoms.push_back(new hw6::LineString(points));
        }
        break;

    case SHPT_POLYGON:
        for (int i = 0; i < pnEntities; ++i)
        {
            points.clear();
            SHPObject *pt = SHPReadObject(file, i);
            for (int j = 0; j < pt->nVertices; ++j)
            {
                if (!showChina)
                    points.push_back(hw6::Point(pt->padfY[j], pt->padfX[j]));
                else
                    points.push_back(hw6::Point(pt->padfX[j], pt->padfY[j]));
            }
            SHPDestroyObject(pt);
            hw6::LineString line(points);
            geoms.push_back(new hw6::Polygon(line));
        }
        break;
    }

    SHPClose(file);
    return geoms;
}

void transformValue(double &res, const char *format = "%.2lf")
{
    char buf[20];
    sprintf(buf, format, res);
    sscanf(buf, "%lf", &res);
}

void wrongMessage(hw6::Envelope e1, hw6::Envelope e2, bool cal)
{
    cout << "Your answer is " << cal << " for test between ";
    e1.print();
    cout << " and ";
    e2.print();
    cout << ", but the answer is " << !cal << endl;
}

void wrongMessage(const hw6::Point &pt1, const hw6::Point &pt2, double dis,
                  double res)
{
    cout << "Your answer is " << dis << " for test ";
    pt1.print();
    cout << " and ";
    pt2.print();
    cout << ", but the answer is " << res << endl;
}

void wrongMessage(hw6::Envelope e1, hw6::Envelope e2, hw6::Envelope cal,
                  hw6::Envelope res)
{
    cout << "Your answer is ";
    cal.print();
    cout << " for test between ";
    e1.print();
    cout << " and ";
    e2.print();
    cout << ", but the answer is ";
    res.print();
    cout << endl;
}

/*
 * print a geom
 */
void printGeom(vector<hw6::Geometry *> &geom)
{
    cout << "Geometry:" << endl;
    for (vector<hw6::Geometry *>::iterator it = geom.begin(); it != geom.end();
         ++it)
    {
        (*it)->print();
    }
}

/*
 * delete a geom
 */
void deleteGeom(vector<hw6::Geometry *> &geom)
{
    for (vector<hw6::Geometry *>::iterator it = geom.begin(); it != geom.end();
         ++it)
    {
        delete *it;
        *it = NULL;
    }
    geom.clear();
}

/*
 * load road data
 */
void loadRoadData()
{
    vector<hw6::Geometry *> geom = readGeom("../data/highway");

    roads.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        roads.push_back(hw6::Feature(to_string(i), geom[i]));

    cout << "road number: " << geom.size() << endl;
    roadTree->setCapacity(20);
    roadTree->constructTree(roads);
}

/*
 * load station data
 */
void loadStationData()
{
    vector<hw6::Geometry *> geom = readGeom("../data/station");
    vector<string> name = readName("../data/station");

    features.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        features.push_back(hw6::Feature(name[i], geom[i]));

    cout << "station number: " << geom.size() << endl;
    pointTree->setCapacity(5);
    pointTree->constructTree(features);
}

/*
 * load taxi data
 */
void loadTaxiData()
{
    vector<hw6::Geometry *> geom = readGeom("../data/taxi");
    vector<string> name = readName("../data/taxi");

    features.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        features.push_back(hw6::Feature(name[i], geom[i]));

    cout << "taxi number: " << geom.size() << endl;
    pointTree->setCapacity(100);
    pointTree->constructTree(features);
}

void loadPolygonData()
{
    vector<hw6::Geometry *> geom = readGeom("../data/polygon");

    polygons.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        polygons.push_back(hw6::Feature(to_string(i), geom[i]));

    cout << "polygon number: " << geom.size() << endl;
    polygonTree->setCapacity(5);
    polygonTree->constructTree(polygons);
}

void loadChinaData()
{
    vector<hw6::Geometry *> geom = readGeom("../data/CN_city");
    vector<string> name = readName("../data/CN_city");

    polygons.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        polygons.push_back(hw6::Feature(name[i], geom[i]));

    cout << "China city number: " << geom.size() << endl;
    polygonTree->setCapacity(5);
    polygonTree->constructTree(polygons);
}

/*
 * range query
 */
void rangeQuery()
{
    vector<hw6::Feature> candidateFeatures;

    // filter step get evenlope intersection
    if (mode == RANGEPOINT)
        pointTree->rangeQuery(selectedRect, candidateFeatures);
    else if (mode == RANGELINE)
        roadTree->rangeQuery(selectedRect, candidateFeatures);
    else if (mode == RANGEPOLYGON)
        polygonTree->rangeQuery(selectedRect, candidateFeatures);

    // refine step define real intersection
    for (size_t i = 0; i < candidateFeatures.size(); ++i)
    {
        if (selectedRect.contain(candidateFeatures[i].getEnvelope()))
            selectedFeatures.push_back(candidateFeatures[i]);
        else if (candidateFeatures[i].getGeom()->intersects(selectedRect))
            selectedFeatures.push_back(candidateFeatures[i]);
    }
}

/*
 * nn Query
 */
void NNQuery(hw6::Point p)
{
    vector<hw6::Feature> candidateFeatures;

    // filter step get candidate use Tree
    if (mode == NNPOINT)
        pointTree->NNQuery(p.getX(), p.getY(), candidateFeatures);
    else if (mode == NNLINE)
        roadTree->NNQuery(p.getX(), p.getY(), candidateFeatures);
    else if (mode == NNPOLYGON)
        polygonTree->NNQuery(p.getX(), p.getY(), candidateFeatures);

    // refine step get real feature
    double distance = candidateFeatures[0].getGeom()->distance(&p);
    nearestFeature = candidateFeatures[0];
    for (size_t i = 1; i < candidateFeatures.size(); ++i)
    {
        if (candidateFeatures[i].getGeom()->distance(&p) < distance)
        {
            distance = candidateFeatures[i].getGeom()->distance(&p);
            nearestFeature = candidateFeatures[i];
        }
    }
}

/***
 * @Description: Spatial Join Between road and station or taxi
 * @Author: jwimd chenjiewei@zju.edu.cn
 * @msg: None
 * @param {double} distance
 * @return {*}
 */
void spatialJoin(double distance)
{
    std::vector<std::pair<Feature, Feature>> joinSet;
    joinSet.clear();

    joinSet = pointTree->spatialJoin(roadTree.get(), distance);

    std::vector<std::pair<Feature, Feature>> resultSet;

    for (size_t i = 0; i < joinSet.size(); ++i)
    {
        const LineString *line = dynamic_cast<const LineString *>(joinSet[i].second.getGeom());
        if (joinSet[i].first.getGeom()->distance(line) <= distance)
            resultSet.push_back(joinSet[i]);
    }
    if (resultSet.empty())
    {
        cout << "empty set" << endl;
        return;
    }

    std::vector<std::string> point, road;
    point.push_back("point");
    road.push_back("road");

    int listLen1 = 6, listLen2 = 5;

    for (size_t i = 0; i < resultSet.size(); ++i)
    {
        point.push_back(resultSet[i].first.getName());
        road.push_back(resultSet[i].second.getName());

        if (resultSet[i].first.getName().length() + 1 > listLen1)
            listLen1 = resultSet[i].first.getName().length() + 1;
        if (resultSet[i].second.getName().length() + 1 > listLen2)
            listLen2 = resultSet[i].second.getName().length() + 1;
    }

    int totalLen = listLen1 + listLen2 + 3;

    cout << left << setw(totalLen) << setfill('-') << '-' << endl;

    for (size_t i = 0; i < resultSet.size() + 1; ++i)
    {
        cout << setw(1) << '|';
        cout << setw(listLen1) << setfill(' ') << point[i];
        cout << setw(1) << '|';
        cout << setw(listLen2) << setfill(' ') << road[i];
        cout << setw(1) << '|' << endl;
        if (i == 0)
            cout << setw(totalLen) << setfill('-') << '-' << endl;
    }

    cout << setw(totalLen) << setfill('-') << '-' << endl;
    cout << "index return " << joinSet.size() << " row" << endl;
    cout << "final return " << resultSet.size() << " row" << endl;
}

/*
 * screen point to geom point
 */
void transfromPt(hw6::Point &pt)
{
    hw6::Envelope bbox = roadTree->getEnvelope();
    double width = bbox.getMaxX() - bbox.getMinX() + 0.002;
    double height = bbox.getMaxY() - bbox.getMinY() + 0.002;

    double x = pt.getX() * width / screenWidth + bbox.getMinX() - 0.001;
    double y = pt.getY() * height / screenHeight + bbox.getMinY() - 0.001;

    if (mode == POLYGON || mode == RANGEPOLYGON || mode == NNPOLYGON)
    {
        bbox = polygonTree->getEnvelope();
        width = bbox.getMaxX() - bbox.getMinX() + 20;
        height = bbox.getMaxY() - bbox.getMinY() + 20;

        x = pt.getX() * width / screenWidth + bbox.getMinX() - 10;
        y = pt.getY() * height / screenHeight + bbox.getMinY() - 10;
    }

    x = max(bbox.getMinX(), x);
    x = min(bbox.getMaxX(), x);
    y = max(bbox.getMinY(), y);
    y = min(bbox.getMaxY(), y);
    pt = hw6::Point(x, y);
}

/*
 * display geom
 */
void display()
{
    // glClearColor(241 / 255.0, 238 / 255.0, 232 / 255.0, 0.0);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    hw6::Envelope bbox;

    if (mode == POLYGON || mode == RANGEPOLYGON || mode == NNPOLYGON)
    {
        bbox = polygonTree->getEnvelope();
        gluOrtho2D(bbox.getMinX() - 10, bbox.getMaxX() + 10,
                   bbox.getMinY() - 10, bbox.getMaxY() + 10);
    }
    else
    {
        bbox = roadTree->getEnvelope();
        gluOrtho2D(bbox.getMinX() - 0.001, bbox.getMaxX() + 0.001,
                   bbox.getMinY() - 0.001, bbox.getMaxY() + 0.001);
    }

    // show road
    if (showRoad)
    {
        glColor3d(252 / 255.0, 214 / 255.0, 164 / 255.0);
        for (size_t i = 0; i < roads.size(); ++i)
            roads[i].draw();
    }

    //  show point
    if (!(mode == RANGELINE || !mode == NNLINE) && !(mode == POLYGON || mode == RANGEPOLYGON || mode == NNPOLYGON))
    {
        glPointSize((float)pointSize);
        glColor3d(0.0, 146 / 255.0, 247 / 255.0);
        for (size_t i = 0; i < features.size(); ++i)
            features[i].draw();
    }

    // show polygon
    if (mode == POLYGON || mode == RANGEPOLYGON || mode == NNPOLYGON)
    {
        glColor3d(123 / 255.0, 200 / 255.0, 200 / 255.0);
        for (size_t i = 0; i < polygons.size(); ++i)
            polygons[i].draw();
    }

    // show tree
    if (showTree)
    {
        glColor3d(0.0, 146 / 255.0, 247 / 255.0);
        if (mode == RANGELINE || mode == NNLINE)
            roadTree->draw();
        else if (mode == POLYGON || mode == RANGEPOLYGON || mode == NNPOLYGON)
            polygonTree->draw();
        else
            pointTree->draw();
    }

    // show nnPoint
    if (mode == NNPOINT)
    {
        glPointSize(5.0);
        glColor3d(0.9, 0.0, 0.0);
        nearestFeature.draw();
    }

    // show nnLine or nnPolygon
    if (mode == NNLINE || mode == NNPOLYGON)
    {
        glLineWidth(3.0);
        glColor3d(0.9, 0.0, 0.0);
        nearestFeature.draw();
        glLineWidth(1.0);
    }

    // show range query result
    if (mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLYGON)
    {
        glColor3d(0.0, 0.0, 0.0);
        selectedRect.draw();
        glColor3d(1.0, 0.0, 0.0);
        for (size_t i = 0; i < selectedFeatures.size(); ++i)
            selectedFeatures[i].draw();
    }

    glFlush();
    glutSwapBuffers();
}

/*
 * mouse event
 */
void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        if (mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLYGON)
        {
            if (firstPoint)
            {
                selectedFeatures.clear();
                corner[0] = hw6::Point(x, screenHeight - y);
                transfromPt(corner[0]);
            }
            else
            {
                corner[1] = hw6::Point(x, screenHeight - y);
                transfromPt(corner[1]);
                selectedRect =
                    hw6::Envelope(min(corner[0].getX(), corner[1].getX()),
                                  max(corner[0].getX(), corner[1].getX()),
                                  min(corner[0].getY(), corner[1].getY()),
                                  max(corner[0].getY(), corner[1].getY()));
                rangeQuery();
            }
            firstPoint = !firstPoint;
            glutPostRedisplay();
        }
    }
}

void passiveMotion(int x, int y)
{
    corner[1] = hw6::Point(x, screenHeight - y);

    if ((mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLYGON) && !firstPoint)
    {
        corner[1] = hw6::Point(x, screenHeight - y);
        transfromPt(corner[1]);
        selectedRect = hw6::Envelope(min(corner[0].getX(), corner[1].getX()),
                                     max(corner[0].getX(), corner[1].getX()),
                                     min(corner[0].getY(), corner[1].getY()),
                                     max(corner[0].getY(), corner[1].getY()));
        rangeQuery();

        glutPostRedisplay();
    }
    else if (mode == NNPOINT || mode == NNLINE || mode == NNPOLYGON)
    {
        hw6::Point p(x, screenHeight - y);
        transfromPt(p);
        NNQuery(p);

        glutPostRedisplay();
    }
}

void changeSize(int w, int h)
{
    screenWidth = w;
    screenHeight = h;
    glViewport(0, 0, w, h);
    glutPostRedisplay();
}

void processNormalKeys(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27:
        exit(0);
        break;
    case 'N':
        mode = NNLINE;
        break;
    case 'n':
        mode = NNPOINT;
        break;
    case 'S':
        mode = RANGELINE;
        firstPoint = true;
        break;
    case 's':
        mode = RANGEPOINT;
        firstPoint = true;
        break;
    case 'J':
    case 'j':
        spatialJoin(0.00001);
    case 'B':
    case 'b':
        loadStationData();
        mode = Default;
        break;
    case 'T':
    case 't':
        loadTaxiData();
        mode = Default;
        break;
    case 'R':
    case 'r':
        showRoad = !showRoad;
        mode = Default;
        break;
    case 'Q':
    case 'q':
        showTree = !showTree;
        break;
    case 'P':
    case 'p':
        showRoad = false;
        showChina = false;
        loadPolygonData();
        mode = POLYGON;
        break;
    case 'C':
    case 'c':
        showRoad = false;
        showChina = true;
        loadChinaData();
        mode = POLYGON;
        break;
    case 'O':
    case 'o':
        mode = RANGEPOLYGON;
        firstPoint = true;
        break;
    case 'I':
    case 'i':
        mode = NNPOLYGON;
        break;
    case '+':
        pointSize *= 1.1;
        break;
    case '-':
        pointSize /= 1.1;
        break;
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
        TreeTy::test(key - '0');
        break;
    default:
        mode = Default;
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[])
{
    cout << "Key Usage:\n"
         << "  S  : range search for roads\n"
         << "  s  : range search for stations\n"
         << "  N  : nearest road search\n"
         << "  n  : nearest station search\n"
         << "  J/j: spatial join between road and station(or taxi), distance = 0.0001\n"
         << "  B/b: Bicycle data\n"
         << "  T/t: Taxi data\n"
         << "  R/r: show Road\n"
         << "  Q/q: show Tree\n"
         << "  P/p: show simple polygon \n"
         << "  C/c: show China city polygon \n"
         << "  O/o: polygon range query \n"
         << "  I/i: polygon nearest search \n"
         << "  +  : increase point size\n"
         << "  -  : decrease point size\n"
         << "  1  : Test Envelope contain, interset and union\n"
         << "  2  : Test distance between Point and LineString\n"
         << "  3  : Test distance between Point and Polygon\n"
         << "  4  : Test tree construction\n"
         << "  5  : Test inner rings inplement of Polygon\n"
         << "  6  : Test (your option here)\n"
         << "  8  : Tree performance analysis\n"
         << "  ESC: quit\n"
         << endl;

    pointTree = make_unique<TreeTy>();
    roadTree = make_unique<TreeTy>();
    polygonTree = make_unique<TreeTy>();

    loadRoadData();

    loadStationData();

    loadPolygonData();

    glutInit(&argc, argv);
    glutInitWindowSize(screenWidth, screenHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutCreateWindow("New York");

    glutMouseFunc(mouse);
    glutDisplayFunc(display);
    glutPassiveMotionFunc(passiveMotion);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(processNormalKeys);

    glutMainLoop();

    return 0;
}
