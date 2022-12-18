#include "QuadTree.h"

#include "float.h"

#include <set>
#include <map>
#include <memory>

namespace hw6
{

    /*
     * QuadNode
     */
    void QuadNode::split(size_t capacity)
    {
        for (int i = 0; i < 4; ++i)
        {
            delete children[i];
            children[i] = nullptr;
        }
        if (features.size() > capacity)
        {
            double midx = (bbox.getMinX() + bbox.getMaxX()) / 2.0;
            double midy = (bbox.getMinY() + bbox.getMaxY()) / 2.0;
            Envelope e0(bbox.getMinX(), midx, midy, bbox.getMaxY());
            children[0] = new QuadNode(e0);
            Envelope e1(midx, bbox.getMaxX(), midy, bbox.getMaxY());
            children[1] = new QuadNode(e1);
            Envelope e2(midx, bbox.getMaxX(), bbox.getMinY(), midy);
            children[2] = new QuadNode(e2);
            Envelope e3(bbox.getMinX(), midx, bbox.getMinY(), midy);
            children[3] = new QuadNode(e3);

            for (Feature f : features)
            {
                for (int i = 0; i < 4; ++i)
                {
                    if (children[i]->getEnvelope().intersect(f.getEnvelope()))
                    {
                        children[i]->add(f);
                    }
                }
            }

            for (int i = 0; i < 4; ++i)
            {
                children[i]->split(capacity);
            }
            features.clear();
        }
    }

    void QuadNode::countNode(int &interiorNum, int &leafNum)
    {
        if (isLeafNode())
        {
            ++leafNum;
        }
        else
        {
            ++interiorNum;
            for (int i = 0; i < 4; ++i)
                children[i]->countNode(interiorNum, leafNum);
        }
    }

    int QuadNode::countHeight(int height)
    {
        ++height;
        if (!isLeafNode())
        {
            int cur = height;
            for (int i = 0; i < 4; ++i)
            {
                height = std::max(height, children[i]->countHeight(cur));
            }
        }
        return height;
    }

    void QuadNode::rangeQuery(const Envelope &rect,
                              std::vector<Feature> &features)
    {
        if (!bbox.intersect(rect))
            return;

        // Task range query
        if (isLeafNode())
        {
            for (auto f : this->features)
            {
                if (f.getEnvelope().intersect(rect))
                {
                    features.push_back(f);
                }
            }
        }
        else
        {
            for (int i = 0; i < 4; ++i)
            {
                children[i]->rangeQuery(rect, features);
            }
        }
    }

    QuadNode *QuadNode::pointInLeafNode(double x, double y)
    {
        if (isLeafNode())
        {
            return this;
        }
        else
        {
            for (int i = 0; i < 4; ++i)
            {
                if (children[i]->bbox.contain(x, y))
                {
                    return children[i]->pointInLeafNode(x, y);
                }
            }
        }
        return nullptr;
    }

    void QuadNode::draw()
    {
        if (isLeafNode())
        {
            bbox.draw();
        }
        else
        {
            for (int i = 0; i < 4; ++i)
                children[i]->draw();
        }
    }

    /*
     * QuadTree
     */
    bool QuadTree::constructTree(const std::vector<Feature> &features)
    {
        if (features.empty())
            return false;
            
        Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
        for (Feature f : features)
        {
            e = e.unionEnvelope(f.getEnvelope());
        }
        root = new QuadNode(e);
        root->add(features);
        root->split(capacity);
        bbox = e;
        return true;
    }

    void QuadTree::countNode(int &interiorNum, int &leafNum)
    {
        interiorNum = 0;
        leafNum = 0;
        if (root)
            root->countNode(interiorNum, leafNum);
    }

    void QuadTree::countHeight(int &height)
    {
        height = 0;
        if (root)
            height = root->countHeight(0);
    }

    void QuadTree::rangeQuery(const Envelope &rect,
                              std::vector<Feature> &features)
    {
        features.clear();
        // Task range query
        std::vector<Feature> features1;
        root->rangeQuery(rect, features1);
        std::set<Feature> featureMap;
        for (auto f : features1)
        {
            if (featureMap.find(f) == featureMap.end())
            {
                if (f.getGeom()->intersects(rect))
                    features.push_back(f);
                featureMap.insert(f);
            }
        }
    }

    bool QuadTree::NNQuery(double x, double y, std::vector<Feature> &features)
    {
        if (!root || !(root->getEnvelope().contain(x, y)))
            return false;
        const Envelope &envelope = root->getEnvelope();
        double minDist = std::max(envelope.getWidth(), envelope.getHeight());
        QuadNode *n = root->pointInLeafNode(x, y);
        for (int i = 0; i < n->getFeatureNum(); ++i)
        {
            minDist = std::min(minDist, n->getFeature(i).maxDistance2Envelope(x, y));
        }
        Envelope rect = Envelope(x - minDist, x + minDist, y - minDist, y + minDist);
        std::vector<Feature> feature_tt;
        rangeQuery(rect, feature_tt);
        // refine step
        double dist;
        for (Feature f : feature_tt)
        {
            std::unique_ptr<Point> p(new Point(x, y));
            dist = f.getGeom()->distance(p.get());
            if (dist <= minDist)
            {
                features.push_back(f);
                minDist = dist;
            }
        }
        return true;
    }

    std::vector<std::pair<Feature, Feature>> QuadTree::spatialJoin(Tree *tree, double distance)
    {
        QuadTree *qtree = dynamic_cast<QuadTree *>(tree); // father class pointer to child class pointer
        std::vector<std::pair<Feature, Feature>> f;
        f.clear();
        if (!qtree)
            return f;
        root->spatialJoin_new(distance, qtree->getRoot(), f);
        return f;
    }

    void QuadNode::spatialJoin_new(double distance, QuadNode *node, std::vector<std::pair<Feature, Feature>> &joinSet)
    {
        if (isLeafNode())
        {
            for (size_t i = 0; i < getFeatureNum(); ++i)
            {
                std::vector<Feature> joinFeature;

                Envelope rect(getFeature(i).getEnvelope().getMinX() - distance, getFeature(i).getEnvelope().getMaxX() + distance, getFeature(i).getEnvelope().getMinY() - distance, getFeature(i).getEnvelope().getMaxY() + distance);

                node->rangeQuery(rect, joinFeature);

                for (size_t j = 0; j < joinFeature.size(); ++j)
                    joinSet.push_back(std::make_pair(getFeature(i), joinFeature[j]));
            }
        }
        else
            for (size_t i = 0; i < 4; i++)
                getChild(i)->spatialJoin_new(distance, node, joinSet);
    }

    void QuadTree::draw()
    {
        if (root)
            root->draw();
    }

} // namespace hw6
