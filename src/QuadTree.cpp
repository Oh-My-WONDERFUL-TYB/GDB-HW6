#include "QuadTree.h"
#include <set>

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

        // Task construction
        // TODO
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
        // TODO
    }

    QuadNode *QuadNode::pointInLeafNode(double x, double y)
    {
        // Task NN query
        // TODO

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

        // Task construction
        // TODO

        bbox = Envelope(-74.1, -73.8, 40.6, 40.8); // ע����д�����Ҫ����Ϊfeatures�İ�Χ�У�����ڵ�İ�Χ��

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
        // TODO

        // filter step (ѡ���ѯ�����뼸�ζ����Χ���ཻ�ļ��ζ���)

        // ע���Ĳ��������ѯ�����غ�ѡ������������hw6��rangeQuery�����
    }

    bool QuadTree::NNQuery(double x, double y, std::vector<Feature> &features)
    {
        if (!root || !(root->getEnvelope().contain(x, y)))
            return false;

        // Task NN query
        // TODO

        // filter step
        // (ʹ��maxDistance2Envelope��������ò�ѯ�㵽���ζ����Χ�е���̵������룬Ȼ�������ѯ��ú�ѡ��)

        const Envelope &envelope = root->getEnvelope();
        double minDist = std::max(envelope.getWidth(), envelope.getHeight());

        // ע���Ĳ����ڽ���ѯ�����غ�ѡ������������hw6��NNQuery�����

        return true;
    }

    std::vector<std::pair<Feature, Feature>> QuadTree::spatialJoin(Tree *tree, double distance)
    {
        QuadTree *qtree = dynamic_cast<QuadTree *>(tree); // father class pointer to child class pointer
        std::vector<std::pair<Feature, Feature>> f;

        // TODO

        // 返回一个vector，这个vector的元素是一个Feature对，这个对满足两元素距离可能小于distance
        // 判断方法是对this的每一个根节点建立一个Envelope——
        // Envelope rect(getFeature(i)->getEnvelope().getMinX() - distance, getFeature(i)->getEnvelope().getMaxX() + distance, getFeature(i)->getEnvelope().getMinY() - distance, getFeature(i)->getEnvelope().getMaxY() + distance);
        // 然后用这个rect去和qtree做rangequery，返回的features和做查询的feature形成一个pair然后返回
        // 可以看r树的实现和最终的输出

        return f;
    }

    void QuadTree::draw()
    {
        if (root)
            root->draw();
    }

} // namespace hw6
