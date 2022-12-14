#ifndef RTREE_SRC_H_INCLUDED
#define RTREE_SRC_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"

#include <algorithm>
#include <array>
#include <queue>
#include <string>
#include <vector>
#include <memory>

namespace hw6
{

    /// <summary>
    /// </summary>
    /// <param name="M">Maximum number of children of each nodes.</param>
    template <uint8_t M>
    class RNode
    {
        constexpr static auto m = M >> 1;
        constexpr static auto M_minus_m = M - m;

    private:
        std::vector<const Feature *> features;
        RNode<M> *parent = nullptr;
        std::array<RNode<M> *, M> children = []()
        {
            decltype(getChildren()) ret;
            for (decltype(M) i = 0; i < M; ++i)
                ret[i] = nullptr;
            return ret;
        }();
        decltype(M) childrenNum = 0;
        Envelope bbox;

    public:
        RNode() = delete;
        RNode(const Envelope &box) : bbox(box) {}

        ~RNode()
        {
            for (decltype(M) i = 0; i < M; ++i)
                children[i] = nullptr;
            parent = nullptr;
            features.clear();
        }

        inline bool isLeafNode() { return childrenNum == 0; }

        inline const Envelope &getEnvelope() { return bbox; }

        inline RNode<M> *getParent() { return parent; }

        inline void setEnvelope(const Envelope &box) { bbox = box; }

        inline RNode<M> *getChild(size_t i)
        {
            return i < childrenNum ? children[i] : nullptr;
        }

        inline const RNode<M> *getChild(size_t i) const
        {
            return i < childrenNum ? children[i] : nullptr;
        }

        inline decltype(M) getChildNum() const { return childrenNum; }

        inline std::array<RNode<M> *, M> getChildren() const { return children; }

        inline size_t getFeatureNum() const
        {
            return features.size();
        }

        inline const Feature *getFeature(size_t i) const { return features[i]; }

        inline const std::vector<const Feature *> &getFeatures() const { return features; }

        inline void add(const Feature *f)
        {
            features.push_back(f);
        }

        inline void add(RNode<M> *child)
        {
            children[childrenNum] = child;
            child->parent = this;
            ++childrenNum;
        }

        void remove(const Feature *f)
        {

            for (auto itr = features.begin(); itr != features.end(); ++itr)
                if ((*itr)->getName() == f->getName())
                {
                    features.erase(itr);
                    break;
                }
            if (features.empty())
                features.shrink_to_fit(); // free memory unused but allocated
        }

        void remove(RNode<M> *child)
        {
            for (decltype(M) i = 0; i < childrenNum; ++i)
                if (children[i] == child)
                {
                    --childrenNum;
                    std::swap(children[i], children[childrenNum]);
                    children[childrenNum] = nullptr;
                    break;
                }
        }

        inline const Feature *popBackFeature()
        {
            const Feature *ret = features.back();
            features.pop_back();
            return ret;
        }

        inline RNode<M> *popBackChild()
        {
            --childrenNum;
            auto ret = children[childrenNum];
            children[childrenNum] = nullptr;
            return ret;
        }

        void countNode(int &interiorNum, int &leafNum)
        {
            if (isLeafNode())
            {
                ++leafNum;
            }
            else
            {
                ++interiorNum;
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    children[i]->countNode(interiorNum, leafNum);
            }
        }

        int countHeight(int height)
        {
            ++height;
            if (!isLeafNode())
            {
                int cur = height;
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    height = std::max(height, children[i]->countHeight(cur));
            }
            return height;
        }

        inline void draw()
        {
            if (isLeafNode())
            {
                bbox.draw();
            }
            else
            {
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    children[i]->draw();
            }
        }

        /***
         * @Description: Range Query for Node(feature's evenlope)
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {Envelope} &rect
         * @param {vector<Feature>} &features
         * @return {*}
         */
        void rangeQuery(const Envelope &rect, std::vector<const Feature *> &features)
        {
            if (!getEnvelope().intersect(rect))
                return;
            if (isLeafNode())
            {
                for (size_t i = 0; i < getFeatureNum(); ++i)
                {
                    if (getFeature(i)->getEnvelope().intersect(rect))
                        features.push_back(getFeature(i));
                }
            }
            else
            {
                for (size_t i = 0; i < getChildNum(); ++i)
                    getChild(i)->rangeQuery(rect, features);
            }
        }

        /***
         * @Description: Find leaf for point or nullptr
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {double} x
         * @param {double} y
         * @param {double} &distance
         * @return {*}
         */
        RNode<M> *pointInLeafNode(double x, double y)
        {
            if (isLeafNode())
                return this;
            else
            {
                for (size_t i = 0; i < getChildNum(); i++)
                    if (getChild(i)->getEnvelope().contain(x, y))
                        return getChild(i)->pointInLeafNode(x, y);
            }
            return nullptr;
        }

        /***
         * @Description: Outer Loop for join
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {double} distance
         * @param {Tree} *node
         * @return {*}
         */
        void spatialJoin(double distance, RNode<M> *node, std::vector<std::pair<Feature, Feature>> &joinSet)
        {
            if (isLeafNode())
            {
                for (size_t i = 0; i < getFeatureNum(); ++i)
                {
                    std::vector<const Feature *> joinFeature;
                    Envelope rect(getFeature(i)->getEnvelope().getMinX() - distance, getFeature(i)->getEnvelope().getMaxX() + distance, getFeature(i)->getEnvelope().getMinY() - distance, getFeature(i)->getEnvelope().getMaxY() + distance);

                    node->rangeQuery(rect, joinFeature);

                    for (size_t j = 0; j < joinFeature.size(); ++j)
                        joinSet.push_back(std::make_pair(*getFeature(i), *joinFeature[j]));
                }
            }
            else
                for (size_t i = 0; i < getChildNum(); i++)
                    getChild(i)->spatialJoin(distance, node, joinSet);
        }

        /***
         * @Description:  Quadratic Split For Rnode
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @return {*}
         */
        std::vector<RNode<M> *>
        split()
        {
            RNode<M> *L1 = nullptr;
            RNode<M> *L2 = nullptr;

            if (this->isLeafNode())
            {

                // Pick Seeds
                double max_area = -1;
                size_t seed1 = -1, seed2 = -1;
                for (size_t i = 0; i < this->getFeatureNum() - 1; i++)
                    for (size_t j = i + 1; j < this->getFeatureNum(); j++)
                    {
                        Envelope split_env(this->getFeature(i)->getEnvelope().unionEnvelope(this->getFeature(j)->getEnvelope()));

                        if (split_env.getArea() > max_area)
                        {
                            max_area = split_env.getArea();
                            seed1 = i;
                            seed2 = j;
                        }
                    }

                L1 = new RNode<M>(this->getFeature(seed1)->getEnvelope());
                L1->add(this->getFeature(seed1));

                L2 = new RNode<M>(this->getFeature(seed2)->getEnvelope());
                L2->add(this->getFeature(seed2));

                this->remove(this->getFeature(seed1));
                if (seed1 < seed2)
                    seed2--;
                this->remove(this->getFeature(seed2));

                // Pick Next

                while (this->getFeatureNum())
                {

                    Envelope env1 = L1->getEnvelope();
                    Envelope env2 = L2->getEnvelope();
                    if (L1->getFeatureNum() == size_t(m))
                    {
                        const Feature *_feature = popBackFeature();
                        L2->setEnvelope(env2.unionEnvelope(_feature->getEnvelope()));
                        L2->add(_feature);
                    }
                    else if (L2->getFeatureNum() == size_t(m))
                    {
                        const Feature *_feature = popBackFeature();
                        L1->setEnvelope(env1.unionEnvelope(_feature->getEnvelope()));
                        L1->add(_feature);
                    }
                    else
                    {
                        const Feature *_feature = nullptr;
                        size_t pick = -1;
                        double max_difference = -1;
                        int flag = -1;
                        for (size_t i = 0; i < this->getFeatureNum(); i++)
                        {
                            _feature = this->getFeature(i);
                            if (env1.unionEnvelope(_feature->getEnvelope()).getArea() - env2.unionEnvelope(_feature->getEnvelope()).getArea() > max_difference)
                            {
                                flag = 0;
                                max_difference = env1.unionEnvelope(_feature->getEnvelope()).getArea() - env2.unionEnvelope(_feature->getEnvelope()).getArea();
                                pick = i;
                            }
                            if (env2.unionEnvelope(_feature->getEnvelope()).getArea() - env1.unionEnvelope(_feature->getEnvelope()).getArea() > max_difference)
                            {
                                flag = 1;
                                max_difference = env2.unionEnvelope(_feature->getEnvelope()).getArea() - env1.unionEnvelope(_feature->getEnvelope()).getArea();
                                pick = i;
                            }
                        }
                        _feature = this->getFeature(pick);
                        if (flag == 1)
                        {
                            L1->setEnvelope(env1.unionEnvelope(_feature->getEnvelope()));
                            L1->add(_feature);
                        }
                        else if (flag == 0)
                        {
                            L2->setEnvelope(env2.unionEnvelope(_feature->getEnvelope()));
                            L2->add(_feature);
                        }
                        this->remove(_feature);
                    }
                }
            }

            else
            {

                // Pick Seeds
                double max_area = 0;
                size_t seed1 = -1, seed2 = -1;
                for (size_t i = 0; i < this->getChildNum() - 1; i++)
                    for (size_t j = i + 1; j < this->getChildNum(); j++)
                    {
                        Envelope split_env(this->getChild(i)->getEnvelope().unionEnvelope(this->getChild(j)->getEnvelope()));

                        if (split_env.getArea() > max_area)
                        {
                            max_area = split_env.getArea();
                            seed1 = i;
                            seed2 = j;
                        }
                    }

                L1 = new RNode<M>(this->getChild(seed1)->getEnvelope());
                L1->add(this->getChild(seed1));

                L2 = new RNode<M>(this->getChild(seed2)->getEnvelope());
                L2->add(this->getChild(seed2));

                if (seed2 == (size_t)M - 1)
                    seed2 = seed1;

                this->remove(this->getChild(seed1));

                this->remove(this->getChild(seed2));

                // Pick Next

                while (this->getChildNum())
                {

                    Envelope env1 = L1->getEnvelope();
                    Envelope env2 = L2->getEnvelope();
                    if (L1->getChildNum() == size_t(m))
                    {
                        RNode<M> *_child = popBackChild();
                        L2->setEnvelope(env2.unionEnvelope(_child->getEnvelope()));
                        L2->add(_child);
                    }
                    else if (L2->getChildNum() == size_t(m))
                    {
                        RNode<M> *_child = popBackChild();
                        L1->setEnvelope(env1.unionEnvelope(_child->getEnvelope()));
                        L1->add(_child);
                    }
                    else
                    {
                        RNode<M> *_child = nullptr;
                        size_t pick = -1;
                        double max_difference = -1;
                        int flag = -1;
                        for (size_t i = 0; i < this->getChildNum(); i++)
                        {
                            _child = this->getChild(i);
                            if (env1.unionEnvelope(_child->getEnvelope()).getArea() - env2.unionEnvelope(_child->getEnvelope()).getArea() > max_difference)
                            {
                                flag = 0;
                                max_difference = env1.unionEnvelope(_child->getEnvelope()).getArea() - env2.unionEnvelope(_child->getEnvelope()).getArea();
                                pick = i;
                            }
                            if (env2.unionEnvelope(_child->getEnvelope()).getArea() - env1.unionEnvelope(_child->getEnvelope()).getArea() > max_difference)
                            {
                                flag = 1;
                                max_difference = env2.unionEnvelope(_child->getEnvelope()).getArea() - env1.unionEnvelope(_child->getEnvelope()).getArea();
                                pick = i;
                            }
                        }
                        _child = this->getChild(pick);
                        if (flag == 1)
                        {
                            L1->setEnvelope(env1.unionEnvelope(_child->getEnvelope()));
                            L1->add(_child);
                        }
                        else if (flag == 0)
                        {
                            L2->setEnvelope(env2.unionEnvelope(_child->getEnvelope()));
                            L2->add(_child);
                        }
                        this->remove(_child);
                    }
                }
            }

            RNode<M> *P = this->getParent();

            // reset parent
            if (P)
            {
                P->remove(this);
                P->setEnvelope(L1->getEnvelope());

                for (size_t i = 0; i < P->getChildNum(); i++)
                    P->setEnvelope(P->getEnvelope().unionEnvelope(P->getChild(i)->getEnvelope()));

                RNode<M> *p = P->getParent();
                RNode<M> *l = P;

                while (p)
                {
                    p->setEnvelope(l->getEnvelope());
                    for (size_t i = 0; i < p->getChildNum(); i++)
                        p->setEnvelope(p->getEnvelope().unionEnvelope(p->getChild(i)->getEnvelope()));
                    l = l->getParent();
                    p = l->getParent();
                }

                P->add(L1);

                p = P->getParent();
                l = P;

                while (p)
                {
                    p->setEnvelope(p->getEnvelope().unionEnvelope(l->getEnvelope()));
                    l = l->getParent();
                    p = l->getParent();
                }
            }

            std::vector<RNode<M> *> ret;
            ret.clear();
            ret.push_back(L1);
            ret.push_back(L2);

            return ret;
        }
    };

    template <uint8_t M>
    class RTree : public Tree
    {
    private:
        // Vars here MAY be useful, but it's ok to delete them {
        constexpr static auto m = M >> 1;
        constexpr static auto M_minus_m = M - m;
        constexpr static double EQ_ERR = 0.0000000005;
        // }

        RNode<M> *root = nullptr;

    public:
        RTree() : Tree(M) { static_assert(M >= 4); }
        ~RTree()
        {
            if (root != nullptr)
                freeTree(root);
            root = nullptr;
        }

        RNode<M> *getRoot() { return root; }

        /***
         * @Description: Check Geom Rule
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {RNode<M>} *node
         * @return {*}
         */
        std::string
        checkGeom(RNode<M> *node)
        {
            if (node->isLeafNode())
            {
                for (size_t i = 0; i < node->getFeatureNum(); ++i)
                    if (!node->getEnvelope().contain(node->getFeature(i)->getEnvelope()))
                        return "leaf dont contain feature\n";
            }
            else
            {
                for (size_t i = 0; i < node->getChildNum(); ++i)
                    if (!node->getEnvelope().contain(node->getChild(i)->getEnvelope()))
                        return "inner dont contain child\n";
                for (size_t i = 0; i < node->getChildNum(); ++i)
                    checkGeom(node->getChild(i));
            }
            return "ok\n";
        }

        /***
         * @Description: Check RTree Structure
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {RNode<M>} *node
         * @return {*}
         */
        std::string checkTree(RNode<M> *node)
        {
            if (root->getParent())
                return "root has parent\n";
            if (!node)
                return "node memory error\n";
            if (node->isLeafNode())
            {
                if (node->getChildNum())
                    return "leaf has child\n";
                if (node->getFeatures().empty())
                    return "leaf feature memory error\n";
                if (node->getFeatureNum() >= M)
                    return "leaf feature num error\n";
                if (node->getFeatureNum() > 1)
                    for (size_t i = 0; i < node->getFeatureNum() - 1; i++)
                        for (size_t j = i + 1; j < node->getFeatureNum(); j++)
                        {
                            if (node->getFeature(i) == node->getFeature(j))
                                return "leaf has two same feature\n";
                        }
            }
            else
            {
                if (node->getFeatureNum())
                    return "inner has feature\n";
                if (node->getChildNum() >= M)
                    return "inner child num error\n";
                if (node->getChildNum() > 1)
                    for (size_t i = 0; i < node->getChildNum() - 1; i++)
                        for (size_t j = i + 1; j < node->getChildNum(); j++)
                        {
                            if (node->getChild(i) == node->getChild(j))
                                return "inner has two same child\n";
                        }
                for (size_t i = 0; i < node->getChildNum(); i++)
                    std::cout << checkTree(node->getChild(i));
            }
            if (node == root)
                return "ok\n";
            if (!node->getParent())
                return "not root node without parent\n";
            bool sign = true;
            for (size_t i = 0; i < node->getParent()->getChildNum(); i++)
                if (node->getParent()->getChild(i) == node)
                {
                    sign = false;
                    break;
                }
            if (sign)
                return "node parent dont have child node\n";
            return "ok\n";
        }

        /***
         * @Description: free malloc for whole tree
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {RNode<M>} *node
         * @return {*}
         */
        void freeTree(RNode<M> *node)
        {
            if (!node->isLeafNode())
                for (size_t i = 0; i < node->getChildNum(); i++)
                    freeTree(node->getChild(i));
            delete node;
        }

        virtual void setCapacity(int capacity) override
        {
            // DO NOTHING, since capacity is immutable in R tree
        }

        /***
         * @Description: Adjust Tree operation for RTree
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {shared_ptr<RNode<M>>} node
         * @return {pointer for root after split}
         */
        std::vector<RNode<M> *> adjustTree(RNode<M> *L1, RNode<M> *L2)
        {
            while (L1 != root)
            {
                RNode<M> *P = L1->getParent();

                P->add(L2);
                P->setEnvelope(P->getEnvelope().unionEnvelope(L2->getEnvelope()));

                RNode<M> *p = P->getParent();
                RNode<M> *l = P;

                while (p)
                {
                    p->setEnvelope(p->getEnvelope().unionEnvelope(l->getEnvelope()));
                    l = l->getParent();
                    p = l->getParent();
                }

                if (P->getChildNum() < size_t(M))
                {
                    L1 = root;
                    L2 = nullptr;
                    break;
                }

                std::vector<RNode<M> *> L = P->split();

                L1 = L[0];
                L2 = L[1];

                if (P == root)
                {
                    root = L1;
                }
                if (P)
                {
                    delete P;
                    P = nullptr;
                }
            }

            std::vector<RNode<M> *> ret;
            ret.clear();
            ret.push_back(L1);
            ret.push_back(L2);

            return ret;
        }

        /***
         * @Description: Choose Leaf operation for RTree
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {Feature} &feature
         * @return {*}
         */
        RNode<M> *chooseLeaf(const Feature *feature) const
        {
            RNode<M> *node = root;
            Envelope feature_env = feature->getEnvelope();
            while (!node->isLeafNode())
            {
                RNode<M> *child = node->getChild(0);
                double min_area_increase = child->getEnvelope().unionEnvelope(feature_env).getArea() - child->getEnvelope().getArea();

                for (decltype(M) i = 1; i < node->getChildNum(); i++)
                {
                    double now_increase = node->getChild(i)->getEnvelope().unionEnvelope(feature_env).getArea() - node->getChild(i)->getEnvelope().getArea();

                    if (now_increase < min_area_increase)
                    {
                        min_area_increase = now_increase;
                        child = node->getChild(i);
                    }
                }

                node = child;
            }

            return node;
        }

        /***
         * @Description: insert operation for RTree
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {Feature} &feature
         * @return {*}
         */
        bool insert(const Feature *feature)
        {

            RNode<M> *node = chooseLeaf(feature);

            node->add(feature);
            node->setEnvelope(node->getEnvelope().unionEnvelope(feature->getEnvelope()));

            RNode<M> *p = node->getParent();
            RNode<M> *l = node;

            while (p)
            {
                p->setEnvelope(p->getEnvelope().unionEnvelope(l->getEnvelope()));
                l = l->getParent();
                p = l->getParent();
            }

            if (node->getFeatureNum() < size_t(M))
                return true;

            std::vector<RNode<M> *> L = node->split();

            RNode<M> *L1 = L[0], *L2 = L[1];

            if (node == root)
                root = L1;
            if (node)
            {
                delete node;
                node = nullptr;
            }

            L = adjustTree(L1, L2);
            L1 = L[0];
            L2 = L[1];

            if (L2)
            {
                root = new RNode<M>(L1->getEnvelope());
                root->add(L1);

                root->setEnvelope(root->getEnvelope().unionEnvelope(L2->getEnvelope()));
                root->add(L2);
            }

            return true;
        }

        /***
         * @Description: Construct RTree By quadratic split
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @param {vector<Feature>} &features
         * @return {*}
         */
        virtual bool constructTree(const std::vector<Feature> &features) override
        {
            if (features.empty())
                return false;

            if (root)
                freeTree(root);
            root = nullptr;

            root = new RNode<M>(features.begin()->getEnvelope());
            root->add(features.begin().base());

            for (auto itr = features.begin() + 1; itr != features.end(); itr++)
            {
                // if (itr->getName() == "10")
                //     std::cout << "get:" << itr->getName() << std::endl;
                if (!insert(itr.base()))
                    return false;
                // std::cout << itr->getName() << std::endl;
                // std::cout << checkGeom(root);
            }

            bbox = root->getEnvelope();
            return true;
        }

        virtual void countNode(int &interiorNum, int &leafNum) override
        {
            interiorNum = leafNum = 0;
            if (root != nullptr)
                root->countNode(interiorNum, leafNum);
        }

        virtual void countHeight(int &height) override
        {
            height = 0;
            if (root != nullptr)
                height = root->countHeight(height);
        }

        /***
         * @Description: Range Query
         * @Author: jwimd chenjiewei@zju.edu.cn
         * @msg: None
         * @return {*}
         */
        virtual void rangeQuery(const Envelope &rect,
                                std::vector<Feature> &features) override
        {
            features.clear();

            std::vector<const Feature *> p_feature;
            p_feature.clear();

            if (root != nullptr)
                root->rangeQuery(rect, p_feature);

            for (size_t i = 0; i < p_feature.size(); ++i)
                features.push_back(*p_feature[i]);
        }

        virtual bool NNQuery(double x, double y,
                             std::vector<Feature> &features) override
        {
            features.clear();
            RNode<M> *node = root->pointInLeafNode(x, y);
            if (node)
            {
                double distance = -1;
                for (size_t i = 0; i < node->getFeatureNum(); ++i)
                    if (distance == -1 || node->getFeature(i)->getEnvelope().maxDistance(x, y) < distance)
                        distance = node->getFeature(i)->getEnvelope().maxDistance(x, y);
                Envelope rect(x - distance, x + distance, y - distance, y + distance);
                rangeQuery(rect, features);
                return true;
            }

            std::unique_ptr<Point>
                p(new Point(x, y));
            std::unique_ptr<Feature> f(new Feature("point", p.get()));
            node = chooseLeaf(f.get());

            double distance = -1;
            for (size_t i = 0; i < node->getFeatureNum(); ++i)
                if (distance == -1 || node->getFeature(i)->getEnvelope().maxDistance(x, y) < distance)
                    distance = node->getFeature(i)->getEnvelope().maxDistance(x, y);
            Envelope rect(x - distance, x + distance, y - distance, y + distance);
            rangeQuery(rect, features);
            return true;
        }

        virtual std::vector<std::pair<Feature, Feature>> spatialJoin(Tree *tree, double distance) override
        {
            RTree<M> *rtree = dynamic_cast<RTree<M> *>(tree);

            std::vector<std::pair<Feature, Feature>> joinSet;
            joinSet.clear();
            if (!rtree)
                return joinSet;
            // int innerCount1 = 0, innerCount2 = 0, leafCount1 = 0, leafCount2 = 0;
            // countNode(innerCount1, leafCount1);
            // rtree->countNode(innerCount2, leafCount2);

            // if (leafCount2 < leafCount1)
            //     rtree->getRoot()->spatialJoin(distance, root, joinSet);
            // else
            root->spatialJoin(distance, rtree->getRoot(), joinSet);

            
            return joinSet;
        }

        std::shared_ptr<RNode<M>> pointInLeafNode(double x, double y)
        {
            if (root != nullptr)
                return root->pointInLeafNode(x, y);
            else
                return nullptr;
        }

        virtual void draw() override
        {
            if (root != nullptr)
                root->draw();
        }

    public:
        static void test(int t);

        static void analyse();
    };

} // namespace hw6

#endif // !RTREE_SRC_H_INCLUDED
