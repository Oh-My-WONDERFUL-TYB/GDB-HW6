#ifndef RTREE_SRC_H_INCLUDED
#define RTREE_SRC_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"

#include <algorithm>
#include <array>
#include <queue>
#include <string>
#include <vector>

#include "CMakeIn.h"

namespace hw6 {

/// <summary>
/// </summary>
/// <param name="M">Maximum number of children of each nodes.</param>
template <uint8_t M> class RNode {
  private:
    RNode<M> *parent = nullptr;
    std::array<RNode<M> *, M> children = []() {
        decltype(children) ret;
        for (decltype(M) i = 0; i < M; ++i)
            ret[i] = nullptr;
        return ret;
    }();
    decltype(M) childrenNum = 0;
    Envelope bbox;
    std::vector<Feature> features;

  public:
    RNode() = delete;
    RNode(const Envelope &box) : bbox(box) {}

    inline bool isLeafNode() { return childrenNum == 0; }

    inline const Envelope &getEnvelope() { return bbox; }

    inline RNode<M> *getParent() { return parent; }

    inline void setEnvelope(const Envelope &box) { bbox = box; }

    inline RNode<M> *getChildNode(size_t i) {
        return i < childrenNum ? children[i] : nullptr;
    }

    inline const RNode<M> *getChildNode(size_t i) const {
        return i < childrenNum ? children[i] : nullptr;
    }

    inline decltype(M) getChildNum() const { return childrenNum; }

    inline size_t getFeatureNum() const { return features.size(); }

    inline const Feature &getFeature(size_t i) const { return features[i]; }

    inline const std::vector<Feature> &getFeatures() const { return features; }

    inline void add(const Feature &f) { features.push_back(f); }

    inline void add(RNode<M> *child) {
        children[childrenNum] = child;
        child->parent = this;
        ++childrenNum;
    }

    inline void remove(const Feature &f) {
        auto where = [&]() {
            for (auto itr = features.begin(); itr != features.end(); ++itr)
                if (itr->getName() == f.getName())
                    return itr;
        }();
        features.erase(where);
        if (features.empty())
            features.shrink_to_fit(); // free memory unused but allocated
    }

    inline void remove(RNode<M> *child) {
        for (decltype(M) i = 0; i < childrenNum; ++i)
            if (children[i] == child) {
                --childrenNum;
                std::swap(children[i], children[childrenNum]);
                children[childrenNum] = nullptr;
                break;
            }
    }

    inline Feature popBackFeature() {
        auto ret = features.back();
        features.pop_back();
        return ret;
    }

    inline RNode<M> *popBackChildNode() {
        --childrenNum;
        auto ret = children[childrenNum];
        children[childrenNum] = nullptr;
        return ret;
    }

    void countNode(int &interiorNum, int &leafNum) {
        if (isLeafNode()) {
            ++leafNum;
        } else {
            ++interiorNum;
            for (decltype(M) i = 0; i < childrenNum; ++i)
                children[i]->countNode(interiorNum, leafNum);
        }
    }

    int countHeight(int height) {
        ++height;
        if (!isLeafNode()) {
            int cur = height;
            for (decltype(M) i = 0; i < childrenNum; ++i)
                height = max(height, children[i]->countHeight(cur));
        }
        return height;
    }

    inline void draw() {
        if (isLeafNode()) {
            bbox.draw();
        } else
            for (decltype(M) i = 0; i < childrenNum; ++i)
                children[i]->draw();
    }

    void rangeQuery(const Envelope &rect, std::vector<Feature> &features) {
        // TODO
    }

    RNode<M> *pointInLeafNode(double x, double y) {
        // TODO
    }
};

template <uint8_t M> class RTree : public Tree {
  private:
    // Vars here MAY be useful, but it's ok to delete them {
    constexpr static auto m = M >> 1;
    constexpr static auto M_minus_m = M - m;
    constexpr static double EQ_ERR = 0.0000000005;
    // }

    RNode<M> *root = nullptr;

  public:
    RTree() : Tree(M) { static_assert(M >= 4); }
    ~RTree() {
        if (root != nullptr)
            delete root;
        root = nullptr;
    }

    virtual void setCapacity(int capacity) override {
        // DO NOTHING, since capacity is immutable in R tree
    }

    virtual bool constructTree(const std::vector<Feature> &features) override {
        // TODO
        return false;
    }

    virtual void countNode(int &interiorNum, int &leafNum) override {
        interiorNum = leafNum = 0;
        if (root != nullptr)
            root->countNode(interiorNum, leafNum);
    }

    virtual void countHeight(int &height) override {
        height = 0;
        if (root != nullptr)
            height = root->countHeight(height);
    }

    virtual void rangeQuery(const Envelope &rect,
                            std::vector<Feature> &features) override {
        features.clear();
        if (root != nullptr)
            root->rangeQuery(rect, features);
    }

    virtual bool NNQuery(double x, double y,
                         std::vector<Feature> &features) override {
        features.clear();
        // TODO
        return false;
    }

    RNode<M> *pointInLeafNode(double x, double y) {
        if (root != nullptr)
            return root->pointInLeafNode(x, y);
        else
            return nullptr;
    }

    virtual void draw() override {
        if (root != nullptr)
            root->draw();
    }

  public:
    static void test(int t);

    static void analyse();
};

} // namespace hw6

#endif // !RTREE_SRC_H_INCLUDED
