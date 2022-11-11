#include <set>
#include "QuadTree.h"

namespace hw6 {

/*
 * QuadNode
 */
void QuadNode::split(size_t capacity)
{
	for (int i = 0; i < 4; ++i) {
		delete nodes[i];
		nodes[i] = NULL;
	}

	// Task construction
	// Write your code here
}

void QuadNode::countNode(int& interiorNum, int& leafNum)
{
	if (isLeafNode()) {
		++leafNum;
	}
	else {
		++interiorNum;
		for (int i = 0; i < 4; ++i)
			nodes[i]->countNode(interiorNum, leafNum);
	}
}

int QuadNode::countHeight(int height)
{
	++height;
	if (!isLeafNode()) {
		int cur = height;
		for (int i = 0; i < 4; ++i) {
			height = max(height, nodes[i]->countHeight(cur));
		}
	}
	return height;
}

void QuadNode::rangeQuery(Envelope& rect, vector<Feature>& features)
{
	if (!bbox.intersect(rect))
		return;

	// Task range query
	// Write your code here

}

QuadNode* QuadNode::pointInLeafNode(double x, double y)
{
	// Task NN query
	// Write your code here

	return NULL;
}

void QuadNode::draw()
{
	if (isLeafNode()) {
		bbox.draw();
	}
	else {
		for (int i = 0; i < 4; ++i)
			nodes[i]->draw();
	}
}

/*
 * QuadTree
 */
bool QuadTree::constructQuadTree(vector<Feature>& features)
{
	if (features.empty())
		return false;

	// Task construction
	// Write your code here

	bbox = Envelope(-74.1, -73.8, 40.6, 40.8); // ע����д�����Ҫ����Ϊfeatures�İ�Χ�У�����ڵ�İ�Χ��

	return true;
}

void QuadTree::countQuadNode(int& interiorNum, int& leafNum)
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

void QuadTree::rangeQuery(Envelope& rect, vector<Feature>& features) 
{ 
	features.clear();

	// Task range query
	// Write your code here

	// filter step (ѡ���ѯ�����뼸�ζ����Χ���ཻ�ļ��ζ���)

	// ע���Ĳ��������ѯ�����غ�ѡ������������hw6��rangeQuery�����
}

bool QuadTree::NNQuery(double x, double y, vector<Feature>& features)
{
	if (!root || !(root->getEnvelope().contain(x, y)))
		return false;

	// Task NN query
	// Write your code here

	// filter step (ʹ��maxDistance2Envelope��������ò�ѯ�㵽���ζ����Χ�е���̵������룬Ȼ�������ѯ��ú�ѡ��)

	const Envelope& envelope = root->getEnvelope();
	double minDist = max(envelope.getWidth(), envelope.getHeight());
	

	// ע���Ĳ����ڽ���ѯ�����غ�ѡ������������hw6��NNQuery�����

	return true;
}

void QuadTree::draw()
{
	if (root)
		root->draw();
}

}
