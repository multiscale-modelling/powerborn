/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef TERNARYNET_H_
#define TERNARYNET_H_

#include "Typedefs.h"
#include "ObjectPool.h"

namespace powerborn
{

class TernaryNode;
class TernaryNet;
struct VertexData;

typedef TernaryNet PowerVertices;
typedef TernaryNode Vertex;

struct VertexData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    UInt4 generators; // first entry is VertexData id, 1-3 are generator ids
    UInt2 commonGenerators(const VertexData& other) const
    {
        // this function assumes that there are at least two equal common generators
        const UInt4 gen1 = generators;
        const UInt4 gen2 = other.generators;
        UInt2 tmp2 = UInt2::Zero();
        for(unsigned int i=1; i<4; ++i)
        {
            tmp2[0] |=  gen1[1] == gen2[i];
            tmp2[1] |=  gen1[3] == gen2[i];
        }
        tmp2[0] = tmp2[0] ? gen1[1] : gen1[2];
        tmp2[1] = tmp2[1] ? gen1[3] : gen1[2];
        return tmp2;
    };
    VertexData(): generators(UInt4::Zero())
    {
    }
    VertexData(const VertexData& other): generators(other.generators)
    {
    }
};

class TernaryNode
{
public:
    struct Connects
    {
        TernaryNode* con_[3];

        Connects() {}
        Connects(TernaryNode* p1, TernaryNode* p2, TernaryNode* p3)
        {
            con_[0] = p1;
            con_[1] = p2;
            con_[2] = p3;
        }
        TernaryNode*& operator[](unsigned int i) {return con_[i];}
        const TernaryNode* operator[](unsigned int i) const {return con_[i];}
    };
	//typedef Eigen::Array<TernaryNode*, 3, 1> Connects;
private:
	friend class TernaryNet;
    VertexData data_;
	Connects connects_;
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	unsigned int id() const {return data_.generators[0];}
	void setId(unsigned int i) {data_.generators[0] = i;}
	Connects& connects() {return connects_;}
	TernaryNode*& connects(unsigned int i) {return connects_[i];}
	const TernaryNode* connects(unsigned int i) const {return connects_[i];}
	const Connects& connects() const {return connects_;}
	//void reset() {connects_.setConstant(NULL);}
	TernaryNode() {}
	~TernaryNode(){}
	VertexData& data() {return data_;}
	const VertexData& data() const {return data_;}
};

class TernaryNet
{
public:
	typedef TernaryNode NodeType;
	typedef AlignedVector<NodeType*>::Type Container;
	typedef AlignedDeque<NodeType>::Type NodeContainer;
	typedef Container::iterator Iterator;
	typedef Container::const_iterator ConstIterator;

private:
	Pool<NodeType> nodes_;
	Container net_;
	unsigned int size_;

	bool checkConnections(const NodeType* const node) const
	{
		bool ok = 1;
		for(int i=0; i<3; ++i)
		{
			bool nb_ok = 0;
			for(int j=0; j<3; ++j)
			{
				nb_ok |= node->connects_[i]->connects_[j] == node;
			}
			ok &= nb_ok;
		}
		return ok;
	}
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	TernaryNet(): size_(0) {}
	void reconnect(NodeType::Connects& connects, const NodeType* op, NodeType* np)
	{
		// replace all connections to op with connections to np
		for(int i=0; i<3; ++i)
		{
			connects[i] = op == connects[i] ? np : connects[i];
		}
	}
	void reset()
	{
		size_ = 0;
		nodes_.reset();
	}
	NodeType& getNewNode() {return nodes_.getObject();}
	NodeType*& node(unsigned int i) {return net_[i];}
	const NodeType* node(unsigned int i) const {return net_[i];}
	Container& nodes() {return net_;}
    const Container& nodes() const {return net_;}

	Iterator begin() {return net_.begin();}
	ConstIterator begin() const {return net_.begin();}
	Iterator end() {return net_.end();}
	ConstIterator end() const {return net_.end();}

	unsigned int size() const {return size_;}
	void setSize(unsigned int i) {size_ = i;}
	void setZero() {size_ = 0;}
	bool check() const
	{
		// consistency check for connections
		bool all_ok = 1;
		for(unsigned int i=0; i<this->size(); ++i)
		{
			bool ok = this->checkConnections(this->node(i));
			if(!ok)
			{
				std::cout << "Conn Check not ok for node id " << this->node(i)->id() << std::endl;
				this->printConnections(this->node(i));
			}
			all_ok &= ok;
		}
		if(!all_ok)
		{
			std::cout << "Connection check failed " << std::endl;
		}
		return all_ok;
	}
	static void printConnections(const NodeType* const node)
	{
		std::cout << "connections for " << node->id() << ": ";
		for(int i=0; i<3; ++i)
		{
			std::cout << node->connects_[i]->id() << " ";
		}
		std:: cout << "| ";
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				std::cout << node->connects_[i]->connects_[j]->id() << " ";
			}
			std::cout << "| ";
		}
		std::cout << std::endl;
	}
};

}

#endif /* TERNARYNET_H_ */
