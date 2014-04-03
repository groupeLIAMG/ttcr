//
//  Node.h
//  ttcr
//
//  Created by Bernard Giroux on 12-08-14.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#ifndef __NODE_H__
#define __NODE_H__


template<typename T>
class Node {
public:
	virtual T getTT(const size_t n) const = 0;
	virtual int getDimension() const = 0;
	virtual T getX() const = 0;
	virtual T getY() const = 0;
	virtual T getZ() const = 0;
};



template<typename T>
class CompareNodePtr {
	// Overloaded operator for the priority queue, compare the "n"th traveltimes of two nodes.
private:
	size_t n;
public:
	CompareNodePtr(const size_t nn) : n(nn) {}
    bool operator()(const Node<T>* n1, const Node<T>* n2) const {
        //  The priority_queue must return the minimum time!!!
        return n1->getTT(n) > n2->getTT(n);
    }
};


#endif
