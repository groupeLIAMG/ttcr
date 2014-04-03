//
//  Metric.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2014-02-01.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_Metric_h
#define ttcr_u_Metric_h

#include "Node.h"

template <typename T>
class Metric {
public:
	virtual T l(const Node<T>&, const sxz<T>&) const = 0;
	virtual T l(const Node<T>&, const sxyz<T>&) const = 0;
	virtual ~Metric() {}
};


template <typename T>
class Metric1 : public Metric<T> {
public:
	T l(const Node<T>& n, const sxz<T>& s) const {
		return fabs(n.getX()-s.x) + fabs(n.getZ()-s.z);
	}
	T l(const Node<T>& n, const sxyz<T>& s) const {
		return fabs(n.getX()-s.x) + fabs(n.getY()-s.y) + fabs(n.getZ()-s.z);
	}
	~Metric1() {}
};

template <typename T>
class Metric2 : public Metric<T> {
public:
	T l(const Node<T>& n, const sxz<T>& s) const {
		return sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
					 (n.getZ()-s.z)*(n.getZ()-s.z) );
	}
	T l(const Node<T>& n, const sxyz<T>& s) const {
		return sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
					(n.getY()-s.y)*(n.getY()-s.y)  +
					(n.getZ()-s.z)*(n.getZ()-s.z) );
	}
	~Metric2() {}
};

#endif
