//
//  SrcRcv.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2013-01-19.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_Interface_h
#define ttcr_u_Interface_h

#include <string>
#include <vector>

#include "ttcr_t.h"

template<typename P, class T>
class Interface {
public:
	Interface() {}
	Interface(const std::vector<P> &c, const std::vector<T> &t) :
	coord(c), tt(t) {}
	
	const std::vector<P>& get_coord() const { return coord; }
	
    const std::vector<T>& get_tt() const { return tt; }
    std::vector<T>& get_tt() { return tt; }
	
	const std::vector<std::vector<P>>& get_r_data() const { return r_data; }
	std::vector<std::vector<P>>& get_r_data() { return r_data; }
	
	void addPoint(const P& pt) {
		coord.push_back( pt );
		tt.push_back( 0 );
	}

private:
    std::vector<P> coord;
    std::vector<T> tt;
	std::vector<std::vector<P>> r_data;
};

#endif
