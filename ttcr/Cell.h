//
//  Cell.h
//  ttcr
//
//  Created by Bernard Giroux on 16-02-27.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#ifndef Cell_h
#define Cell_h

#include <vector>

template <typename T, typename NODE, typename S>
class Cell {
public:
    Cell(const size_t n) : slowness(std::vector<T>(n)) { }
    
    const T getSlowness(const size_t n) const { return slowness[n]; }
    
    void setSlowness(const T s) {
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s;
        }
    }
    
    int setSlowness(const T *s, const size_t ns) {
        if ( slowness.size() != ns ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    int setSlowness(const std::vector<T>& s) {
        if ( slowness.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    T computeDt(const NODE& source, const S& node,
                const size_t cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
    T computeDt(const NODE& source, const NODE& node,
                const size_t cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
private:
    std::vector<T> slowness;
};

#endif /* Cell_h */
