//
//  Cell.h
//  ttcr
//
//  Created by Bernard Giroux on 16-02-27.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef Cell_h
#define Cell_h

#include <cmath>
#include <iostream>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {
    
    template <typename T, typename NODE, typename S>
    class Cell {
    public:
        Cell(const size_t n) : slowness(std::vector<T>(n)) { }
        
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
        
        int setXi(const std::vector<T>& s) {
            std::cerr << "Error: xi not defined for Cell.";
            return 1;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for Cell.";
            return 1;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for Cell.";
            return 1;
        }
        
        int setVs0(const std::vector<T>& s) {
            std::cerr << "Error: Vs0 not defined for Cell.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for Cell.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for Cell.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for Cell.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            return slowness[cellNo] * source.getDistance( node );
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            return slowness[cellNo] * source.getDistance( node );
        }
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v = source.getDistance( node );
        }
        
    private:
        std::vector<T> slowness;
    };
    
    
    // Elliptical anisotropy in 2D (Y Dimension ignored)
    
    template <typename T, typename NODE, typename S>
    class CellElliptical {
    public:
        CellElliptical(const size_t n) :
        slowness(std::vector<T>(n)),
        xi(std::vector<T>(n,0.0)) {
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
        
        int setXi(const std::vector<T>& s) {
            if ( xi.size() != s.size() ) {
                std::cerr << "Error: xi vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<xi.size(); ++n ) {
                xi[n] = s[n]*s[n];
            }
            return 0;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellElliptical.";
            return 1;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for CellElliptical.";
            return 1;
        }
        
        int setVs0(const std::vector<T>& s) {
            std::cerr << "Error: Vs0 not defined for CellElliptical.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for CellElliptical.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for CellElliptical.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for CellElliptical.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T lz = node.z - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T lz = node.getZ() - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = fabs(node.x - source.getX());
            cell.v2 = fabs(node.z - source.getZ());
        }
        
    private:
        std::vector<T> slowness;
        std::vector<T> xi;        // anisotropy ratio, xi = sz / sx, *** squared ***
    };
    
    
    // Elliptical anisotropy in 2D (Y Dimension ignored)
    
    template <typename T, typename NODE, typename S>
    class CellTiltedElliptical {
    public:
        CellTiltedElliptical(const size_t n) :
        slowness(std::vector<T>(n)),
        xi(std::vector<T>(n,0.0)),
        tAngle(std::vector<T>(n,0.0)),
        ca(std::vector<T>(n,1.0)),
        sa(std::vector<T>(n,0.0)) {
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
        
        int setXi(const std::vector<T>& s) {
            if ( xi.size() != s.size() ) {
                std::cerr << "Error: xi vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<xi.size(); ++n ) {
                xi[n] = s[n]*s[n];
            }
            return 0;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            if ( tAngle.size() != s.size() ) {
                std::cerr << "Error: angle vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<tAngle.size(); ++n ) {
                tAngle[n] = s[n];
                ca[n] = std::cos(s[n]);
                sa[n] = std::sin(s[n]);
            }
            return 0;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for CellTiltedElliptical.";
            return 1;
        }
        
        int setVs0(const std::vector<T>& s) {
            std::cerr << "Error: Vs0 not defined for CellTiltedElliptical.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for CellTiltedElliptical.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for CellTiltedElliptical.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for CellTiltedElliptical.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T lz = node.z - source.getZ();
            T t1 = lx * ca[cellNo] + lz * sa[cellNo];
            T t2 = lz * ca[cellNo] - lx * sa[cellNo];
            
            return slowness[cellNo] * std::sqrt( t1*t1 + xi[cellNo]*xi[cellNo]*t2*t2 );
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T lz = node.getZ() - source.getZ();
            T t1 = lx * ca[cellNo] + lz * sa[cellNo];
            T t2 = lz * ca[cellNo] - lx * sa[cellNo];
            
            return slowness[cellNo] * std::sqrt( t1*t1 + xi[cellNo]*xi[cellNo]*t2*t2 );
        }
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = fabs(node.x - source.getX());
            cell.v2 = fabs(node.z - source.getZ());
        }
        
    private:
        std::vector<T> slowness;
        std::vector<T> xi;        // anisotropy ratio, xi = sz / sx, *** squared ***
        std::vector<T> tAngle;    // column-wise (z axis) anisotropy angle of the cells, in radians
        std::vector<T> ca;        // cosine of tAngle
        std::vector<T> sa;        // sine of tAngle
    };
    
    
    
    
    //  VTI anisotropy, P or SV phase, in 2D (Y Dimension ignored)
    template <typename T, typename NODE, typename S>
    class CellVTI_PSV {
    public:
        CellVTI_PSV(const size_t n) :
        sign(1.0),
        Vp0(std::vector<T>(n)),
        Vs0(std::vector<T>(n)),
        epsilon(std::vector<T>(n)),
        delta(std::vector<T>(n)) {
        }
        
        int setVp0(const std::vector<T>& s) {
            if ( Vp0.size() != s.size() ) {
                std::cerr << "Error: Vp0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vp0.size(); ++n ) {
                Vp0[n] = s[n];
            }
            return 0;
        }
        
        int setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                std::cerr << "Error: Vs0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vs0.size(); ++n ) {
                Vs0[n] = s[n];
            }
            return 0;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            if ( epsilon.size() != s.size() ) {
                std::cerr << "Error: epsilon vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<epsilon.size(); ++n ) {
                epsilon[n] = s[n];
            }
            return 0;
        }
        
        int setDelta(const std::vector<T>& s) {
            if ( delta.size() != s.size() ) {
                std::cerr << "Error: delta vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<delta.size(); ++n ) {
                delta[n] = s[n];
            }
            return 0;
        }
        
        void setPhase(const int p) {
            if ( p==1 ) sign = 1.;  // P wave
            else sign = -1.;        // SV wave
        }
        
        int setXi(const std::vector<T>& s) {
            std::cerr << "Error: xi not defined for CellVTI_PSV.";
            return 1;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellVTI_PSV.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for CellVTI_PSV.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T theta = atan2(node.x - source.getX(), node.z - source.getZ());
            T f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);
            
            T tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;
            
            tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
            sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );
            
            T v = Vp0[cellNo] * sqrt( tmp );
            return source.getDistance( node ) / v;
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T theta = atan2(node.getX() - source.getX(), node.getZ() - source.getZ());
            T f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);
            
            T tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;
            
            tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
            sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );
            
            T v = Vp0[cellNo] * sqrt( tmp );
            return source.getDistance( node ) / v;
        }
        
    private:
        T sign;
        std::vector<T> Vp0;
        std::vector<T> Vs0;
        std::vector<T> epsilon;
        std::vector<T> delta;
    };
    
    
    
    //  VTI anisotropy, SH phase, in 2D (Y Dimension ignored)
    template <typename T, typename NODE, typename S>
    class CellVTI_SH {
    public:
        CellVTI_SH(const size_t n) :
        Vs0(std::vector<T>(n)),
        gamma(std::vector<T>(n)) {
        }
        
        int setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                std::cerr << "Error: Vs0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vs0.size(); ++n ) {
                Vs0[n] = s[n];
            }
            return 0;
        }
        
        int setGamma(const std::vector<T>& s) {
            if ( gamma.size() != s.size() ) {
                std::cerr << "Error: gamma vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<gamma.size(); ++n ) {
                gamma[n] = s[n];
            }
            return 0;
        }
        
        int setXi(const std::vector<T>& s) {
            std::cerr << "Error: xi not defined for CellVTI_SH.";
            return 1;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellVTI_SH.";
            return 1;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for CellVTI_SH.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for CellVTI_SH.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for CellVTI_SH.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T theta = atan2(node.x - source.getX(), node.z - source.getZ());
            T v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
            return source.getDistance( node ) / v;
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T theta = atan2(node.getX() - source.getX(), node.getZ() - source.getZ());
            T v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
            return source.getDistance( node ) / v;
        }
        
    private:
        std::vector<T> Vs0;
        std::vector<T> gamma;
    };
    
    
    
    
    
    
    
    // Elliptical anisotropy in 3D
    
    template <typename T, typename NODE, typename S>
    class CellElliptical3D {
    public:
        CellElliptical3D(const size_t n) :
        slowness(std::vector<T>(n)),
        xi(std::vector<T>(n)) {
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
        
        int setXi(const std::vector<T>& s) {
            if ( xi.size() != s.size() ) {
                std::cerr << "Error: xi vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<xi.size(); ++n ) {
                xi[n] = s[n]*s[n];
            }
            return 0;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellElliptical3D.";
            return 1;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for CellElliptical3D.";
            return 1;
        }
        
        int setVs0(const std::vector<T>& s) {
            std::cerr << "Error: Vs0 not defined for CellElliptical3D.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for CellElliptical3D.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for CellElliptical3D.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for CellElliptical3D.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T lz = node.z - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T lz = node.getZ() - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }
        
    private:
        std::vector<T> slowness;
        std::vector<T> xi;        // anisotropy ratio, xi = sz / sx, *** squared ***
    };
    
    
    
    //  VTI anisotropy, P or SV phase, in 2D (Y Dimension ignored)
    template <typename T, typename NODE, typename S>
    class CellVTI_PSV3D {
    public:
        CellVTI_PSV3D(const size_t n) :
        sign(1.0),
        Vp0(std::vector<T>(n)),
        Vs0(std::vector<T>(n)),
        epsilon(std::vector<T>(n)),
        delta(std::vector<T>(n)) {
        }
        
        int setVp0(const std::vector<T>& s) {
            if ( Vp0.size() != s.size() ) {
                std::cerr << "Error: Vp0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vp0.size(); ++n ) {
                Vp0[n] = s[n];
            }
            return 0;
        }
        
        int setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                std::cerr << "Error: Vs0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vs0.size(); ++n ) {
                Vs0[n] = s[n];
            }
            return 0;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            if ( epsilon.size() != s.size() ) {
                std::cerr << "Error: epsilon vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<epsilon.size(); ++n ) {
                epsilon[n] = s[n];
            }
            return 0;
        }
        
        int setDelta(const std::vector<T>& s) {
            if ( delta.size() != s.size() ) {
                std::cerr << "Error: delta vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<delta.size(); ++n ) {
                delta[n] = s[n];
            }
            return 0;
        }
        
        void setPhase(const int p) {
            if ( p==1 ) sign = 1.;  // P wave
            else sign = -1.;        // SV wave
        }
        
        int setXi(const std::vector<T>& s) {
            std::cerr << "Error: xi not defined for CellVTI_PSV3D.";
            return 1;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellVTI_PSV3D.";
            return 1;
        }
        
        int setGamma(const std::vector<T>& s) {
            std::cerr << "Error: gamma not defined for CellVTI_PSV3D.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T theta = atan2(lx, node.z - source.getZ());
            T f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);
            
            T tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;
            
            tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
            sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );
            
            T v = Vp0[cellNo] * sqrt( tmp );
            return source.getDistance( node ) / v;
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T theta = atan2(lx, node.getZ() - source.getZ());
            T f = 1. - (Vs0[cellNo]*Vs0[cellNo]) / (Vp0[cellNo]*Vp0[cellNo]);
            
            T tmp = 1. + (2.*epsilon[cellNo]*sin(theta)*sin(theta)) / f;
            
            tmp = 1. + epsilon[cellNo]*sin(theta)*sin(theta) - f/2. +
            sign*f/2.*sqrt( tmp*tmp - (2.*(epsilon[cellNo]-delta[cellNo])*sin(2.*theta)*sin(2.*theta))/f );
            
            T v = Vp0[cellNo] * sqrt( tmp );
            return source.getDistance( node ) / v;
        }
        
    private:
        T sign;
        std::vector<T> Vp0;
        std::vector<T> Vs0;
        std::vector<T> epsilon;
        std::vector<T> delta;
    };
    
    
    
    //  VTI anisotropy, SH phase, in 2D (Y Dimension ignored)
    template <typename T, typename NODE, typename S>
    class CellVTI_SH3D {
    public:
        CellVTI_SH3D(const size_t n) :
        Vs0(std::vector<T>(n)),
        gamma(std::vector<T>(n)) {
        }
        
        int setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                std::cerr << "Error: Vs0 vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<Vs0.size(); ++n ) {
                Vs0[n] = s[n];
            }
            return 0;
        }
        
        int setGamma(const std::vector<T>& s) {
            if ( gamma.size() != s.size() ) {
                std::cerr << "Error: gamma vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<gamma.size(); ++n ) {
                gamma[n] = s[n];
            }
            return 0;
        }
        
        int setXi(const std::vector<T>& s) {
            std::cerr << "Error: xi not defined for CellVTI_SH3D.";
            return 1;
        }
        
        int setTiltAngle(const std::vector<T>& s) {
            std::cerr << "Error: TiltAngle not defined for CellVTI_SH3D.";
            return 1;
        }
        
        int setVp0(const std::vector<T>& s) {
            std::cerr << "Error: Vp0 not defined for CellVTI_SH3D.";
            return 1;
        }
        
        int setDelta(const std::vector<T>& s) {
            std::cerr << "Error: delta not defined for CellVTI_SH3D.";
            return 1;
        }
        
        int setEpsilon(const std::vector<T>& s) {
            std::cerr << "Error: epsilon not defined for CellVTI_SH3D.";
            return 1;
        }
        
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T theta = atan2(lx, node.z - source.getZ());
            T v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
            return source.getDistance( node ) / v;
        }
        
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T theta = atan2(lx, node.getZ() - source.getZ());
            T v = Vs0[cellNo] * sqrt(1. + 2.*gamma[cellNo]*sin(theta)*sin(theta));
            return source.getDistance( node ) / v;
        }
        
    private:
        std::vector<T> Vs0;
        std::vector<T> gamma;
    };
    
}

#endif /* Cell_h */
