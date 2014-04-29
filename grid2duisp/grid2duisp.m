%GRID2DUISP class to perform raytracing on undulated surfaces
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = grid2duisp(nodes, triangles, nsecondary)
%    clear g
%
%   Input for instantiation
%    nodes: coordonates of mesh nodes (nNodes by 3)
%             1st column contains X coordinates, 2nd contains Y coordinates,
%             3rd contains Z coordinates
%    triangles: indices of nodes making mesh triangles (nCells by 3)
%    nsecondary: number of secondary nodes (optional: 5 if not set)
%
%  Raytracing
%    [tt] = raytrace(g, s, Tx, Rx, t0)
%    [tt, rays] = raytrace(g, s, Tx, Rx, t0)
%
%   Input
%    g: grid instance
%    s: slowness vector ( nSlowness by 1 )
%    Tx: source coordinates, nTx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates,
%          3rd contains Z coordinates
%    Rx: receiver coordinates, nRx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates,
%          3rd contains Z coordinates
%    t0: source epoch, nTx by 1
%          t0 is optional (0 if not given)
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%    *** nSlowness must equal the number of nodes in g ***
%    *** Indexing of slowness values must follow node order ***
%
%
%   Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 3
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2014-04-26


classdef grid2duisp < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = grid2duisp(varargin)
            this.objectHandle = grid2duisp_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            grid2duisp_mex('delete', this.objectHandle);
        end

        %% setSlowness
        function varargout = setSlowness(this, varargin)
            [varargout{1:nargout}] = grid2duisp_mex('setSlowness', this.objectHandle, varargin{:});
        end

        %% raytrace
        function varargout = raytrace(this, varargin)
            [varargout{1:nargout}] = grid2duisp_mex('raytrace', this.objectHandle, varargin{:});
        end
    end
end