%GRID3DUNFS class to perform raytracing in 3D with the fast sweeping method
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = grid3dunfs(nodes, tetrahedra, nthreads)
%    clear g
%
%   Input for instantiation
%    nodes: coordonates of mesh nodes (nNodes by 3)
%             1st column contains X coordinates, 2nd contains Y coordinates,
%             3rd contains Z coordinates
%    tetrahedra: indices of nodes making mesh tetrahedra (nCells by 4)
%    nthreads: number of threads (optional, default = 1)
%
%  Raytracing
%    [tt] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays, M] = g.raytrace(s, Tx, Rx, t0)
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
%    *** IMPORTANT: Tx or Rx should _not_ lie on (or close to) an external
%                   face of the grid when rays are needed ***
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%    *** nSlowness must equal the number of nodes in g ***
%    *** Indexing of slowness values must follow node order ***
%
%
%   Output
%     tt: vector of traveltimes, nRx by 1
%     rays: cell object containing the matrices of coordinates of the ray
%           paths, nRx by 1.  Each matrix is nPts by 3
%     M:    cell object containing the matrices of partial derivative dt/dV and
%           dt/dsc along rays (t is time, V is velocity ans sc is static
%           correction)
%       *** Important ***
%           M contains as many cells as unique Tx coordinates.  Ordering of Tx
%           and Rx coordinates is important to retrieve the correct results: Tx
%           & Rx data should be ordered so that redundant Tx should be
%           contiguous, i.e. for 2 Tx and 2 Rx we should have
%
%           Tx = [Tx_x1 Tx_y1 Tx_z1
%                 Tx_x1 Tx_y1 Tx_z1
%                 Tx_x2 Tx_y2 Tx_z2
%                 Tx_x2 Tx_y2 Tx_z2];
%           Rx = [Rx_x1 Rx_y1 Rx_z1
%                 Rx_x2 Rx_y2 Rx_z2
%                 Rx_x1 Rx_y1 Rx_z1
%                 Rx_x2 Rx_y2 Rx_z2];
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2015-03-04


classdef grid3dunfs < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = grid3dunfs(varargin)
            this.objectHandle = grid3dunfs_mex('new', varargin{:});
        end

        %% Destructor - Destroy the C++ class instance
        function delete(this)
            grid3dunfs_mex('delete', this.objectHandle);
        end

        %% setSlowness
        function varargout = setSlowness(this, varargin)
            [varargout{1:nargout}] = grid3dunfs_mex('setSlowness', this.objectHandle, varargin{:});
        end

        %% raytrace
        function varargout = raytrace(this, varargin)
            [varargout{1:nargout}] = grid3dunfs_mex('raytrace', this.objectHandle, varargin{:});
        end
    end
end
