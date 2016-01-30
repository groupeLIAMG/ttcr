%GRID3DRCFS class to perform raytracing in 3D with the fast sweeping method
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = grid3drcfs(par, nthreads)
%    clear g
%
%   Input for instantiation
%    par: struct variable for grid definition with the following fields
%          xmin: origin in X
%          ymin: origin in Y
%          zmin: origin in Z
%          dx: cell size in X
%          dy: cell size in Y
%          dz: cell size in Z
%          nx: number of cells in X
%          ny: number of cells in Y
%          nz: number of cells in Z
%    nthreads: number of threads (optional, default = 1)
%
%  Raytracing
%    [tt] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays, L] = g.raytrace(s, Tx, Rx, t0)
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
%    *** nSlowness must equal g.nx*g.ny*g.nz ***
%    *** Indexing of slowness values is done by "vectorizing" a 3D array,
%        i.e. if slowness field s is of size (nx,ny,nz), enter s(:) as
%        first argument
%
%
%   Output
%    tt:   vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 3
%    L:    data kernel matrix (tt = L*s)
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2015-03-04


classdef grid3drcfs < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = grid3drcfs(varargin)
            this.objectHandle = grid3drcfs_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            grid3drcfs_mex('delete', this.objectHandle);
        end
        
        %% setSlowness
        function varargout = setSlowness(this, varargin)
            [varargout{1:nargout}] = grid3drcfs_mex('setSlowness', this.objectHandle, varargin{:});
        end
        
        %% raytrace
        function varargout = raytrace(this, varargin)
            [varargout{1:nargout}] = grid3drcfs_mex('raytrace', this.objectHandle, varargin{:});
        end
    end
end
