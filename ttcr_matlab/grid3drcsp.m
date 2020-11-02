%GRID3DRCSP class to perform raytracing in 3D with the shortest-path method
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = grid3drcsp(par, nthreads)
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
%          nsx: number of secondary nodes per cell in X
%          nsy: number of secondary nodes per cell in Y
%          nsz: number of secondary nodes per cell in Z
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


classdef grid3drcsp < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        % Constructor - Create a new C++ class instance
        function this = grid3drcsp(varargin)
            this.objectHandle = grid3drcsp_mex('new', varargin{:});
        end

        % Destructor - Destroy the C++ class instance
        function delete(this)
            grid3drcsp_mex('delete', this.objectHandle);
        end

        % setSlowness
        function varargout = setSlowness(this, varargin)
            [varargout{1:nargout}] = grid3drcsp_mex('setSlowness', this.objectHandle, varargin{:});
        end

        % raytrace
        function varargout = raytrace(this, varargin)
            [varargout{1:nargout}] = grid3drcsp_mex('raytrace', this.objectHandle, varargin{:});
        end

        % for saving in mat-files
        function s = saveobj(obj)
            s.xmin = grid2drcfs_mex('get_xmin', obj.objectHandle);
            s.ymin = grid2drcfs_mex('get_ymin', obj.objectHandle);
            s.zmin = grid2drcfs_mex('get_zmin', obj.objectHandle);
            s.dx = grid2drcfs_mex('get_dx', obj.objectHandle);
            s.dy = grid2drcfs_mex('get_dy', obj.objectHandle);
            s.dz = grid2drcfs_mex('get_dz', obj.objectHandle);
            s.nx = grid2drcfs_mex('get_nx', obj.objectHandle);
            s.ny = grid2drcfs_mex('get_ny', obj.objectHandle);
            s.nz = grid2drcfs_mex('get_nz', obj.objectHandle);
            s.nsx = grid2drcfs_mex('get_nsx', obj.objectHandle);
            s.nsy = grid2drcfs_mex('get_nsy', obj.objectHandle);
            s.nsz = grid2drcfs_mex('get_nsz', obj.objectHandle);
            s.nthreads = grid2drcfs_mex('get_nthreads', obj.objectHandle);
        end
    end
    methods(Static)
        % for loading from mat-files
        function obj = loadobj(s)
            if isstruct(s)
                obj = grid2drcfs(s, s.nthreads);
            else
                error('Wrong input arguments')
            end
        end
    end
end
