%GRID2DRCFS class to perform raytracing in 2D with the fast sweeping method
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = grid2drcfs(par, nthreads)
%    clear g
%
%   Input for instantiation
%    par: struct variable for grid definition with the following fields
%          xmin: origin in X
%          zmin: origin in Z
%          dx: cell size in X
%          dz: cell size in Z
%          nx: number of cells in X
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
%    Tx: source coordinates, nTx by 2
%          1st column contains X coordinates,
%          3rd contains Z coordinates
%    Rx: receiver coordinates, nRx by 2
%          1st column contains X coordinates,
%          3rd contains Z coordinates
%    t0: source epoch, nTx by 1
%          t0 is optional (0 if not given)
%
%    *** IMPORTANT: Tx or Rx should _not_ lie on (or close to) an external
%                   face of the grid when rays are needed ***
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%    *** nSlowness must equal g.nx*g.nz ***
%
%
%   Output
%    tt:   vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 2
%    L:    data kernel matrix (tt = L*s)
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2015-03-04


classdef grid2drcfs < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        % Constructor - Create a new C++ class instance
        function this = grid2drcfs(varargin)
            this.objectHandle = grid2drcfs_mex('new', varargin{:});
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            grid2drcfs_mex('delete', this.objectHandle);
        end
        
        % setSlowness
        function varargout = setSlowness(this, varargin)
            [varargout{1:nargout}] = grid2drcfs_mex('setSlowness', this.objectHandle, varargin{:});
        end
        
        % raytrace
        function varargout = raytrace(this, varargin)
            [varargout{1:nargout}] = grid2drcfs_mex('raytrace', this.objectHandle, varargin{:});
        end
        
        
        function s = saveobj(obj)
            s.xmin = grid2drcfs_mex('get_xmin', obj.objectHandle);
            s.zmin = grid2drcfs_mex('get_zmin', obj.objectHandle);
            s.dx = grid2drcfs_mex('get_dx', obj.objectHandle);
            s.dz = grid2drcfs_mex('get_dz', obj.objectHandle);
            s.nx = grid2drcfs_mex('get_nx', obj.objectHandle);
            s.nz = grid2drcfs_mex('get_nz', obj.objectHandle);
            s.nthreads = grid2drcfs_mex('get_nthreads', obj.objectHandle);
        end
    end
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                obj = grid2drcfs(s, s.nthreads);
            else
                error('Wrong input arguments')
            end
        end
    end
end
