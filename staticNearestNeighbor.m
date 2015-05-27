function [v] = staticNearestNeighbor(a,d)
% performs static nearest-neighbor (mid-tread) quantization
% function [v] = staticNearestNeighbor(a,d)
% inputs: a is the argument to be quantized, d is quantization interval
% output: v is quantized version of a
% rounds towards -Inf on both edges
% future additions will include vector quantization using nearestneighbor()

% BR, 2/13/2013

% for 1D:

% mid-tread quantizer:
if( (a<=(d/2)) && (a>(-d/2)) )
    v = 0;
elseif( mod(a,d/2)==0 )  % round twds -Inf on bin edges
    v = d*floor(a/d);
else
    v = d*round(a/d);
end


% modify so if d is a vector, discrete set V^m (a lattice) is created...

% create discrete set V:

% call nearestneighbour search



%function [idx, tri] = nearestneighbour(varargin)
%NEARESTNEIGHBOUR    find nearest neighbours
%   IDX = NEARESTNEIGHBOUR(X) finds the nearest neighbour by Euclidean
%   distance to each point (column) in X from X. X is a matrix with points
%   as columns. IDX is a vector of indices into X, such that X(:, IDX) are
%   the nearest neighbours to X. e.g. the nearest neighbour to X(:, 2) is
%   X(:, IDX(2))
%
%   IDX = NEARESTNEIGHBOUR(P, X) finds the nearest neighbour by Euclidean
%   distance to each point in P from X. P and X are both matrices with the
%   same number of rows, and points are the columns of the matrices. Output
%   is a vector of indices into X such that X(:, IDX) are the nearest
%   neighbours to P