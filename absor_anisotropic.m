function [tform, s, ErrorStats] = absor_anisotropic(A,B,w)
%ABSOR_ANISOTROPIC   Horn's method with anisotropic scaling.
%
% Function for mapping a collection of source points to a collectio of
% target points using Horn's method with anisotropic scaling. See absor.m 
% for more details on Horn's method.
%
% SYNTAX:
%   [tform, s, ErrorStats] = absor_anisotropic(A,B,w)
%
% INPUT:
%   A   - A 3-by-n matrix with source points.
%   B   - A 3-by-n matrix with target points.
%   w   - (optional) A 1-by-n vector with weights.
%
% OUTPUT:
%   tform      - Affine transformation on augmented matrix form.
%   s          - Vector with scale factors.
%   ErrorStats - Structure with error statistics.

% Weights:
if nargin==2
    w = ones(size(A,1),1);
end
% Find best anisotropic scaling:
f = @(s)findscaling(s,A,B,w); % Anonymous function used for parameterization.
s = fminunc(f, [1 1 1]);
S = diag(s);
% Horn's method:
[regParams,Bfit,ErrorStats] = absor(S*A', B', 'doTrans', true, 'doScale', true, 'weights', w/sum(w));
% Create affine transformation with anisotropic scaling:
T1 = regParams.M;
tform = T1*[S [0;0;0]; 0 0 0 1];

function err = findscaling(s, X, Y, w)
S = diag(s);
[~,~,ErrorStats] = absor(S*X', Y', ...
    'doTrans', true, 'doScale', true, 'weights', w/sum(w));
err = ErrorStats.errlsq;