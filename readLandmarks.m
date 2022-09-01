function [coordinates, labels, type] = readLandmarks(fname)
%READLANDMARKS   Function to read a landmark text file.
%
% SYNTAX:
%   [coordinates, labels, type] = readLandmarks(fname)
%
% INPUT:
%   fname       - Filename.
%
% OUTPUT:
%   coordinates - Landmark coordinates defined as ML,AP,DV with positive 
%                 directions towards right, anterior, superior, 
%                 respectively.
%   labels      - Cell array with landmark names.
%   type        - Cell array with landmark types. Possible values are
%                   M = single landmark on the median plane
%                   L = left landmark in a bilateral pair
%                   R = right landmark in a bilateral pair

dummy = readcell(fname, 'NumHeaderLines', 1, ...
    'ExpectedNumVariables', 5, 'Delimiter', ';');
coordinates = cell2mat(dummy(:,1:3));
if isnumeric(dummy{1,4})
    labels = cellfun(@num2str,dummy(:,4),'UniformOutput',false);
else
    labels = dummy(:,4);
end
if size(dummy,2)==5
    type = dummy(:,5);
    if ~all(ismember(type, {'M' 'L' 'R'}))
        error('Unrecognized landmark type.')
    end
    if ~all(contains(lower(labels(type=="L")), 'left'))
        error('All left markers must have the word "left" in their label.')
    end
    if ~all(contains(lower(labels(type=="R")), 'right'))
        error('All right markers must have the word "right" in their label.')
    end
else
    type = [];
end