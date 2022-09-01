function writeLandmarks(fname, coordinates, labels)
%WRITELANDMARKS   Function to write landmarks to a text file.
%
% SYNTAX:
%   writeLandmarks(fname, coordinates, labels)
%
% INPUT:
%   fname       - Filename.
%   coordinates - Landmark coordinates. Coordinates are slice numbers 
%                 defined as ML,AP,DV with positive directions towards 
%                 right, anterior, superior, respectively.
%   labels      - Cell array with landmark names.

fid = fopen(fname,'w');
fprintf(fid, '%s\n', '% Coordinates are slice numbers defined as ML;AP;DV with positive directions towards anterior, right, superior.');
for iL = 1:numel(labels)
    fprintf(fid, '%d;%d;%d;%s\n', coordinates(iL,:), labels{iL});
end
fclose(fid);