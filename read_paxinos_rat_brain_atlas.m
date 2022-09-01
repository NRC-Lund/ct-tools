function S = read_paxinos_rat_brain_atlas(ml,ap,dv,edition)
% SYNTAX:
%   S = read_paxinos_rat_brain_atlas(ml,ap,dv)
%   S = read_paxinos_rat_brain_atlas(ml,ap,dv,edition)
%
% INPUT:
%   ml      - Mediolateral coordinate in mm. 0 is the midline, positive is
%             right.
%   ap      - Anteroposterior coordinate in mm. 0 is bregma, positive is
%             anterior.
%   dv      - Dorsoventral coordinate in mm. 0 is bregma, positive is
%             dorsal.
%   edition - (optional scalar, default=6). Select edition 5 or 6. Edition
%             5 is read from the web, edition 6 is read from file.
%
% OUTPUT:
%   S       - Structure with the fields 'coronal', 'sagittal' and
%             'horizontal'. (For the 6th edition only the coronal sections 
%             are included.) Each of the section fields contains an image
%             and a reference object RA. The reference object can be used
%             to convert between atlas coordinates and image coordinates.
%
% EXAMPLE:
%   To get the image coordinate of ML=5, DV=4, do
%   [x,y] = S.coronal.RA.worldToSubscript(5,4)
%
% REFERENCE:
%   - 5th edition:
%     Matt Gaidica (2020). mattgaidica/RatBrainAtlasAPI (https://www.github.com/mattgaidica/RatBrainAtlasAPI), GitHub. Retrieved February 13, 2020.
%     Atlas Source: Paxinos, George, and Charles Watson. The rat brain in stereotaxic coordinates: hard cover edition. Access Online via Elsevier, 2006.
%
%   - 6th edition:
%﻿    Paxinos, G., and Watson, C. (2007). The rat brain in stereotaxic coordinates (6th ed.). Amsterdam ; Boston ; Academic Press/Elsevier.


if nargin==3 || isempty(edition)
    edition = 6;
end
    
switch edition
    case 5
        url  = atlasUrl(ml,ap,dv);
        try
            S = webread(url);
        catch ME
            error('Unable to complete web request.');
            S = [];
        end
        S = jsondecode(S);
        S.coronal.image = webread(S.coronal.image_url,'ContentType','image');
        S.sagittal.image = webread(S.sagittal.image_url,'ContentType','image');
        S.horizontal.image = webread(S.horizontal.image_url,'ContentType','image');
        % supply TRUE to make marking
        if license('test', 'image_toolbox')
            S.coronal.image_marked = insertShape(S.coronal.image,'FilledCircle',[S.coronal.left S.coronal.top 10],'Color','r','Opacity',1);
            S.sagittal.image_marked = insertShape(S.sagittal.image,'FilledCircle',[S.sagittal.left S.sagittal.top 10],'Color','r','Opacity',1);
            S.horizontal.image_marked = insertShape(S.horizontal.image,'FilledCircle',[S.horizontal.left S.horizontal.top 10],'Color','r','Opacity',1);
        end
    case 6
        pathname = fileparts(mfilename('fullpath'));
        coronal = load(fullfile(pathname, 'Atlases', 'Paxinos_and_Watson_2007_The_Rat_Brain_6th_edition', ...
            'Paxinos_coronal_6th_edition.mat'));
        ix = interp1(coronal.T.AP, 1:height(coronal.T), ap, 'nearest');
        if isnan(ix) || ix<1 || ix>height(coronal.T)
            S = [];
            return
        end
        S.coronal.image = coronal.Img(:,:,ix);
        S.coronal.RA = coronal.RA(ix);
        S.coronal.image = flip(S.coronal.image,1);
        S.coronal.RA.YWorldLimits = -flip(S.coronal.RA.YWorldLimits);
        [S.coronal.top, S.coronal.left] = S.coronal.RA.worldToSubscript(ml,-dv);
        if isnan(S.coronal.top)
            S = [];
            return
        end
        S.coronal.image_marked = insertShape(S.coronal.image,'FilledCircle',[S.coronal.left S.coronal.top 10],'Color','r','Opacity',1);
end



function url  = atlasUrl(ml,ap,dv)
api = 'http://labs.gaidi.ca/rat-brain-atlas/api.php?';
url = [api 'ml=' num2str(ml) '&ap=' num2str(ap) '&dv=' num2str(dv)];