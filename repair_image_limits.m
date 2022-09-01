function ImgOut = repair_image_limits(ImgIn)
%REPAIR_IMAGE_LIMITS Fix limits of float images
%
% Fix limits of float images (they are not supposed to be outside [0 1]).
%
% SYNTAX:
%   ImgOut = repair_image_limits(ImgIn)
%
ImgOut = ImgIn;
if isfloat(ImgOut) % Floaty images must be [0 1].
    MinLim = min(ImgOut(:));
    if MinLim<0
        ImgOut = ImgOut - MinLim;
        MaxLim = max(ImgOut(:));
        if MaxLim>1
            ImgOut = ImgOut / MaxLim;
        end
    end
end