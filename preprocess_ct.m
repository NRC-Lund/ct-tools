classdef preprocess_ct < matlab.apps.AppBase
%PREPROCESS_CT   A GUI for cropping and re-orienting CT scans.
%
% Use this GUI to crop and re-orient your CT scans to comply with the
% coordinate system used in this toolbox (Rigth, Anterior, Superior). 
% Follow the steps outlined in the GUI, check that the pink labels are 
% correct and save to a new image. All operations are "rigid" in the sense
% that they do not mirror the image (with the exception of "Flip ML axis").
%
% No resampling of the data occurs in the processing steps, so image
% quality is unaffected. 
%
% SYNTAX:
%   preprocess_ct
%   preprocess_ct(img)
%   obj = preprocess_ct
%
% INPUT:
%   img     - (optional) 3-dimensional grayscale image matrix. See
%             documentation for orthosliceViewer for compatible data types.
%
% OUTPUT:
%   obj     - Returns the GUI object (mostly useful for development and
%             debugging).


    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        GridLayout              matlab.ui.container.GridLayout
        LeftPanel               matlab.ui.container.Panel
        DefineaxesButton        matlab.ui.control.Button
        Step1Label              matlab.ui.control.Label
        DefinedirectionsButton  matlab.ui.control.Button
        Step2Label              matlab.ui.control.Label
        HighSliderLabel         matlab.ui.control.Label
        HighSlider              matlab.ui.control.Slider
        LowSliderLabel          matlab.ui.control.Label
        LowSlider               matlab.ui.control.Slider
        IntensityrangeLabel     matlab.ui.control.Label
        Step3Label              matlab.ui.control.Label
        CropXYButton            matlab.ui.control.Button
        CropXZButton            matlab.ui.control.Button
        Step4Label               matlab.ui.control.Label
        FlipButton              matlab.ui.control.Button
        Step5Label              matlab.ui.control.Label
        SaveButton              matlab.ui.control.Button
        LoadimageButton         matlab.ui.control.Button
        Step0Label              matlab.ui.control.Label
        SubsamplingLabel        matlab.ui.control.Label
        SubsamplingSpinner      matlab.ui.control.Spinner
        RightPanel              matlab.ui.container.Panel
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = public)
        Img;            % 3D image
        h_ortho;        % Handle to orthosliceViewer.
        cfg = struct('X', [], 'Y', [], 'Z', []); % Configuration of axes.
        State = 0;      % Which step are we on.
        DisplayRange = [0 100]; % Intensity range in percent.
        ImgMinMax;      % Min and max intensity values of image.
        Subsampling;    % Subsampling factor.
        SubsampledSize; % Size of subsampled image.
        LabelXY;        % Text label for the XY plane.
        LabelZY;        % Text label for the ZY plane.
        LabelXZ;        % Text label for the XZ plane.
        h_progress      % Handle to progress bar.
    end
    
    methods (Access = private)
        
        function UpdateView(app)
            if ~isempty(app.h_ortho) % Delete old viewer if it exists...
                delete(app.h_ortho)  % For some reason this deletes the panel as well...
                % Recreate RightPanel
                app.RightPanel = uipanel(app.GridLayout, "BackgroundColor", app.LeftPanel.BackgroundColor);
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
            % Start orthosliceViewer:
            app.h_ortho = orthosliceViewer(app.Img(1:app.Subsampling:end, 1:app.Subsampling:end, 1:app.Subsampling:end), 'Parent', app.RightPanel);
            % Included to try to mitigate panning problem
            [tmp1, tmp2, tmp3] = app.h_ortho.getAxesHandles;
            disableDefaultInteractivity(tmp1);
            disableDefaultInteractivity(tmp2);
            disableDefaultInteractivity(tmp3);
            % Use our own controls for intensity:
            app.h_ortho.DisplayRangeInteraction = "off";
            app.h_ortho.DisplayRange = app.ImgMinMax.*app.DisplayRange/100;
            % Change text colors:
            h = findall(app.RightPanel, 'Type', 'text');
            [h.Color] = deal([1 1 1]);
            h = findobj(h, 'Tag', 'SliceNumber');
            h.Color = app.h_ortho.CrosshairColor;
            % Add labels:
            [ax1, ax2, ax3] = app.h_ortho.getAxesHandles;
            h_ax = [ax1 ax2 ax3];
            if app.State>1 % If we have assigned labels to dimensions:
                app.LabelXY = text(h_ax(1), h_ax(1).XLim(2), 1, "Axial");
                app.LabelZY = text(h_ax(2), h_ax(2).XLim(2), 1, "Coronal");
                app.LabelXZ = text(h_ax(3), h_ax(3).XLim(2), 1, "Sagital");
                h = [app.LabelXY app.LabelZY app.LabelXZ];
                set(h, 'Color', 'm');
                set(h, 'FontWeight', "bold")
                set(h, 'VerticalAlignment', "top")
                set(h, 'HorizontalAlignment', "right")
                set(h, 'FontSize', 12)
            end
            if app.State>2 % If we have assigned labels to directions:
                for iAx = 1:3
                    arrowlength = min(diff(h_ax(iAx).XLim), diff(h_ax(iAx).YLim)) * 0.1;
                    offset = arrowlength/2;
                    headsize = arrowlength/4;
                    h_line1 = line(h_ax(iAx), 'XData', [h_ax(iAx).XLim(1) h_ax(iAx).XLim(1)] + offset, ...
                        'YData', [h_ax(iAx).YLim(1) h_ax(iAx).YLim(1)+arrowlength] + offset, ...
                        'Color', 'm', 'LineWidth', 1);
                    arrowhead = [h_line1.XData(2) h_line1.YData(2); h_line1.XData(2)+headsize/2 h_line1.YData(2)-headsize; h_line1.XData(2)-headsize/2 h_line1.YData(2)-headsize];
                    h_head1 = patch(h_ax(iAx), arrowhead(:,1), arrowhead(:,2), 'm');
                    h_txt1 = text(h_ax(iAx), h_line1.XData(2), h_line1.YData(2), '', "Color", 'm', "VerticalAlignment","top", "HorizontalAlignment","center");
                    h_line2 = line(h_ax(iAx), 'XData', [h_ax(iAx).XLim(1) h_ax(iAx).XLim(1)+arrowlength] + offset, ...
                        'YData', [h_ax(iAx).YLim(1) h_ax(iAx).YLim(1)] + offset, ...
                        'Color', 'm', 'LineWidth', 1);
                    arrowhead = [h_line2.XData(2) h_line2.YData(2); h_line2.XData(2)-headsize h_line2.YData(2)-headsize/2; h_line2.XData(2)-headsize h_line2.YData(2)+headsize/2];
                    h_head2 = patch(h_ax(iAx), arrowhead(:,1), arrowhead(:,2), 'm');
                    h_txt2 = text(h_ax(iAx), h_line2.XData(2), h_line2.YData(2), '', "VerticalAlignment","middle", "HorizontalAlignment","left");
                    set([h_txt1 h_txt2], "Color", 'm', "FontSize", 12)
                    switch iAx
                        case 1
                            h_txt1.String = "Right";
                            h_txt2.String = "Anterior";
                        case 2
                            h_txt1.String = "Right";
                            h_txt2.String = "Dorsal";
                        case 3
                            h_txt1.String = "Dorsal";
                            h_txt2.String = "Anterior";
                    end
                end
            end
        end
        
        function NewImage(app)
            % Subsample for better responsiveness:
            app.Subsampling = ceil(max(size(app.Img))/500);
            app.SubsamplingSpinner.Value = app.Subsampling;
            app.SubsampledSize = floor(size(app.Img)/app.Subsampling);
            % Limits:
            app.ImgMinMax = double([min(app.Img(:)) max(app.Img(:))]);
            app.DisplayRange = [0 100];
            app.LowSlider.Value = app.DisplayRange(1);
            app.HighSlider.Value = app.DisplayRange(2);
            % Update view:
            UpdateState(app, 1);
            UpdateView(app);
        end
        
            function UpdateState(app, NewState)
            % Enable/disable controls depending on where in the process we are:
            app.State = NewState;
            set([app.CropXYButton app.CropXZButton ...
                app.DefineaxesButton app.DefinedirectionsButton ...
                app.FlipButton ...
                app.SaveButton app.HighSlider app.LowSlider], ...
                'Enable', 'off');
            if app.State>0 % We have loaded an image...
                app.DefineaxesButton.Enable = "on";
                app.CropXYButton.Enable = "on";
                app.CropXZButton.Enable = "on";
                app.HighSlider.Enable = "on";
                app.LowSlider.Enable = "on";
            end
            if app.State>1 % We have assigned dimensions...
                app.DefinedirectionsButton.Enable = "on";
            end
            if app.State>2 % We have assigned directions...
                app.FlipButton.Enable = "on";
                app.SaveButton.Enable = "on";
            end
        end
    end
    
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            narginchk(1,2) % Check that we have the correct number of arguments.
            if ~isempty(varargin) % If we pass an image as an input argument...
                if isnumeric(varargin{1})
                    app.Img = squeeze(varargin{1});
                    if ndims(app.Img)~=3
                        error('The image must have 3 non-singleton dimensions.')
                    end
                elseif ischar(varargin{1})
                    if exist(varargin{1}, 'dir')
                        LoadDicom(app, varargin{1})
                    elseif exist(varargin{1}, 'file')
                        LoadNifti(app, varargin{1})
                    else
                        error('Could not find the data.')
                    end
                else
                    error('Wrong input type.')
                end
                UpdateState(app,1);
                NewImage(app)
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {660, 660};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {147, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Button pushed function: DefineaxesButton
        function DefineaxesButtonPushed(app, event)
            % Select the rostro-caudal axis:
            Unselected = {'X' 'Y' 'Z'};
            answer = uiconfirm(app.UIFigure, 'Please select the anterior-posterior axis:', 'Select axes',...
                "Options", {Unselected{1}, Unselected{2}, Unselected{3} 'Cancel'}, "DefaultOption", 1, "CancelOption", 4);
            if answer=="Cancel"
                return
            end
            app.cfg.(answer) = 'Anterior-posterior';
            % Select the medio-lateral axis:
            Unselected = setdiff(Unselected, answer);
            answer = uiconfirm(app.UIFigure, 'Please select the medial-lateral axis:', 'Select axes',...
                "Options", {Unselected{1}, Unselected{2} 'Cancel'}, "DefaultOption", 1, "CancelOption", 3);
            if answer=="Cancel"
                return
            end
            app.cfg.(answer) = 'Medial-lateral';
            % Set dorso-ventral axis:
            Unselected = setdiff(Unselected, answer);
            app.cfg.(Unselected{1}) = 'Dorsal-ventral';
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Permuting image...', ...
                    'Indeterminate','on');
            % Permute:
            standard = {'Medial-lateral' 'Anterior-posterior' 'Dorsal-ventral'}; % Y,X,Z!
            [~,ix] = ismember(standard,{app.cfg.Y app.cfg.X app.cfg.Z});
            app.Img = permute(app.Img, ix);
            % Is this permutation a reflection?
            neg = [1 3 2; 2 1 3; 3 2 1]; % All permutations that will result in a negative coordinate system.
            if ismember(ix, neg, 'rows')
                disp('Flipping...   ')
                app.Img = flip(app.Img,1); % The permutation resulted in a negative system, which means that the image is mirrored. Let's flip the ML axis.
            end
            % We are now done with step 1:
            UpdateState(app,2);
            % Update view:
            UpdateView(app)
            delete(app.h_progress)
        end

        % Button pushed function: DefinedirectionsButton
        function DefinedirectionsButtonPushed(app, event)
            % Where is rostral?
            answer1 = uiconfirm(app.UIFigure, 'In the sagital view, anterior is to the...', 'Select orientation',...
                "Options", {'Left' 'Right' 'Cancel'}, "DefaultOption", 2, "CancelOption", 3);
            if answer1=="Cancel"
                return
            end
            % Where is dorsal?
            answer2 = uiconfirm(app.UIFigure, 'In the sagital view, dorsal is ...', 'Select orientation',...
                "Options", {'Up' 'Down' 'Cancel'}, "DefaultOption", 2, "CancelOption", 3);
            if answer2=="Cancel"
                return
            end
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Flipping image...', ...
                    'Indeterminate','on');
            % Perform a 180� rotation around the dorso-ventral axis if necessary:
            if answer1=="Left"
                app.Img = flip(app.Img,2); % Flip the x-axis.
                app.Img = flip(app.Img,1); % Flip the y-axis to keep handedness.
            end
            % Perform a 180� rotation around the rostro-caudal axis if necessary:
            if answer2=="Up"
                app.Img = flip(app.Img,3); % Flip the z-axis.
                app.Img = flip(app.Img,1); % Flip the y-axis to keep handedness.
            end
            % We are now done with step 2:
            UpdateState(app,3);
            % Update view:
            UpdateView(app)
            delete(app.h_progress)
        end
        
        % Button pushed function: DefinedirectionsButton
        function FlipButtonPushed(app, event)
            answer = uiconfirm(app.UIFigure, 'Are you sure that you want to flip the image? This will result in a mirror image.', 'Flip?',...
                "Options", {'Yes' 'Cancel'}, "DefaultOption", 2, "CancelOption", 2);
            if answer=="Yes"
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Flipping image...', ...
                    'Indeterminate','on');
                app.Img = flip(app.Img,1);
                % Update view:
                UpdateView(app)
                delete(app.h_progress)
            end
        end

        % Value changed function: HighSlider, LowSlider
        function LimitSliderValueChanged(app, event)
            app.DisplayRange = sort([app.LowSlider.Value app.HighSlider.Value]);
            app.LowSlider.Value = app.DisplayRange(1); app.HighSlider.Value = app.DisplayRange(2);
            app.h_ortho.DisplayRange = app.ImgMinMax(2)*app.DisplayRange/100;
        end

        % Button pushed function: CropXYButton
        function CropXYButtonPushed(app, event)
            app.h_ortho.CrosshairEnable = 'off';
            [ax1, ax2, ax3] = app.h_ortho.getAxesHandles;
            roi = drawrectangle(ax1);
            pos = round(roi.Position);
            ix_Y = [pos(2) pos(2)+pos(4)];  % Position in subsampled image...
            ix_X = [pos(1) pos(1)+pos(3)];
            ix_Y2 = ix_Y * app.Subsampling; % Position in original image...
            ix_X2 = ix_X * app.Subsampling;
            ix_X2 = [max(ix_X2(1),1) min(ix_X2(2),size(app.Img,2))]; % Assert range
            ix_Y2 = [max(ix_Y2(1),1) min(ix_Y2(2),size(app.Img,1))]; % Assert range
            delete(roi)
            answer = uiconfirm(app.UIFigure, 'Are you sure that you want to crop the image?', 'Crop?',...
                "Options", {'Yes' 'Cancel'}, "DefaultOption", 2, "CancelOption", 2);
            if answer=="Yes"
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Cropping image...', ...
                    'Indeterminate','on');
                app.Img = app.Img(ix_Y2(1):ix_Y2(2),ix_X2(1):ix_X2(2),:);
                % Update view:
                UpdateView(app)
                delete(app.h_progress)
            else
                app.h_ortho.CrosshairEnable = 'on';
            end
        end

        % Button pushed function: CropXZButton
        function CropXZButtonPushed(app, event)
            app.h_ortho.CrosshairEnable = 'off';
            [ax1, ax2, ax3] = app.h_ortho.getAxesHandles;
            roi = drawrectangle(ax3);
            pos = round(roi.Position);
            ix_Z = [pos(2) pos(2)+pos(4)];  % Position in subsampled image...
            ix_X = [pos(1) pos(1)+pos(3)];
            ix_Z2 = ix_Z * app.Subsampling; % Position in original image...
            ix_X2 = ix_X * app.Subsampling;
            ix_Z2 = [max(ix_Z2(1),1) min(ix_Z2(2),size(app.Img,3))]; % Assert range
            ix_X2 = [max(ix_X2(1),1) min(ix_X2(2),size(app.Img,2))]; % Assert range
            delete(roi)
            answer = uiconfirm(app.UIFigure, 'Are you sure that you want to crop the image?', 'Crop?',...
                "Options", {'Yes' 'Cancel'}, "DefaultOption", 2, "CancelOption", 2);
            if answer=="Yes"
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Cropping image...', ...
                    'Indeterminate','on');
                app.Img = app.Img(:,ix_X2(1):ix_X2(2),ix_Z2(1):ix_Z2(2),:);
                % Update view:
                UpdateView(app)
                delete(app.h_progress)
            else
                app.h_ortho.CrosshairEnable = 'on';
            end
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            [file, path] = uiputfile({'*.nii'});
            if isequal(file,0)
                return
            else
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Saving image...', ...
                    'Indeterminate','on');
                niftiwrite(app.Img, fullfile(path, file))
                disp('File saved!')
                delete(app.h_progress)
            end
        end

        % Button pushed function: LoadimageButton
        function LoadimageButtonPushed(app, event)
            answer = uiconfirm(app.UIFigure, 'Please select image format:', 'Load image',...
                "Options", {'DICOM' 'NIfTI' 'Cancel'}, "DefaultOption", 1, "CancelOption", 3);
            switch answer
                case 'DICOM'
                    dname = uigetdir('', 'Select DICOM folder');
                    if isequal(dname,0)
                        return
                    end
                    LoadDicom(app, dname)
                case 'NIfTI'
                    [file,path] = uigetfile({'*.nii'}, 'Select NIfTI file');
                    if isequal(file,0)
                        return
                    end
                    LoadNifti(app, fullfile(path, file))
                case 'Cancel'
                    return
            end
            NewImage(app)
        end
        
        function LoadDicom(app, dname)
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Loading image...', ...
                'Indeterminate','on');
            app.Img = squeeze(dicomreadVolume(dname));
            app.Img = permute(app.Img, [2 1 3]); % This flip was introduced after looking at a reference object scanned in Umeå.
            delete(app.h_progress)
        end
        
        function LoadNifti(app, fname)
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Loading image...', ...
                'Indeterminate','on');
            info = niftiinfo(fname);
            app.Img = squeeze(niftiread(fname));
            if info.Transform.T(1,1)<0
                disp('Flipping 1st dimension due to negative scaling factor...')
                app.Img = flip(app.Img,1);
            end
            if info.Transform.T(2,2)<0
                disp('Flipping 2nd dimension due to negative scaling factor...')
                app.Img = flip(app.Img,2);
            end
            if info.Transform.T(3,3)<0
                disp('Flipping 3rd dimension due to negative scaling factor...')
                app.Img = flip(app.Img,3);
            end
            delete(app.h_progress)
        end

        % Value changed function: SubsamplingSpinner
        function SubsamplingSpinnerValueChanged(app, event)
            app.Subsampling = app.SubsamplingSpinner.Value;
            app.SubsampledSize = floor(size(app.Img)/app.Subsampling);
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Changing display resolution...', ...
                    'Indeterminate','on');
            UpdateView(app)
            delete(app.h_progress)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.149 0.149 0.149];
            app.UIFigure.Position = [100 100 1078 660];
            app.UIFigure.Name = 'preprocess_ct';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {147, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.ForegroundColor = [0.651 0.651 0.651];
            app.LeftPanel.BackgroundColor = [0.149 0.149 0.149];
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create Step0Label
            app.Step0Label = uilabel(app.LeftPanel);
            app.Step0Label.FontColor = [1 1 1];
            app.Step0Label.Position = [24 621 100 22];
            app.Step0Label.Text = 'Step 0:';
            
            % Create LoadimageButton
            app.LoadimageButton = uibutton(app.LeftPanel, 'push');
            app.LoadimageButton.ButtonPushedFcn = createCallbackFcn(app, @LoadimageButtonPushed, true);
            app.LoadimageButton.BackgroundColor = [0.502 0.502 0.502];
            app.LoadimageButton.FontColor = [1 1 1];
            app.LoadimageButton.Tooltip = {'Load NIfTI image from file.'};
            app.LoadimageButton.Position = [24 600 100 22];
            app.LoadimageButton.Text = 'Load image';
            
            % Create Step1Label
            app.Step1Label = uilabel(app.LeftPanel);
            app.Step1Label.FontColor = [1 1 1];
            app.Step1Label.Position = [24 564 100 22];
            app.Step1Label.Text = 'Step 1 (optional):';

            % Create CropXYButton
            app.CropXYButton = uibutton(app.LeftPanel, 'push');
            app.CropXYButton.ButtonPushedFcn = createCallbackFcn(app, @CropXYButtonPushed, true);
            app.CropXYButton.BackgroundColor = [0.502 0.502 0.502];
            app.CropXYButton.FontColor = [1 1 1];
            app.CropXYButton.Enable = 'off';
            app.CropXYButton.Tooltip = {'Crop image in the XY plane.'};
            app.CropXYButton.Position = [24 543 99 22];
            app.CropXYButton.Text = 'Crop XY';

            % Create CropXZButton
            app.CropXZButton = uibutton(app.LeftPanel, 'push');
            app.CropXZButton.ButtonPushedFcn = createCallbackFcn(app, @CropXZButtonPushed, true);
            app.CropXZButton.BackgroundColor = [0.502 0.502 0.502];
            app.CropXZButton.FontColor = [1 1 1];
            app.CropXZButton.Enable = 'off';
            app.CropXZButton.Tooltip = {'Crop image in the XZ plane.'};
            app.CropXZButton.Position = [24 518 99 22];
            app.CropXZButton.Text = 'Crop XZ';
            
            % Create Step2Label
            app.Step2Label = uilabel(app.LeftPanel);
            app.Step2Label.FontColor = [1 1 1];
            app.Step2Label.Position = [24 483 100 22];
            app.Step2Label.Text = 'Step 2:';
            
            % Create DefineaxesButton
            app.DefineaxesButton = uibutton(app.LeftPanel, 'push');
            app.DefineaxesButton.ButtonPushedFcn = createCallbackFcn(app, @DefineaxesButtonPushed, true);
            app.DefineaxesButton.BackgroundColor = [0.502 0.502 0.502];
            app.DefineaxesButton.FontColor = [1 1 1];
            app.DefineaxesButton.Enable = 'off';
            app.DefineaxesButton.Tooltip = {'Assign image axes to anatomical axes. The resulting transormation can be factorized into 90-degree rotations (i.e. no reflection or resampling occurs).'};
            app.DefineaxesButton.Position = [23 462 100 22];
            app.DefineaxesButton.Text = 'Define axes';
            
            % Create Step3Label
            app.Step3Label = uilabel(app.LeftPanel);
            app.Step3Label.FontColor = [1 1 1];
            app.Step3Label.Position = [24 426 100 22];
            app.Step3Label.Text = 'Step 3:';
            
            % Create DefinedirectionsButton
            app.DefinedirectionsButton = uibutton(app.LeftPanel, 'push');
            app.DefinedirectionsButton.ButtonPushedFcn = createCallbackFcn(app, @DefinedirectionsButtonPushed, true);
            app.DefinedirectionsButton.BackgroundColor = [0.502 0.502 0.502];
            app.DefinedirectionsButton.FontColor = [1 1 1];
            app.DefinedirectionsButton.Enable = 'off';
            app.DefinedirectionsButton.Tooltip = {'Define image directions so that Anterior/Right/Dorsal points in the positive direction.'};
            app.DefinedirectionsButton.Position = [21 405 106 22];
            app.DefinedirectionsButton.Text = 'Define directions';

            % Create Step4Label
            app.Step4Label = uilabel(app.LeftPanel);
            app.Step4Label.FontColor = [1 1 1];
            app.Step4Label.Position = [24 369 110 22];
            app.Step4Label.Text = 'Mirror (if necessary):';
            
            % Create FlipButton
            app.FlipButton = uibutton(app.LeftPanel, 'push');
            app.FlipButton.ButtonPushedFcn = createCallbackFcn(app, @FlipButtonPushed, true);
            app.FlipButton.BackgroundColor = [0.502 0.502 0.502];
            app.FlipButton.FontColor = [1 1 1];
            app.FlipButton.Enable = 'off';
            app.FlipButton.Tooltip = {'Flip the medial-lateral axis. WARNING: This will create a mirror image!'};
            app.FlipButton.Position = [21 348 106 22];
            app.FlipButton.Text = 'Flip ML axis';
            
            % Create Step5Label
            app.Step5Label = uilabel(app.LeftPanel);
            app.Step5Label.FontColor = [1 1 1];
            app.Step5Label.Position = [24 312 100 22];
            app.Step5Label.Text = 'Step 4:';

            % Create SaveButton
            app.SaveButton = uibutton(app.LeftPanel, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.BackgroundColor = [0.502 0.502 0.502];
            app.SaveButton.FontColor = [1 1 1];
            app.SaveButton.Enable = 'off';
            app.SaveButton.Tooltip = {'Save image to NIfTI file. First check that the pink labels indicating anatomical planes and directions are correct.'};
            app.SaveButton.Position = [24 291 99 22];
            app.SaveButton.Text = 'Save';

            % Create SubsamplingLabel
            app.SubsamplingLabel = uilabel(app.LeftPanel);
            app.SubsamplingLabel.HorizontalAlignment = 'right';
            app.SubsamplingLabel.FontColor = [1 1 1];
            app.SubsamplingLabel.Position = [10 121 75 22];
            app.SubsamplingLabel.Text = 'Subsampling:';

            % Create SubsamplingSpinner
            app.SubsamplingSpinner = uispinner(app.LeftPanel);
            app.SubsamplingSpinner.RoundFractionalValues = 'on';
            app.SubsamplingSpinner.ValueChangedFcn = createCallbackFcn(app, @SubsamplingSpinnerValueChanged, true);
            app.SubsamplingSpinner.FontColor = [1 1 1];
            app.SubsamplingSpinner.BackgroundColor = [0.502 0.502 0.502];
            app.SubsamplingSpinner.Tooltip = {'Set subsampling factor of displayed image.'};
            app.SubsamplingSpinner.Position = [92 121 48 22];

            % Create IntensityrangeLabel
            app.IntensityrangeLabel = uilabel(app.LeftPanel);
            app.IntensityrangeLabel.FontColor = [1 1 1];
            app.IntensityrangeLabel.Position = [11 80 88 22];
            app.IntensityrangeLabel.Text = 'Intensity range:';
            
            % Create HighSliderLabel
            app.HighSliderLabel = uilabel(app.LeftPanel);
            app.HighSliderLabel.HorizontalAlignment = 'right';
            app.HighSliderLabel.FontColor = [1 1 1];
            app.HighSliderLabel.Position = [7 50 30 22];
            app.HighSliderLabel.Text = 'High';

            % Create HighSlider
            app.HighSlider = uislider(app.LeftPanel);
            app.HighSlider.MajorTicks = [];
            app.HighSlider.ValueChangedFcn = createCallbackFcn(app, @LimitSliderValueChanged, true);
            app.HighSlider.MinorTicks = [];
            app.HighSlider.Enable = 'off';
            app.HighSlider.Tooltip = {'Set intensity range.'};
            app.HighSlider.FontColor = [1 1 1];
            app.HighSlider.Position = [58 59 73 3];

            % Create LowSliderLabel
            app.LowSliderLabel = uilabel(app.LeftPanel);
            app.LowSliderLabel.HorizontalAlignment = 'right';
            app.LowSliderLabel.FontColor = [1 1 1];
            app.LowSliderLabel.Position = [8 15 29 22];
            app.LowSliderLabel.Text = 'Low';

            % Create LowSlider
            app.LowSlider = uislider(app.LeftPanel);
            app.LowSlider.MajorTicks = [];
            app.LowSlider.ValueChangedFcn = createCallbackFcn(app, @LimitSliderValueChanged, true);
            app.LowSlider.MinorTicks = [];
            app.LowSlider.Enable = 'off';
            app.LowSlider.Tooltip = {'Set intensity range.'};
            app.LowSlider.FontColor = [1 1 1];
            app.LowSlider.Position = [58 24 73 3];
            
            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.BackgroundColor = [0.149 0.149 0.149];
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = preprocess_ct(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end