classdef register_landmarks < matlab.apps.AppBase
%REGISTER_LANDMARKS   A GUI for registering landmarks in a CT image.
%
% Use this GUI to register predefined landmarks in a CT image. The default
% landmarks are defined in Atlas/atlas-landmarks.txt. The CT image should
% first have been processed with the preprocess_ct GUI to assure correct
% coordinate system orientation. The landmarks can be saved to a text file.
% The estimated affine transform from image to atlas is stored in the
% property app.tform.
%
% SYNTAX:
%   register_landmarks
%   register_landmarks(img)
%   obj = register_landmarks
%
% INPUT:
%   img     - (optional) 3-dimensional grayscale image matrix conforming
%             with the preprocess_ct conventions (RAS).
%
% OUTPUT:
%   obj     - Returns the GUI object (mostly useful for development and
%             debugging).

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        LoadimageMenu                  matlab.ui.container.Menu
        LoadlandmarksMenu              matlab.ui.container.Menu
        SavelandmarksMenu              matlab.ui.container.Menu
        LoadtransformationMenu         matlab.ui.container.Menu
        SavetransformationMenu         matlab.ui.container.Menu
        TransformMenu                  matlab.ui.container.Menu
        RegistrationTargetMenu         matlab.ui.container.Menu
        RegistrationTargetSubMenu      matlab.ui.container.Menu
        RegistrationMethodMenu         matlab.ui.container.Menu
        ModifiedHornsMethodMenu        matlab.ui.container.Menu
        StandardHornsMethodMenu        matlab.ui.container.Menu
        SVDMethodMenu                  matlab.ui.container.Menu
        LambdaBregmaMethodMenu         matlab.ui.container.Menu
        CalculateLandmarksMenu         matlab.ui.container.Menu
        GridLayout                     matlab.ui.container.GridLayout
        LeftPanel                      matlab.ui.container.Panel
        TransformationPanel            matlab.ui.container.Panel
        PreviewregistrationButton      matlab.ui.control.Button
        CalculatetransformationButton  matlab.ui.control.Button
        ViewfusedimageButton           matlab.ui.control.Button
        AdjusttransformationButton     matlab.ui.control.Button
        ViewPanel                      matlab.ui.container.Panel
        CursorLabel                    matlab.ui.control.Label
        UITable                        matlab.ui.control.Table
        PaxinosButton                  matlab.ui.control.Button
        SubsamplingLabel               matlab.ui.control.Label
        SubsamplingSpinner             matlab.ui.control.Spinner
        IntensityrangeLabel            matlab.ui.control.Label
        HighSliderLabel                matlab.ui.control.Label
        HighSlider                     matlab.ui.control.Slider
        LowSliderLabel                 matlab.ui.control.Label
        LowSlider                      matlab.ui.control.Slider
        LandmarksPanel                 matlab.ui.container.Panel
        JumptolandmarkButton           matlab.ui.control.Button
        UpdatelandmarkButton           matlab.ui.control.Button
        APLabel                        matlab.ui.control.Label
        APEditField_2                  matlab.ui.control.NumericEditField
        MLEditField_2                  matlab.ui.control.NumericEditField
        DVEditField_2                  matlab.ui.control.NumericEditField
        DeletelandmarkButton           matlab.ui.control.Button
        MLLabel                        matlab.ui.control.Label
        DVLabel                        matlab.ui.control.Label
        SelectlandmarkDropDownLabel    matlab.ui.control.Label
        SelectlandmarkDropDown         matlab.ui.control.DropDown
        RightPanel                     matlab.ui.container.Panel
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = public)
        Img;            % Native image. dim1=ML, dim2=AP, dim3=DV. Right-Anterior-Superior.
        h_ortho;        % Handle to orthosliceViewer. Note that h_ortho.SliceNumber is [AP,ML,DV].
        State_ImageExists = false;
        State_LandmarksExist = false;
        State_TransformExists = false;
        State_WaxholmT2Loaded = false;
        DisplayRange = [0 100]; % Intensity range in percent.
        ImgMinMax;      % Min and max intensity values of image.
        Subsampling;    % Subsampling factor.
        SubsampledSize; % Size of subsampled image.
        h_progress      % Handle to progress bar.
        PreferredDir = cd;   % Preferred directory.
        Filename_T2 = 'WHS_SD_rat_T2star_v1.01.nii'; % REMOVE
        Pathname_WDS = fullfile(fileparts(mfilename('fullpath')), 'Atlases', 'WHS_SD_rat_v1.01'); % REMOVE
        Pathname_Pax = fullfile(fileparts(mfilename('fullpath')), 'Atlases', 'Paxinos_and_Watson_2007_The_Rat_Brain_6th_edition'); % REMOVE
        WaxholmImg        % REMOVE! Atlas T2 image (WHS SDr v1.01). dim1=ML, dim2=AP, dim3=DV. Right-Anterior-Superior.
        LandmarkLabels = {}; % Landmark labels.
        Landmarks;        % Landmark slice numbers in the native image. dim1=ML, dim2=AP, dim3=DV. Right-Anterior-Superior.
        LandmarkTypes = {};
        AtlasLandmarks;
        WaxholmSize = [512 1024 512]; % REMOVE!
        tform;            % affine3d object defining the affine transform to atlas space (slice numbers in ML,AP,DV).
        Atlases;
        RegistrationTarget = '';
        RegistrationMethod = 'Modified Horn';
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
            % Calculate image and start orthoSliceViewer:
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
            % Add listener:
            addlistener(app.h_ortho,'CrosshairMoving',@app.CrosshairMoved);
            addlistener(app.h_ortho,'CrosshairMoved',@app.CrosshairMoved);
            CrosshairMoved(app)
            % Add labels:
            [ax1, ax2, ax3] = app.h_ortho.getAxesHandles;
            h_ax = [ax1 ax2 ax3];
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
            app.State_ImageExists = true;
            UpdateState(app);
            UpdateView(app);
        end
        
        function LoadAtlasDefinitions(app)
            % Select atlas:
            if isempty(app.RegistrationTarget)
                ix = 1;
                app.RegistrationTarget = app.Atlases(ix).name;
            else
                ix = ismember({app.Atlases.name}, app.RegistrationTarget);
            end
            AtlasDir = fullfile(fileparts(mfilename('fullpath')), 'Atlases', app.Atlases(ix).dir);
            % Load landmarks:
            fname = fullfile(AtlasDir,app.Atlases(ix).landmarks);
            fprintf('Loading landmarks from %s...', fname)
            [Landmarks, LandmarkLabels] = readLandmarks(fname);
            app.AtlasLandmarks = NaN(numel(app.LandmarkLabels),3);
            for iCo = 1:size(Landmarks,1)
                b = ismember(app.LandmarkLabels, LandmarkLabels{iCo});
                if any(b)
                    app.AtlasLandmarks(b,:) = Landmarks(iCo,:);
                end
            end
            fprintf('done.\n')
            % Update menu:
            set(app.RegistrationTargetSubMenu, 'Checked', 'off')
            app.RegistrationTargetSubMenu(...
                ismember({app.RegistrationTargetSubMenu.Text}, app.RegistrationTarget)).Checked = 'on';
            if app.RegistrationTarget=="Paxinos and Watson: The Rat Brain (6th edition)" || ...
               app.RegistrationTarget=="Franklin and Paxinos: The Mouse Brain (3rd edition)" % FIX
                app.LambdaBregmaMethodMenu.Enable = 'on';
            else
                app.LambdaBregmaMethodMenu.Enable = 'off';
            end
        end
        
        function LoadWaxholm(app)
            fprintf('Loading Waxholm atlas...')
            % T2:
            fname = fullfile(app.Pathname_WDS, app.Filename_T2);
            if exist(fname,"file")
                T2_exists = true;
            else
                answer = uiconfirm(app.UIFigure, ['Missing file ' fname], 'Missing T2 image',...
                    "Options", {'Browse' 'Skip'}, "DefaultOption", 2);
                switch answer
                    case 'Browse'
                        [file,path] = uigetfile({fname}, 'Select NIfTI file', ...
                            fileparts(mfilename('fullpath')));
                        if isequal(file,0) || ~exist(fullfile(path,file),"file")
                            T2_exists = false;
                        else
                            app.Pathname_WDS = path;
                            app.Filename_T2 = file;
                            T2_exists = true;
                        end
                    case 'Skip'
                        T2_exists = false;
                end
            end
            if T2_exists
                app.WaxholmImg = niftiread(fullfile(app.Pathname_WDS, app.Filename_T2));
                fprintf('done.\n')
                app.State_WaxholmT2Loaded = true;
            else
                fprintf('Failed to load T2 image.\n')
            end
        end
        
        function UpdateState(app)
            % Enable/disable controls depending on the state.
            
            % We have loaded an image:
            if app.State_ImageExists 
                app.HighSlider.Enable = "on";
                app.LowSlider.Enable = "on";
                app.UITable.Enable = "on";
                app.SubsamplingSpinner.Enable = "on";
                app.UpdatelandmarkButton.Enable = "on";
            else
                app.HighSlider.Enable = "off";
                app.LowSlider.Enable = "off";
                app.UITable.Enable = "off";
                app.SubsamplingSpinner.Enable = "off";
                app.UpdatelandmarkButton.Enable = "off";
            end
            
            % We have a transformation:
            if app.State_TransformExists 
                app.SavetransformationMenu.Enable = "on";
                app.CalculateLandmarksMenu.Enable = "on";
                app.AdjusttransformationButton.Enable = "on";
            else
                app.SavetransformationMenu.Enable = "off";
                app.CalculateLandmarksMenu.Enable = "off";
                app.AdjusttransformationButton.Enable = "off";
            end
            
            % We have landmarks:
            if app.State_LandmarksExist
                app.SavelandmarksMenu.Enable = "on";
                app.CalculatetransformationButton.Enable = "on";
            else
                app.SavelandmarksMenu.Enable = "off";
                app.CalculatetransformationButton.Enable = "off";
            end
            
            % We have image and transformation:
            if app.State_ImageExists && app.State_TransformExists
                app.ViewfusedimageButton.Enable = "on";
                app.PaxinosButton.Enable = "on";
            else
                app.ViewfusedimageButton.Enable = "off";
                app.PaxinosButton.Enable = "off";
            end
            
            % We have image and landmarks:
            if app.State_ImageExists && app.State_LandmarksExist
                app.JumptolandmarkButton.Enable = "on";
            else
                app.JumptolandmarkButton.Enable = "off";
            end
            
            % We have transformation and landmarks:
            if app.State_TransformExists && app.State_LandmarksExist
                app.PreviewregistrationButton.Enable = "on";
            else
                app.PreviewregistrationButton.Enable = "off";
            end
            
        end
        
        function UpdateCursorTable(app, ix)
            if isempty(app.h_ortho)
                return
            end
            % ix is the slice number in the full image (not subsampled)
            if nargin==1
                ix = (app.h_ortho.SliceNumbers-1) * app.Subsampling + 1; % AP,ML,DV
                ix = [ix(2) ix(1) ix(3)]; % ML,AP,DV
            end
            if app.State_TransformExists % If we have a tranform...
                co_atlas = app.tform.transformPointsForward(ix);
            else
                co_atlas = NaN(1,3);
            end
            % Table:
            DataCell = num2cell([ix; co_atlas]);
            DataCell(1,:) = cellfun(@(x) num2str(x,'%d'), DataCell(1,:), 'UniformOutput', false);
            DataCell(2,:) = cellfun(@(x) num2str(x,'%0.2f'), DataCell(2,:), 'UniformOutput', false);
            app.UITable.Data = DataCell;
        end
        
        function UpdateLandmarkTable(app)
            pos = app.Landmarks(strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels),:); % ML,AP,DV
            bNan = isnan(pos);
            if any(bNan)
                app.JumptolandmarkButton.Enable = "off";
                app.DeletelandmarkButton.Enable = "off";
            else
                app.JumptolandmarkButton.Enable = "on";
                app.DeletelandmarkButton.Enable = "on";
            end
            h = [app.MLEditField_2 app.APEditField_2 app.DVEditField_2];
            for ih = 1:3
                if bNan(ih)
                    h(ih).Value = 0;
                    h(ih).FontColor = [1 0 0];
                else
                    h(ih).Value = pos(ih);
                    h(ih).FontColor = [1 1 1];
                end
            end
        end
        
        function UpdateSelectLandmarkDropDown(app)
            dummy = repmat({' '}, numel(app.LandmarkLabels), 1);
            dummy(all(~isnan(app.Landmarks),2)) = {'*'};
            Items = cellfun(@(x,y) [x ' ' y], dummy, app.LandmarkLabels, 'UniformOutput', false);
            if isempty(app.SelectlandmarkDropDown.Value)
                PreviousValue = Items{1}(3:end);
            else
                PreviousValue = app.SelectlandmarkDropDown.Value(3:end);
            end
            app.SelectlandmarkDropDown.Items = Items;
            app.SelectlandmarkDropDown.Value = ...
                app.SelectlandmarkDropDown.Items(ismember(app.LandmarkLabels,PreviousValue));
        end
        
        function UpdateTransform(app)
            switch app.RegistrationMethod
                case 'Modified Horn'
                    bNan = any(isnan(app.Landmarks),2) | any(isnan(app.AtlasLandmarks),2);
                    A = app.Landmarks(~bNan,:);
                    B = app.AtlasLandmarks(~bNan,:);
                    if size(A,1)<2
                        h_error = errordlg('Not enough landmarks', ...
                            'Error', 'modal');
                        uiwait(h_error)
                        return
                    end
                    [tform, s, ErrorStats] = absor_anisotropic(A,B);
                    fprintf('Scaling: ML=%0.4f, AP=%0.4f, DV=%0.4f\n', s(1), s(2), s(3))
                    fprintf('Error: errlsq=%0.2f, errmax=%0.2f\n', ErrorStats.errlsq, ErrorStats.errmax)
                    app.tform = affine3d(tform');
                case 'Standard Horn'
                    bNan = any(isnan(app.Landmarks),2) | any(isnan(app.AtlasLandmarks),2);
                    A = app.Landmarks(~bNan,:);
                    B = app.AtlasLandmarks(~bNan,:);
                    if size(A,1)<2
                        h_error = errordlg('Not enough landmarks', ...
                            'Error', 'modal');
                        uiwait(h_error)
                        return
                    end
                    % Horn's method:
                    w = ones(size(A,1),1);
                    [regParams,Bfit,ErrorStats] = absor(A', B', 'doTrans', true, 'doScale', true, 'weights', w/sum(w));
                    fprintf('Error: errlsq=%0.2f, errmax=%0.2f\n', ErrorStats.errlsq, ErrorStats.errmax)
                    % Create affine transformation with anisotropic scaling:
                    app.tform = affine3d(regParams.M');
                case 'SVD'
                    bNan = any(isnan(app.Landmarks),2) | any(isnan(app.AtlasLandmarks),2);
                    A = app.Landmarks(~bNan,:);
                    B = app.AtlasLandmarks(~bNan,:);
                    if size(A,1)<2
                        h_error = errordlg('Not enough landmarks', ...
                            'Error', 'modal');
                        uiwait(h_error)
                        return
                    end
                    [~,~,~,~,T] = svd_transform(A, B);
                    app.tform = affine3d(T');
                case 'Bregma-Lambda'
                    Labels = lower(app.LandmarkLabels);
                    Bregma = app.Landmarks(Labels=="bregma",:);
                    Lambda = app.Landmarks(Labels=="lambda",:);
                    % Midline points:
                    bMidline = app.LandmarkTypes=="M";
                    A_mid = app.Landmarks(bMidline,:);
                    bNan = any(isnan(A_mid),2);
                    A_mid = A_mid(~bNan,:);
                    % Bilateral points:
                    bLeft = app.LandmarkTypes=="L";
                    A_left = app.Landmarks(bLeft,:);
                    A_right = NaN(size(A_left));
                    tmp = Labels(bLeft);
                    for iPoint = 1:numel(tmp)
                        bSel = ismember(Labels, strrep(tmp{iPoint},'left', 'right')) & ...
                            app.LandmarkTypes=="R";
                        if sum(bSel)==1
                            A_right(iPoint,:) = app.Landmarks(bSel,:);
                        end
                    end
                    bNan = any(isnan([A_left A_right]),2);
                    A_left = A_left(~bNan,:);
                    A_right = A_right(~bNan,:);
                    % Check:
                    if any(isnan(Bregma)) || any(isnan(Lambda)) || isempty(A_mid) || isempty(A_left) || isempty(A_right)
                        h_error = errordlg('Not enough landmarks', ...
                            'Error', 'modal');
                        uiwait(h_error)
                        return
                    end
                    % Align:
                    this_atlas = app.Atlases( ismember({app.Atlases.name}, app.RegistrationTarget));
                    [tform_scale, tform_roll, tform_align, tform_trans] = ...
                        align_bregma_lambda(Bregma, Lambda, A_left, A_right, A_mid, this_atlas.lambda_bregma);
                    app.tform = affine3d((tform_scale*tform_roll*tform_align*tform_trans)');
            end
            % Update state:
            app.State_TransformExists = true;
            UpdateState(app)
        end
        
    end
    
    
    methods (Access = public)
        
        function results = CrosshairMoved(app, src, evt)
            ix = (app.h_ortho.SliceNumbers-1) * app.Subsampling + 1; % AP,ML,DV
            ix = [ix(2) ix(1) ix(3)]; % ML,AP,DV
            UpdateCursorTable(app, ix)
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            % Default landmark labels:
            [~, app.LandmarkLabels, app.LandmarkTypes] = readLandmarks('default_landmarks.txt');
            app.Landmarks = NaN(numel(app.LandmarkLabels),3);
            % Load atlas definitions:
            app.Atlases = atlas_definitions;
            % Create RegistrationTargetSubMenu:
            for iAtlas = 1:numel(app.Atlases)
                app.RegistrationTargetSubMenu(iAtlas) = uimenu(app.RegistrationTargetMenu);
                app.RegistrationTargetSubMenu(iAtlas).MenuSelectedFcn = createCallbackFcn(app, @RegistrationTargetMenuSelected, true);
                app.RegistrationTargetSubMenu(iAtlas).Checked = 'off';
                app.RegistrationTargetSubMenu(iAtlas).Text = app.Atlases(iAtlas).name;
            end
            % Load atlas:
            app.LoadAtlasDefinitions;
            % Pass image argument:
            if numel(varargin)>0 % If we pass an image as an input argument...
                app.Img = squeeze(varargin{1});
                if ndims(app.Img)~=3
                    error('The image must have 3 non-singleton dimensions.')
                else
                    app.State_ImageExists = true;
                    UpdateState(app);
                    NewImage(app)
                end
            end
            
            % Update GUI state:
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
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
                app.GridLayout.ColumnWidth = {198, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Value changed function: HighSlider, LowSlider
        function LimitSliderValueChanged(app, event)
            app.DisplayRange = sort([app.LowSlider.Value app.HighSlider.Value]);
            app.LowSlider.Value = app.DisplayRange(1); app.HighSlider.Value = app.DisplayRange(2);
            app.h_ortho.DisplayRange = app.ImgMinMax(2)*app.DisplayRange/100;
        end

        % Value changed function: SubsamplingSpinner
        function SubsamplingSpinnerValueChanged(app, event)
            app.Subsampling = app.SubsamplingSpinner.Value;
            app.SubsampledSize = floor(size(app.Img)/app.Subsampling);
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Changing display resolution...', ...
                'Indeterminate','on');
            UpdateView(app)
            JumptolandmarkButtonPushed(app)
            delete(app.h_progress)
        end

        % Value changed function: SelectlandmarkDropDown
        function SelectlandmarkDropDownValueChanged(app, event)
            pos = app.Landmarks(strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels),:);
            if any(isnan(pos))
                app.JumptolandmarkButton.Enable = "off";
                app.DeletelandmarkButton.Enable = "off";
                set([app.MLEditField_2 app.APEditField_2 app.DVEditField_2], 'FontColor', [1 0 0])
            else
                set([app.MLEditField_2 app.APEditField_2 app.DVEditField_2], 'FontColor', [1 1 1])
                app.JumptolandmarkButton.Enable = "on";
                app.DeletelandmarkButton.Enable = "on";
            end
            UpdateLandmarkTable(app)
        end

        % Button pushed function: JumptolandmarkButton
        function JumptolandmarkButtonPushed(app, event)
            pos = app.Landmarks(strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels),:); % ML,AP,DV
            if any(isnan(pos))
                return
            end
            pos = round(1 + (pos-1)/app.Subsampling);
            if any(pos>app.SubsampledSize)
                uialert(app.UIFigure, 'Cannot jump beause the landmark is outside the image.', 'Out of bounds');
            else
                app.h_ortho.SliceNumbers = [pos(2) pos(1) pos(3)];
                ix = (app.h_ortho.SliceNumbers-1) * app.Subsampling + 1; % AP,ML,DV
                ix = [ix(2) ix(1) ix(3)]; % ML,AP,DV
                UpdateCursorTable(app, ix)
            end
            
        end

        % Menu selected function: LoadimageMenu
        function Loadimage(app, event)
            [file,path] = uigetfile({'*.nii'}, 'Select NIfTI file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Loading image...', ...
                'Indeterminate','on');
            app.Img = squeeze(niftiread(fullfile(path,file)));
            app.Img = repair_image_limits(app.Img);
            NewImage(app)
            delete(app.h_progress)
        end

        % Button pushed function: UpdatelandmarkButton
        function UpdatelandmarkButtonPushed(app, event)
            pos = cellfun(@str2num, app.UITable.Data(1,:));
            app.Landmarks(strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels),:) = pos;
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end

        % Value changed function: APEditField_2, DVEditField_2, 
        % MLEditField_2
        function LandmarkEditFieldValueChanged(app, event)
            bLandmark = strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels);
            pos = app.Landmarks(bLandmark,:);
            pos(ismember([app.MLEditField_2 app.APEditField_2 app.DVEditField_2], event.Source)) = event.Value;
            app.Landmarks(bLandmark,:) = pos;
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end

        % Menu selected function: SavelandmarksMenu
        function Savelandmarks(app, event)
            [file, path] = uiputfile({'*.txt'}, 'Save landmarks', app.PreferredDir);
            if isequal(file,0)
                return
            else
                app.PreferredDir = path;
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Saving landmarks...', ...
                    'Indeterminate','on');
                writeLandmarks(fullfile(path, file), app.Landmarks, app.LandmarkLabels);
                disp('File saved!')
                delete(app.h_progress)
            end
        end

        % Menu selected function: LoadlandmarksMenu
        function Loadlandmarks(app, event)
            [file,path] = uigetfile({'*.txt'}, 'Select landmark file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            [Landmarks, LandmarkLabels, LandmarkTypes] = ...
                readLandmarks(fullfile(path,file));
            app.Landmarks = NaN(numel(app.LandmarkLabels),3);
            for iCo = 1:size(Landmarks,1)
                b = ismember(app.LandmarkLabels, LandmarkLabels{iCo});
                if any(b)
                    app.Landmarks(b,:) = Landmarks(iCo,:);
                end
            end
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end

        % Button pushed function: PreviewregistrationButton
        function PreviewregistrationButtonPushed(app, event)
            A = app.Landmarks; % Image landmarks in native coordinates.
            B = app.AtlasLandmarks;
            % Plot landmarks:
            bNanNative = any(isnan(A),2);
            bNanAtlas = any(isnan(B),2);
            bNanAny = bNanNative | bNanAtlas;
            B2 = app.tform.transformPointsForward(app.Landmarks);   % Image landmarks in atlas coordinates.
            figure;
            hold on
            h1 = plot3(B2(:,1), B2(:,2), B2(:,3), 'o');
            h2 = plot3(B(:,1), B(:,2), B(:,3), '*');
            plot3([B(~bNanAny,1) B2(~bNanAny,1)]', ...
                [B(~bNanAny,2) B2(~bNanAny,2)]', ...
                [B(~bNanAny,3) B2(~bNanAny,3)]', 'k')
            text(B(~bNanAtlas,1), B(~bNanAtlas,2), B(~bNanAtlas,3), app.LandmarkLabels(~bNanAtlas))
            text(B2(bNanAtlas&~bNanNative,1), B2(bNanAtlas&~bNanNative,2), B2(bNanAtlas&~bNanNative,3), app.LandmarkLabels(bNanAtlas&~bNanNative))
            hold off
            grid on
            legend([h1 h2], 'Image', 'Atlas')
            xlabel('ML'); ylabel('AP'); zlabel('DV')
        end

        % Button pushed function: DeletelandmarkButton
        function DeletelandmarkButtonPushed(app, event)
            app.Landmarks(strcmp(app.SelectlandmarkDropDown.Value(3:end), app.LandmarkLabels),:) = [NaN NaN NaN];
            h = [app.MLEditField_2 app.APEditField_2 app.DVEditField_2];
            for ih = 1:3
                h(ih).Value = 0;
                h(ih).FontColor = [1 0 0];
            end
%             app.JumptolandmarkButton.Enable = "off";
%             app.DeletelandmarkButton.Enable = "off";
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end

        % Button pushed function: CalculatetransformationButton
        function CalculatetransformationButtonPushed(app, event)
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', ...
                'Message', 'Calculating transform...', 'Indeterminate','on');
            UpdateTransform(app)
            close(app.h_progress)
            UpdateCursorTable(app)
        end

        % Button pushed function: ViewfusedimageButton
        function ViewfusedimageButtonPushed(app, event)
            % moving this to the start to remove wait between clicks
            Atlas = app.Atlases(ismember({app.Atlases.name}, app.RegistrationTarget));
            if Atlas.format == "slices"
                % Get input:
                answer = inputdlg(...
                    {'Coronal step size (number of atlas slices to skip)', ...
                    'Sagittal step size (number of atlas slices to skip)'}, ...
                    'Step size', 1, ...
                    {'11' '2'});
                if isempty(answer)
                    return
                end
                step_coronal = str2num(answer{1});
                step_sagittal = str2num(answer{2});
            end
 
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Scaling image...', ...
                'Indeterminate', 'on', 'Cancelable', 'on');
            % Scale intensity:
            img = app.Img;
            if isinteger(img)
                MinImg = double(min(img(:)));
                MaxImg = double(max(img(:)));
                MinClass = double(intmin(class(img)));
                MaxClass = double(intmax(class(img)));
                MinLim = (MinImg-MinClass)/(MaxClass-MinClass);
                MaxLim = (MaxImg-MinClass)/(MaxClass-MinClass);
                for iSlice = 1:size(img,3)
                    img(:,:,iSlice) = imadjust(img(:,:,iSlice), [MinLim MaxLim]);
                    img(:,:,iSlice) = imadjust(img(:,:,iSlice), app.DisplayRange/100);
                end
            else
                img = img - min(img(:));
                img = img / max(img(:));
                for iSlice = 1:size(img,3)
                    img(:,:,iSlice) = imadjust(img(:,:,iSlice), app.DisplayRange/100);
                end
            end
            if app.h_progress.CancelRequested
                return
            end
            
            % Fuse:
            switch Atlas.format
                case 'volume' % FIX
                    if ~app.State_WaxholmT2Loaded
                        LoadWaxholm(app)
                    end
                    app.h_progress.Message = 'Transforming image...';
                    % Warp image:
                    img = imwarp(permute(img, [2 1 3]), app.tform, 'OutputView',imref3d([app.WaxholmSize(2) app.WaxholmSize(1) app.WaxholmSize(3)]));
                    img = permute(img, [2 1 3]);
                    app.h_progress.Message = 'Fusing image with atlas...';
                    % Fuse image:
                    ImgFuse = zeros(app.WaxholmSize(1), app.WaxholmSize(2), app.WaxholmSize(3), 3, 'uint8');
                    app.h_progress.Indeterminate = 'off';
                    for iZ = 1:app.WaxholmSize(3)
                        app.h_progress.Value = iZ/app.WaxholmSize(3);
                        ImgFuse(:,:,iZ,:) = imfuse(app.WaxholmImg(:,:,iZ), img(:,:,iZ));
                        if app.h_progress.CancelRequested
                            break
                        end
                    end
                    close(app.h_progress)
                    % Plot image:
                    orthosliceViewer(ImgFuse);
                case 'slices'
%                     % Get input:
%                     answer = inputdlg(...
%                         {'Coronal step size (number of atlas slices to skip)', ...
%                         'Sagittal step size (number of atlas slices to skip)'}, ...
%                         'Step size', 1, ...
%                         {'11' '2'});
%                     if isempty(answer)
%                         return
%                     end
%                     step_coronal = str2num(answer{1});
%                     step_sagittal = str2num(answer{2});
                    % Load:
                    app.h_progress.Message = 'Loading atlas...';
                    fprintf('Loading atlas...')
                    AtlasDir = fullfile(fileparts(mfilename('fullpath')), 'Atlases', Atlas.dir);
                    for iFile = 1:numel(Atlas.files)
                        L(iFile) = load(fullfile(AtlasDir, Atlas.files(iFile).name));
                    end
                    fprintf('done.\n')
                    if app.h_progress.CancelRequested
                        return
                    end
                    % Fuse:
                    app.h_progress.Message = 'Fusing image with atlas...';
                    if numel(L)==1
                        [coronal, sagittal] = fuse_image_with_atlas(img, app.tform, ...
                            Atlas.files(1).type, L(1), ...
                            'apstep', step_coronal, ...
                            'mlstep', step_sagittal);
                    else
                        [coronal, sagittal] = fuse_image_with_atlas(img, app.tform, ...
                            Atlas.files(1).type, L(1), ...
                            Atlas.files(2).type, L(2), ...
                            'apstep', step_coronal, ...
                            'mlstep', step_sagittal);
                    end
                    if app.h_progress.CancelRequested
                        return
                    end
                    % Plot:
                    if ~isempty(coronal)
                        figure;
                        clf
                        n = size(coronal.Img,4);
                        n1 = ceil(sqrt(n));
                        n2 = ceil(n./n1);
                        imshow(imtile(flip(coronal.Img), 'GridSize', [n2 n1]))
                    end
                    if ~isempty(sagittal)
                        figure;
                        clf
                        n = size(sagittal.Img,4);
                        n1 = ceil(sqrt(n));
                        n2 = ceil(n./n1);
                        imshow(imtile(flip(flip(sagittal.Img),2), 'GridSize', [n1 n2]))
                    end
                    % Unbusy:
                    close(app.h_progress)
            end
        end

        function PaxinosButtonPushed(app, event) % FIX
            % Get busy:
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Accessing Paxinos atlas...', ...
                'Indeterminate','on');
            % Get Paxinos:
            Paxinos_mm = cellfun(@str2num, app.UITable.Data(4,:)); % mm in Paxinos space. ML,AP,DV
            S = read_paxinos_rat_brain_atlas(Paxinos_mm(1),Paxinos_mm(2),Paxinos_mm(3),6);
            if isempty(S)
                h_error = errordlg('Coordinates are outside the atlas', ...
                    'Error', 'modal');
                uiwait(h_error)
            else 
                % Plot:
                figure
                clf
                imshow(S.coronal.image_marked);
            end
            % Unbusy:
            delete(app.h_progress)
        end
        
        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            NewData = str2num(event.NewData);
            if isempty(NewData) % If not numerical...
                ix = NaN(1,3);
            elseif event.Indices(1)>1 && ~app.State_TransformExists % If we have no transform...
                ix = NaN(1,3);
            else
                switch event.Indices(1)
                    case 1
                        ix = cellfun(@str2num, app.UITable.Data(1,:)); % ML,AP,DV
                    case 2
                        co_atlas = cellfun(@str2num, app.UITable.Data(2,:)); % Atlas coordinate [ML,AP,DV].
                        ix = round(app.tform.transformPointsInverse(co_atlas)); % Transforming from WHS to native.
                end
            end
            ix2 = round(1 + (ix-1)/app.Subsampling); % Slice in subsampled image. ML,AP,DV
            if any(ix2>app.SubsampledSize) || ~all(ix2>0)
                uialert(app.UIFigure, 'Value is out of range.', 'Out of bounds');
                ix = (app.h_ortho.SliceNumbers-1) * app.Subsampling + 1; % AP,ML,DV
                ix = [ix(2) ix(1) ix(3)]; % ML,AP,DV
                UpdateCursorTable(app,ix)
            else
                app.h_ortho.SliceNumbers = [ix2(2) ix2(1) ix2(3)];
                UpdateCursorTable(app, ix)
            end
        end

        % Menu selected function: LoadtransformationMenu
        function LoadtransformationMenuSelected(app, event)
            [file,path] = uigetfile({'*.mat'}, 'Select transformation file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            L = load(fullfile(path, file));
            % Transform:
            if isfield(L, 'tform')
                app.tform = L.tform;
                app.State_TransformExists = true;
                UpdateState(app)
            else
                uialert(app.UIFigure, 'Could not find a tform variable in this file.', 'Wrong file format');
            end
            % Registration target:
            if isfield(L, 'RegistrationTarget')
                if L.RegistrationTarget=="Paxinos (6th edition)" % Legacy...
                   L.RegistrationTarget = 'Paxinos and Watson: The Rat Brain (6th edition)';
                end
                app.RegistrationTarget = L.RegistrationTarget;
            else
                app.RegistrationTarget = 'Waxholm SD v1.01'; % Assume legacy target.
            end
            % Update menu:
            set(app.RegistrationTargetSubMenu, 'Checked', 'off')
            app.RegistrationTargetSubMenu(...
                ismember({app.RegistrationTargetSubMenu.Text}, app.RegistrationTarget)).Checked = 'on';
            % Registration method:
            if isfield(L, 'RegistrationMethod')
                app.RegistrationMethod = L.RegistrationMethod;
            else
                app.RegistrationMethod = 'Modified Horn'; % Assume legacy method.
            end
            set([app.ModifiedHornsMethodMenu ...
                 app.StandardHornsMethodMenu ...
                 app.SVDMethodMenu ...
                 app.LambdaBregmaMethodMenu], ...
                 'Checked', 'off')
            switch app.RegistrationMethod
                case "Modified Horn"
                    app.ModifiedHornsMethodMenu.Checked = "on";
                case "Standard Horn"
                    app.StandardHornsMethodMenu.Checked = "on";
                case "SVD"
                    app.SVDMethodMenu.Checked = "on";
                case "Bregma-Lambda"
                    app.LambdaBregmaMethodMenu.Checked = "on";
            end
            % Update:
            app.State_TransformExists = true;
            UpdateState(app)
        end
        
        % Menu selected function: SavetransformationMenu
        function SavetransformationMenuSelected(app, event)
            [file, path] = uiputfile({'*.mat'}, 'Save transform', app.PreferredDir);
            if isequal(file,0)
                return
            else
                app.PreferredDir = path;
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Saving transformation...', ...
                    'Indeterminate','on');
                tform = app.tform;
                RegistrationTarget = app.RegistrationTarget;
                RegistrationMethod = app.RegistrationMethod;
                save(fullfile(path, file), 'tform', 'RegistrationTarget', 'RegistrationMethod');
                disp('File saved!')
                delete(app.h_progress)
            end
        end
        
        % RegistrationTargetMenu selected function: 
        function RegistrationTargetMenuSelected(app, event)
            RegistrationTargetChanged(app,event.Source.Text)
        end

        function RegistrationTargetChanged(app, NewTarget)
            if isequal(NewTarget,app.RegistrationTarget)
                return
            end
            if ~isempty(app.tform)
%                 switch app.RegistrationTarget
%                     case {'Waxholm SD v1.01' 'Paxinos and Watson: The Rat Brain (6th edition)'} % We have methods to translate between these atlases.
%                         str = sprintf('What do you want to do with the current transformation? Delete, leave unchanged, or convert to the new target atlas (using the %s method)?', app.TranslationMethod);
%                         selection = uiconfirm(app.UIFigure,str,'Transformation', ...
%                             'Options', {'Delete' 'Leave unchanged' 'Convert'}, ...
%                             'DefaultOption', 'Convert');
%                         switch selection
%                             case 'Delete'
%                                 app.tform = [];
%                                 app.State_TransformExists = false;
%                             case 'Convert'
%                                 tform_slice2vox = [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
%                                 tform_vox2mm = whs_voxel2mm;
%                                 tform_wax2pax = get_whs2paxinos_tform(app.TranslationMethod);
%                                 switch app.RegistrationTarget
%                                     case 'Waxholm SD v1.01'
%                                         tform_wax2pax = get_whs2paxinos_tform(app.TranslationMethod);
%                                         app.tform.T = (tform_wax2pax * tform_vox2mm * tform_slice2vox * app.tform.T')';
%                                     case 'Paxinos and Watson: The Rat Brain (6th edition)'
%                                         app.tform.T = (inv(tform_wax2pax) * inv(tform_vox2mm) * inv(tform_slice2vox) * app.tform.T')';
%                                 end
%                         end
%                     otherwise
                        str = 'What do you want to do with the current transformation? Delete or leave unchanged?';
                        selection = uiconfirm(app.UIFigure,str,'Transformation', ...
                            'Options', {'Delete' 'Leave unchanged'}, ...
                            'DefaultOption', 'Delete');
                        switch selection
                            case 'Delete'
                                app.tform = [];
                                app.State_TransformExists = false;
                        end
%                 end
            end
            app.RegistrationTarget = NewTarget;
            % Load atlas:
            app.LoadAtlasDefinitions;
            % Update menu:
            set(app.RegistrationTargetSubMenu, 'Checked', 'off')
            app.RegistrationTargetSubMenu(...
                ismember({app.RegistrationTargetSubMenu.Text}, app.RegistrationTarget)).Checked = 'on';
            % Update state:
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end
        
        % RegistrationMethodMenu selected function: 
        function RegistrationMethodMenuSelected(app, event)
            set([app.ModifiedHornsMethodMenu ...
                 app.StandardHornsMethodMenu ...
                 app.SVDMethodMenu ...
                 app.LambdaBregmaMethodMenu], ...
                 'Checked', 'off')
            switch event.Source.Text
                case "Modified Horn's method (anisotropic)"
                    app.RegistrationMethod = 'Modified Horn';
                    app.ModifiedHornsMethodMenu.Checked = "on";
                case "Standard Horn's method (isotropic)"
                    app.RegistrationMethod = 'Standard Horn';
                    app.StandardHornsMethodMenu.Checked = "on";
                case "SVD method"
                    app.RegistrationMethod = 'SVD';
                    app.SVDMethodMenu.Checked = "on";
                case "Lambda-bregma method"
                    app.RegistrationMethod = 'Bregma-Lambda';
                    app.LambdaBregmaMethodMenu.Checked = "on";
            end
        end
        
        function CalculateLandmarksMenuSelected(app, event)
            [ML, AP, DV] = app.tform.transformPointsInverse(...
                app.AtlasLandmarks(:,1),...
                app.AtlasLandmarks(:,2),...
                app.AtlasLandmarks(:,3));
            app.Landmarks = [ML AP DV];
            app.State_LandmarksExist = ~all(any(isnan(app.Landmarks),2));
            UpdateState(app)
            UpdateSelectLandmarkDropDown(app)
            UpdateLandmarkTable(app)
        end
        
        function AdjusttransformationButtonPushed(app, event)
            % Get input:
            answer = inputdlg(...
                {'ML rotation (degrees, + is pitch up)', ...
                'AP rotation (degrees, + is roll right)', ...
                'DV rotation (degrees, + is yaw left)', ...
                'ML translation (+ is right)', ...
                'AP translation (+ is anterior)', ...
                'DV translation (+ is dorsal)', ...
                'ML scaling', ...
                'AP scaling', ...
                'DV scaling'}, ...
                'Adjust transform', 1, ...
                {'0' '0' '0' '0' '0' '0' '1' '1' '1'});
            if isempty(answer)
                return
            end
            v_ML = str2num(answer{1})/180*pi;
            v_AP = str2num(answer{2})/180*pi;
            v_DV = str2num(answer{3})/180*pi;
            t_ML = str2num(answer{4});
            t_AP = str2num(answer{5});
            t_DV = str2num(answer{6});
            s_ML = str2num(answer{7});
            s_AP = str2num(answer{8});
            s_DV = str2num(answer{9});
            % Define matrices for rotation, translation and scaling:
            rot = @(vx,vy,vz) [cos(vz)*cos(vy) cos(vz)*sin(vy)*sin(vx)-sin(vz)*cos(vx) cos(vz)*sin(vy)*cos(vx)+sin(vz)*sin(vx) 0;...
                sin(vz)*cos(vy) sin(vz)*sin(vy)*sin(vx)+cos(vz)*cos(vx) sin(vz)*sin(vy)*cos(vx)-cos(vz)*sin(vx) 0;...
                -sin(vy) cos(vy)*sin(vx) cos(vy)*cos(vx) 0;...
                0 0 0 1];
            trans = @(x,y,z) [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
            scaling = @(x,y,z) [x 0 0 0; 0 y 0 0; 0 0 z 0; 0 0 0 1];
            % Find origo (we want rotations around origo of the target space):
            [Orig(1), Orig(2), Orig(3)] = app.tform.transformPointsInverse(0,0,0);
            % Calculate transformation:
%             T_rot = trans(Orig(1),Orig(2),Orig(3)) * rot(v_ML,v_AP,v_DV) * trans(-Orig(1),-Orig(2),-Orig(3));
            T = scaling(s_ML,s_AP,s_DV) * trans(t_ML,t_AP,t_DV) * rot(v_ML,v_AP,v_DV) * app.tform.T';
            app.tform.T = T';
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
            app.UIFigure.Name = 'register_landmarks';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadimageMenu
            app.LoadimageMenu = uimenu(app.FileMenu);
            app.LoadimageMenu.MenuSelectedFcn = createCallbackFcn(app, @Loadimage, true);
            app.LoadimageMenu.Text = 'Load image';

            % Create LoadlandmarksMenu
            app.LoadlandmarksMenu = uimenu(app.FileMenu);
            app.LoadlandmarksMenu.MenuSelectedFcn = createCallbackFcn(app, @Loadlandmarks, true);
            app.LoadlandmarksMenu.Separator = 'on';
            app.LoadlandmarksMenu.Text = 'Load landmarks';

            % Create SavelandmarksMenu
            app.SavelandmarksMenu = uimenu(app.FileMenu);
            app.SavelandmarksMenu.MenuSelectedFcn = createCallbackFcn(app, @Savelandmarks, true);
            app.SavelandmarksMenu.Text = 'Save landmarks';

            % Create LoadtransformationMenu
            app.LoadtransformationMenu = uimenu(app.FileMenu);
            app.LoadtransformationMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadtransformationMenuSelected, true);
            app.LoadtransformationMenu.Separator = 'on';
            app.LoadtransformationMenu.Enable = 'on';
            app.LoadtransformationMenu.Text = 'Load transformation';
            
            % Create SavetransformationMenu
            app.SavetransformationMenu = uimenu(app.FileMenu);
            app.SavetransformationMenu.MenuSelectedFcn = createCallbackFcn(app, @SavetransformationMenuSelected, true);
            app.SavetransformationMenu.Enable = 'off';
            app.SavetransformationMenu.Text = 'Save transformation';

            % Create TransformMenu
            app.TransformMenu = uimenu(app.UIFigure);
            app.TransformMenu.Text = 'Transform';
            
            % Create RegistrationTargetMenu
            app.RegistrationTargetMenu = uimenu(app.TransformMenu);
            app.RegistrationTargetMenu.Text = 'Registration target';
            
            % Create RegistrationMethodMenu
            app.RegistrationMethodMenu = uimenu(app.TransformMenu);
            app.RegistrationMethodMenu.Text = 'Registration method';

            % Create ModifiedHornsMethodMenu
            app.ModifiedHornsMethodMenu = uimenu(app.RegistrationMethodMenu);
            app.ModifiedHornsMethodMenu.MenuSelectedFcn = createCallbackFcn(app, @RegistrationMethodMenuSelected, true);
            app.ModifiedHornsMethodMenu.Checked = 'on';
            app.ModifiedHornsMethodMenu.Text = 'Modified Horn''s method (anisotropic)';

            % Create StandardHornsMethodMenu
            app.StandardHornsMethodMenu = uimenu(app.RegistrationMethodMenu);
            app.StandardHornsMethodMenu.MenuSelectedFcn = createCallbackFcn(app, @RegistrationMethodMenuSelected, true);
            app.StandardHornsMethodMenu.Text = 'Standard Horn''s method (isotropic)';

            % Create SVDMethodMenu
            app.SVDMethodMenu = uimenu(app.RegistrationMethodMenu);
            app.SVDMethodMenu.MenuSelectedFcn = createCallbackFcn(app, @RegistrationMethodMenuSelected, true);
            app.SVDMethodMenu.Text = 'SVD method';
            
            % Create LambdaBregmamethodMenu
            app.LambdaBregmaMethodMenu = uimenu(app.RegistrationMethodMenu);
            app.LambdaBregmaMethodMenu.MenuSelectedFcn = createCallbackFcn(app, @RegistrationMethodMenuSelected, true);
            app.LambdaBregmaMethodMenu.Text = 'Lambda-bregma method';
            app.LambdaBregmaMethodMenu.Tooltip = 'Fix bregma at [0,0,0] and lambda at [0,-8.7772,0] and use other landmarks only for rotation around the AP axis.';

            % Create CalculateLandmarksMenu
            app.CalculateLandmarksMenu = uimenu(app.TransformMenu);
            app.CalculateLandmarksMenu.MenuSelectedFcn = createCallbackFcn(app, @CalculateLandmarksMenuSelected, true);
            app.CalculateLandmarksMenu.Separator = 'on';
            app.CalculateLandmarksMenu.Text = 'Calculate landmarks with inverse transform';
            app.CalculateLandmarksMenu.Tooltip = 'Calculate landmarks by using the inverse transform on the atlas landmarks.';
            app.CalculateLandmarksMenu.Enable = 'off';
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {215, '1x'};
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

            % Create TransformationPanel
            app.TransformationPanel = uipanel(app.LeftPanel);
            app.TransformationPanel.ForegroundColor = [1 1 1];
            app.TransformationPanel.BorderType = 'none';
            app.TransformationPanel.Title = 'Transformation';
            app.TransformationPanel.BackgroundColor = [0.149 0.149 0.149];
            app.TransformationPanel.FontWeight = 'bold';
            app.TransformationPanel.Position = [1 308 191 121];
            w = 148;
            x = round((app.TransformationPanel.Position(3)-w)/2);

            % Create CalculatetransformationButton
            app.CalculatetransformationButton = uibutton(app.TransformationPanel, 'push');
            app.CalculatetransformationButton.ButtonPushedFcn = createCallbackFcn(app, @CalculatetransformationButtonPushed, true);
            app.CalculatetransformationButton.BackgroundColor = [0.302 0.302 0.302];
            app.CalculatetransformationButton.FontColor = [1 1 1];
            app.CalculatetransformationButton.Enable = 'off';
            app.CalculatetransformationButton.Tooltip = {'Calculate affine transformation using current landmark definitions.'};
            app.CalculatetransformationButton.Position = [x 95 w 22];
            app.CalculatetransformationButton.Text = 'Calculate transformation';

            % Create PreviewregistrationButton
            app.PreviewregistrationButton = uibutton(app.TransformationPanel, 'push');
            app.PreviewregistrationButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewregistrationButtonPushed, true);
            app.PreviewregistrationButton.BackgroundColor = [0.302 0.302 0.302];
            app.PreviewregistrationButton.FontColor = [1 1 1];
            app.PreviewregistrationButton.Enable = 'off';
            app.PreviewregistrationButton.Tooltip = {'View how native landmarks relate to atlas landmarks after transformation.'};
            app.PreviewregistrationButton.Position = [x 66 w 22];
            app.PreviewregistrationButton.Text = 'Preview registration';
            
            % Create ViewfusedimageButton
            app.ViewfusedimageButton = uibutton(app.TransformationPanel, 'push');
            app.ViewfusedimageButton.ButtonPushedFcn = createCallbackFcn(app, @ViewfusedimageButtonPushed, true);
            app.ViewfusedimageButton.BackgroundColor = [0.302 0.302 0.302];
            app.ViewfusedimageButton.FontColor = [1 1 1];
            app.ViewfusedimageButton.Enable = 'off';
            app.ViewfusedimageButton.Tooltip = {'Fuse the image with the WHS atlas T2 image and view the result.'};
            app.ViewfusedimageButton.Position = [x 37 w 22];
            app.ViewfusedimageButton.Text = 'View fused image';
            
            % Create AdjusttransformationButton
            app.AdjusttransformationButton = uibutton(app.TransformationPanel, 'push');
            app.AdjusttransformationButton.ButtonPushedFcn = createCallbackFcn(app, @AdjusttransformationButtonPushed, true);
            app.AdjusttransformationButton.BackgroundColor = [0.302 0.302 0.302];
            app.AdjusttransformationButton.FontColor = [1 1 1];
            app.AdjusttransformationButton.Enable = 'off';
            app.AdjusttransformationButton.Tooltip = {'Positive angles are defined as ML=pitch up, AP=roll right, DV=yaw left. Positive translations are defined as ML=right, AP=anterior, DV=dorsal in units of the target system.'};
            app.AdjusttransformationButton.Position = [x 8 w 22];
            app.AdjusttransformationButton.Text = 'Adjust transformation';

            % Create ViewPanel
            app.ViewPanel = uipanel(app.LeftPanel);
            app.ViewPanel.ForegroundColor = [1 1 1];
            app.ViewPanel.BorderType = 'none';
            app.ViewPanel.Title = 'View';
            app.ViewPanel.BackgroundColor = [0.149 0.149 0.149];
            app.ViewPanel.FontWeight = 'bold';
            app.ViewPanel.Position = [1 6 191 290-25];

            % Create CursorLabel
            app.CursorLabel = uilabel(app.ViewPanel);
            app.CursorLabel.FontColor = [1 1 1];
            app.CursorLabel.Position = [14 220 45 22];
            app.CursorLabel.Text = 'Cursor:';

            % Create UITable
            app.UITable = uitable(app.ViewPanel);
            app.UITable.ColumnName = {'ML'; 'AP'; 'DV'};
            app.UITable.ColumnWidth = {39, 39, 39};
            app.UITable.RowName = {'Slice'; 'Atlas'};
            app.UITable.ColumnEditable = true;
            app.UITable.RowStriping = 'off';
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.Tooltip = {'Cursor position:'; '- Slice: Native slice number.'; '- Atlas: Estimated atlas coordinate.'};
            app.UITable.Enable = 'off';
            app.UITable.FontSize = 11;
            app.UITable.Position = [13 141 168 80];

            % Create PaxinosButton
            app.PaxinosButton = uibutton(app.ViewPanel, 'push');
            app.PaxinosButton.ButtonPushedFcn = createCallbackFcn(app, @PaxinosButtonPushed, true);
            app.PaxinosButton.BackgroundColor = [0.302 0.302 0.302];
            app.PaxinosButton.FontColor = [1 1 1];
            app.PaxinosButton.Enable = 'off';
            app.PaxinosButton.Position = [68 141 113 22];
            app.PaxinosButton.Text = 'View in Paxinos';
            app.PaxinosButton.Visible = 'off'; % FIX
            
            % Create SubsamplingLabel
            app.SubsamplingLabel = uilabel(app.ViewPanel);
            app.SubsamplingLabel.FontColor = [1 1 1];
            app.SubsamplingLabel.Position = [14 104 75 22];
            app.SubsamplingLabel.Text = 'Subsampling:';

            % Create SubsamplingSpinner
            app.SubsamplingSpinner = uispinner(app.ViewPanel);
            app.SubsamplingSpinner.RoundFractionalValues = 'on';
            app.SubsamplingSpinner.ValueChangedFcn = createCallbackFcn(app, @SubsamplingSpinnerValueChanged, true);
            app.SubsamplingSpinner.FontColor = [1 1 1];
            app.SubsamplingSpinner.BackgroundColor = [0.302 0.302 0.302];
            app.SubsamplingSpinner.Tooltip = {'Set subsampling factor of displayed image.'};
            app.SubsamplingSpinner.Position = [96 104 48 22];

            % Create IntensityrangeLabel
            app.IntensityrangeLabel = uilabel(app.ViewPanel);
            app.IntensityrangeLabel.FontColor = [1 1 1];
            app.IntensityrangeLabel.Position = [12 70 88 22];
            app.IntensityrangeLabel.Text = 'Intensity range:';

            % Create HighSliderLabel
            app.HighSliderLabel = uilabel(app.ViewPanel);
            app.HighSliderLabel.HorizontalAlignment = 'right';
            app.HighSliderLabel.FontColor = [1 1 1];
            app.HighSliderLabel.Position = [25 40 30 22];
            app.HighSliderLabel.Text = 'High';

            % Create HighSlider
            app.HighSlider = uislider(app.ViewPanel);
            app.HighSlider.MajorTicks = [];
            app.HighSlider.ValueChangedFcn = createCallbackFcn(app, @LimitSliderValueChanged, true);
            app.HighSlider.MinorTicks = [];
            app.HighSlider.Enable = 'off';
            app.HighSlider.Tooltip = {'Set intensity range.'};
            app.HighSlider.FontColor = [1 1 1];
            app.HighSlider.Position = [76 49 73 3];

            % Create LowSliderLabel
            app.LowSliderLabel = uilabel(app.ViewPanel);
            app.LowSliderLabel.HorizontalAlignment = 'right';
            app.LowSliderLabel.FontColor = [1 1 1];
            app.LowSliderLabel.Position = [26 5 29 22];
            app.LowSliderLabel.Text = 'Low';

            % Create LowSlider
            app.LowSlider = uislider(app.ViewPanel);
            app.LowSlider.MajorTicks = [];
            app.LowSlider.ValueChangedFcn = createCallbackFcn(app, @LimitSliderValueChanged, true);
            app.LowSlider.MinorTicks = [];
            app.LowSlider.Enable = 'off';
            app.LowSlider.Tooltip = {'Set intensity range.'};
            app.LowSlider.FontColor = [1 1 1];
            app.LowSlider.Position = [76 14 73 3];

            % Create LandmarksPanel
            app.LandmarksPanel = uipanel(app.LeftPanel);
            app.LandmarksPanel.ForegroundColor = [1 1 1];
            app.LandmarksPanel.BorderType = 'none';
            app.LandmarksPanel.Title = 'Landmarks';
            app.LandmarksPanel.BackgroundColor = [0.149 0.149 0.149];
            app.LandmarksPanel.FontWeight = 'bold';
            app.LandmarksPanel.Position = [1 431 191 221];

            % Create JumptolandmarkButton
            app.JumptolandmarkButton = uibutton(app.LandmarksPanel, 'push');
            app.JumptolandmarkButton.ButtonPushedFcn = createCallbackFcn(app, @JumptolandmarkButtonPushed, true);
            app.JumptolandmarkButton.BackgroundColor = [0.302 0.302 0.302];
            app.JumptolandmarkButton.FontColor = [1 1 1];
            app.JumptolandmarkButton.Tooltip = {'Move the cursor to the selected landmark.'};
            app.JumptolandmarkButton.Position = [34 21 133 22];
            app.JumptolandmarkButton.Text = 'Jump to landmark';

            % Create UpdatelandmarkButton
            app.UpdatelandmarkButton = uibutton(app.LandmarksPanel, 'push');
            app.UpdatelandmarkButton.ButtonPushedFcn = createCallbackFcn(app, @UpdatelandmarkButtonPushed, true);
            app.UpdatelandmarkButton.BackgroundColor = [0.302 0.302 0.302];
            app.UpdatelandmarkButton.FontColor = [1 1 1];
            app.UpdatelandmarkButton.Tooltip = {'Use the current cursor position as landmark coordinates.'};
            app.UpdatelandmarkButton.Position = [34 77 133 22];
            app.UpdatelandmarkButton.Text = 'Update landmark';

            % Create APLabel
            app.APLabel = uilabel(app.LandmarksPanel);
            app.APLabel.HorizontalAlignment = 'center';
            app.APLabel.FontColor = [1 1 1];
            app.APLabel.Position = [89 128 25 22];
            app.APLabel.Text = 'AP';

            % Create APEditField_2
            app.APEditField_2 = uieditfield(app.LandmarksPanel, 'numeric');
            app.APEditField_2.RoundFractionalValues = 'on';
            app.APEditField_2.ValueChangedFcn = createCallbackFcn(app, @LandmarkEditFieldValueChanged, true);
            app.APEditField_2.FontColor = [1 1 1];
            app.APEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.APEditField_2.Tooltip = {'Anterior-posterior slice number of currently selected landmark'};
            app.APEditField_2.Position = [78 107 45 22];

            % Create MLEditField_2
            app.MLEditField_2 = uieditfield(app.LandmarksPanel, 'numeric');
            app.MLEditField_2.RoundFractionalValues = 'on';
            app.MLEditField_2.ValueChangedFcn = createCallbackFcn(app, @LandmarkEditFieldValueChanged, true);
            app.MLEditField_2.FontColor = [1 1 1];
            app.MLEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.MLEditField_2.Tooltip = {'Medial-lateral slice number of currently selected landmark'};
            app.MLEditField_2.Position = [34 107 45 22];

            % Create DVEditField_2
            app.DVEditField_2 = uieditfield(app.LandmarksPanel, 'numeric');
            app.DVEditField_2.RoundFractionalValues = 'on';
            app.DVEditField_2.ValueChangedFcn = createCallbackFcn(app, @LandmarkEditFieldValueChanged, true);
            app.DVEditField_2.FontColor = [1 1 1];
            app.DVEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.DVEditField_2.Tooltip = {'Dorsal-ventral slice number of currently selected landmark'};
            app.DVEditField_2.Position = [122 107 45 22];

            % Create DeletelandmarkButton
            app.DeletelandmarkButton = uibutton(app.LandmarksPanel, 'push');
            app.DeletelandmarkButton.ButtonPushedFcn = createCallbackFcn(app, @DeletelandmarkButtonPushed, true);
            app.DeletelandmarkButton.BackgroundColor = [0.302 0.302 0.302];
            app.DeletelandmarkButton.FontColor = [1 1 1];
            app.DeletelandmarkButton.Tooltip = {'Clear coordinates of the selected landmark.'};
            app.DeletelandmarkButton.Position = [34 49 133 22];
            app.DeletelandmarkButton.Text = 'Delete landmark';

            % Create MLLabel
            app.MLLabel = uilabel(app.LandmarksPanel);
            app.MLLabel.HorizontalAlignment = 'center';
            app.MLLabel.FontColor = [1 1 1];
            app.MLLabel.Position = [45 128 25 22];
            app.MLLabel.Text = 'ML';

            % Create DVLabel
            app.DVLabel = uilabel(app.LandmarksPanel);
            app.DVLabel.HorizontalAlignment = 'center';
            app.DVLabel.FontColor = [1 1 1];
            app.DVLabel.Position = [133 128 25 22];
            app.DVLabel.Text = 'DV';

            % Create SelectlandmarkDropDownLabel
            app.SelectlandmarkDropDownLabel = uilabel(app.LandmarksPanel);
            app.SelectlandmarkDropDownLabel.FontColor = [1 1 1];
            app.SelectlandmarkDropDownLabel.Position = [14 177 96 22];
            app.SelectlandmarkDropDownLabel.Text = 'Select landmark:';

            % Create SelectlandmarkDropDown
            app.SelectlandmarkDropDown = uidropdown(app.LandmarksPanel);
            app.SelectlandmarkDropDown.Items = {''};
            app.SelectlandmarkDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectlandmarkDropDownValueChanged, true);
            app.SelectlandmarkDropDown.Tooltip = {'Select landmark'};
            app.SelectlandmarkDropDown.FontColor = [1 1 1];
            app.SelectlandmarkDropDown.BackgroundColor = [0.302 0.302 0.302];
            app.SelectlandmarkDropDown.Position = [14 155 171 22];
%             UpdateSelectLandmarkDropDown(app)

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
        function app = register_landmarks(varargin)

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