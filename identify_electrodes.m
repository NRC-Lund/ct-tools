classdef identify_electrodes < matlab.apps.AppBase
%IDENTIFY_ELECTRODES   A GUI for identifying electrodes in a CT image.
%
% Use this GUI to identify electrodes in a CT image. You need a NIfTI
% image (for example created with preprocess_ct), an affine3d tranformation
% object (for example created with register_landmarks), and electrode
% coordinates (can be loaded from the database or from an ElectrodeInfo.xls
% file).
%
% SYNTAX:
%   identify_electrodes
%   identify_electrodes(img)
%   identify_electrodes(img, tform)
%   obj = identify_electrodes(...)
%
% INPUT:
%   img     - (optional) 2- or 3-dimensional grayscale image matrix 
%             conforming with the preprocess_ct conventions (RAS).
%   tform   - (optional) affine3d object defining the transformation from
%             CT voxels to WHS voxels.
%
% OUTPUT:
%   obj     - Contains the GUI object (mostly useful for development and
%             debugging).

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        LoadimageMenu                  matlab.ui.container.Menu
        LoadtransformationMenu         matlab.ui.container.Menu
        LoadWiresXlsMenu               matlab.ui.container.Menu
        LoadWiresDatabaseMenu          matlab.ui.container.Menu
        LoadMarkersMatMenu             matlab.ui.container.Menu
        SaveMarkersMenu                matlab.ui.container.Menu
        SaveToDatabaseMenu             matlab.ui.container.Menu
        SegmentationMenu               matlab.ui.container.Menu
        CreateMaskMenu                 matlab.ui.container.Menu
        InitializeSegmentationMenu     matlab.ui.container.Menu
        PlotSegmentationLabelsMenu     matlab.ui.container.Menu
        ExtendSegmentionMenu           matlab.ui.container.Menu
        PlotSegmentationTipsMenu       matlab.ui.container.Menu
        PlotSegmentationOverlapMenu    matlab.ui.container.Menu
        MarkerMenu                     matlab.ui.container.Menu
        MarkerDeleteMenu               matlab.ui.container.Menu
        MarkerAddMenu                  matlab.ui.container.Menu
        SetDvMenu                      matlab.ui.container.Menu
        SetGaussianFilterSizeMenu      matlab.ui.container.Menu
        ShowMarkerLabelsMenu           matlab.ui.container.Menu
        GridLayout                     matlab.ui.container.GridLayout
        LeftPanel                      matlab.ui.container.Panel
        SelectMarkersButton            matlab.ui.control.Button
        ViewPanel                      matlab.ui.container.Panel
        CursorLabel                    matlab.ui.control.Label
        UITable                        matlab.ui.control.Table
        PaxinosButton                  matlab.ui.control.Button
        DvSliceLabel                   matlab.ui.control.Label
        DvSliceSpinner                 matlab.ui.control.Spinner
        DvSliceStepLabel               matlab.ui.control.Label
        DvSliceStepSpinner             matlab.ui.control.Spinner
        ShowMaskCheckBox               matlab.ui.control.CheckBox
        ShowMarkersCheckBox            matlab.ui.control.CheckBox
        IntensityrangeLabel            matlab.ui.control.Label
        HighSliderLabel                matlab.ui.control.Label
        HighSlider                     matlab.ui.control.Slider
        LowSliderLabel                 matlab.ui.control.Label
        LowSlider                      matlab.ui.control.Slider
        MarkersPanel                   matlab.ui.container.Panel
        AdjustMarkersButton            matlab.ui.control.Button
        APLabel                        matlab.ui.control.Label
        APEditField_2                  matlab.ui.control.NumericEditField
        MLEditField_2                  matlab.ui.control.NumericEditField
        DVEditField_2                  matlab.ui.control.NumericEditField
        MLLabel                        matlab.ui.control.Label
        DVLabel                        matlab.ui.control.Label
        SelectMarkerDropDownLabel      matlab.ui.control.Label
        SelectMarkerDropDown           matlab.ui.control.DropDown
        SelectMarkerTypeDropDownLabel  matlab.ui.control.Label
        SelectMarkerTypeDropDown       matlab.ui.control.DropDown
        RightPanel                     matlab.ui.container.Panel
        UIAxes                         matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = public, SetObservable)
        Img;            % Native image. dim1=ML, dim2=AP, dim3=DV. Right-Anterior-Superior.
        h_image;        % Handle to image object. Note that coordinates are [AP,ML].
        State_Image = false;
        State_Transformation = false;
        State_WiresTarget = false;
        State_WiresEnd = false;
        State_Mask = false;
        State_SegmentationInit = false;
        State_SegmentationExtend = false;
        DisplayRange = [0 100]; % Intensity range in percent.
        ImgMinMax;      % Min and max intensity values of image.
        bShowMarkers = true;
        bShowMask = true;
        bShowMarkerLabels = false;
        DvSlice = 1;
        DvPax = 0;
        GaussianFilterSize = 10;
        h_progress      % Handle to progress bar.
        PreferredDir = cd;   % Preferred directory.
        WireLabels;     % Marker labels.
        StartMarkers;   % Marker slice numbers in the native image. dim1=ML, dim2=AP, dim3=DV. Right-Anterior-Superior.
        TargetMarkers;
        EndMarkers;
        tform;          % affine3d object defining the affine transform to atlas space (slice numbers in ML,AP,DV).
        points;         % ROI objects.
        electroslice_nRoiSlices = 6; % Slices used to initialize electroslice.
        electroslice_limits = [];
        electroslice_bundle;
        electroslice_MaskSize = 2;
        electroslice_min_dist = 5;
        electroslice_mask = [];
        listener_points_moving;
        listener_points_moved;
        listener_points_clicked;
        RegistrationTarget = 'Paxinos and Watson: The Rat Brain (6th edition)';
    end
    
    methods (Access = public)
        function mask = CreateMask(app, MaskSize)
            if nargin==1
                MaskSize = 2;
            end
            % DV slice:
            if numel(unique(app.StartMarkers(:,3)))>1
                uialert(app.UIFigure, 'Start point markers must have the same DV coordinate.', 'Different DV');
                return
            end
            DvSlice = app.StartMarkers(1,3);
            % Find coordinates:;
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', sprintf('Tracing wires in slice %d...', 1));
            dv_frame_ix = fliplr(DvSlice-app.electroslice_nRoiSlices+1:DvSlice);
            nMarkers = size(app.StartMarkers,1);
            pos = NaN(nMarkers,2,app.electroslice_nRoiSlices);
            pos(:,:,1) = app.StartMarkers(:,[1 2]);
            %for iSlice = 1:app.electroslice_nRoiSlices-1
            for iSlice = 1:app.electroslice_nRoiSlices % attempted fix
                %app.h_progress.Value = (iSlice+1) / app.electroslice_nRoiSlices;
                app.h_progress.Value = iSlice / app.electroslice_nRoiSlices; % fix
                %app.h_progress.Message = sprintf('Tracing wires in slice %d...', iSlice+1);
                app.h_progress.Message = sprintf('Tracing wires in slice %d...', iSlice); % fix
                %Img = app.Img(:,:,DvSlice-iSlice);
                Img = app.Img(:,:,dv_frame_ix(iSlice)); % fix
                Img = im2double(Img);
                dummy = zeros(size(Img));
                for iMarker = 1:size(pos,1)
                    dummy(round(pos(iMarker,1,iSlice)),round(pos(iMarker,2,iSlice))) = 1;
                    dummy = imgaussfilt(dummy,app.GaussianFilterSize);
                    dummy = immultiply(dummy,Img);
                    [~,ix] = max(dummy(:));
                    [pos(iMarker,1,iSlice+1), pos(iMarker,2,iSlice+1)] = ind2sub(size(dummy),ix);
                end
            end
            pos = round(pos);
            % Create mask:
            mask = zeros(size(app.Img), 'uint8');
            SE = strel('disk',MaskSize);
            for iSlice = 1:app.electroslice_nRoiSlices
                for iMarker = 1:nMarkers
                    mask(pos(iMarker,1,iSlice),pos(iMarker,2,iSlice), dv_frame_ix(iSlice)) = 1;
                end
                mask(:,:,dv_frame_ix(iSlice)) = imdilate(mask(:,:,dv_frame_ix(iSlice)), SE);
            end
            delete(app.h_progress)
        end
        
        function bundle = ElectroSliceBundle(app)
            nMarkers = size(app.StartMarkers,1);
            Img = app.Img;
            bundle = electroslice.Bundle(nMarkers, app.electroslice_nRoiSlices, ...
                1, app.electroslice_mask, Img, true, app.electroslice_limits, app.electroslice_min_dist); % fix: delegate adjustment to electroslice
            if isempty(bundle.bundle_index)
                uialert(app.UIFigure, 'Segmentation failed to initialize.', 'Segmentation');
                return
            end
            co = app.StartMarkers(:,[1 2]);
            for iWire = 1:nMarkers
                [~,ix] = min(sum((bundle.wires(iWire).manual_init(1,1:2) - co).^2,2));
                bundle.wires(iWire).wire_label = str2num(app.WireLabels{ix});
            end
            % TODO: add check for repeated wire labels before plotting
            % if there are repeated wire labls, initialization should be
            % run again (probably with higher minimum distance)
            bundle.plot_labels(); % it's better to plot after labels have been corrected
        end
        
        function markers = ElectroSliceExtend(app, bundle)
            Img = app.Img;
            bundle.extend_all( Img, true, app.electroslice_limits); % attempted fix
            nMarkers = size(app.StartMarkers,1);
            tips = NaN(nMarkers,3);
            co = app.StartMarkers(:,[1 2]);
            for iWire = 1:nMarkers
                [~,ix] = min(sum((bundle.wires(iWire).manual_init(1,1:2) - co).^2,2));
                tips(ix,:) = bundle.wires(iWire).get_tip();
            end
            markers = table(cellfun(@str2num,app.WireLabels), ...
                tips(:,1), tips(:,2), tips(:,3), ...
                'VariableNames', {'dbid' 'ML' 'AP' 'DV'});
            if any(isnan(markers.ML))
                ix = find(isnan(markers.ML));
                for iWire = 1:numel(ix)
                    fprintf('Wire %s has an undefined end point and was replaced with the starting point.\n', app.WireLabels{ix(iWire)})
                    markers.ML(ix(iWire)) = app.StartMarkers(ix(iWire),1);
                    markers.AP(ix(iWire)) = app.StartMarkers(ix(iWire),2);
                    markers.DV(ix(iWire)) = app.StartMarkers(ix(iWire),3);
                end
            end
        end
        
        function RepairImageLimits(app)
            if isfloat(app.Img) % Floaty images must be [0 1].
                MinLim = min(app.Img(:));
                if MinLim<0
                    app.Img = app.Img - MinLim;
                    MaxLim = max(app.Img(:));
                    if MaxLim>1
                        app.Img = app.Img / MaxLim;
                    end
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function NewMarkers(app)
            if app.State_WiresTarget || app.State_WiresEnd
                delete(app.points)
                cmap = lines(5);
                switch app.SelectMarkerTypeDropDown.Value
                    case 'Target'
                        col = cmap(4,:);
                        co = app.TargetMarkers;
                    case 'Starting point'
                        col = cmap(1,:);
                        co = app.StartMarkers;
                    case 'End point'
                        col = cmap(5,:);
                        co = app.EndMarkers;
                end
                for iPoint = 1:size(co,1)
                    points(iPoint) = drawpoint('Parent', app.UIAxes, ...
                        'Position',[co(iPoint,2) co(iPoint,1)], ...
                        'Color', col); % AP, ML
                    points(iPoint).SelectedColor = cmap(2,:);
                end
                [points.UIContextMenu] = deal([]);
                if app.bShowMarkerLabels
                    [points.Label] = deal(app.WireLabels{:});
                    [points.LabelAlpha] = deal(0.7);
                end
                app.points = points;
                delete(app.listener_points_moving)
                delete(app.listener_points_moved)
                delete(app.listener_points_clicked)
                app.listener_points_moving = addlistener(points,'MovingROI',createCallbackFcn(app, @PointEvent, true));
                app.listener_points_moved = addlistener(points,'ROIMoved',createCallbackFcn(app, @PointEvent, true));
                app.listener_points_clicked = addlistener(points,'ROIClicked',createCallbackFcn(app, @PointEvent, true));
            end
            UpdateMarkers(app); % Updates visibility.
        end
        
        function UpdateMarkers(app)
            % Set marker visibility:
            if ~isempty(app.points)
                switch app.SelectMarkerTypeDropDown.Value
                    case 'Target'
                        pos = app.TargetMarkers; % ML,AP,DV
                    case 'Starting point'
                        pos = app.StartMarkers; % ML,AP,DV
                    case 'End point'
                        pos = app.EndMarkers; % ML,AP,DV
                end
                bSel = round(pos(:,3))==app.DvSlice;
                [app.points.Visible] = deal(false);
                if any(bSel) && app.bShowMarkers
                    [app.points(bSel).Visible] = deal(true);
                end
            end
        end
        
        function NewImage(app)
            % DV slice:
            app.DvSlice = 1;
            app.DvSliceSpinner.Value = 1;
%             app.DVEditField_2.Value = app.DvSlice;
            if ndims(app.Img)==3
                app.DvSliceSpinner.Limits = [1 size(app.Img,3)];
            end
            % Limits:
            app.ImgMinMax = double([min(app.Img(:)) max(app.Img(:))]);
            app.DisplayRange = [0 100];
            app.LowSlider.Value = app.DisplayRange(1);
            app.HighSlider.Value = app.DisplayRange(2);
            % Delete any transformation and markers:
            app.StartMarkers = [];
            app.EndMarkers = [];
            app.TargetMarkers = [];
            app.WireLabels = [];
            if ~isempty(app.points)
                delete(app.points);
            end
            % Update state:
            app.State_Image = true;
            app.State_WiresTarget = false;
            app.State_WiresEnd = false;
            app.State_Mask = false;
            app.State_Mask = false;
            app.State_SegmentationInit = false;
            app.State_SegmentationExtend = false;
            UpdateState(app);
            % Update image:
            app.h_image = imshow(app.Img(:,:,app.DvSlice), 'Parent', app.UIAxes);
            app.UIAxes.CLim = app.ImgMinMax(2)*app.DisplayRange/100;
            % Arrows:
            arrowlength = min(diff(app.UIAxes.XLim), diff(app.UIAxes.YLim)) * 0.1;
            offset = arrowlength/2;
            headsize = arrowlength/4;
            h_line1 = line(app.UIAxes, 'XData', [app.UIAxes.XLim(1) app.UIAxes.XLim(1)] + offset, ...
                'YData', [app.UIAxes.YLim(1) app.UIAxes.YLim(1)+arrowlength] + offset, ...
                'Color', 'm', 'LineWidth', 1);
            arrowhead = [h_line1.XData(2) h_line1.YData(2); h_line1.XData(2)+headsize/2 h_line1.YData(2)-headsize; h_line1.XData(2)-headsize/2 h_line1.YData(2)-headsize];
            h_head1 = patch(app.UIAxes, arrowhead(:,1), arrowhead(:,2), 'm');
            h_txt1 = text(app.UIAxes, h_line1.XData(2), h_line1.YData(2), '', "Color", 'm', "VerticalAlignment","top", "HorizontalAlignment","center");
            h_line2 = line(app.UIAxes, 'XData', [app.UIAxes.XLim(1) app.UIAxes.XLim(1)+arrowlength] + offset, ...
                'YData', [app.UIAxes.YLim(1) app.UIAxes.YLim(1)] + offset, ...
                'Color', 'm', 'LineWidth', 1);
            arrowhead = [h_line2.XData(2) h_line2.YData(2); h_line2.XData(2)-headsize h_line2.YData(2)-headsize/2; h_line2.XData(2)-headsize h_line2.YData(2)+headsize/2];
            h_head2 = patch(app.UIAxes, arrowhead(:,1), arrowhead(:,2), 'm');
            h_txt2 = text(app.UIAxes, h_line2.XData(2), h_line2.YData(2), '', "VerticalAlignment","middle", "HorizontalAlignment","left");
            set([h_txt1 h_txt2], "Color", 'm', "FontSize", 12)
            h_txt1.String = "Right";
            h_txt2.String = "Anterior";
            % Update markers:
            NewMarkers(app);
        end
        
        function UpdateImage(app)
            if isempty(app.electroslice_mask) || ~app.bShowMask
                app.h_image.CData = app.Img(:,:,app.DvSlice);
            else
                mask = double(logical(app.electroslice_mask(:,:,app.DvSlice)));
                mask = cat(3, zeros(size(mask,1),size(mask,2),2),mask);
                app.h_image.CData = imfuse(app.Img(:,:,app.DvSlice), mask, 'blend');
            end
        end
        
        function NewWires(app, co_start, co_target, co_end, labels, dv_slice)
            if ~isempty(co_target) && ~any(isnan(co_target(:)))
                app.StartMarkers = co_start; % ML, AP, DV.
                app.TargetMarkers = co_target; % ML, AP, DV.
                app.State_WiresTarget = true;
            else
                app.State_WiresTarget = false;
            end
            if ~isempty(co_end) && ~any(isnan(co_end(:)))
                app.EndMarkers = co_end;
                app.State_WiresEnd = true;
            else
                app.State_WiresEnd = false;
            end
            app.WireLabels = labels;
            if ~isempty(app.points)
                delete(app.points);
            end
            UpdateState(app);
            NewMarkers(app);
            app.SelectMarkerDropDown.Items = [{'-'}; app.WireLabels];
            UpdateMarkerTable(app);
            if nargin==6
                app.DvSlice = dv_slice;
                app.DvSliceSpinner.Value = app.DvSlice;
                UpdateImage(app)
            end
        end
        
        function co_start = CalculateStart(app, co_target, DvSlice)
            co_zero = app.tform.transformPointsForward(co_target);
            co_zero(:,3) = 0;
            co_zero = app.tform.transformPointsInverse(co_zero);
            v_line = co_target-co_zero;
            v_plane = [0 0 1];   % Plane normal.
            p_plane = [0 0 DvSlice];   % Point on plane.
            co_start = co_target;
            for iWire = 1:size(co_target,1)
                co_start(iWire,:) = intersect_line_plane(...
                    v_line(iWire,:), co_target(iWire,:), v_plane, p_plane);
            end
        end
        
        function UpdateState(app)
            % Enable/disable controls depending on where in the process we are:
            set([app.HighSlider app.LowSlider app.UITable ...
                app.SelectMarkerDropDown ...
                app.SelectMarkerTypeDropDown ...
                app.APEditField_2 ...
                app.MLEditField_2 ...
                app.DVEditField_2 ...
                app.AdjustMarkersButton ...
                app.SelectMarkersButton ...
                app.DvSliceSpinner ...
                app.ShowMaskCheckBox ...
                app.ShowMarkersCheckBox ...
                app.LoadtransformationMenu ...
                app.LoadWiresDatabaseMenu ...
                app.LoadMarkersMatMenu ...
                app.LoadWiresXlsMenu ...
                app.SaveMarkersMenu ...
                app.SaveToDatabaseMenu ...
                app.CreateMaskMenu ...
                app.InitializeSegmentationMenu ...
                app.PlotSegmentationLabelsMenu ...
                app.ExtendSegmentionMenu ...
                app.PlotSegmentationTipsMenu ...
                app.PlotSegmentationOverlapMenu ...
                app.MarkerDeleteMenu ...
                app.MarkerAddMenu ...
                app.SetDvMenu ...
                app.PaxinosButton], ...
                'Enable', 'off'); % PlotSegLabels, PlotSegTips and PlotSegOverlap are NEW!
            if app.State_Image % We have loaded an image...
                app.HighSlider.Enable = "on";
                app.LowSlider.Enable = "on";
                app.UITable.Enable = "on";
                app.LoadtransformationMenu.Enable = "on";
                app.LoadMarkersMatMenu.Enable = "on";
                if ndims(app.Img)==3
                    app.DvSliceSpinner.Enable = "on";
                end
            end
            if app.State_Transformation % We have loaded a transform...
                app.LoadWiresXlsMenu.Enable = "on";
                app.LoadWiresDatabaseMenu.Enable = "on";
            end
            if app.State_WiresTarget % We have loaded wire coordinates...
                app.SelectMarkerDropDown.Enable = "on";
                app.SelectMarkerTypeDropDown.Enable = "on";
                app.APEditField_2.Enable = "on";
                app.MLEditField_2.Enable = "on";
                app.DVEditField_2.Enable = "on";
                app.AdjustMarkersButton.Enable = "on";
                app.SelectMarkersButton.Enable = "on";
                app.SaveMarkersMenu.Enable = "on";
                app.SelectMarkerTypeDropDown.Items = {'Target' 'Starting point'};
                app.SelectMarkerTypeDropDown.Value = 'Starting point';
                app.CreateMaskMenu.Enable = "on";
                %app.PaxinosButton.Enable = "on"; % Disabled functionality at the moment.
                app.ShowMarkersCheckBox.Enable = "on";
                app.MarkerDeleteMenu.Enable = "on";
                app.MarkerAddMenu.Enable = "on";
                app.SetDvMenu.Enable = "on";
            end
            if app.State_Mask % We have a mask... 
                app.InitializeSegmentationMenu.Enable = "on";
                app.ShowMaskCheckBox.Enable = "on";
                if app.bShowMask==1
                    app.LowSlider.Enable = "off";
                    app.HighSlider.Enable = "off";
                else
                    app.LowSlider.Enable = "on";
                    app.HighSlider.Enable = "on";
                end
            end
            if app.State_SegmentationInit % We have initialized segmentation...
                app.PlotSegmentationLabelsMenu.Enable = "on"; % NEW!
                app.ExtendSegmentionMenu.Enable = "on";
            end
            if app.State_WiresEnd % We have wire end points...
                app.SelectMarkerTypeDropDown.Items = {'Target' 'Starting point' 'End point'};
                app.SaveToDatabaseMenu.Enable = "on";
            end
            if app.State_SegmentationExtend
                app.PlotSegmentationTipsMenu.Enable = "on"; % NEW!
                app.PlotSegmentationOverlapMenu.Enable = "on"; % NEW!
            end
        end
        
        function UpdateCursorTable(app, ix)
            if isempty(ix) % If we don't have a coordinate...
                ix = NaN(1,3);
            end
            if app.State_Transformation % If we have a tranform...
                co_atlas = app.tform.transformPointsForward(ix);
            else
                co_atlas = NaN(1,3);
            end
            % Table:
            DataCell = num2cell([round(ix); co_atlas]);
            DataCell(1,:) = cellfun(@(x) num2str(x,'%d'), DataCell(1,:), 'UniformOutput', false);
            DataCell(2,:) = cellfun(@(x) num2str(x,'%0.2f'), DataCell(2,:), 'UniformOutput', false);
            app.UITable.Data = DataCell;
        end
        
        function UpdateMarkerTable(app)
            pos = vertcat(app.points.Position);
            switch app.SelectMarkerTypeDropDown.Value
                case 'Target'
                    app.TargetMarkers(:,[1 2]) = pos(:,[2 1]);
                case 'Starting point'
                    app.StartMarkers(:,[1 2]) = pos(:,[2 1]);
                case 'End point'
                    app.EndMarkers(:,[1 2]) = pos(:,[2 1]);
            end
            h = [app.MLEditField_2 app.APEditField_2 app.DVEditField_2];
            bSel = [app.points.Selected];
            if ~any(bSel)
                [h.Enable] = deal('off');
                [h.Value] = deal(0);
                app.SelectMarkerDropDown.Value = '-';
                app.SelectMarkerDropDown.BackgroundColor = [0.302 0.302 0.302];
                app.AdjustMarkersButton.Enable = "off";
                UpdateCursorTable(app, [])
            else
                if sum(bSel)==1
                    [h.Enable] = deal('on');
                    app.SelectMarkerDropDown.Value = app.WireLabels(bSel);
                    app.SelectMarkerDropDown.BackgroundColor = [0.8500 0.3250 0.0980];
                    bSel = strcmp(app.SelectMarkerDropDown.Value, app.WireLabels);
                    switch app.SelectMarkerTypeDropDown.Value
                        case 'Target'
                            pos = app.TargetMarkers(bSel,:); % ML,AP,DV
                        case 'Starting point'
                            pos = app.StartMarkers(bSel,:); % ML,AP,DV
                        case 'End point'
                            pos = app.EndMarkers(bSel,:); % ML,AP,DV
                    end
                    app.MLEditField_2.Value = pos(1);
                    app.APEditField_2.Value = pos(2);
                    app.DVEditField_2.Value = pos(3);
                    UpdateCursorTable(app, pos)
                else
                    [h.Enable] = deal('off');
                    [h.Value] = deal(0);
                    app.SelectMarkerDropDown.Value = '-';
                    app.SelectMarkerDropDown.BackgroundColor = [0.302 0.302 0.302];
                    UpdateCursorTable(app, [])
                end
                app.AdjustMarkersButton.Enable = "on";
            end
        end
        
        function DvChanged(app, src, evnt)
            app.DvSliceSpinner.Value = app.DvSlice;
            UpdateImage(app)
            UpdateMarkers(app)
        end
        
        function RegistrationTargetOut = CheckRegistrationTarget(app, RegistrationTargetIn)
            Atlases = atlas_definitions;
            names = {Atlases.name};
            if ismember(RegistrationTargetIn, names)
                RegistrationTargetOut = RegistrationTargetIn;
            else
                if isempty(RegistrationTargetIn)
                    RegistrationTargetIn = '<EMPTY>';
                end
                [sel,ok] = listdlg('PromptString',...
                    {'This file specifies an unknown registration target:', ...
                    RegistrationTargetIn, ...
                    '', 'Please select a valid registration target:'},...
                    'SelectionMode','single',...
                    'ListString', names, 'ListSize', [500 200]);
                if ok
                    RegistrationTargetOut = names{sel};
                else
                    RegistrationTargetOut = '';
                end
            end
        end
        
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            % Pass image argument:
            if numel(varargin)>0 % If we pass an image as an input argument...
                app.Img = squeeze(varargin{1});
                app.State_Image = true;
                NewImage(app)
            end
            if numel(varargin)>1 % If we pass a tform...
                app.tform = squeeze(varargin{2});
                app.State_Transformation = true;
                UpdateState(app)
                NewMarkers(app)
            end
            addlistener(app,'DvSlice','PostSet',@app.DvChanged);
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
            app.LowSlider.Value = app.DisplayRange(1);
            app.HighSlider.Value = app.DisplayRange(2);
            app.LowSlider.Tooltip = {sprintf('%0.2f%%', app.LowSlider.Value)};
            app.HighSlider.Tooltip = {sprintf('%0.2f%%', app.HighSlider.Value)};
            app.UIAxes.CLim = app.ImgMinMax(2)*app.DisplayRange/100;
        end

        % Value changed function: DvSliceSpinner
        function DvSliceSpinnerValueChanged(app, event)
            app.DvSlice = app.DvSliceSpinner.Value;
        end
        
        % Value changed function: DvSliceStepSpinner
        function DvSliceStepSpinnerValueChanged(app, event)
            app.DvSliceSpinner.Step = event.Value;
        end
        
        % Value changed function: ShowMarkersCheckBoxValueChanged
        function ShowMarkersCheckBoxValueChanged(app, event)
            app.bShowMarkers = event.Value;
            UpdateMarkers(app);
        end
        
        % Value changed function: ShowMarkersCheckBoxValueChanged
        function ShowMaskCheckBoxValueChanged(app, event)
            app.bShowMask = event.Value;
            UpdateImage(app)
            if app.bShowMask
                app.LowSlider.Enable = 'off';
                app.HighSlider.Enable = 'off';
            else
                app.LowSlider.Enable = 'on';
                app.HighSlider.Enable = 'on';
            end
        end
        
        % Value changed function: SelectMarkerDropDown
        function SelectMarkerDropDownValueChanged(app, ~)
            bSel = ismember(app.WireLabels, app.SelectMarkerDropDown.Value);
            [app.points.Selected] = deal(false);
            if any(bSel)
                app.points(bSel).Selected = true;
                switch app.SelectMarkerTypeDropDown.Value
                    case 'Target'
                        app.DvSlice = round(app.TargetMarkers(bSel,3));
                    case 'Starting point'
                        app.DvSlice = round(app.StartMarkers(bSel,3));
                    case 'End point'
                        app.DvSlice = round(app.EndMarkers(bSel,3));
                end
            end
            UpdateMarkerTable(app)
        end
        
        % Value changed function: SelectMarkerDropDown
        function SelectMarkerTypeDropDownValueChanged(app, event)
            if ~isempty(app.points)
                delete(app.points);
            end
            cmap = lines(5);
            switch app.SelectMarkerTypeDropDown.Value
                case 'Target'
                    app.SelectMarkerTypeDropDown.BackgroundColor=cmap(4,:);
                case 'Starting point'
                    app.SelectMarkerTypeDropDown.BackgroundColor=cmap(1,:);
                case 'End point'
                    app.SelectMarkerTypeDropDown.BackgroundColor=cmap(5,:);
            end
            NewMarkers(app)
            SelectMarkerDropDownValueChanged(app)
        end

        % Button pushed function: AdjustMarkersButton
        function AdjustMarkersButtonPushed(app, event)
            PreviousPositions = vertcat(app.points.Position);
            Img = app.Img(:,:,app.DvSlice);
            Img = im2double(Img);
            mask = zeros(size(Img));
            bSel = [app.points.Selected];
            if any(bSel)
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Adjusting markers...');
                PointIx = find(bSel);
                for iPoint = 1:numel(PointIx)
                    mask(round(app.points(PointIx(iPoint)).Position(2)),round(app.points(PointIx(iPoint)).Position(1))) = 1;
                    mask = imgaussfilt(mask,app.GaussianFilterSize);
                    mask = immultiply(mask,Img);
                    [~,ix] = max(mask(:));
                    [x, y] = ind2sub(size(mask),ix);
                    app.points(PointIx(iPoint)).Position = [y x];
                    app.h_progress.Value = iPoint/numel(PointIx);
                end
            end
            UpdateMarkerTable(app)
            delete(app.h_progress)
            answer = uiconfirm(app.UIFigure, 'Keep changes?', 'Adjust markers',...
                    "Options", {'OK' 'Cancel'}, "DefaultOption", 1);
            if answer=="Cancel"
                for iPoint = 1:numel(app.points)
                    app.points(iPoint).Position = PreviousPositions(iPoint,:);
                end
            end
        end

        function PointEvent(app,event)
            switch event.EventName
                case{'MovingROI'}
                    d = event.CurrentPosition-event.PreviousPosition;
                    SelRoi = findobj(event.Source.Parent, 'Type', 'images.roi.point', 'Selected', true);
                    SelRoi = setdiff(SelRoi, event.Source);
                    for iRoi = 1:numel(SelRoi)
                        SelRoi(iRoi).Position = SelRoi(iRoi).Position + d;
                    end
                    event.Source.Selected = true;
                    
                case 'ROIClicked'
                    event.Source.Selected = ~event.PreviousSelected;
            end
            UpdateMarkerTable(app)
        end
        
        function AxesButtonDownFcn(app, event)
%             [app.points.Selected] = deal(false);
        end
        
        function ScrollWheelFcn( app, event)
          app.DvSlice = app.DvSlice - event.VerticalScrollCount; % negative sign to keep directions intuitive
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
            RepairImageLimits(app);
            NewImage(app)
            delete(app.h_progress)
        end

        % Value changed function: APEditField_2, DVEditField_2, 
        % MLEditField_2
        function MarkerEditFieldValueChanged(app, event)
            bMarker = strcmp(app.SelectMarkerDropDown.Value, app.WireLabels);
            switch app.SelectMarkerTypeDropDown.Value
                case 'Target'
                    name = 'TargetMarkers';
                case 'Starting point'
                    name = 'StartMarkers';
                case 'End point'
                    name = 'EndMarkers';
            end
            pos = app.(name)(bMarker,:);
            pos(ismember([app.MLEditField_2 app.APEditField_2 app.DVEditField_2], event.Source)) = event.Value;
            app.(name)(bMarker,:) = pos;
            app.points(bMarker).Position = [pos(2) pos(1)];
            UpdateMarkerTable(app)
            UpdateMarkers(app)
        end

        % Button pushed function: SelectMarkersButton
        function SelectMarkersButtonPushed(app, event)
            cmap = lines(2);
            rect = drawrectangle(app.UIAxes, 'Color', 'w', 'StripeColor', 'k', 'FaceAlpha', 0, ...
                'LineWidth', 2);
            bSel = arrayfun(@(x) inROI(rect, x.Position(1), x.Position(2)), app.points);
            [app.points(bSel).Selected] = deal(true);
            [app.points(~bSel).Selected] = deal(false);
            delete(rect)
            UpdateMarkerTable(app)
        end
        
        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
        end

        function PaxinosButtonPushed(app, event)
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
        
        % Menu selected function: LoadtransformationMenu
        function LoadtransformationMenuSelected(app, event)
            [file,path] = uigetfile({'*.mat'}, 'Select transformation file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            L = load(fullfile(path, file));
            if isfield(L, 'tform')
                % Registration target:
                if isfield(L, 'RegistrationTarget')
                    RegistrationTarget = CheckRegistrationTarget(app, L.RegistrationTarget);
                else
                    RegistrationTarget = CheckRegistrationTarget(app, '');
                end
                if isempty(RegistrationTarget)
                    return
                else
                    app.RegistrationTarget = RegistrationTarget;
                end
                % Transform:
                app.tform = L.tform;
                app.State_Transformation = true;
                UpdateState(app)
                uialert(app.UIFigure, 'Done!', 'Transformation loaded', 'Icon', 'success');
            else
                uialert(app.UIFigure, 'Could not find a tform variable in this file.', 'Wrong file format');
            end
        end
        
        function ClearWires(app)
            app.WireLabels = [];
            app.TargetMarkers = [];
            app.StartMarkers = [];
            app.EndMarkers = [];
            app.State_WiresTarget = false;
            app.State_WiresEnd = false;
            delete(app.listener_points_clicked)
            delete(app.listener_points_moved)
            delete(app.listener_points_moving)
            if ~isempty(app.points)
                delete(app.points);
            end
            UpdateState(app)
        end
        
        function LoadWiresXlsMenuSelected(app, event)
            % Selet DV plane:
            answer = inputdlg('Please set the DV plane to use as start position:', ...
                'DV start position', 1, {num2str(app.DvSlice)});
            if isempty(answer)
                return
            end
            DvSlice = str2num(answer{1});
            % Select file:
            [file,path] = uigetfile({'*.xls;*.xlsx'}, 'Select ElectrodeInfo file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            % Read Excel file:
            try
                [~,txt,raw]=xlsread(fullfile(path,file), 1);
            catch
                h_error = errordlg(['Could not read the file ' fullfile(path,file)], ...
                    'Error', 'modal');
                uiwait(h_error)
                return
            end
            % Find start tag:
            DiveIn = find(contains(txt(:,1),'<MATLAB>'), 1) + 1;
            if isempty(DiveIn)
                DiveIn = 1;
            end
            % Check column names:
            Vars = {'Electrode connector', 'Row on electrode connector', ...
                'Pin on electrode connector', 'Type', 'Structure', ...
                'Hemisphere', 'AP', 'ML', 'DV', 'Comments'};
            if ~isequal(Vars, txt(DiveIn,:)) || ~isequal(Vars, raw(DiveIn,:))
                h_error = errordlg(['The file has the wrong format (variables does not match).'], ...
                    'Error', 'modal');
                uiwait(h_error)
                return
            end
            % Find end tag:
            DiveOut = find(contains(txt(:,1),'</MATLAB>'), 1) - 1;
            if isempty(DiveOut)
                DiveOut = min(size(raw,1), size(txt,1));
            end
            % Crop:
            txt = txt(DiveIn:DiveOut,:);
            raw = raw(DiveIn:DiveOut,:);
            % Pin IDs:
            PinID_Electrode = [txt(2:end,1:2) raw(2:end,3)];    % Index is pin, i.e. table row.
            % Wire info:
            WireInfo = [txt(2:end,4:6) raw(2:end,7:9)];
            % Table:
            VarNames = strrep(txt(1,:), ' ', '_');
            wires = cell2table([PinID_Electrode WireInfo txt(2:end,10)], ...
                'VariableNames', VarNames);
            labels = cell(height(wires),1);
            for iWire = 1:numel(labels)
                labels{iWire} = sprintf('%s_%s_%d', ...
                    wires.Electrode_connector{iWire}, ...
                    wires.Row_on_electrode_connector{iWire}, ...
                    wires.Pin_on_electrode_connector(iWire));
            end
            % Define coordinate system:
            systems = {'Paxinos and Watson: The Rat Brain (6th edition)'; 'Waxholm SD v1.01'};
            [sel,ok] = listdlg('PromptString','Database coordinate system:',...
                'SelectionMode','single',...
                'ListString', systems, 'ListSize', [300 200]);
            if ~ok
                return
            end
            co_system = systems{sel};
            % Target coordinates:
            target_org = [wires.ML wires.AP wires.DV]; % These coordinates should be RAS.
            % Valid?
            bNan = any(isnan(target_org),2);
            if all(bNan)
                uialert(app.UIFigure, 'Found no wires with valid coordinates.', 'No wires');
                return
            end
            % Transform to native:
            switch app.RegistrationTarget
                case 'Waxholm SD v1.01'
                    % Target:
                    [target_nat(:,1), target_nat(:,2), target_nat(:,3)] = whs_mm2voxel(...
                        target_org(:,1), target_org(:,2), target_org(:,3)); % Now in WHS voxel.
                    target_nat = target_nat+1; % Voxel to slice.
                    target_nat = round(app.tform.transformPointsInverse(target_nat)); % The Waxholm transform is always defined between native slice and WHS slice.
                case 'Paxinos and Watson: The Rat Brain (6th edition)'
                    % Target:
                    target_nat = round(app.tform.transformPointsInverse(target_org)); % The Paxinos transform is always defined between native slice and Paxinos mm.
            end
            % Calculate start positions:
            co_start = CalculateStart(app, target_nat, DvSlice);
            % Update GUI:
            NewWires(app, co_start(~bNan,:), target_nat(~bNan,:), [], labels(~bNan));
            app.DvSlice = DvSlice;
        end
        
        function LoadWiresDatabaseMenuSelected(app, event)
            % Selet DV plane:
            answer = inputdlg('Please set the DV plane to use as start position:', ...
                'DV start position', 1, {num2str(app.DvSlice)});
            if isempty(answer)
                return
            end
            DvSlice = str2num(answer{1});
            % Get database names:
            dbnames = dbtools.getDatabases;
            [sel,ok] = listdlg('PromptString','Select database:',...
                'SelectionMode','single',...
                'ListString', dbnames);
            if ~ok
                return
            end
            dbname = dbnames{sel};
            % Get implants:
            implants = dbtools.getImplant(dbname);
            [sel,ok] = listdlg('PromptString','Select implant:',...
                'SelectionMode','single',...
                'ListString', {implants.name}');
            if ~ok
                return
            end
            implant = implants(sel);
            % Get wires:
            wires = dbtools.getWire(dbname, 'implant_id', implant.dbid);
            wires = struct2table(wires);
            % Define coordinate system:
            systems = {'Paxinos and Watson: The Rat Brain (6th edition)'; 'Waxholm SD v1.01'; 'Franklin and Paxinos: The Mouse Brain (3rd edition)'};
            [sel,ok] = listdlg('PromptString','Database coordinate system:',...
                'SelectionMode','single',...
                'ListString', systems, 'ListSize', [300 200]);
            if ~ok
                return
            end
            co_system = systems{sel};
            % Target coordinates:
            target_org = [wires.target_ml wires.target_ap wires.target_dv]; % These coordinates should be RAS.
            if iscell(target_org)
                target_org(cellfun(@isempty, target_org)) = {NaN};
                target_org = cell2mat(target_org);
            end
            % CT coordinates:
            ct_org = [wires.ct_ml wires.ct_ap wires.ct_dv]; % These coordinates should be RAS.
            if iscell(ct_org)
                ct_org(cellfun(@isempty, ct_org)) = {NaN};
                ct_org = cell2mat(ct_org);
            end
            % Valid?
            bNan = any(isnan(target_org),2);
            bNanCt = any(isnan(ct_org),2);
            if all(bNan)
                uialert(app.UIFigure, 'Found no wires with valid coordinates.', 'No wires');
                return
            end
            % Transform to native:
            switch app.RegistrationTarget
                case 'Waxholm SD v1.01'
                    % Target:
                    [target_nat(:,1), target_nat(:,2), target_nat(:,3)] = whs_mm2voxel(...
                        target_org(:,1), target_org(:,2), target_org(:,3)); % Now in WHS voxel.
                    target_nat = target_nat+1; % Voxel to slice.
                    target_nat = round(app.tform.transformPointsInverse(target_nat)); % The Waxholm transform is always defined between native slice and WHS slice.
                    % End:
                    if isequal(bNan,bNanCt)
                        [ct_nat(:,1), ct_nat(:,2), ct_nat(:,3)] = whs_mm2voxel(...
                            ct_org(:,1), ct_org(:,2), ct_org(:,3)); % Now in WHS voxel.
                        ct_nat = ct_nat+1; % Voxel to slice.
                        ct_nat = round(app.tform.transformPointsInverse(ct_nat)); % The Waxholm transform is always defined between native slice and WHS slice.
                    else
                        ct_nat = NaN(size(target_nat));
                    end
                case {'Paxinos and Watson: The Rat Brain (6th edition)', 'Franklin and Paxinos: The Mouse Brain (3rd edition)'}
                    % Target:
                    target_nat = round(app.tform.transformPointsInverse(target_org)); % The Paxinos transform is always defined between native slice and Paxinos mm.
                    % End:
                    if isequal(bNan,bNanCt)
                        ct_nat = round(app.tform.transformPointsInverse(ct_org));
                    else
                        ct_nat = NaN(size(target_nat));
                    end
            end
            % Calculate start positions:
            co_start = CalculateStart(app, target_nat, DvSlice);
            % Update GUI:
            NewWires(app, co_start(~bNan,:), target_nat(~bNan,:), ...
                ct_nat(~bNan,:), cellstr(num2str(wires(~bNan,:).dbid)));
            app.DvSlice = DvSlice;
        end
        
        function LoadMarkersMatMenuSelected(app, event)
            [file,path] = uigetfile({'*.mat'}, 'Select MAT file', app.PreferredDir);
            if isequal(file,0)
                return
            end
            app.PreferredDir = path;
            L = load(fullfile(path,file));
            if isfield(L, 'markers')
                % Do we need to calculate the starting points?
                co_start = [L.markers.Start_ML L.markers.Start_AP L.markers.Start_DV];
                if isempty(co_start) || any(isnan(co_start(:))) % We need to calculate the start if it does not exist...
                    % Selet DV plane:
                    answer = inputdlg('Please set the DV plane to use as start position:', ...
                        'DV start position', 1, {num2str(app.DvSlice)});
                    if isempty(answer)
                        return
                    end
                    DvSlice = str2num(answer{1});
                    co_start = app.CalculateStart(co_target, DvSlice);
                end
                % Valid:
                if ismember('exist',L.markers.Properties.VariableNames)
                    valid = L.markers.exist;
                else
                    valid = true(height(L.markers),1);
                end
                % Update GUI:
                NewWires(app, co_start, ...
                    [L.markers.Target_ML L.markers.Target_AP L.markers.Target_DV], ...
                    [L.markers.End_ML L.markers.End_AP L.markers.End_DV], ...
                    cellstr(num2str(L.markers.dbid)), ...
                    round(median(L.markers.Start_DV)));
            else
                uialert(app.UIFigure, 'Could not find a ''markers'' variable in this file.', 'Wrong file format');
                return
            end
            if isfield(L, 'tform')
                answer = uiconfirm(app.UIFigure, 'This file contains a transformation matrix. Do you want to use it?', 'Load transform',...
                    "Options", {'Yes' 'No'}, "DefaultOption", 1);
                if answer=="Yes"
                    % Registration target:
                    if isfield(L, 'RegistrationTarget')
                        RegistrationTarget = CheckRegistrationTarget(app, L.RegistrationTarget);
                    else
                        RegistrationTarget = CheckRegistrationTarget(app, '');
                    end
                    if isempty(RegistrationTarget)
                        return
                    else
                        app.RegistrationTarget = RegistrationTarget;
                    end
                    % Transform:
                    app.tform = L.tform;
                    app.State_Transformation = true;
                    UpdateState(app)
                end
            end
        end
        
        function SaveMarkersMenuSelected(app, event)
            [file, path] = uiputfile({'*.mat'}, 'Save markers', app.PreferredDir);
            if isequal(file,0)
                return
            else
                app.PreferredDir = path;
                app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Saving markers...', ...
                    'Indeterminate','on');
                if isempty(app.EndMarkers)
                    EndMarkers = NaN(size(app.StartMarkers));
                else
                    EndMarkers = app.EndMarkers;
                end
                markers = table(cellfun(@str2num,app.WireLabels), ...
                    app.TargetMarkers(:,1), ...
                    app.TargetMarkers(:,2), ...
                    app.TargetMarkers(:,3), ...
                    app.StartMarkers(:,1), ...
                    app.StartMarkers(:,2), ...
                    app.StartMarkers(:,3), ...
                    EndMarkers(:,1), ...
                    EndMarkers(:,2), ...
                    EndMarkers(:,3), ...
                    'VariableNames', {'dbid' 'Target_ML' 'Target_AP' 'Target_DV' ...
                                             'Start_ML' 'Start_AP' 'Start_DV' ...
                                             'End_ML' 'End_AP' 'End_DV'});
                tform = app.tform;
                RegistrationTarget = app.RegistrationTarget;
                save(fullfile(path, file), 'markers',  'tform', 'RegistrationTarget');
                disp('File saved!')
                delete(app.h_progress)
            end
        end
        
        function SaveToDatabaseMenuSelected(app, event)
            % Get database names:
            dbnames = dbtools.getDatabases;
            [sel,ok] = listdlg('PromptString','Select database:',...
                'SelectionMode','single',...
                'ListString', dbnames);
            if ~ok
                return
            end
            dbname = dbnames{sel};
            % Transformation:
            co = app.tform.transformPointsForward(app.EndMarkers);
            ML = co(:,1);
            AP = co(:,2);
            DV = co(:,3);
            % Write:
            sides = containers.Map( {-1,1}, {'left', 'right'});
            app.h_progress = uiprogressdlg(app.UIFigure,'Title','Please wait...', 'Message', 'Writing to database...');
            for iWire = 1:numel(app.WireLabels)
                if isempty(dbtools.getWire(dbname, str2num(app.WireLabels{iWire})))
                    warning(sprintf('Could not find %s in the database.', app.WireLabels{iWire}))
                else
                    dbtools.set(dbname, 'wires', app.WireLabels{iWire}, ...
                        'ct_ml', sprintf('%0.2f', ML(iWire)));
                    dbtools.set(dbname, 'wires', app.WireLabels{iWire}, ...
                        'ct_ap', sprintf('%0.2f', AP(iWire)));
                    dbtools.set(dbname, 'wires', app.WireLabels{iWire}, ...
                        'ct_dv', sprintf('%0.2f', DV(iWire)));
                    dbtools.set(dbname, 'wires', app.WireLabels{iWire}, ...
                        'side', sides( sign( ML(iWire))) );
                end
                app.h_progress.Value = iWire/numel(app.WireLabels);
            end
            delete(app.h_progress)
        end
        
        function CreateMaskMenuSelected(app, event)
            % Set mask size:
            answer = inputdlg(...
                {'Enter mask size (approximately the radius in pixels):', 'Enter filter radius (decrease if wires "steal" from each other):'}, ...
                'Mask settings', 1, {num2str(app.electroslice_MaskSize) num2str(app.GaussianFilterSize)});
            if isempty(answer)
                return
            end
            app.electroslice_MaskSize = str2num(answer{1});
            app.GaussianFilterSize = str2num(answer{2});
            % Create mask:
            app.electroslice_mask = app.CreateMask(app.electroslice_MaskSize);
            app.bShowMask = true;
            app.ShowMaskCheckBox.Value = true;
            app.State_Mask = true;
            app.State_SegmentationInit = false;
            app.State_SegmentationExtend = false;
            UpdateState(app);
            UpdateImage(app);
        end
        
        function InitializeSegmentationMenuSelected(app, event)
            % Set limits:
            DvSlice = app.StartMarkers(1,3);
%             if isempty(app.electroslice_limits)
%                 app.electroslice_limits = stretchlim(app.Img(:,:,DvSlice));
                 app.electroslice_limits = app.DisplayRange/100; % attempted fix
%             end
            answer = inputdlg(...
                {'Low intensity limit [0 1]:', 'High intensity limit [0 1]:', 'Minimum wire distance:'}, ...
                'Settings', 1, {num2str(app.electroslice_limits(1)) ...
                num2str(app.electroslice_limits(2)) num2str(app.electroslice_min_dist)});
            if isempty(answer)
                return
            end
            app.electroslice_limits(1) = str2num(answer{1});
            app.electroslice_limits(2) = str2num(answer{2});
            app.electroslice_min_dist = str2num(answer{3});
            % Create Bundle:
            app.electroslice_bundle = app.ElectroSliceBundle;
            app.State_SegmentationInit = ~isempty(app.electroslice_bundle.bundle_index); % False if initialization failed.
            app.State_SegmentationExtend = false;
            UpdateState(app);
        end
        
        function PlotSegmentationLabelsMenuSelected(app, event)
          app.electroslice_bundle.plot_labels();
        end
                
        function ExtendSegmentationMenuSelected(app, event)
            markers = app.ElectroSliceExtend(app.electroslice_bundle);
            app.State_SegmentationExtend = true;
            app.EndMarkers = [markers.ML markers.AP markers.DV];
            app.State_WiresEnd = true;
            UpdateState(app);
        end
        
        function PlotSegmentationTipsMenuSelected(app, event)
          % update wire tip positions reading from markers
          num_labels = cellfun( @str2num, app.WireLabels);
          for this_wire = [app.electroslice_bundle.wires]
            % get corresponding marker position
            pos = app.EndMarkers( num_labels == this_wire.wire_label, :);
            % recreate cube
            this_wire.tip_coords = this_wire.get_tip_cube( pos);
          end
          % now plot
          app.electroslice_bundle.plot_tips();
        end
        
        function PlotSegmentationOverlapMenuSelected(app, event)
          % update wire tip positions reading from markers
          num_labels = cellfun( @str2num, app.WireLabels);
          for this_wire = [app.electroslice_bundle.wires]
            % get corresponding marker position
            pos = app.EndMarkers( num_labels == this_wire.wire_label, :);
            % recreate cube
            this_wire.tip_coords = this_wire.get_tip_cube( pos);
          end
          % now plot
          wires = app.electroslice_bundle.check_overlap();
          if numel(wires)==1 && wires==0
              uialert(app.UIFigure, 'No ovelaps to plot.', 'Plot overlaps', 'icon', 'info');
          end
        end
        
        function DeleteMarkerMenuSelected(app, event)
            bSel = [app.points.Selected];
            if any(bSel)
                app.WireLabels = app.WireLabels(~bSel);
                app.StartMarkers = app.StartMarkers(~bSel,:);
                app.TargetMarkers = app.TargetMarkers(~bSel,:);
                if ~isempty(app.EndMarkers)
                    app.EndMarkers = app.EndMarkers(~bSel,:);
                end
                NewMarkers(app)
            end
        end
        
        function AddMarkerMenuSelected(app, event)
            h_fig = uifigure('Position',[100 100 437 317]);
            g = uigridlayout(h_fig, [7 4]);
            
            ImgSize = size(app.Img);
            
            % Label:
            NewLabel = min(min(str2double(app.WireLabels),0))-1;
            h = uilabel(g);
            h.Layout.Row = 1;
            h.Layout.Column = 1;
            h.Text = 'Label:';
            h.HorizontalAlignment = 'right';
            h_label = uieditfield(g,'text', 'Value', num2str(NewLabel));
            h_label.Layout.Row = 1;
            h_label.Layout.Column = 2;
            h_label.Enable = 'off';
            
            % ML label:
            h = uilabel(g);
            h.Layout.Row = 2;
            h.Layout.Column = 2;
            h.Text = 'ML';
            h.VerticalAlignment = 'bottom';
            h.HorizontalAlignment = 'center';
            
            % AP label:
            h = uilabel(g);
            h.Layout.Row = 2;
            h.Layout.Column = 3;
            h.Text = 'AP';
            h.VerticalAlignment = 'bottom';
            h.HorizontalAlignment = 'center';
            
            % DV label:
            h = uilabel(g);
            h.Layout.Row = 2;
            h.Layout.Column = 4;
            h.Text = 'DV';
            h.VerticalAlignment = 'bottom';
            h.HorizontalAlignment = 'center';
            
            % Target label:
            h = uilabel(g);
            h.Layout.Row = 3;
            h.Layout.Column = 1;
            h.Text = 'Target:';
            h.VerticalAlignment = 'center';
            h.HorizontalAlignment = 'right';
            
            % Start point label:
            h = uilabel(g);
            h.Layout.Row = 4;
            h.Layout.Column = 1;
            h.Text = 'Start point:';
            h.VerticalAlignment = 'center';
            h.HorizontalAlignment = 'right';
            
            % End point label:
            if app.State_WiresEnd
                h = uilabel(g);
                h.Layout.Row = 5;
                h.Layout.Column = 1;
                h.Text = 'End point:';
                h.VerticalAlignment = 'center';
                h.HorizontalAlignment = 'right';
            end
            
            % Target
            iRow = 3;
            h_target_ML = uieditfield(g,'numeric', 'Limits', [1 ImgSize(1)], 'Value', 1);
            h_target_ML.Layout.Row = iRow;
            h_target_ML.Layout.Column = 2;
            h_target_AP = uieditfield(g,'numeric', 'Limits', [1 ImgSize(2)], 'Value', 1);
            h_target_AP.Layout.Row = iRow;
            h_target_AP.Layout.Column = 3;
            h_target_DV = uieditfield(g,'numeric', 'Limits', [1 ImgSize(3)], 'Value', 1);
            h_target_DV.Layout.Row = iRow;
            h_target_DV.Layout.Column = 4;
            
            % Start
            iRow = 4;
            h_start_ML = uieditfield(g,'numeric', 'Limits', [1 ImgSize(1)], 'Value', 1);
            h_start_ML.Layout.Row = iRow;
            h_start_ML.Layout.Column = 2;
            h_start_AP = uieditfield(g,'numeric', 'Limits', [1 ImgSize(2)], 'Value', 1);
            h_start_AP.Layout.Row = iRow;
            h_start_AP.Layout.Column = 3;
            h_start_DV = uieditfield(g,'numeric', 'Limits', [1 ImgSize(3)], 'Value', 1);
            h_start_DV.Layout.Row = iRow;
            h_start_DV.Layout.Column = 4;
            
            % End
            if app.State_WiresEnd
                iRow = 5;
                h_end_ML = uieditfield(g,'numeric', 'Limits', [1 ImgSize(1)], 'Value', 1);
                h_end_ML.Layout.Row = iRow;
                h_end_ML.Layout.Column = 2;
                h_end_AP = uieditfield(g,'numeric', 'Limits', [1 ImgSize(2)], 'Value', 1);
                h_end_AP.Layout.Row = iRow;
                h_end_AP.Layout.Column = 3;
                h_end_DV = uieditfield(g,'numeric', 'Limits', [1 ImgSize(3)], 'Value', 1);
                h_end_DV.Layout.Row = iRow;
                h_end_DV.Layout.Column = 4;
                h_end_ML.Enable = 'off';
                h_end_AP.Enable = 'off';
                h_end_DV.Enable = 'off';
            end
            
            % Buttons:
            h_Ok = uibutton(g);
            h_Ok.Layout.Row = 7;
            h_Ok.Layout.Column = 4;
            h_Ok.Text = 'OK';
            h_Ok.ButtonPushedFcn = @ButtonCallback;
            
            h_Cancel = uibutton(g);
            h_Cancel.Layout.Row = 7;
            h_Cancel.Layout.Column = 3;
            h_Cancel.Text = 'Cancel';
            h_Cancel.ButtonPushedFcn = @ButtonCallback;
            
            uiwait(h_fig)
            
            function ButtonCallback(src, event)
                uiresume(h_fig)
                if src.Text=="OK"
                    app.WireLabels = [app.WireLabels; h_label.Value];
                    app.TargetMarkers = [app.TargetMarkers; ...
                        h_target_ML.Value h_target_AP.Value h_target_DV.Value];
                    app.StartMarkers = [app.StartMarkers; ...
                        h_start_ML.Value h_start_AP.Value h_start_DV.Value];
                    if app.State_WiresEnd
                        app.EndMarkers = [app.EndMarkers; ...
                            h_end_ML.Value h_end_AP.Value h_end_DV.Value];
                    end
                    NewMarkers(app)
                    app.SelectMarkerDropDown.Items = [{'-'}; app.WireLabels];
%                     UpdateMarkerTable(app);
                end
                close(h_fig)
            end
            
        end
        
        function SetDvMenuSelected(app, event)
            DvSlice = app.StartMarkers(1,3);
            answer = inputdlg({'New DV coordinate:'}, 'Settings', 1, {num2str(DvSlice)});
            if isempty(answer)
                return
            end
            app.StartMarkers(:,3) = str2num(answer{1});
            NewMarkers(app)
            app.DvSlice = str2num(answer{1});
        end
        
        function SetGaussianFilterSizeMenuSelected(app, event)
            answer = inputdlg({'Enter filter radius (decrease if wires "steal" from each other):'}, ...
                'Gaussian filter', 1, {num2str(app.GaussianFilterSize)});
            if isempty(answer)
                return
            end
            app.GaussianFilterSize = str2num(answer{1});
        end
        
        function ShowMarkerLabelsMenuSelected(app, event)
            app.bShowMarkerLabels = ~app.bShowMarkerLabels;
            if app.bShowMarkerLabels
                [app.points.Label] = deal(app.WireLabels{:});
                [app.points.LabelAlpha] = deal(0.7);
                event.Source.Checked = 'on';
            else
                [app.points.Label] = deal('');
                event.Source.Checked = 'off';
            end
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
            app.UIFigure.Name = 'identify_electrodes';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @AxesButtonDownFcn, true);
            app.UIFigure.WindowScrollWheelFcn = createCallbackFcn( app, @ScrollWheelFcn, true); % NEW!

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadimageMenu
            app.LoadimageMenu = uimenu(app.FileMenu);
            app.LoadimageMenu.MenuSelectedFcn = createCallbackFcn(app, @Loadimage, true);
            app.LoadimageMenu.Text = 'Load image';

            % Create LoadtransformationMenu
            app.LoadtransformationMenu = uimenu(app.FileMenu);
            app.LoadtransformationMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadtransformationMenuSelected, true);
            app.LoadtransformationMenu.Separator = 'off';
            app.LoadtransformationMenu.Enable = 'off';
            app.LoadtransformationMenu.Text = 'Load transformation';
            
            % Create LoadWiresXlsMenu
            app.LoadWiresXlsMenu = uimenu(app.FileMenu);
            app.LoadWiresXlsMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadWiresXlsMenuSelected, true);
            app.LoadWiresXlsMenu.Separator = 'on';
            app.LoadWiresXlsMenu.Enable = 'off';
            app.LoadWiresXlsMenu.Text = 'Load wires from .xls file';
            
            % Create LoadWiresDatabaseMenu
            app.LoadWiresDatabaseMenu = uimenu(app.FileMenu);
            app.LoadWiresDatabaseMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadWiresDatabaseMenuSelected, true);
            app.LoadWiresDatabaseMenu.Enable = 'off';
            app.LoadWiresDatabaseMenu.Text = 'Load wires from database';
            
            % Create LoadMarkersMatMenu
            app.LoadMarkersMatMenu = uimenu(app.FileMenu);
            app.LoadMarkersMatMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadMarkersMatMenuSelected, true);
            app.LoadMarkersMatMenu.Enable = 'off';
            app.LoadMarkersMatMenu.Text = 'Load wires from .mat file';
            app.LoadMarkersMatMenu.Separator = 'off';
            
            % Create SaveMarkersMenu
            app.SaveMarkersMenu = uimenu(app.FileMenu);
            app.SaveMarkersMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveMarkersMenuSelected, true);
            app.SaveMarkersMenu.Enable = 'off';
            app.SaveMarkersMenu.Text = 'Save wires to .mat file';
            app.SaveMarkersMenu.Separator = 'on';

            % Create SaveToDatabaseMenu
            app.SaveToDatabaseMenu = uimenu(app.FileMenu);
            app.SaveToDatabaseMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveToDatabaseMenuSelected, true);
            app.SaveToDatabaseMenu.Enable = 'off';
            app.SaveToDatabaseMenu.Text = 'Save end points to database';
            
            % Create SegmentationMenu
            app.SegmentationMenu = uimenu(app.UIFigure);
            app.SegmentationMenu.Text = 'Segmentation';
            
            % Create CreateMaskMenu
            app.CreateMaskMenu = uimenu(app.SegmentationMenu);
            app.CreateMaskMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @CreateMaskMenuSelected, true);
            app.CreateMaskMenu.Enable = 'off';
            app.CreateMaskMenu.Text = 'Step 1: Create mask';
            
            % Create InitializeSegmentationMenu
            app.InitializeSegmentationMenu = uimenu(app.SegmentationMenu);
            app.InitializeSegmentationMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @InitializeSegmentationMenuSelected, true);
            app.InitializeSegmentationMenu.Enable = 'off';
            app.InitializeSegmentationMenu.Text = 'Step 2: Initialize';
            
            % Create PlotSegmentationLabelsMenu
            app.PlotSegmentationLabelsMenu = uimenu(app.SegmentationMenu);
            app.PlotSegmentationLabelsMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @PlotSegmentationLabelsMenuSelected, true);
            app.PlotSegmentationLabelsMenu.Enable = 'off';
            app.PlotSegmentationLabelsMenu.Text = 'Plot wire labels';
            
            % Create ExtendSegmentationMenu
            app.ExtendSegmentionMenu = uimenu(app.SegmentationMenu);
            app.ExtendSegmentionMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @ExtendSegmentationMenuSelected, true);
            app.ExtendSegmentionMenu.Enable = 'off';
            app.ExtendSegmentionMenu.Text = 'Step 3: Extend';
            
            % Create PlotSegmentationTipsMenu
            app.PlotSegmentationTipsMenu = uimenu(app.SegmentationMenu);
            app.PlotSegmentationTipsMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @PlotSegmentationTipsMenuSelected, true);
            app.PlotSegmentationTipsMenu.Enable = 'off';
            app.PlotSegmentationTipsMenu.Text = 'Plot wire tips';
            
            % Create PlotSegmentationOverlapMenu
            app.PlotSegmentationOverlapMenu = uimenu(app.SegmentationMenu);
            app.PlotSegmentationOverlapMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @PlotSegmentationOverlapMenuSelected, true);
            app.PlotSegmentationOverlapMenu.Enable = 'off';
            app.PlotSegmentationOverlapMenu.Text = 'Plot wire tip overlaps';
            
            % Create MarkerMenu
            app.MarkerMenu = uimenu(app.UIFigure);
            app.MarkerMenu.Text = 'Marker';
            
            % Create MarkerDeleteMenu
            app.MarkerDeleteMenu = uimenu(app.MarkerMenu);
            app.MarkerDeleteMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @DeleteMarkerMenuSelected, true);
            app.MarkerDeleteMenu.Enable = 'off';
            app.MarkerDeleteMenu.Text = 'Delete selected markers';
            
            % Create MarkerAddMenu
            app.MarkerAddMenu = uimenu(app.MarkerMenu);
            app.MarkerAddMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @AddMarkerMenuSelected, true);
            app.MarkerAddMenu.Enable = 'off';
            app.MarkerAddMenu.Text = 'Add marker';
            
            % Create SetDvMenu
            app.SetDvMenu = uimenu(app.MarkerMenu);
            app.SetDvMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @SetDvMenuSelected, true);
            app.SetDvMenu.Enable = 'off';
            app.SetDvMenu.Text = 'Change DV coordinate for all starting points';
            
            % Create SetGaussianFilterSizeMenu
            app.SetGaussianFilterSizeMenu = uimenu(app.MarkerMenu);
            app.SetGaussianFilterSizeMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @SetGaussianFilterSizeMenuSelected, true);
            app.SetGaussianFilterSizeMenu.Enable = 'on';
            app.SetGaussianFilterSizeMenu.Text = 'Set Gaussian filter size (used for marker adjustment)';

            % Create ShowMarkerLabelsMenu
            app.ShowMarkerLabelsMenu = uimenu(app.MarkerMenu);
            app.ShowMarkerLabelsMenu.MenuSelectedFcn = ...
                createCallbackFcn(app, @ShowMarkerLabelsMenuSelected, true);
            app.ShowMarkerLabelsMenu.Enable = 'on';
            app.ShowMarkerLabelsMenu.Checked = 'off';
            app.ShowMarkerLabelsMenu.Text = 'Show marker labels';
            app.ShowMarkerLabelsMenu.Separator = 'on';
            
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

            
            y = 30;
            
            % Create ViewPanel
            app.ViewPanel = uipanel(app.LeftPanel);
            app.ViewPanel.ForegroundColor = [1 1 1];
            app.ViewPanel.BorderType = 'none';
            app.ViewPanel.Title = 'View';
            app.ViewPanel.BackgroundColor = [0.149 0.149 0.149];
            app.ViewPanel.FontWeight = 'bold';
            app.ViewPanel.Position = [1 6 191 340+2*y];

            % Create CursorLabel
            app.CursorLabel = uilabel(app.ViewPanel);
            app.CursorLabel.FontColor = [1 1 1];
            app.CursorLabel.Position = [14 104+80+5.5*y 120 22];
            app.CursorLabel.Text = 'Selected marker:';

            % Create UITable
            app.UITable = uitable(app.ViewPanel);
            app.UITable.ColumnName = {'ML'; 'AP'; 'DV'};
            app.UITable.ColumnWidth = {39, 39, 39};
            app.UITable.RowName = {'Slice'; 'Atlas'};
            app.UITable.ColumnEditable = false;
            app.UITable.RowStriping = 'off';
            app.UITable.Tooltip = {'Cursor position:'; '- Slice: Native slice number.'; '- Atlas: Estimated atlas coordinate.'};
            app.UITable.Enable = 'off';
            app.UITable.FontSize = 11;
            app.UITable.Position = [13 104+5.5*y 168 80];

            % Create PaxinosButton
            app.PaxinosButton = uibutton(app.ViewPanel, 'push');
            app.PaxinosButton.ButtonPushedFcn = createCallbackFcn(app, @PaxinosButtonPushed, true);
            app.PaxinosButton.BackgroundColor = [0.302 0.302 0.302];
            app.PaxinosButton.FontColor = [1 1 1];
            app.PaxinosButton.Enable = 'off';
            app.PaxinosButton.Position = [91 104+4.5*y 90 22];
            app.PaxinosButton.Text = 'View in atlas';
            
            % Create DvSliceLabel
            app.DvSliceLabel = uilabel(app.ViewPanel);
            app.DvSliceLabel.FontColor = [1 1 1];
            app.DvSliceLabel.Position = [14 104+3*y 70 22];
            app.DvSliceLabel.Text = 'DV slice:';
            
            % Create DvSliceSpinner
            app.DvSliceSpinner = uispinner(app.ViewPanel);
            app.DvSliceSpinner.RoundFractionalValues = 'on';
            app.DvSliceSpinner.ValueChangedFcn = createCallbackFcn(app, @DvSliceSpinnerValueChanged, true);
            app.DvSliceSpinner.FontColor = [1 1 1];
            app.DvSliceSpinner.BackgroundColor = [0.302 0.302 0.302];
            app.DvSliceSpinner.Tooltip = {'Set DV slice to view.'};
            app.DvSliceSpinner.Position = [91 104+3*y 68 22];
            app.DvSliceSpinner.Enable = 'off';
            
            % Create DvSliceStepLabel
            app.DvSliceStepLabel = uilabel(app.ViewPanel);
            app.DvSliceStepLabel.FontColor = [1 1 1];
            app.DvSliceStepLabel.Position = [14 104+2*y 70 22];
            app.DvSliceStepLabel.Text = 'DV step:';
            
            % Create DvSliceStepSpinner
            app.DvSliceStepSpinner = uispinner(app.ViewPanel);
            app.DvSliceStepSpinner.RoundFractionalValues = 'on';
            app.DvSliceStepSpinner.ValueChangedFcn = createCallbackFcn(app, @DvSliceStepSpinnerValueChanged, true);
            app.DvSliceStepSpinner.FontColor = [1 1 1];
            app.DvSliceStepSpinner.BackgroundColor = [0.302 0.302 0.302];
            app.DvSliceStepSpinner.Tooltip = {'Set DV step.'};
            app.DvSliceStepSpinner.Position = [91 104+2*y 68 22];
            app.DvSliceStepSpinner.Value = 1;
            app.DvSliceStepSpinner.Limits = [1 Inf];
            
            % Create ShowMarkersCheckBox
            app.ShowMarkersCheckBox = uicheckbox(app.ViewPanel);
            app.ShowMarkersCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowMarkersCheckBoxValueChanged, true);
            app.ShowMarkersCheckBox.FontColor = [1 1 1];
            app.ShowMarkersCheckBox.Enable = 'off';
            app.ShowMarkersCheckBox.Position = [14 104+y 150 22];
            app.ShowMarkersCheckBox.Text = 'Show markers';
            app.ShowMarkersCheckBox.Value = app.bShowMarkers;
            
            % Create ShowMaskCheckBox
            app.ShowMaskCheckBox = uicheckbox(app.ViewPanel);
            app.ShowMaskCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowMaskCheckBoxValueChanged, true);
            app.ShowMaskCheckBox.FontColor = [1 1 1];
            app.ShowMaskCheckBox.Position = [14 104 150 22];
            app.ShowMaskCheckBox.Enable = 'off';
            app.ShowMaskCheckBox.Text = 'Show mask';
            app.ShowMaskCheckBox.Value = app.bShowMask;
            
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

            % Create MarkersPanel
            app.MarkersPanel = uipanel(app.LeftPanel);
            app.MarkersPanel.ForegroundColor = [1 1 1];
            app.MarkersPanel.BorderType = 'none';
            app.MarkersPanel.Title = 'Wires';
            app.MarkersPanel.BackgroundColor = [0.149 0.149 0.149];
            app.MarkersPanel.FontWeight = 'bold';
            app.MarkersPanel.Position = [1 420 191 232];

            y = 40;
            
            % Create SelectMarkersButton
            app.SelectMarkersButton = uibutton(app.MarkersPanel, 'push');
            app.SelectMarkersButton.ButtonPushedFcn = createCallbackFcn(app, @SelectMarkersButtonPushed, true);
            app.SelectMarkersButton.BackgroundColor = [0.302 0.302 0.302];
            app.SelectMarkersButton.FontColor = [1 1 1];
            app.SelectMarkersButton.Enable = 'off';
            app.SelectMarkersButton.Tooltip = {'Select markers by drawing a square.'};
            app.SelectMarkersButton.Position = [34 77-y 133 22];
            app.SelectMarkersButton.Text = 'Select markers';
            
            % Create AdjustMarkersButton
            app.AdjustMarkersButton = uibutton(app.MarkersPanel, 'push');
            app.AdjustMarkersButton.ButtonPushedFcn = createCallbackFcn(app, @AdjustMarkersButtonPushed, true);
            app.AdjustMarkersButton.BackgroundColor = [0.302 0.302 0.302];
            app.AdjustMarkersButton.FontColor = [1 1 1];
            app.AdjustMarkersButton.Tooltip = {'Adjust the positions of selected markers to the nearest bright spot.'};
            app.AdjustMarkersButton.Position = [34 49-y 133 22];
            app.AdjustMarkersButton.Text = 'Adjust markers';
            app.AdjustMarkersButton.Enable = 'off';

            % Create APLabel
            app.APLabel = uilabel(app.MarkersPanel);
            app.APLabel.HorizontalAlignment = 'center';
            app.APLabel.FontColor = [1 1 1];
            app.APLabel.Position = [89 128-y 25 22];
            app.APLabel.Text = 'AP';

            % Create APEditField_2
            app.APEditField_2 = uieditfield(app.MarkersPanel, 'numeric');
            app.APEditField_2.RoundFractionalValues = 'on';
            app.APEditField_2.ValueChangedFcn = createCallbackFcn(app, @MarkerEditFieldValueChanged, true);
            app.APEditField_2.FontColor = [1 1 1];
            app.APEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.APEditField_2.Tooltip = {'Anterior-posterior slice number of currently selected marker'};
            app.APEditField_2.Position = [78 107-y 45 22];
            app.APEditField_2.Enable = 'off';

            % Create MLEditField_2
            app.MLEditField_2 = uieditfield(app.MarkersPanel, 'numeric');
            app.MLEditField_2.RoundFractionalValues = 'on';
            app.MLEditField_2.ValueChangedFcn = createCallbackFcn(app, @MarkerEditFieldValueChanged, true);
            app.MLEditField_2.FontColor = [1 1 1];
            app.MLEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.MLEditField_2.Tooltip = {'Medial-lateral slice number of currently selected marker'};
            app.MLEditField_2.Position = [34 107-y 45 22];
            app.MLEditField_2.Enable = 'off';

            % Create DVEditField_2
            app.DVEditField_2 = uieditfield(app.MarkersPanel, 'numeric');
            app.DVEditField_2.RoundFractionalValues = 'on';
            app.DVEditField_2.ValueChangedFcn = createCallbackFcn(app, @MarkerEditFieldValueChanged, true);
            app.DVEditField_2.FontColor = [1 1 1];
            app.DVEditField_2.BackgroundColor = [0.302 0.302 0.302];
            app.DVEditField_2.Tooltip = {'Dorsal-ventral slice number of currently selected marker'};
            app.DVEditField_2.Position = [122 107-y 45 22];
            app.DVEditField_2.Enable = 'off';

            % Create MLLabel
            app.MLLabel = uilabel(app.MarkersPanel);
            app.MLLabel.HorizontalAlignment = 'center';
            app.MLLabel.FontColor = [1 1 1];
            app.MLLabel.Position = [45 128-y 25 22];
            app.MLLabel.Text = 'ML';

            % Create DVLabel
            app.DVLabel = uilabel(app.MarkersPanel);
            app.DVLabel.HorizontalAlignment = 'center';
            app.DVLabel.FontColor = [1 1 1];
            app.DVLabel.Position = [133 128-y 25 22];
            app.DVLabel.Text = 'DV';

            % Create SelectMarkerDropDownLabel
            app.SelectMarkerDropDownLabel = uilabel(app.MarkersPanel);
            app.SelectMarkerDropDownLabel.FontColor = [1 1 1];
            app.SelectMarkerDropDownLabel.Position = [14 177+10 110 22];
            app.SelectMarkerDropDownLabel.Text = 'Select wire:';

            % Create SelectMarkerDropDown
            app.SelectMarkerDropDown = uidropdown(app.MarkersPanel);
            app.SelectMarkerDropDown.Items = {'-'};
            app.SelectMarkerDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectMarkerDropDownValueChanged, true);
            app.SelectMarkerDropDown.Tooltip = {'Select wire'};
            app.SelectMarkerDropDown.FontColor = [1 1 1];
            app.SelectMarkerDropDown.BackgroundColor = [0.302 0.302 0.302];
            app.SelectMarkerDropDown.Position = [14 155+10 171 22];
            app.SelectMarkerDropDown.Value = '-';
            app.SelectMarkerDropDown.Enable = 'off';
            
            % Create SelectMarkerTypeDropDownLabel
            app.SelectMarkerTypeDropDownLabel = uilabel(app.MarkersPanel);
            app.SelectMarkerTypeDropDownLabel.FontColor = [1 1 1];
            app.SelectMarkerTypeDropDownLabel.Position = [14 177-y 110 22];
            app.SelectMarkerTypeDropDownLabel.Text = 'Select marker type:';

            % Create SelectMarkerTypeDropDown
            app.SelectMarkerTypeDropDown = uidropdown(app.MarkersPanel);
            app.SelectMarkerTypeDropDown.Items = {'-'};
            app.SelectMarkerTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectMarkerTypeDropDownValueChanged, true);
            app.SelectMarkerTypeDropDown.Tooltip = {'Select marker type'};
            app.SelectMarkerTypeDropDown.FontColor = [1 1 1];
            app.SelectMarkerTypeDropDown.BackgroundColor = [0.302 0.302 0.302];
            app.SelectMarkerTypeDropDown.Position = [14 155-y 171 22];
            app.SelectMarkerTypeDropDown.Value = '-';
            app.SelectMarkerTypeDropDown.Enable = 'off';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.BackgroundColor = [0.149 0.149 0.149];
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            
            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            app.UIAxes.BackgroundColor = [0.149 0.149 0.149];
            app.UIAxes.Color = app.UIAxes.BackgroundColor;
            app.UIAxes.Position = [10 10 860 618];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = identify_electrodes(varargin)

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