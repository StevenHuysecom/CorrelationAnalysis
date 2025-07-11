classdef DDMMovie < Core.MPMovie
properties (SetAccess = 'private')
    
    IsEnoughRam
    IsCudaDevice

end
properties (SetAccess ='public')
    DDMOutput
    AllFrames
    DDMInfo
end
methods 
    
%%  Instantiation
        function obj = DDMMovie(raw, MPCal, info,DDMInfo)
             obj  = obj@Core.MPMovie(raw,MPCal,info); 

             assert(~isempty(DDMInfo), 'No MovieInfo provided, check parameter input'); 
             
             obj.DDMInfo = DDMInfo;

             try
                 gpuDevice;
                 obj.IsCudaDevice = true;         
             catch
                 obj.IsCudaDevice = false;
             end

        end
        
%%  Analysis
        function  [RadialValueInQSpace,ValidIndeces]=Get3DGrid(obj,sizes,critangle)
            
        % Prepares the grid for 3d averaging, with proper values for the
        % wavectors.
            if length(obj.calibrated{1}.oRelZPos)>1
              zpixel = abs(mean(diff(obj.calibrated{1}.oRelZPos)));     
            else
              zpixel = 1;   
            end

            for i=1:length(sizes)
                x{i} = 2*pi*(-round(sizes(i)/2)+1:1:round(sizes(i)/2)-1)*1/obj.DDMInfo.PixelSize*1/sizes(i);
                %y{i} = 2*pi./([1:1:round(sizes(i)./2), round(sizes(i)./2)-1:-1:1]*obj.DDMInfo.PixelSize);
            end
            if sizes(end)~=1
                z = 2*pi*(-round(sizes(end)/2):1:round(sizes(end)/2)-1)*(1/zpixel)*1/sizes(end);
            else
                z = 0;   
            end
                
            [X,Y,Z] = meshgrid(x{2},x{1},z);
            RadialValueInQSpace = sqrt(X.^2 + Y.^2 + Z.^2);
            ValidIndeces =  Z./sqrt(X.^2 + Y.^2 + Z.^2)<abs(cosd(critangle));
        end
        
        function [ROI] = DefineROI(obj)
            ROIFile = append(obj.raw.movInfo.Path, filesep, 'ROI.mat');
            
            if obj.DDMInfo.MultiPlane == 1
                if exist(ROIFile)
                    ROI = load(ROIFile);
                    ROI = ROI.ROI;
                else
                    for i = 1:8
                        if i < 5
                            frame = obj.getFrame(1);
                            x = 1;
                        else
                            frame = obj.getFrame(2);
                            x = 2;
                        end
                        % h = questdlg('Do you want to use a ROI?','Question to user','Yes','No', 'Yes');
                        frame = obj.getFrame(1);
                        if strcmp(obj.DDMInfo.useROI  ,'on')
                            f = figure
            
                            imagesc(frame(:,:,1))
                            title(append('Camera ', num2str(x)))
                            hold on
                            if and(i ~= 1, i ~=5)
                                if x == 1
                                    for z = 1:size(ROI, 1)
                                        rectangle('Position',ROI{z,1});
                                    end
                                elseif x == 2
                                    for z = 5:size(ROI, 1)
                                        rectangle('Position',ROI{z,1});
                                    end
                                end
                            end
                            test = drawrectangle();
                            ROI{i,1}  = round(test.Position);
                            close(f)
            
                        else
                            ROI{i,1} = [1,1,size(frame, 1)-1, size(frame,2)-1];  
                        end
                    end
                end
                
                Filename = append(obj.raw.movInfo.Path, filesep, 'ROI.mat');
                save(Filename, 'ROI');
            else
                if exist(ROIFile)
                    ROI = load(ROIFile);
                    ROI = ROI.ROI;
                else
                    % h = questdlg('Do you want to use a ROI?','Question to user','Yes','No', 'Yes');
                    frame = obj.getFrame(1);
                    if strcmp(obj.DDMInfo.useROI  ,'on')
                        figure
        
                        imagesc(frame(:,:,1))
                        test = drawrectangle();
                        ROI  = round(test.Position);
        
                    else
                        ROI = [1,1,size(frame, 1)-1, size(frame,2)-1];  
                    end
        
                    Filename = append(obj.raw.movInfo.Path, filesep, 'ROI.mat');
                    save(Filename, 'ROI');
                end
            end
        end

        function [FullMask, SegmentedMask, SegmSelect] = SelectSegmentation(obj)
            for c = 1:size(obj.AllFrames, 1)
                if strcmp(obj.DDMInfo.Segmentation, 'off')
                    [Mask] = ones(size(obj.AllFrames{1, 1}));
                    FullMask = Mask;
                    SegmentedMask = NaN(size(FullMask));
                    SegmSelect = [1, 0];
                    mkdir(append(obj.DDMInfo.file.path, filesep, 'Results_NoSegmentation'))
                elseif strcmp(obj.DDMInfo.Segmentation, 'on')
                    if strcmp(obj.DDMInfo.RerunSegmentation, 'off')
                        if exist(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'))
                            Rerun = 0;
                        elseif ~exist(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'))
                            Rerun = 1;
                        end
                    elseif strcmp(obj.DDMInfo.RerunSegmentation, 'on')
                        Rerun = 1;
                    end
                    if Rerun == 1
                        if strcmp(obj.DDMInfo.Ilastik, 'on')
                            if contains(obj.raw.movInfo.Path, 'Actin')
                                [Mask2] = A_SarahV.IlastikFirstMask_Actin(obj.raw.movInfo.Path, obj.DDMInfo.setMapping, obj.DDMInfo.ManualRegionRemoval);
                            else
                                [Mask2] = A_SarahV.IlastikFirstMask(obj.raw.movInfo.Path, obj.DDMInfo.setMapping, obj.DDMInfo.ManualRegionRemoval);
                            end
                        else
                            [Mask2] = obj.FrameSegmentation;
                        end
                        SegmentedMask = Mask2;
                        FullMask = NaN(size(SegmentedMask));
                        SegmSelect = [0, 1];
                        mkdir(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation'));
                        save(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c)),'Mask2');
                    elseif Rerun == 0
                        Mask = load(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'));
                        try
                            Mask = Mask.Mask;
                        catch
                            Mask = Mask.Mask2;
                        end
                        SegmentedMask = Mask;
                        FullMask = NaN(size(SegmentedMask));
                        SegmSelect = [0, 1];
                    end
                elseif strcmp(obj.DDMInfo.Segmentation, 'both')
                    [Mask1] = ones(obj.DDMInfo.Dimensions, obj.DDMInfo.Dimensions);
                    FullMask = Mask1;
                    mkdir(append(obj.DDMInfo.file.path, filesep, 'Results_NoSegmentation'))
                    if strcmp(obj.DDMInfo.RerunSegmentation, 'off')
                        if exist(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'))
                            Rerun = 0;
                        elseif ~exist(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'))
                            Rerun = 1;
                        end
                    elseif strcmp(obj.DDMInfo.RerunSegmentation, 'on')
                        Rerun = 1;
                    end
                    if Rerun == 1
                        if strcmp(obj.DDMInfo.Ilastik, 'on')
                            if contains(obj.raw.movInfo.Path, 'Actin')
                                [Mask2] = A_SarahV.IlastikFirstMask_Actin(obj.raw.movInfo.Path, obj.DDMInfo.setMapping, obj.DDMInfo.ManualRegionRemoval);
                            else
                                [Mask2] = A_SarahV.IlastikFirstMask(obj.raw.movInfo.Path, obj.DDMInfo.setMapping, obj.DDMInfo.ManualRegionRemoval);
                            end
                        else
                            [Mask2] = obj.FrameSegmentation;
                        end
                        SegmentedMask = Mask2;
                        SegmSelect = [1, 1];
                        mkdir(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation'));
                        save(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c)),'Mask2');
                    elseif Rerun == 0
                        Mask2 = load(append(obj.DDMInfo.file.path, filesep, 'Results_Segmentation', filesep, 'SegmentMask', num2str(c), '.mat'));
                        try
                            Mask2 = Mask2.Mask2;
                        catch
                            Mask2 = Mask2.Mask;
                        end
                        SegmentedMask = Mask2;
                        SegmSelect = [1, 1];
                    end
                end
                fullmask{c} = FullMask; 
                segmentedmask{c} = SegmentedMask;
                segmselect{c} = SegmSelect;
            end
        end

        function [bw] = FrameSegmentation(obj)  
            obj.AllFrames = cell(1,obj.DDMInfo.nFrames);
            f = waitbar(0,'Loading frames');
            
            for i = 1:obj.DDMInfo.nFrames 
                try 
                    waitbar(i./obj.DDMInfo.nFrames, f, 'Loading frames');
                    user = memory;
                    if user.MemAvailableAllArrays > 1e+09
                        try
                            ZeroFrame = single(obj.getFrame(i));
                        catch
                            ZeroFrame = single(obj.getFrame(i).Cam1);
                        end
                        close all
        
                        for k = 1:length(size(ZeroFrame))
                            if mod(size(ZeroFrame,1), 2) == 1 
                                ZeroFrame = ZeroFrame(2:end,:,:);
                            end
                            ZeroFrame = shiftdim(ZeroFrame,1);        
                        end
        
                        obj.AllFrames{i} = ZeroFrame;  
                    else
                        obj.DDMInfo.nFrames = i-1;
                        warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) ' frames. Loaded only ' num2str(i-1)] )
                        break;
                    end
                catch
                end
            end 
            close(f)
            obj.AllFrames = obj.AllFrames(~cellfun('isempty',obj.AllFrames));
        
            for i = 1:size(obj.AllFrames, 2)
                if i == 1
                    SumImage = double(obj.AllFrames{1,i});
                else
                    SumImage = SumImage + obj.AllFrames{1, i};  
                end
            end
            MeanImage = SumImage ./ size(obj.AllFrames, 2);
        
            ZeroFrame = SumImage;
            if contains(obj.raw.fullPath, 'Vinculin')
                minpix = round(5./((obj.DDMInfo.PixelSize).^2));
                maxpix = round(10./((obj.DDMInfo.PixelSize).^2));
            elseif contains(obj.raw.fullPath, 'Cadherin')
                minpix = round(5./((obj.DDMInfo.PixelSize).^2));
                maxpix = round(1000./((obj.DDMInfo.PixelSize).^2));
            end
            
            % Get threshold interactively
            if strcmp(obj.DDMInfo.ManualThresholdSelection, 'on')
                threshold = interactive_threshold_selection();
            else
                threshold = obj.DDMInfo.Threshold;
            end

            % Final processing
            [ShowImage, bw] = process_image(threshold);
            if strcmp(obj.DDMInfo.ManualRegionRemoval, 'on')
               bw = interactive_region_removal(ZeroFrame, bw);
            end
            Mask = bw; 

            f = figure();
            subplot(1, 2, 1);
            imagesc(ZeroFrame);
            title('Raw data (SumImage)');
            subplot(1, 2, 2);
            ShowImage(bw == 0) = 0;
            imagesc(ShowImage);
            title(['Thr = ', num2str(threshold)]);

            ImagePath = append(obj.DDMInfo.file.path, filesep, 'SegmentedImage.png');
            saveas(f,ImagePath)
        
            function threshold = interactive_threshold_selection()
                threshold = 0.10; 
                [ShowImage, ~] = process_image(threshold);

                fig = figure('Name', 'Threshold Selection', 'NumberTitle', 'off');
                subplot(1, 2, 1);
                imagesc(ZeroFrame);
                title('Raw Data');
                axis image;
                subplot(1, 2, 2);
                hImg = imagesc(ShowImage);
                title(['Threshold = ', num2str(threshold)]);
                axis image;
        
                uicontrol('Style', 'text', 'String', 'Threshold', 'Position', [50, 10, 80, 20]);
                hSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, 'Value', threshold, ...
                                    'Position', [140, 10, 300, 20], 'Callback', @update_image,...
                                    'SliderStep', [0.005, 0.1]);
                uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
                          'Position', [460, 10, 80, 20], 'Callback', @confirm_threshold);
                uiwait(fig);
        
                function update_image(hObject, ~)
                    threshold = get(hObject, 'Value');
                    [ShowImage, ~] = process_image(threshold);
                    hImg.CData = ShowImage;
                    title(['Threshold = ', num2str(threshold)]);
                end
        
                function confirm_threshold(~, ~)
                    uiresume(gcbf);
                    close(gcbf);
                end
            end
        
            function [ShowImage, bw] = process_image(threshold)
                ImgFiltered = imgaussfilt(ZeroFrame, 2);
                se = strel('disk', 7);
                bg = imopen(ImgFiltered, se);
                ImgCorrected = ImgFiltered - bg;
                ImgCorrected = ImgCorrected ./ max(ImgCorrected);
                ImgCorrected = imadjust(ImgCorrected);
        
                bw = imbinarize(ImgCorrected, threshold);
                % bw = imopen(bw, strel('disk', 3));
                % bw = imclose(bw, strel('disk', 5));
                % bw = bwareafilt(bw, [minpix maxpix]);
                % bw = imfill(bw, 'holes');
        
                ShowImage = ZeroFrame;
                ShowImage(bw == 0) = 0;
            end

            function bw = interactive_region_removal(ZeroFrame, bw)
                fig = figure('Name', 'Region Removal', 'NumberTitle', 'off');
                subplot(1,2,1);
                imagesc(ZeroFrame); title('Raw Data'); axis image;
                subplot(1,2,2);
                ShowImage = ZeroFrame;
                ShowImage(bw == 0) = 0;
                hMask = imagesc(ShowImage); title('Select Regions to Remove');
                axis image;
                
                uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
                          'Position', [250, 10, 100, 30], 'Callback', @confirm_removal);
                hold on;
                while true
                    h = drawfreehand(gca, 'Color', 'r');
                    if isempty(h) || ~isvalid(h)
                        break;
                    end
                    mask = createMask(h);
                    bw(mask) = 0;
                    ShowImage = ZeroFrame;
                    ShowImage(bw == 0) = 0;
                    if isvalid(hMask) % ✅ Check if hMask is still valid
                        set(hMask, 'CData', ShowImage);
                    end
                end
                hold off;
        
                function confirm_removal(~, ~)
                    uiresume(gcbf);
                    close(gcbf);
                end
            end
        end
        
      

        
        function  [Qmin, Qmax] = LoadAllFramesWithROI(obj,driftCorr, ROI)
             
            obj.AllFrames = cell(1,obj.DDMInfo.nFrames);
            f = waitbar(0,'Loading frames');
            for i=1:obj.DDMInfo.nFrames 
                try 
                    waitbar(i./obj.DDMInfo.nFrames,f,'Loading frames');
                    user = memory;
                    if  user.MemAvailableAllArrays>1e+09
                        try
                        ZeroFrame= single(obj.getFrame(i));
                        catch
                        ZeroFrame= single(obj.getFrame(i).Cam1);
                        end
                        close all
    
                        
                        for k=1:length(size(ZeroFrame))
                           if mod(size(ZeroFrame,1),2) == 1 
                               ZeroFrame = ZeroFrame(2:end,:,:);
                               
                           end
                           ZeroFrame = shiftdim(ZeroFrame,1);        
                        end
                                    
                        obj.AllFrames{i}= ZeroFrame;  
                        
                    else
                        obj.DDMInfo.nFrames = i-1;
                        warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) 'frames. Loaded only ' num2str(i-1)] )
                        break;
                    end
                catch
                end
            end 
            close(f);
            obj.AllFrames = obj.AllFrames(~cellfun('isempty',obj.AllFrames));
            obj.DDMInfo.nFrames = size(obj.AllFrames, 2);
            correlationInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
            %correlation function
            correlationInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
            %for driftCalculation ==> 1 mean that drift is calculated for each frame
            scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
            %objects of interest

            if driftCorr
                [corrData,~] = PreProcess.CorrelationDrift(obj.AllFrames,scalingFactor,correlationInfo);
            else
                corrData = obj.AllFrames;
            end
            
            if isempty(ROI)
               ROI = [1 1 size(corrData{1},2),size(corrData{1},1)];
            end
            
            if ~iscell(ROI)
                ROIcell = {ROI};
                MaxROI = 1;
            else
                MaxROI = size(ROI, 1);
                ROIcell = ROI;
            end
            
            obj.AllFrames = {};
            for c = 1:MaxROI
                ROI = ROIcell{c,1};
                %Fix ROI
                if mod(ROI(3),2) ~= 0
                    ROI(3) = ROI(3)-1;
                end
                
                if mod(ROI(4),2) ~= 0
                    ROI(4) = ROI(4)-1;
                end
                
                dim = size(corrData);
                %apply ROI
                if isempty(obj.DDMInfo.Qmin)
                    Qmin(c) = (2*pi)./((min(ROI(1,3:4)))*obj.DDMInfo.PixelSize);
                else
                    if size(obj.DDMInfo.Qmin, 2) > 1
                        Qmin(c) = obj.DDMInfo.Qmin(c);
                    else
                        Qmin(c) = obj.DDMInfo.Qmin;
                    end
                end
    
                if isempty(obj.DDMInfo.Qmax)
                    DiffLimit = obj.DDMInfo.Wavelength./(2*obj.DDMInfo.NA);
                    if obj.DDMInfo.ParticleSize < DiffLimit
                        Qmax(c) = (2*pi)./DiffLimit;
                    else
                        Qmax(c) = (2*pi)./obj.DDMInfo.ParticleSize;
                    end
                else
                    if size(obj.DDMInfo.Qmax, 2) > 1
                        Qmax(c) = obj.DDMInfo.Qmax(c);
                    else
                        Qmax(c) = obj.DDMInfo.Qmax;
                    end
                end
                
                for i = 1:dim(end)
                   corrDataC{c,i} = corrData{i}(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1) + ROI(3),:); 
                end
                obj.AllFrames = corrDataC;
                obj.DDMInfo.Qmax = min(Qmax);
                obj.DDMInfo.Qmin = max(Qmin);
            end
            
        end


        function  [Qmin, Qmax] = LoadAllFramesWithMask(obj,driftCorr, ROI, Mask, s)
                   
            obj.AllFrames = cell(1,obj.DDMInfo.nFrames);
            f = waitbar(0,'Loading frames');
            for i=1:obj.DDMInfo.nFrames 
                % try 
                    waitbar(i./obj.DDMInfo.nFrames,f,'Loading frames');
                    user = memory;
                    if  user.MemAvailableAllArrays>1e+09
                        try
                        ZeroFrame= single(obj.getFrame(i));
                        catch
                        ZeroFrame= single(obj.getFrame(i).Cam1);
                        end
                        close all
    
                        
                        for k=1:length(size(ZeroFrame))
                           if mod(size(ZeroFrame,1),2) == 1 
                               ZeroFrame = ZeroFrame(2:end,:,:);
                               
                           end
                           ZeroFrame = shiftdim(ZeroFrame,1);        
                        end
                                    
                        obj.AllFrames{i}= ZeroFrame;  
                        
                    else
                        obj.DDMInfo.nFrames = i-1;
                        warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) 'frames. Loaded only ' num2str(i-1)] )
                        break;
                    end
                % catch
                % end
            end 
            close(f);
            obj.AllFrames = obj.AllFrames(~cellfun('isempty',obj.AllFrames));
            if contains(obj.raw.movInfo.Path, 'Vinculin')
                se = strel('disk', 5);
            else
                se = strel('disk', 15);
            end
            if isstruct(Mask)
                MaxLoop = 3;
            else
                MaxLoop = 1;
            end
            FramesStore = obj.AllFrames;
            obj.AllFrames = {};
            for c = 1:MaxLoop
                if isstruct(Mask)
                    if c == 1
                        MaskM = Mask.Mask_Ruff;
                    elseif c == 2
                        MaskM = Mask.Mask_Trans;
                    elseif c == 3
                        MaskM = Mask.Mask_Cort;
                    end
                else
                    MaskM = Mask;
                end
                for i = 1:size(FramesStore, 2)
                    Frame = FramesStore{1,i};
                    if strcmp(obj.DDMInfo.GlobalBgSubstraction, 'on')
                        bg = imerode(Frame, strel('disk', 3));
                        NewFrame = Frame - bg;
                    else
                        NewFrame = Frame;
                    end

                    if s == 2
                        NewFrame = NewFrame .* MaskM;%+ RefFrame .* (1 - MaskM);
                    end
                    obj.AllFrames{c,i} = NewFrame;
                end
            end
            obj.DDMInfo.nFrames = size(obj.AllFrames, 2);
            correlationInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
            %correlation function
            correlationInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
            %for driftCalculation ==> 1 mean that drift is calculated for each frame
            scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
            %objects of interest

            if driftCorr
                [corrData,~] = PreProcess.CorrelationDrift(obj.AllFrames,scalingFactor,correlationInfo);
            else
                corrData = obj.AllFrames;
            end
            
            if isempty(ROI)
               ROI = [1 1 size(corrData{1},2),size(corrData{1},1)];
            end

            if ~iscell(ROI)
                if size(obj.AllFrames, 1) ~= 1
                    MaxROI = size(obj.AllFrames, 1);
                    for c = 1:MaxROI
                        ROIcell{c,1} = ROI;
                    end
                else
                    ROIcell = {ROI};
                    MaxROI = 1;
                end
            else
                MaxROI = size(ROI, 1);
                ROIcell = ROI;
            end

            obj.AllFrames = {};
            for c = 1:MaxROI
                ROI = ROIcell{c,1};
            
                %Fix ROI
                if mod(ROI(3),2) ~=0
                    ROI(3) = ROI(3)-1;
                end
                
                if mod(ROI(4),2) ~=0
                    ROI(4) = ROI(4)-1;
                end
                
                dim = size(corrData);
                %apply ROI
                if isempty(obj.DDMInfo.Qmin)
                    Qmin(c) = (2*pi)./((min(ROI(1,3:4)))*obj.DDMInfo.PixelSize);
                else
                    if size(obj.DDMInfo.Qmin, 2) > 1
                        Qmin(c) = obj.DDMInfo.Qmin(c);
                    else
                        Qmin(c) = obj.DDMInfo.Qmin;
                    end
                end
    
                if isempty(obj.DDMInfo.Qmax)
                    DiffLimit = obj.DDMInfo.Wavelength./(2*obj.DDMInfo.NA);
                    if obj.DDMInfo.ParticleSize < DiffLimit
                        Qmax(c) = (2*pi)./DiffLimit;
                    else
                        Qmax(c) = (2*pi)./obj.DDMInfo.ParticleSize;
                    end
                else
                    if size(obj.DDMInfo.Qmax, 2) > 1
                        Qmax(c) = obj.DDMInfo.Qmax(c);
                    else
                        Qmax(c) = obj.DDMInfo.Qmax;
                    end
                end
                
                for i = 1:dim(end)
                   corrDataC{c,i} = corrData{1,i}(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1) + ROI(3),:); 
                end
                obj.AllFrames = corrDataC;
                obj.DDMInfo.Qmax = min(Qmax);
                obj.DDMInfo.Qmin = max(Qmin);
            end
            
        end
        
        function  [Qmin, Qmax] = LoadAllFrames(obj,driftCorr)
        % Loads all data. Stops if it almost runs out of memory ,- when the
        % amount of memory left is less then 2 GB, in order to leave some
        % of it for other operations. 
            
            h = questdlg('Do you want to use a ROI?','Question to user','Yes','No', 'Yes');
            frame = obj.getFrame(1);
            if strcmp(h,'Yes')
                figure

                imagesc(frame(:,:,1))
                test = drawrectangle();
                ROI  = round(test.Position);    
            else
                ROI = [];

            end
                   
            obj.AllFrames = cell(1,obj.DDMInfo.nFrames);
            f = waitbar(0,'Loading frames');
            for i=1:obj.DDMInfo.nFrames 
                waitbar(i./obj.DDMInfo.nFrames,f,'Loading frames');
                user = memory;
                if  user.MemAvailableAllArrays>1e+09
                    try
                    ZeroFrame= single(obj.getFrame(i));
                    catch
                    ZeroFrame= single(obj.getFrame(i).Cam1);
                    end
                    close all

                    
                    for k=1:length(size(ZeroFrame))
                       if mod(size(ZeroFrame,1),2) == 1 
                           ZeroFrame = ZeroFrame(2:end,:,:);
                           
                       end
                       ZeroFrame = shiftdim(ZeroFrame,1);        
                    end
                    
                    obj.AllFrames{i}= ZeroFrame;  
                    
                else
                    obj.DDMInfo.nFrames = i-1;
                    warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) 'frames. Loaded only ' num2str(i-1)] )
                    break;
                end
            end 
            close(f);
            correlationInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
            %correlation function
            correlationInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
            %for driftCalculation ==> 1 mean that drift is calculated for each frame
            scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
            %objects of interest
            

            if driftCorr
                [corrData,~] = PreProcess.CorrelationDrift(obj.AllFrames,scalingFactor,correlationInfo);
            else
                corrData = obj.AllFrames;
            end
            
            if isempty(ROI)
               ROI = [1 1 size(corrData{1},2)-2,size(corrData{1},1)-2];
            end
            
            %Fix ROI
            if mod(ROI(3),2) ~=0
                ROI(3) = ROI(3)-1;
            end

            if mod(ROI(4),2) ~=0
                ROI(4) = ROI(4)-1;
            end
            
            dim = size(corrData);
            %apply ROI
            if isempty(obj.DDMInfo.Qmin)
                Qmin = (2*pi)./((min(ROI(1,3:4)))*obj.DDMInfo.PixelSize);
                obj.DDMInfo.Qmin = Qmin;
            else
                Qmin = obj.DDMInfo.Qmin;
            end

            if isempty(obj.DDMInfo.Qmax)
                DiffLimit = obj.DDMInfo.Wavelength./(2*obj.DDMInfo.NA);
                if obj.DDMInfo.ParticleSize < DiffLimit
                    Qmax = (2*pi)./DiffLimit;
                    obj.DDMInfo.Qmax = Qmax;
                else
                    Qmax = (2*pi)./obj.DDMInfo.ParticleSize;
                    obj.DDMInfo.Qmax = Qmax;
                end
            else
                Qmax = obj.DDMInfo.Qmax;
            end
            
            for i = 1:dim(end)
               corrData{i} = corrData{i}(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1) + ROI(3),:); 
            end
            obj.AllFrames = corrData;
            
        end

        function  LoadAllFrames2(obj,driftCorr)
        % Loads all data. Stops if it almost runs out of memory ,- when the
        % amount of memory left is less then 2 GB, in order to leave some
        % of it for other operations.
                   
            if isempty(obj.DDMInfo.TotTime)
                obj.DDMInfo.nFrames = obj.raw.maxFrame;
            end

            obj.AllFrames = cell(1,obj.DDMInfo.nFrames{1});
            f = waitbar(0,'Loading frames');
            for i=1:obj.DDMInfo.nFrames{1} 
                waitbar(i./obj.DDMInfo.nFrames{1},f,'Loading frames');
                user = memory;
                if  user.MemAvailableAllArrays>1e+09
                    try
                    ZeroFrame= single(obj.getFrame(i));
                    catch
                    ZeroFrame= single(obj.getFrame(i).Cam1);
                    end
                    close all

                    
                    for k=1:length(size(ZeroFrame))
                       if mod(size(ZeroFrame,1),2) == 1 
                           ZeroFrame = ZeroFrame(2:end,:,:);
                           
                       end
                       ZeroFrame = shiftdim(ZeroFrame,1);        
                    end
                    
                    obj.AllFrames{i}= ZeroFrame;  
                    
                else
                    obj.DDMInfo.nFrames = i-1;
                    warning(['Not enough RAM to load all of the '  num2str(obj.DDM.nFrames) 'frames. Loaded only ' num2str(i-1)] )
                    break;
                end
            end 
            close(f);
            DDMInfo.corrSz = 100; %in px. Radius of the ROI used for correlation
            %correlation function
            DDMInfo.driftPeriod = 1; %in Frame, Number of frame that are averaged
            %for driftCalculation ==> 1 mean that drift is calculated for each frame
            scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
            %objects of interest
            

            if driftCorr
                [corrData,~] = PreProcess.CorrelationDrift(obj.AllFrames,scalingFactor,DDMInfo);
            else
                corrData = obj.AllFrames;
            end

            ROIFile = append(obj.raw.movInfo{1}.Path, filesep, 'ROI.mat');
            if exist(ROIFile)
                ROI = load(ROIFile);
                ROI = ROI.ROI;
            else
                % h = questdlg('Do you want to use a ROI?','Question to user','Yes','No', 'Yes');
                frame = obj.getFrame(1);
                if strcmp(obj.DDMInfo.useROI  ,'on')
                    figure
    
                    imagesc(frame(:,:,1))
                    test = drawrectangle();
                    ROI  = round(test.Position);
    
                else
                    ROI = [1,1,size(frame, 1)-1, size(frame,2)-1];  
                end
    
                Filename = append(obj.raw.movInfo{1}.Path, filesep, 'ROI.mat');
                save(Filename, 'ROI');
            end
            
            % ROI = [1 1 size(corrData{1},2)-2,size(corrData{1},1)-2];
            
            %Fix ROI
            if mod(ROI(3),2) ~=0
                ROI(3) = ROI(3)-1;
            end

            if mod(ROI(4),2) ~=0
                ROI(4) = ROI(4)-1;
            end
            
            dim = size(corrData);
            %apply ROI            
            for i = 1:dim(end)
               corrData{i} = corrData{i}(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1) + ROI(3),:); 
            end
            obj.AllFrames = corrData;
            
        end
        
        
        function [AvgFFT, SingleFrameFFT] =  CalculateDelta(obj,AvgFFT,FrameSize,dt,ROI,c)
            nFrames = obj.DDMInfo.nFrames{1} - dt;
            SingleFrameFFT = AvgFFT;
            DiffFrames = cell(2, nFrames);
            counts = 0;

            xRange = ROI(1,1):ROI(1,2);
            yRange = ROI(2,1):ROI(2,2);
            zRange = ROI(3,1):ROI(3,2);

            useGPU = (obj.IsCudaDevice == 1);

            for t=1:nFrames

                Frame1 = obj.AllFrames{c, t}(xRange, yRange, zRange);
                Frame2 = obj.AllFrames{c, t + dt}(xRange, yRange, zRange);

                if useGPU == 1
                    FrameDelta = gpuArray(Frame2 - Frame1);
                else
                    FrameDelta = Frame2 - Frame1;
                end

                FrameDeltaMean = mean(FrameDelta, 'all');
                FrameDeltaZeroMean = FrameDelta - FrameDeltaMean;

                FrameDeltaFFT = abs(fftshift(fftn(FrameDeltaZeroMean, [FrameSize(1), FrameSize(2), FrameSize(3)]))).^2;
                AvgFFT = AvgFFT + FrameDeltaFFT;
                counts = counts + 1;
        
                SingleFrame = Frame1;
                SingleFrameZeroMean = SingleFrame - mean(SingleFrame, 'all');
                SingleFrameFFT = SingleFrameFFT + abs(fftshift(fftn(SingleFrameZeroMean))).^2;
            end
            AvgFFT = AvgFFT./counts;
            SingleFrameFFT = SingleFrameFFT./counts;
        end 
        
        function IntroduceShift(obj, shift)
            TakeFrame = obj.AllFrames;
            NumberOfFrames = size(obj.AllFrames,2);
            obj.AllFrames = struct([]);
            for i = 1:round(NumberOfFrames/shift)
                try
                    Cutted = TakeFrame{1,i}(:, ((end-round((size(TakeFrame{1,i}, 2)/2)))-i*shift: end-i*shift));
                    obj.AllFrames{1,i} = Cutted;
                catch
                end
            end
            obj.DDMInfo.nFrames = size(obj.AllFrames, 2);
        end         
        
        function  DDMOutput = main(obj,Info,file,varargin)
            for c = 1:size(obj.AllFrames, 1)
                DDMOutputfile = append(file.path, filesep, 'DDMOutput' , num2str(c), '.mat');
                if ~exist(DDMOutputfile)
                    obj.info.runMethod  = 'run';
                end

            
                if strcmp(obj.info.runMethod, 'run')
                    %Parse inputs conditionally
                    p = inputParser;  
                    addOptional(p, 'ROI',  [1, size(obj.AllFrames{c,1},1), size(obj.AllFrames{c,1},1);  %Default ROI as image size
                                            1, size(obj.AllFrames{c,1},2), size(obj.AllFrames{c,1},2);
                                            1 ,size(obj.AllFrames{c,1},3) ,size(obj.AllFrames{c,1},3)]);
                    addOptional(p, 'Padsize',zeros(1, 3));
                    addOptional(p, 'NumBins', 200);
                    addOptional(p, 'CriticalAngle',0);
                    parse(p,varargin{:});
                    FrameSize = p.Results.ROI(:,3)'+p.Results.Padsize;
                    ExpTime = Info.ExpTime; 
                    DDMOutput = [];
                    
                    
                    %Generate struct grid for n-D averaging
                    f = waitbar(0,'Calculating AvgFFT frame by frame');
                    if isempty(obj.DDMInfo.MinTimeLag)
                        dtMin = 1;
                    else
                        dtMin = obj.DDMInfo.MinTimeLag;
                    end
                    
                    DDMOutput = [];
                    AnisotropyOutput = [];
                    Asymptote = [];
                    Baseline = [];
                    for dt=dtMin:obj.DDMInfo.nFrames{c}-dtMin
                        obj.IsCudaDevice = 0;
                        if obj.IsCudaDevice==1
                            AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3),'gpuArray');
                        else
                            AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3));
                        end
                            
                        [AvgFFT, SingleFrameFFT] = obj.CalculateDelta(AvgFFT,FrameSize,dt, p.Results.ROI, c); 
                        
                        [RadiallyAveragedDDMSignal, AnisotropyValues, A, B]=  obj.AverageRadialy3D(AvgFFT,FrameSize, p.Results.CriticalAngle, SingleFrameFFT);
                        % QImages{dt, :} = reconstructions;
                        DDMOutput(:,1) =[NaN ; RadiallyAveragedDDMSignal(:,1)];
                        DDMOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)];
    
                        AnisotropyOutput(:,1) = [NaN ; RadiallyAveragedDDMSignal(:,1)];
                        AnisotropyOutput(:,end+1) = [dt ; RadiallyAveragedDDMSignal(:,2)]; 
                        waitbar(dt./(obj.DDMInfo.nFrames{c}-dtMin),f,'Calculating AvgFFT frame by frame');

                        Asymptote = [Asymptote,A];
                        Baseline = [Baseline, B];
                    end
                    close(f)
                    close all
                    DDMOutput =  obj.ConvertOutput(DDMOutput, ExpTime);
                    AnisotropyOutput = obj.ConvertOutput(AnisotropyOutput, ExpTime);                  
                    obj.DDMOutput{c,1} = DDMOutput;
        
                    filename = append(file.path, filesep, 'DDMOutput', num2str(c), '.mat');
                    save(filename, 'DDMOutput')
                    filename = append(file.path, filesep, 'AnisotropyOutput', num2str(c), '.mat');
                    save(filename, 'AnisotropyOutput')
                else
                    disp('Found DDMOuput file - loading it');
                    DDMOutput = load(DDMOutputfile);
                    DDMOutput = DDMOutput.DDMOutput;
                    obj.DDMOutput{c,1} = DDMOutput;
                    disp('Found DDMOuput file - Done');
                end
    
                %%% Calculate average DDMSignal over interested QRange
                Start = find(DDMOutput.QVector > obj.DDMInfo.Qmin, 1);
                End = find(DDMOutput.QVector < obj.DDMInfo.Qmax & DDMOutput.QVector ~= 0, 1, 'last');
                for i = Start:1:End
                    DDMSignal(i,:) = DDMOutput.DDMSignalValue{i,1};  
                end
                DDMSignalAv = nanmean(DDMSignal,1);
                DDMSignalvar = nanvar(DDMSignal,1);
                AverageDDMSignal = table(DDMSignalAv.', DDMSignalvar.','VariableNames',["Av Signal","Var signal"]);
    
                filename = append(file.path, 'AverageDDMSignal', num2str(c), '.mat');
                save(filename, 'AverageDDMSignal')
            end

        end

        function  DDMOutput = FFTReconstructions(obj,Info,file, MaskStruct, varargin)
            %Parse inputs conditionally
            if contains(obj.raw.movInfo.Path, 'Actin')
                maxLoop  = 3;
            else
                maxLoop = 1;
            end
            for c = 1:maxLoop
                if contains(obj.raw.movInfo.Path, 'Actin')
                    if c == 1
                        Mask = MaskStruct.Mask_Ruff;
                    elseif c == 2
                        Mask = MaskStruct.Mask_Trans;
                    elseif c == 3
                        Mask = MaskStruct.Mask_Cort;
                    end
                else
                    Mask = MaskStruct;
                end
                p = inputParser;  
                addOptional(p, 'ROI',  [1, size(obj.AllFrames{1},1), size(obj.AllFrames{1},1);  %Default ROI as image size
                                        1, size(obj.AllFrames{1},2), size(obj.AllFrames{1},2);
                                        1 ,size(obj.AllFrames{1},3) ,size(obj.AllFrames{1},3)]);
                addOptional(p, 'Padsize',zeros(1, 3));
                addOptional(p, 'NumBins', 200);
                addOptional(p, 'CriticalAngle',0);
                parse(p,varargin{:});
                FrameSize = p.Results.ROI(:,3)'+p.Results.Padsize;
                ExpTime = Info.ExpTime; 
                DDMOutput = [];
                
                
                %Generate struct grid for n-D averaging
                f = waitbar(0,'Calculating AvgFFT frame by frame');
                if isempty(obj.DDMInfo.MinTimeLag)
                    dtMin = 1;
                else
                    dtMin = obj.DDMInfo.MinTimeLag;
                end
    
                [RadialValueInQSpace, ValidRange]=obj.Get3DGrid(FrameSize,p.Results.CriticalAngle);
                SimilarityMatrix = {};
                for dt=dtMin:5:dtMin+50
                    waitbar(dt./(obj.DDMInfo.nFrames-1), f, 'Calculating Diff Images & FFT')
                    obj.IsCudaDevice = 0;
                    if obj.IsCudaDevice==1
                        AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3),'gpuArray');
                    else
                        AvgFFT = zeros(FrameSize(1),FrameSize(2),FrameSize(3));
                    end
                    NumBins = round((obj.DDMInfo.Qmax - obj.DDMInfo.Qmin)./(4*pi./(min(size(RadialValueInQSpace))*obj.DDMInfo.PixelSize))-1);
                    for t=obj.DDMInfo.nFrames-dt-50:obj.DDMInfo.nFrames-dt
                        ROI = p.Results.ROI;
                        if obj.IsCudaDevice==1
                            FrameDelta =gpuArray(obj.AllFrames{t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))); 
                        else
                            FrameDelta =(obj.AllFrames{t+dt}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2))-obj.AllFrames{t}(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),ROI(3,1):ROI(3,2)));  
                        end
                        DiffFrame = FrameDelta;
                        AvgFFT = fftshift(fft2(FrameDelta)); 
                        Mask = Mask(1:size(DiffFrame, 1), 1:size(DiffFrame,2));
                        DiffFrame(Mask ~= 1) = 0;
    
                        MinQTracked = min(RadialValueInQSpace(:));
                        Min = max(MinQTracked, obj.DDMInfo.Qmin);
                        BinSize = (obj.DDMInfo.Qmax - Min)./NumBins;
                        FoundRadii = 0;
                        R = Min;
                        next = 1;
        
                        %%% Calculate Average
                        Reconstructions = {};
                        SimilarityIdx = [];
                        while ~isempty(FoundRadii)
                            R = R+BinSize;
                            if R <= obj.DDMInfo.Qmax
                                if size(AvgFFT) == size(RadialValueInQSpace)    
                                    QMask = RadialValueInQSpace;
                                    QMask(QMask>=R-BinSize & RadialValueInQSpace<R) = 1;
                                    QMask(QMask ~= 1) = 0;
                                    FourierIm = AvgFFT;
                                    FourierIm(QMask ~= 1) = 0;
                                    ReconstructedIm = ifft2(ifftshift(FourierIm));
                                    Reconstructions{end+1,1} = ReconstructedIm;
                                    Reconstructions{end, 2} = R;
                                    ReconstructedIm(Mask ~= 1) = 0;
                                    SimilarityIdx(end+1) = ssim(ReconstructedIm, DiffFrame);
                                elseif size(AvgFFT)
                                    error('Differential image and reference frame q vector do not have same dimensions')
                                end
                            else
                                break
                            end
                        end
                        if t == obj.DDMInfo.nFrames-dt-50
                            Similarityt = SimilarityIdx;
                        else
                            Similarityt(1,:) = Similarityt(1,:) + SimilarityIdx;
                        end
                        Similarityt(2,:) = cell2mat(Reconstructions(:,2))';
                    end
                    SimilarityMatrix{end+1,1} = dt;
                    SimilarityMatrix{end, 2} = Similarityt;
    
                end
                close(f)
                SimilarityMatrix = table(SimilarityMatrix(:,1), SimilarityMatrix(:,2), 'VariableNames', {'TimeLag', 'SimilarityIndex'});
                if contains(obj.raw.movInfo.Path, 'Actin')
                    if c == 1
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Ruff.mat');
                    elseif c == 2
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Trans.mat');
                    elseif c == 3
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Cort.mat');
                    end
                else
                    Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors.mat');
                end
                
                save(Filename, "SimilarityMatrix");
    
                cc = colormap(jet(size(SimilarityMatrix,1)));
                Fig = figure()
                for i = 1:size(SimilarityMatrix,1)
                    plot(SimilarityMatrix.SimilarityIndex{i,1}(2,:), SimilarityMatrix.SimilarityIndex{i,1}(1,:), 'Color', cc(i,:), 'LineWidth', 1') ;
                    hold on
                end
                colormap(cc)
                d = colorbar;
                caxis([dtMin dtMin+50]);
                title(d, 'TimeLag')
                ylabel('Similarity')
                xlabel('QVector')
                
                if contains(obj.raw.movInfo.Path, 'Actin')
                    if c == 1
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Ruff.png');
                    elseif c == 2
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Trans.png');
                    elseif c == 3
                        Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors_Cort.png');
                    end
                else
                    Filename = append(obj.raw.movInfo.Path, filesep, 'SimilarityIdxQVectors.png');
                end
                
                saveas(Fig,Filename);
            end
        end
        
        
        
        function [AverageDDMValueAtR, anisotropy_values, A, B] = AverageRadialy3D(obj, AvgFFT,FrameSize,critangle,SingleFrame)
            
        % Averages the scattering function in 3d radially. Returns the
        % average values in function of time and q-vector.

                Reconstructions = [];
                [RadialValueInQSpace, ValidRange]=obj.Get3DGrid(FrameSize,critangle);
                % BinSize = max(RadialValueInQSpace,[],'all')/NumBins;
                NumBins = round((obj.DDMInfo.Qmax - obj.DDMInfo.Qmin)./(4*pi./(min(size(RadialValueInQSpace))*obj.DDMInfo.PixelSize))-1);
                MinQTracked = min(RadialValueInQSpace(:));
                Min = max(MinQTracked, obj.DDMInfo.Qmin);
                BinSize = (obj.DDMInfo.Qmax - Min)./NumBins;
                FoundRadii = 0;
                R = Min;
                next = 1;
                AverageDDMValueAtR = nan(NumBins,2);
                AB = nan(NumBins,1);
                AvgFFT = gather(AvgFFT);

                [N, M] = size(AvgFFT);
                cx = floor(N/2) + 1;
                cy = floor(M/2) + 1;
                num_theta = 360;
                theta_vals = linspace(0, 2*pi, num_theta);

                %%% Calculate Average
                while ~isempty(FoundRadii)
                    R = R+BinSize;
                    if R <= obj.DDMInfo.Qmax
                        if size(AvgFFT) == size(RadialValueInQSpace)
                            for t = 1:num_theta
                                x = round(cx + R * cos(theta_vals(t)));
                                y = round(cy + R * sin(theta_vals(t)));
                                if x >= 1 && x <= N && y >= 1 && y <= M
                                    intensity_profile(t) = AvgFFT(x,y);
                                end
                            end

                            I_max = max(intensity_profile);
                            I_min = min(intensity_profile);
                            anisotropy_values(next,:) = [R, (I_max - I_min) / (I_max + I_min)];

                            FoundRadii = AvgFFT(RadialValueInQSpace>=R-BinSize & RadialValueInQSpace<R &  ValidRange );
                            AverageDDMValueAtR(next,:) = [R ,  nanmean(FoundRadii, 'all')];
                            next = next+1;

                            Mask = RadialValueInQSpace;
                            Mask(Mask>=R-BinSize & RadialValueInQSpace<R) = 1;
                            Mask(Mask ~= 1) = 0;

                            AB(next, 1) = 2*mean(SingleFrame(Mask == 1));

                        elseif size(AvgFFT)
                            error('Differential image and reference frame q vector do not have same dimensions')
                        end
                    else
                        break
                    end
                end
                B = AB(end);
                A = AB - B;
        end

        function [Quality] = CheckISFQuality(obj)
          
             for i = 1:size(obj.DDMOutput, 1)
                isf = obj.DDMOutput.DDMSignalValue{i,1};
                tau = obj.DDMOutput.Time{i,1};
            
                isf = isf(:);
                tau = tau(:);
            
                % Normalize ISF to range [0, 1] for analysis
                isf_min = min(isf);
                isf_max = max(isf);
                normalized_isf = (isf - isf_min) / (isf_max - isf_min);
            
                isf_derivative = diff(normalized_isf) ./ diff(tau); % Numerical derivative
                smoothness_measure(i,1) = mean(isf_derivative); % Average increase rate
            
                if smoothness_measure(i,1) > 0.0020
                    isf_quality(i,1) = 1; 
                else
                    isf_quality(i,1) = 0;
                end     
             end

             if any(isf_quality )
                 Quality = 'Good';
                 obj.DDMOutput.Quality = isf_quality;
                
                 diff_array = [0; diff(isf_quality)];
                 starts = find(diff_array == 1); % Start of a sequence of 1's
                 ends = find(diff_array == -1) - 1; % End of a sequence of 1's
                
                 if isf_quality(end) == 1
                     ends = [ends; length(isf_quality)];
                 end

                 if size(ends, 1) > size(starts, 2)
                     if isf_quality(1) == 1
                        starts = [1; starts];
                     end
                 end
                 lengths = ends - starts + 1;
                
                 [~, max_idx] = max(lengths);
                
                 if isempty(max_idx)
                     % No sequences of 1's found
                     start_idx = [];
                     end_idx = [];
                 else
                     start_idx = starts(max_idx);
                     end_idx = ends(max_idx);
                 end
                
                 if isempty(start_idx)
                     start_idx = 1;
                 end
                 if isempty(end_idx)
                     end_idx = size(obj.DDMOutput, 1);
                 end
                 obj.DDMInfo.Qmin = obj.DDMOutput.QVector(start_idx);
                 obj.DDMInfo.Qmax = obj.DDMOutput.QVector(end_idx-1);
             else
                 Quality = 'Bad';
             end
        end

     
end 
methods(Static)
    

        function DDMOutput = ConvertOutput(DDMOutList,ExpTime)
            DDMOutput = table(0,{[]},{[]},'VariableNames',{'QVector','Time','DDMSignalValue'});
            warning('off','all')
            for i=2:size(DDMOutList,1)
                DDMOutput.QVector(i-1) = DDMOutList(i,1);
                DDMOutput.Time(i-1) = {DDMOutList(1,2:end)*ExpTime};
                DDMOutput.DDMSignalValue(i-1) ={DDMOutList(i,2:end)};
    
            end
            warning('on','all')
        end
     



end  
end