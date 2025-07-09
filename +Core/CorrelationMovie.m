classdef CorrelationMovie < Core.MPMovie
properties (SetAccess = 'private')
    
    IsEnoughRam
    IsCudaDevice

end
properties (SetAccess ='public')
    CorrelationOutput
    AllFrames
    CorrelationInfo
end
    
    methods 
    
%%  Instantiation
        function obj = CorrelationMovie(raw, MPCal, info,CorrelationInfo)
             obj  = obj@Core.MPMovie(raw,MPCal,info); 

             assert(~isempty(CorrelationInfo), 'No MovieInfo provided, check parameter input'); 
             
             obj.CorrelationInfo = CorrelationInfo;

             try
                 gpuDevice;
                 obj.IsCudaDevice = true;         
             catch
                 obj.IsCudaDevice = false;
             end

        end
        
%%  Analysis

        function  LoadAllFrames(obj,driftCorr)
        % Loads all data. Stops if it almost runs out of memory ,- when the
        % amount of memory left is less then 2 GB, in order to leave some
        % of it for other operations.
                   
            if isempty(obj.CorrelationInfo.TotTime)
                obj.CorrelationInfo.nFrames = obj.raw.maxFrame;
            end

            obj.AllFrames = cell(1,obj.CorrelationInfo.nFrames); % TOOK {1} out (for TotTime not 0))
            f = waitbar(0,'Loading frames');
            for i=1:obj.CorrelationInfo.nFrames % TOOK {1} out
                waitbar(i./obj.CorrelationInfo.nFrames,f,'Loading frames'); % TOOK {1} out
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
                    obj.CorrelationInfo.nFrames = i-1;
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
    
                Filename = append(obj.raw.movInfo.Path, filesep, 'ROI.mat');
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
            
        
        function  [CorrelationOutput] = main(obj,Info,file, varargin)
                Corr1_cell = cell(obj.CorrelationInfo.nFrames-1, 1); % Each dt gets its own cell % TOOK {1} out

                for dt = 1:obj.CorrelationInfo.nFrames-1 % TOOK {1} out
                    stp = dt;
                    cnt = 1;
                    tempCorr = [];
                
                    while cnt <= stp && cnt + stp <= obj.CorrelationInfo.nFrames % TOOK {1} out
                        idx = cnt:stp:obj.CorrelationInfo.nFrames; % TOOK {1} out
                        for c = 1:numel(idx) - 1
                            Frame1 = obj.AllFrames{1,idx(c)};
                            Frame2 = obj.AllFrames{1,idx(c+1)};
                            if strcmp(obj.raw.movToLoad, 'Mask')
                                Frame1(Frame1 == 0) = 1;
                                Frame1(Frame1 ~= 1) = 0;
                                Frame2(Frame2 == 0) = 1;
                                Frame2(Frame2 ~= 1) = 0;
                            end
                            tempCorr(end+1,1) = multissim(Frame1, Frame2);
                        end
                        cnt = cnt + 1;
                    end
                
                    Corr1_cell{dt} = tempCorr; % Safe!
                    fprintf('Processing dt = %d / %d\n', dt, obj.CorrelationInfo.nFrames - 1); % TOOK {1} out
                end
                
                % After the parfor loop, convert cell array to matrix (if needed)
                Corr1 = NaN(obj.CorrelationInfo.nFrames-1, obj.CorrelationInfo.nFrames-1); % TOOK {1} out (2x)
                for dt = 1:obj.CorrelationInfo.nFrames-1 % TOOK {1} out
                    nVals = numel(Corr1_cell{dt});
                    Corr1(1:nVals, dt) = Corr1_cell{dt};
                end
                
                if nargin == 3
                    r = 1;
                    e = 1;
                else
                    r = varargin{1};
                    e = varargin{2};
                end

                Time = [1:size(Corr1, 2)]*obj.CorrelationInfo.ExpTime;
                CorrelationOutput = struct('Time', Time, 'Correlation', Corr1);
                obj.CorrelationOutput{r, e} = CorrelationOutput;
        end

        function PlotCorrelation(obj, MeanCorr, Title, SavePath)
            
            Fig1 = figure(1);
            fieldNames = fieldnames(MeanCorr);
            count = 0;
            for i = 1:numel(fieldnames(MeanCorr))               
                FieldName = fieldNames{i};
                if strcmp(FieldName, 'CCM')
                    Colour = 'r';
                else
                    Colour = 'b';
                end
                Matrix = [];
                for j = 1:size(MeanCorr.(FieldName), 1)
                    count = count + 1;
                    plot(MeanCorr.(FieldName).Time{j,1}./60, nanmean(MeanCorr.(FieldName).Corr{j,1}, 1),...
                        'Color', Colour, 'LineWidth', 1);
                    hold on

                    if j == 1
                       LegendColors{1,count} = FieldName;
                    else
                       LegendColors{1,count} = '';
                    end

                    Matrix =  [Matrix; MeanCorr.(FieldName).Corr{j,1}];
                end

                CorrMatrix.(FieldName) = Matrix;
            end

            legend(LegendColors)
            xlabel('TimeLag (min)')
            ylabel('Multiscale structural similarity')
            ylim([0 1])
            title(append('Correlation: ', Title(3:end), ' - all samples'))

            Fig2 = figure(2);
            LegendColors = {};

            for i = 1:numel(fieldnames(CorrMatrix)) 
                StdError = [];
                FieldName = fieldNames{i};
                LegendColors{1,2*i-1} = '';
                LegendColors{1,2*i} = FieldName;

                if strcmp(FieldName, 'CCM')
                    ColourError = "#F267C5";
                    ColourPlot = "#A2152F";
                else
                    ColourError = "#268CDD";
                    ColourPlot = "#0072BD";
                end

                for j = 1:size(CorrMatrix.(FieldName), 2) 
                    StdError(1,j) = nanstd(CorrMatrix.(FieldName)(j, :))./(sqrt(nnz(~isnan(CorrMatrix.(FieldName)(j, :)))));
                end

                errorbar(MeanCorr.(FieldName).Time{1,1}./60, nanmean(CorrMatrix.(FieldName), 1), StdError,...
                    'CapSize', 0, 'Color', ColourError, 'LineWidth', 3); 
                hold on
                plot(MeanCorr.(FieldName).Time{1,1}./60, nanmean(CorrMatrix.(FieldName), 1), 'Color', ColourPlot, 'LineWidth', 2);
                hold on
            end
            legend(LegendColors)
            xlabel('TimeLag (min)')
            ylabel('Multiscale structural similarity')
            title(append('Correlation: ', Title(3:end), ' - all samples together + error'))
            ylim([0 1])

            if isempty(obj.raw.movToLoad)
                Figure1Path = append(SavePath, filesep, 'CorrelationSamplesSeparate_NoSegmentation.png');
                Figure2Path = append(SavePath, filesep, 'CorrelationSamplesError_NoSegmentation.png');
            else
                Figure1Path = append(SavePath, filesep, 'CorrelationSamplesSeparate_', obj.raw.movToLoad, '.png');
                Figure2Path = append(SavePath, filesep, 'CorrelationSamplesError_', obj.raw.movToLoad, '.png');
            end

            saveas(Fig1, Figure1Path);
            disp("Correlation figure 1: Saved")
            saveas(Fig2, Figure2Path);
            disp("Correlation figure 2: Saved")

            disp(append("Done with ", Title, ' - Moving to the next one...'))
        end
    end
end