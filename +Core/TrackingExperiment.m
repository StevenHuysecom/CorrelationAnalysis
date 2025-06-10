classdef TrackingExperiment < handle

    properties
    
    path
    ext
    trackMovies
    cal2D
    info
    ZCal
    SRCal
    SRCal2
    traces3D
    traces3Dcommon
        
    end

     methods
        function obj = TrackingExperiment(folder2Data,cal2D,info,SRCalPath,zCalPath)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            obj.cal2D = cal2D;
            obj.info = info;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath; 
            obj.SRCal2 = SRCalPath;
        end

        %Set function
        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            folderContent = dir(path);
            %Get how many folder are in the main folder
            idx = sum(cellfun(@sum,{folderContent.isdir}));
            %Matlab always store ., .. as folder for relative path so we
            %want to find more than 2 folder in folderContent.
            assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            
            obj.path = path;
            
        end

        function set.cal2D(obj,cal2D)
            if isempty(cal2D)
                obj.cal2D = [];
            else
                assert(ischar(cal2D), 'Path should be given as a string');
                assert(isfolder(cal2D), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')

                [file2Analyze] = Core.Movie.getFileInPath(cal2D,'2DCal.mat');

                if isempty(file2Analyze)
                    error('No 2D calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(and(isstruct(cal), and(isfield(cal,'camConfig'),isfield(cal,'file'))),...
                        '2D calibration is supposed to be a struct with 4 fields');
                    obj.cal2D = cal;

                end
            end

        end

        function set.SRCal(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal.cal = SRCal;
                obj.SRCal.path = SRCal;
            else
                assert(isfolder(SRCal), 'The given path is not a folder');
               

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration.mat');

                if isempty(file2Analyze)
                    error('No SR calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                        'SR calibration is supposed to be a struct with 2 fields');

                    obj.SRCal.cal = cal;
                    obj.SRCal.path = SRCal;
                end
            end
            
        end

        function set.SRCal2(obj,SRCal)
            if isempty(SRCal)
                obj.SRCal2.cal = SRCal;
                obj.SRCal2.path = SRCal;
            else
                assert(isfolder(SRCal), 'The given path is not a folder');

                if isfield(obj.info, 'multiModal')
                    if obj.info.multiModal == true
                        MultiModal = 1;
                    end
                elseif isfield(obj.info, 'multiTracking')
                    MultiModal = 1;
                end

                if MultiModal == true
                     [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration2.mat');

                    if isempty(file2Analyze)
                        error('No SR calibration file found in the given folder');
                    else
                        fileName = [file2Analyze.folder filesep file2Analyze.name];
                        cal = load(fileName);
                        field = fieldnames(cal);
                        cal = cal.(field{1});
                        assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                            'SR calibration is supposed to be a struct with 2 fields');
    
                        obj.SRCal2.cal = cal;
                        obj.SRCal2.path = SRCal;
                    end
                end
            end
        end

         function set.ZCal(obj,zCal)
            if isempty(zCal)
                obj.ZCal.cal = zCal;
                obj.ZCal.path = zCal;
            else
                assert(isfolder(zCal), 'The given path is not a folder');

                %Check Given path
                [file2Analyze] = Core.Movie.getFileInPath(zCal,'zCalibration.mat');

                if isempty(file2Analyze)
                    error('No z calibration file found in the given folder');
                else
                    fileName = [file2Analyze.folder filesep file2Analyze.name];
                    cal = load(fileName);
                    field = fieldnames(cal);
                    cal = cal.(field{1});
                    assert(isstruct(cal),'zCalibration is supposed to be in cells format');
                    assert(and(isfield(cal,'fitZParam'),isfield(cal,'calib')),...
                        'Something is wrong in the fields of your Z calibration');

                    obj.ZCal.cal = cal;
                    obj.ZCal.path = zCal;
                end
            end    
         end

        %get 3D traces
        
        function [traces3D] = getTraces3D(obj)
            traces3D = obj.traces3D;
        end


        function retrieveMovies(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory

            %Load the movies for planes 1-8
            for i = 3:size(folder2Mov,1)
                %Check if the directory
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);
               
                if ~isempty(file2Analyze)
                    
                    count = 0;
                    count = count+1;
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    if isfield(obj.info.file, "MovToLoad")
                        file.MovToLoad = obj.info.file.MovToLoad;
                    end
                    tmp = Core.MPTrackingMovie(file , obj.cal2D, obj.info, obj.SRCal.path, obj.ZCal.path);
                    
                    if count == 1
                        tmp.giveInfo;
                    else
                        tmp.info = obj.trackMovies.(['mov' num2str(1)]).getInfo; 
                    end
                    tmp.calibrate;
                    obj.trackMovies.(['mov' num2str(((i-2)*2)-1)]) = tmp;

                else
                    
                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ' obj.ext ' file and is therefore ignored']);
                    
                end
                
            end                
            
            if isempty(obj.trackMovies)              
               error(['No %s was found for planes 9-16. Please check:\n',...
                   '1) that you gave the correct file extension.\n',...
                   '2) that you gave the path of a folder containing folders containing movies with the given extension'],obj.ext);        
                
            end
            disp('=======> DONE ! <========')
        end

        function retrieveTrackData(obj,detectParam, trackParam)
            %Checking user input
            assert(nargin==3, 'retrieveZCalData expects 2 inputs, 1)detection Parameters, tracking parameter');
            %assert(and(isstruct(detectParam),and(isfield(detectParam,'chi2'),isfield(detectParam,'delta'))),'Detection parameter is expected to be a struct with 2 fields : "chi2"(~threshold for detection) and "delta"(size of window for test)');
            assert(and(isfield(trackParam,'radius'),isfield(trackParam,'memory')),...
                'Tracking parameter is expected to be a struct with two field "radius" and "memory"')
            fieldsN = fieldnames(obj.trackMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            allTraces = [];
            for i = 1: nfields
                
                disp(['Retrieving data from tracking file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                currentTrackMov = obj.trackMovies.(fieldsN{i});
                
                %Molecule detection
                currentTrackMov.findCandidatePos(detectParam);

                %Convert candidates to particles
                currentTrackMov.ConvCandToPart;
                
                %tracking occurs here
                currentTrackMov.trackParticle(trackParam);
                
                [traces] = currentTrackMov.getTraces;
                for q = 1:length(traces)
                    allTraces = [];
                    fileN = cell(length(traces{q,1}),1);
                    fileN(:,1) = {i};

                    [xStep,xMotor] = currentTrackMov.getXPosMotor;
                    [yStep,yMotor] = currentTrackMov.getYPosMotor;
                    [zSt,zMotor]   = currentTrackMov.getZPosMotor;

                    colMot = cell(length(traces{q,1}),1);
                    colMot(:,1) = {xMotor};
                    colStep = cell(length(traces{q,1}),1);
                    colStep(:,1) = {xStep};

                    rowMot = cell(length(traces{q,1}),1);
                    rowMot(:,1) = {yMotor};
                    rowStep = cell(length(traces{q,1}),1);
                    rowStep(:,1) = {yStep};

                    zMot = cell(length(traces{q,1}),1);
                    zMot(:,1) = {zMotor};
                    zStep = cell(length(traces{q,1}),1);
                    zStep(:,1) = {zSt};

                    allTraces = [allTraces; traces{q,1}(:), fileN,colStep,colMot,rowStep,rowMot,zStep,zMot ];
                    obj.traces3D{q,1} = allTraces;
                    obj.traces3D{q,2} = q; 
                end
            end
            
            
%             filename = [obj.path filesep 'traces3D.mat'];
%             save(filename,'allTraces');
            
            
            disp('=================> DONE <===================');
        end

        function ConsolidateChannels3(obj)
            channel1 = obj.traces3D{1,1};
            channel2 = obj.traces3D{2,1};

            %%% make sure channel1 has lowest number of traces
            numTraces1 = size(channel1, 1);
            numTraces2 = size(channel2, 1);
            if numTraces2 > numTraces1
                Swapchannel = channel1;
                channel1 = channel2;
                channel2 = Swapchannel;
                numTraces1 = size(channel1, 1);
                numTraces2 = size(channel2, 1);
            else
            end
        
            distances = zeros(numTraces1, numTraces2);
        
            f = waitbar(0,'Calculating distances between channels');
            for i = 1:numTraces1
                waitbar(i./numTraces1,f,'Calculating distances between channels');
                trace1 = channel1{i,1};
                for j = 1:numTraces2
                    trace2 = channel2{j,1}; 
                    
                    coords1 = table2array(trace1(:, 1:3)); 
                    time1 = table2array(trace1(:, 10));   
                    coords2 = table2array(trace2(:, 1:3));
                    time2 = table2array(trace2(:, 10)); 
                
                    common_time = intersect(time1, time2);
                    
                    if isempty(common_time)
                        avg_distance = Inf;
                    else
                        coords1_common = coords1(ismember(time1, common_time), :);
                        coords2_common = coords2(ismember(time2, common_time), :);
                        distance = sqrt(sum((coords1_common - coords2_common).^2, 2));
                        avg_distance = mean(distance);
                    end

                    distances(i, j) = avg_distance;
                end
            end
            close(f)
            
            [~, TracesIdx] = max(size(distances));
            [closest, closest_indices] = min(distances, [], TracesIdx);

            unique_vals = unique(closest_indices);
            n = 0;

            % while numel(unique_vals) < numel(closest_indices)
            %     n = n+1;
            %     duplicates = unique_vals(histc(closest_indices, unique_vals) > 1);
            % 
            %     for i = 1:numel(duplicates)
            %         duplicate_value = duplicates(i);
            %         duplicate_indices = find(closest_indices == duplicate_value);
            %         [MaxDup, MaxDupIdx] = max(closest(duplicate_indices));
            %         SearchPartner = duplicate_indices(MaxDupIdx);
            % 
            %         PossPartners = distances(SearchPartner, :);
            %         PossPartners(1, duplicate_value) = Inf;
            % 
            %         [closest(SearchPartner, 1), closest_indices(SearchPartner, 1)] = min(PossPartners);
            % 
            % 
            %     end
            % 
            %     MaxIterations = min([size(channel1, 1), size(channel2,1)])+10;
            % 
            %     if n == MaxIterations
            %         break
            %     end
            % 
            %     unique_vals = unique(closest_indices);
            % end

            for i = 1:size(closest_indices, 1)
                traces3Dcommon{i,1} = obj.traces3D{1,1}{i,1};
                traces3Dcommon{i,2} = obj.traces3D{2,1}{closest_indices(i),1};
            end

            obj.traces3Dcommon = traces3Dcommon;
        end
        
       
        % Select particles only when visible at the same time in both
        % channels:
        function ConsolidateChannels2(obj)
            Tresh = obj.info.euDist
            obj.traces3Dcommon = struct([]);

            %%% clean up traces
            TracesCh1 = struct([]);
            for i = 1:size(obj.traces3D{1,1})
                if size(obj.traces3D{1,1}{i,1},1) > 10
                    TracesCh1{end+1,1} = obj.traces3D{1,1}{i,1};
                end
            end
            TracesCh2 = struct([]);
            for i = 1:size(obj.traces3D{2,1})
                if size(obj.traces3D{2,1}{i,1},1) > 10
                    TracesCh2{end+1,1} = obj.traces3D{2,1}{i,1};
                end
            end
            
            mutualIndex = 1;
            for i = 1:size(TracesCh1,1)
                CurrentTrace1 = TracesCh1{i,1};
                for j = 1:size(TracesCh2,1)
                    CurrentTrace2 = TracesCh2{j,1};
                    [commonTimestamps, idx1, idx2] = intersect(CurrentTrace1(:, 10), CurrentTrace2(:, 10));
                    if ~isempty(commonTimestamps)
                        CurrentTrace1Overlap = CurrentTrace1(idx1, :);
                        CurrentTrace2Overlap = CurrentTrace2(idx2, :);
                        distances = sqrt(sum((CurrentTrace1Overlap(:, 1:3) - CurrentTrace2Overlap(:, 1:3)).^2, 2));
                        if all(distances.sum < Tresh)
                            obj.traces3Dcommon{end+1,1} = CurrentTrace1Overlap;
                            obj.traces3Dcommon{end,2} = CurrentTrace2Overlap;
                            mutualIndex = mutualIndex + 1;
                            break
                        end
                    end
                end
            end

            for i = 1:size(obj.traces3Dcommon)
                Int1 = obj.traces3Dcommon{i,1}.intensity;
                Int2 = obj.traces3Dcommon{i,2}.intensity;
                Ratio = Int1./Int2;
                obj.traces3Dcommon{i,3} = Ratio;
                Difference = Int1 - Int2;
                obj.traces3Dcommon{i,4} = Ratio;
            end
        end


        function ConsolidateChannels(obj)
            
            Tresh = obj.info.euDist
            obj.traces3Dcommon = struct([]);
            %Extract traces per channel
            for i = 1:(size(obj.traces3D, 1)./2)
                Traces0 = obj.traces3D{(i*2)-1,1};
                Traces1 = obj.traces3D{(i*2),1};
            end
                                    
            Trajec0 = struct([]);
            Trajec1 = struct([]);
            AlreadyChecked = [];
            %Open traces for channel one and two, check them one by one
            h1 = waitbar(0,'Consolidating Channels - Initializing');
            h2 = waitbar(0,'Consolidating Channels - Initializing');
            pos_w1=get(h1,'position');
            pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4) pos_w1(3) pos_w1(4)];
            set(h2,'position',pos_w2,'doublebuffer','on');
            set(h1,'doublebuffer','on');
            

            for j = 1:size(Traces0, 1)
                CurrentPartCh0 = Traces0{j,1};

                string1 = append('Consolidating Channels - Channel 0 particle ', num2str(j), '/', num2str(max(size(Traces0,1))));
                waitbar(j/max(size(Traces0,1)),h1,string1)
                Stop = 0;

                for k = 1:size(Traces1, 1);
                    if ismember(k, AlreadyChecked)
                        continue
                    end
                    CurrentPartCh1 = Traces1{k,1};

                    string2 = append('Consolidating Channels - Channel 1 particle ', num2str(k), '/', num2str(max(size(Traces1,1))));
                    waitbar(k/max(size(Traces1,1)),h2,string2)

                    % Check if particles have the same Timestamp:
                    SameTime = ismember(CurrentPartCh0.t, CurrentPartCh1.t);
                    Trace0 = [];
                    Trace1 = [];
                    
                    for l = 1:length(SameTime)
                        if SameTime(l) == true
                           Time = CurrentPartCh0.t(l);
                           rowCh0 = CurrentPartCh0.row(l);
                           colCh0 = CurrentPartCh0.col(l);
                           zCh0 = CurrentPartCh0.z(l);

                           idx = find(CurrentPartCh1.t == Time);
                           rowCh1 = CurrentPartCh1.row(idx);
                           colCh1 = CurrentPartCh1.col(idx);
                           zCh1 = CurrentPartCh1.z(idx);

                           EuDist = sqrt((rowCh0 - rowCh1).^2 + (colCh0 - colCh1).^2 + (zCh0 - zCh1).^2);
                           checkRes = EuDist < Tresh;
    
                           if checkRes == true
                               CP0 = table2array(CurrentPartCh0);
                               CP1 = table2array(CurrentPartCh1);

                               Trace0(l,:) = CP0(l, :);
                               Trace1(l,:) = CP1(idx,:); 

                               AlreadyChecked(end+1,1) = k;
                               Stop = 1;
                           end
                       end 
                    end

                   Trace0T = struct([]);
                   Trace1T = struct([]); 
                   VarNames = {'row','col','z', 'rowM', 'colM', 'zM', 'adjR', 'intensity', 'SNR', 't', 'rT', 'NumParticle'};

                   % If traces are not consequtive, split them.
                   % for Trace0
                   if ~isempty(Trace0) 
                       zeroRows = all(Trace0(:,:) == 0, 2);
                       zeroIndices = find(zeroRows);

                       startIdx = 1;

                       for i = 1:length(zeroIndices)
                           endIdx = zeroIndices(i) - 1; % Index before the zero row
                           if startIdx <= endIdx
                               TracesAdd = array2table(Trace0(startIdx:endIdx, :), 'VariableNames', VarNames);
                               Trace0T{end+1} = TracesAdd; % Add the sub-table to cell array
                           end
                           startIdx = zeroIndices(i) + 1; % Start of the next sub-table
                       end

                       if startIdx <= height(Trace0)
                           TracesAdd = array2table(Trace0(startIdx:end, :), 'VariableNames', VarNames);
                           Trace0T{end+1} = TracesAdd;
                       end
                   end

                   % for Trace1
                    if ~isempty(Trace1) 
                       zeroRows = all(Trace1(:,:) == 0, 2);
                       zeroIndices = find(zeroRows);

                       startIdx = 1;

                       for i = 1:length(zeroIndices)
                           endIdx = zeroIndices(i) - 1; % Index before the zero row
                           if startIdx <= endIdx
                              TracesAdd = array2table(Trace1(startIdx:endIdx, :), 'VariableNames', VarNames);
                              Trace1T{end+1} = TracesAdd; % Add the sub-table to cell array
                           end
                           startIdx = zeroIndices(i) + 1; % Start of the next sub-table
                       end

                       if startIdx <= height(Trace1)
                           TracesAdd = array2table(Trace1(startIdx:end, :), 'VariableNames', VarNames);
                           Trace1T{end+1} = TracesAdd;
                       end
                   end

                    if Stop == 1
                        break
                    end
                end
    
                if isempty(Trace0T) == false
                    for i = 1:length(Trace0T)
                        Trajec0{end+1,1} = Trace0T{1,i};
                    end
                    for i = 1:length(Trace1T)
                        Trajec1{end+1,1} = Trace1T{1,i};
                    end
                end

            end
            obj.traces3Dcommon = {Trajec0; Trajec1};
        end

         %Plotting for individual movies
        function showLoc(obj,idx)
             fieldsN = fieldnames(obj.trackMovies);
             maxIdx = length(fieldsN);
             assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
             
             currentMov = obj.trackMovies.(fieldsN{idx});
             
             currentMov.showCorrLoc;
        end
        
        % function showTraces(obj,idx)
        %      fieldsN = fieldnames(obj.trackMovies);
        %      maxIdx = length(fieldsN);
        %      assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
        % 
        %      currentMov = obj.trackMovies.(fieldsN{idx});
        % 
        %      currentMov.showTraces;
        % end

        function showTraces(obj)
            figure()
            for j = 1:length(obj.traces3Dcommon)
                traces = obj.traces3Dcommon{j,1};
                
                subplot(round(j./2), 2, j)
                hold on 
                for i = 1: length(traces)
                    currentTrace = traces{i};
                    plot3(currentTrace.col, currentTrace.row, currentTrace.z)
                    
                end
                hold on
            end
        end
        
        function evalAccuracy(obj,dim,idx)
            
            fieldsN = fieldnames(obj.trackMovies);
            maxIdx = length(fieldsN);
            assert(idx <= maxIdx,['Requested index to Movie is too large, only ' num2str(maxIdx) ' movies']);
            
            currentMov = obj.trackMovies.(fieldsN{idx});
            
            currentMov.evalAccuracy(dim);
            
        end
        
        function [int,SNR] = getAvgIntensity(obj)
            if ~isempty(obj.traces3Dcommon)
                assert(~isempty(obj.traces3Dcommon),'You need to extract 3D traces before extracting average intensity');
                for j = 1:length(obj.traces3Dcommon)
                    traces = obj.traces3Dcommon{j,1};
                    nTraces = length(traces);
                    int = zeros(nTraces,1);
                    SNR = zeros(nTraces,1);
                    for i = 1: length(traces)
                        currentTrace = traces{i};
                        intensity(i) = mean(currentTrace.intensity);
                        SNRatio(i) = mean(currentTrace.SNR);
                        
                    end
                    
                    int(j,1) = mean(intensity);
                    SNR(j,1) = mean(SNRatio);
                end
            elseif ~isempty(obj.traces3D)
                 for j = 1:length(obj.traces3D)
                    traces = obj.traces3D{j,1};
                    nTraces = length(traces);
                    int = zeros(nTraces,1);
                    SNR = zeros(nTraces,1);
                    for i = 1: length(traces)
                        currentTrace = traces{i};
                        intensity(i) = mean(currentTrace.intensity);
                        SNRatio(i) = mean(currentTrace.SNR);
                        
                    end
                    
                    int(j,1) = mean(intensity);
                    SNR(j,1) = mean(SNRatio);
                end
            else
                assert(~isempty(obj.traces3Dcommon),'You need to extract 3D traces before extracting average intensity');
            end
            
        end
        
        function [msd,traces] = getMSD(obj,dimension)
            
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting RMSD');
            
            traces = obj.traces3D;
            
            switch nargin
                case 1
                    
                    dimension = '3D';
                    
                case 2
                    
                otherwise
                    
                    error('too many input arguments');
                    
            end
            msd = cell(traces);
            MSDmat = zeros(obj.trackMovies.('mov1').raw.movInfo.maxFrame(1),length(traces));
            for i = 1 : size(traces,1)
                
                currentTrace = traces{i,:};
                
                if size(traces{i},1)>1
                    coord = [currentTrace.col,currentTrace.row,currentTrace.z];
                    [msdTmp,~] = MSD.calc(coord,dimension);
                    currentTrace.MSD = zeros(size(currentTrace.row,1),1);
                    currentTrace.MSD(1:end-1) = msdTmp;
                    traces{i} = currentTrace;
                    msd{i} = msdTmp;
                    MSDmat(1:length(msdTmp),i) = msdTmp(:);
                end
               
            end 
            sizes = cellfun(@size,traces,'UniformOutput',false);
            idxMat   = cellfun(@(x) x==11,sizes(:,1), 'UniformOutput', 0);
            idx = cellfun(@sum,idxMat,'UniformOutput',1);
            %delete traces where no MSD was calculated 
            traces(logical(~idx),:) = [];
            msd(logical(~idx),:) = [];
            obj.MSD = msd;
            obj.traces3D = traces;    
            
        end
        
        function saveData(obj)
            
            trackRes = struct; 
            disp('Saving Data');
            
            if ~isempty(obj.traces3Dcommon)
                
                trackData = obj.traces3Dcommon;
                MSDs = obj.MSD;
                    
                if ~isempty(MSDs)
                    trackRes.MSD = MSDs;
                end
                
                trackRes.traces = trackData; 
                trackRes.info = obj.info;
                trackRes.path = obj.path;
                filename = [obj.path filesep 'trackResultsCommonCh.mat'];
                save(filename,'trackRes');
                disp('Common channeldata was succesfully saved');
            end

            if ~isempty(obj.traces3D)
                for i = 1:length(obj.traces3D)
                    trackData = obj.traces3D{i,1};
                    % MSDs = obj.MSD;
                    % 
                    % if ~isempty(MSDs)
                    %     trackRes.MSD = MSDs;
                    % end
                    
                    trackRes.traces = trackData; 
                    trackRes.info = obj.info;
                    trackRes.path = obj.path;
                    name = append('trackResults',num2str(i),'.mat');
                    filename = [obj.path filesep name];
                    save(filename,'trackRes');
                    disp('Data were succesfully saved');
                end    
            else
                
                warning('No Data was saved because no traces or MSD could be found, please make sure you ran the analysis first');     
            end   
        end
        
        function showMSD(obj)
            MSD = obj.MSD;
            
            figure()
            hold on
            
            for i = 1:length(MSD)
                currentMSD = MSD{i};
                plot(currentMSD)
            end
                       
        end

        function MakeMovie(obj, sizeParticles, minSize, trailing, frameRate)
            %% Make Awesome movie
            filename = append(obj.info.file.path, filesep, 'AwesomeTraceMovie.gif');

            radius = 1000;
            yLimit = [0 obj.trackMovies.mov1.calibrated{1, 1}.Height];
            xLimit = [0 obj.trackMovies.mov1.calibrated{1, 1}.Width];

            Fig = figure;
            xlim(xLimit);
            ylim(yLimit);
            xlim manual;
            ylim manual;
            gcf;
            hold on
            
            [x,y,z] = sphere(32);
            x = x*sizeParticles/2;
            y = y*sizeParticles/2;
            z = z*sizeParticles/2;
            
            camlight
            lighting('gouraud');
            traces = obj.traces3D{1, 1};
            
            for i = 1 :obj.trackMovies.mov1.raw.maxFrame{1, 1}
                [frameA] = obj.trackMovies.mov1.getFrame(i, 1);
                imagesc(frameA);
                hold on
                for j = 1:size(traces,1)
                    currTrace = traces{j,1};
                    if height(currTrace) > minSize
                        idx2Frame = currTrace.t==i;
                        idx = i-trailing:i;
                        idx = ismember(currTrace.t,idx);
                        if and(~all(idx==0), ~all(idx2Frame ==0))
            
                            
                            data2Plot = currTrace(idx,:);
                            axis image
                            xlim(xLimit);
                            ylim(yLimit);
                            xlim manual;
                            ylim manual;
            
                            gcf;
                            hold on
            
                            plot((data2Plot.col)./(obj.info.PxSize),(data2Plot.row)./(obj.info.PxSize),'color',[0 0 0])
            
                            X = x+currTrace.col(idx2Frame);
                            Y = y+currTrace.row(idx2Frame);  

                            title(append(num2str(i*obj.info.ExpTime), 'sec out of ', num2str(obj.trackMovies.mov1.raw.maxFrame{1, 1}*obj.info.ExpTime)));
            
                        end
                    end
                end
                drawnow;
                frame = getframe(Fig);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
            
                if i == 1
            
                    imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);
            
                else
            
                    imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');
            
                end
                clf;
            
            end
            
            
        end
    end
end
