classdef MPParticleMovie < Core.MPMovie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        
        candidatePos
        unCorrLocPos
        corrLocPos
        particles
        intData
        BgcorrectionCh1Ch2
        
    end
    
    methods
        function obj = MPParticleMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParamFull, frames)
            %Method to perform localization on each plane for each frame
            %Check if some candidate exists already in the folder (previously saved)
            for q = 1:obj.info.multiModal+1
                if iscell(detectParamFull)
                    if q == 1
                        detectParam = detectParamFull{1};
                    elseif q == 2
                        detectParam = detectParamFull{2};
                    end
                else
                    detectParam = detectParamFull;
                end

                switch nargin
                        case 2
                            
                            frames = 1: obj.calibrated{1,q}.nFrames;
                            disp('Running detection on every frame');
         
                            
                        case 3
                            
                            [frames] = obj.checkFrame(frames,obj.calibrated{1,q}.nFrames);
                            
                        otherwise
                            
                            error('too many inputs');
                            
                 end

                 if isfield(obj.raw, "movToLoad")
                     folder = append('calibrated', obj.raw.movToLoad);
                 else
                     folder = append('calibrated', num2str(q));
                 end

                 path = append(obj.raw.movInfo{1,1}.Path, filesep, folder);
                 [run, candidate] = obj.existCandidate(obj.raw.movInfo{1,1}.Path, '.mat');
                
                %if we only ask 1 frame we always run
                if length(frames) == 1
                    run = true;
                end
                if run
                   
                   
                    %Localization occurs here
                    assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first or load previous data');
                    assert(nargin>1,'not enough input argument or accept loading of previous data (if possible)');
                    BgCorrFactor = {};

                    [candidate] = obj.detectCandidate(detectParam,frames,q);
                    
                elseif ~isempty(candidate)
                else
                    %help message
                    disp('getCandidatePos is a function that detects features in a movie');
                    disp('To work, it needs to receive a structure containing 2 detection parameter:');
                    disp('delta which is the radius around which detection is performed usually 6 pixels');
                    disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                    disp('Typically between 20 and 200');
                    
                end
                %if we only ask 1 frame we do not save
                if length(frames) >1
                    %save the data
                    if isfield(obj.raw, "movToLoad")
                        folder = append('calibrated', obj.raw.movToLoad);
                    else
                        folder = append('calibrated', num2str(q));
                    end
                      
                    fileName = sprintf('%s%s%s%scandidatePos.mat',obj.raw.movInfo{1,1}.Path,'\',folder,'\');
                    save(fileName,'candidate');
                else
                end

                candidatePos{q,1} = candidate;
                obj.candidatePos{q,1} = candidate;
                
            end
            obj.candidatePos = candidatePos;
            obj.info.detectParam = detectParam;
        end

        function ConvCandToPart(obj)
            f = waitbar(0,'Initializing...');
            SRList = [];
            for q = 1:size(obj.candidatePos, 1)
                
                for frame = 1:size(obj.candidatePos{q, 1},1)
                    waitbar(frame./size(obj.candidatePos{q, 1},1),f,append('Constructing particles - Channel ', num2str(q), ' out of ', num2str(size(obj.candidatePos, 1))));
                    NumPart = size(obj.candidatePos{q, 1}{frame, 1}, 1);
                    SRListFrame = [obj.candidatePos{q, 1}{frame, 1}.row, obj.candidatePos{q, 1}{frame, 1}.col, zeros(NumPart, 1), obj.candidatePos{q, 1}{frame, 1}.Intensity,...
                                        obj.candidatePos{q, 1}{frame, 1}.AreaPx, zeros(NumPart,1)+frame, obj.candidatePos{q, 1}{frame, 1}.Eccentricity, obj.candidatePos{q, 1}{frame, 1}.Orientation];
                    SRList = [SRList; SRListFrame];
                    nParticles(frame) = NumPart;
                    idx2TP(frame) = frame;
                    List{frame} = obj.candidatePos{q,1}{frame, 1};
                end

                Traces = [];
                nTraces = [];

                SRList = array2table(SRList, "VariableNames", {'row', 'col', 'z', 'Intensity', 'Area', 't', 'Eccentricity', 'Orientation'});

                Particle.List = List; 
                Particle.nParticles = nParticles;
                Particle.idx2TP = idx2TP;
                Particle.Traces = [];
                Particle.nTraces = [];
                Particle.SRList = SRList;

                obj.particles{q,1} = Particle;
            end
            close(f)
        end
        
        function getROIs(obj)
            for q = 1:obj.info.multiModal+1
                T = obj.candidatePos{q,1}{1,1};
                max_distance = obj.info.detectParam.consThresh;
                particles = struct([]);
                cleanedT = table();
                valid_rows = false(height(T), 1);

                for i = 1:height(T)
                    match_count = 1; 
                    current_row = T.row(i);
                    current_col = T.col(i);
                    current_plane = T.plane(i);
                    for j = i+1:height(T)
                        next_row = T.row(j);
                        next_col = T.col(j);
                        next_plane = T.plane(j);
                        if next_plane ~= current_plane + 1
                            break;
                        end
                        dist = sqrt((next_row - current_row)^2 + (next_col - current_col)^2);
                        if dist <= max_distance
                            match_count = match_count + 1;
                            current_row = next_row; 
                            current_col = next_col;
                            current_plane = next_plane;
                        end
                        if match_count >= 3
                            valid_rows(i:j) = true;
                        end
                    end
                end
                for i = 1:size(obj.candidatePos{q,1},1)
                    %obj.candidatePos{q,1}{i,1} = T(valid_rows, :); 
                    obj.ROI{q,1} = max_distance;
                end 
            end      
        end
            
            function [candidate] = getCandidatePos(obj, frames, q)
                %Extract the position of the candidate of a given frame
                [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                candidate = obj.candidatePos{q,1}{idx};
                
                if isempty(candidate)
                    
                    warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                    
                end
            end  
            
            function SRLocalizeCandidate(obj,detectParam,frames)
                for q = 1:obj.info.multiModal+1
                    if iscell(detectParam)
                        roiSize = detectParam{1,q}.delta;
                    else
                        roiSize = detectParam.delta;
                    end
                    assert(~isempty(obj.calibrated{1,q}),'Data should be calibrated to consolidate');
                    assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
                    assert(~isempty(obj.candidatePos{q,1}), 'No candidate found, please run findCandidatePos before consolidation');
                    folder = append('calibrated', num2str(q));

                    path = append(obj.raw.movInfo.Path, filesep, folder);
                    [run,locPos] = obj.existLocPos(path,'.mat');
                    % run = 1;
                    
                    if run
                        switch nargin
        
                            case 1
                                roiSize = 6;
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running SRLocalization on every frame with ROI of 6 pixel radius');
        
                            case 2
        
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running SRLocalization on every frame');
        
                            case 3
        
                                [frames] = obj.checkFrame(frames,obj.calibrated{1,q}.nFrames);
        
                            otherwise
        
                                error('too many inputs');
        
                        end
                                               
                        locPos = cell(size(obj.candidatePos{q,1}));
                        h = waitbar(0,'Fitting candidates ...');
                        nFrames = length(frames);
                        %Localization occurs here
                        for i = 1 : 1:nFrames
                            disp(['Fitting candidates: frame ' num2str(i) ' / ' num2str(nFrames)]);
                            idx = frames(i);
                            %#1 Extract Candidate Position for specific frame
                            % if obj.info.rotationalCalib == 1
                            %     data =  obj.calibrated{2, q};  
                            % else
                                [data] = obj.getFrame(idx, q);
                            % end


                            [frameCandidate] = obj.getCandidatePos(idx, q);
                            
                            if isempty(frameCandidate)
                                
                                warning('Frame %d did not contain any candidate',idx);
                                locPos{i} = [];
                                
                            else
                                
                                
                            locPos{idx} = obj.superResLocFit(data,frameCandidate,roiSize);
                                
                            end
                            waitbar(i/nFrames,h,['Fitting candidates: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                        end
                        close(h);
                    else
                    end
                        %save the data

                    folder = append('calibrated',num2str(q));

                    if run == 1
                        if nFrames > 1
                            fileName = sprintf('%s%s%s%sSRLocPos.mat',obj.raw.movInfo.Path,'\', folder, '\');
                            save(fileName,'locPos');
                        end
                    end
                    
                    %store in the object
                    obj.unCorrLocPos{q,1} = locPos;
                    obj.corrLocPos{q,1}   = locPos;
                end
            end
            
            function [locPos] = getLocPos(obj, frames, q, z)
                 %Extract the position of the candidate of a given frame
                [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                if isnan(z)
                    locPos = obj.corrLocPos{q,1}{idx};
                else
                    locPos = obj.corrLocPos{q,1}{z,1}{idx};
                end
               
                if isempty(locPos)
                    
                    warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                    
                end
            end
            
            function consolidatePlanes(obj,frames,detectParam)
                for q = 1: obj.info.multiModal +1
                    if iscell(detectParam)
                        consThresh = detectParam{1,q}.consThresh;
                    else
                        consThresh = detectParam.consThresh;
                    end
                    %Consolidation refers to connect molecules that were localized
                    %at similar position in different plane on a single frame.
                    assert(~isempty(obj.calibrated{1,q}),'Data should be calibrated to consolidate');
                    assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
                    assert(~isempty(obj.candidatePos{q,1}), 'No candidate found, please run findCandidatePos before consolidation');
                    assert(~isempty(obj.unCorrLocPos{q,1}),'Localization needs to be performed before consolidation');
                   
                    %Check if some particles were saved already.
                    folder = append('calibrated', num2str(q));

                    path = append(obj.raw.movInfo.Path, filesep, folder);
                    [run, particle] = obj.existParticles(path, '.mat');

                    
                    if run
                        %Check the number of function input
                        switch nargin
                            case 1
                                
                                frames = 1: obj.calibrated{1,q}.nFrames;
                                disp('Running consolidation on every frame with roi of 6 pixel');
                                consThresh = 4;
                            case 2
                                [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                                consThresh = 4;                       
                            case 3
                                [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                                assert(isnumeric(consThresh),'Consolidation threshold should be numeric');
                            otherwise
                                
                                error('Something wrong with number of input');
                                
                        end
                        
                        
                        nFrames = length(frames);
                        %allocate for storage
                        particleList = cell(1,obj.raw.maxFrame(1));
                        nParticles = zeros(1,obj.raw.maxFrame(1));
                        idx2TP = zeros(1,obj.raw.maxFrame(1));
                        h = waitbar(0,'Consolidating candidate ...');
                        
                        %Consolidation occurs here
                        for i = 1 :nFrames
                            disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                            idx = frames(i);
                            %#1 Extract localized Position for specific frame
                            z = NaN;
        
                            [fCandMet] = obj.getLocPos(idx, q, z);
                            
                            if isempty(fCandMet)
                                
                                warning('Frame %d did not contain any localized positions',idx);
                                particleList{idx} = [];
                                %particleList{idx}{1} = nan(5);
                                nParticles(idx) = 0;
                                
                            else
                                  %#2 Consolidate the position of the given frame
                                  %across plane
                                  %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
                                  switch obj.info.zMethod
                                      case 'Intensity'
                                         %we do not do anythin at the moment.
                                         focusMetric = fCandMet.magX+fCandMet.magY;
                                         
                                      case '3DFit'
                                         %we do not do anythin at the moment.
                                         focusMetric = fCandMet.magX+fCandMet.magY;
                                      case 'PSFE' 
                                        [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet.ellip,fCandMet.fMetric);
                                  end
                                    %reformating to keep the same format as how the data is saved
                                    %later
                                    fCandMet.fMetric = focusMetric;
                                    
                                    %focusMetric((1-corrEllip)>0.3) = NaN;
                                    
                                    %Plane Consolidation occur here
                                    [part] = obj.planeConsolidation(fCandMet,focusMetric,consThresh, q);
        
                                    %we delete empty cells from the array
                                    idx2Empty = cellfun(@isempty,part);
                                    part(idx2Empty(:,1),:) = [];
                           
                                    particleList{idx} = part;
                                    nParticles(idx) = length(part);
                                    
                                    if ~isempty(part)
                                        idx2TP(idx) = idx;
                                    end
                                
                            end
                            waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
                        end
                        close(h);
                        
                        %#3 Storing List
                        particle.List       = particleList;
                        particle.nParticles = nParticles;
                        particle.tPoint     = nFrames;
                        particle.idx2TP     = nonzeros(idx2TP);
                        particle.Traces     = [];
                        particle.nTraces    = [];
                        
                        folder = append('calibrated', num2str(q));

                        fileName = sprintf('%s%s%s%sparticle.mat',obj.raw.movInfo.Path,'\', folder, '\');
                        save(fileName,'particle');
                        obj.particles{q,1} = particle;
                        
                    elseif run == 0
                        obj.particles{q,1} = particle;
                    end

                end

                %%% Pass the particles that are not detected in the other
                %%% channel. When particles are added in the new channel,
                %%% recalculate their intensity. 
                if obj.info.rotational == 1
                    obj.PartChannelConsolidation;
                end
            end


            function PartChannelConsolidation(obj)
                %%% loop through frames

                tform = load(append(obj.raw.movInfo.Path, filesep, 'ChannelTransformations.mat'));
                tform = tform.transformation;
                f = waitbar(0, 'Initizalizing');
                for i = 1:size(obj.particles{1, 1}.List, 2)
                    D = [];
                    matches = [];
                    waitbar(i./size(obj.particles{1, 1}.List, 2), f, 'Consolidating particles in both channels')
                    ParticlesCh1 = obj.particles{1, 1}.List{1,i};
                    ParticlesCh2 = obj.particles{2, 1}.List{1,i};

                    %%% Delete particles that are visible in both channels.
                    %%% Those will not be passed to the other channel
                    Coord1 = [];
                    Coord2 = [];
                    for j = 1:size(ParticlesCh1,1)
                        %Coord1(j,:) = table2array(nanmean(ParticlesCh1{j, 1} (:,{'row', 'col'}), 1));
                        Coord1(j,:) = [table2array(nanmean(ParticlesCh1{j, 1} (:,{'col'}), 1)), table2array(nanmean(ParticlesCh1{j, 1} (:,{'row'}), 1))];
                    end
                    for j = 1:size(ParticlesCh2,1)
                        %Coord2(j,:) = table2array(nanmean(ParticlesCh2{j, 1} (:,{'row', 'col'}), 1));
                        Coord2(j,:) = [table2array(nanmean(ParticlesCh2{j, 1} (:,{'col'}), 1)), table2array(nanmean(ParticlesCh2{j, 1} (:,{'row'}), 1))];
                    end
                    if ~or(isempty(Coord2), isempty(Coord1))
                        % Coord2New = Transformation.Coords2toCoords1.b*Coord2*Transformation.Coords2toCoords1.T + Transformation.Coords2toCoords1.c;
                        Coord2New = transformPointsForward(tform, Coord2);
                        D = pdist2(Coord1, Coord2New);
                        [matches, costs] = matchpairs(D, 10);

                        PlaneCheck = [];
                        %%%Check if the matched particles have at least a plane in common & if the distance is small:
                        for k = 1:size(matches, 1)
                            CommonPlanes = intersect(ParticlesCh1{matches(k,1)}.plane, ParticlesCh2{matches(k,2)}.plane);
                            if numel(CommonPlanes) >= 1
                                PlaneCheck(k,1) = 1;
                            else
                                PlaneCheck(k,1) = 0;
                            end
                        end
                        matches(PlaneCheck == 0, :) = [];

                        %%%Delete particles that have a partner from the lists
                        ToRemove1 = zeros(size(ParticlesCh1, 1),1);
                        ToRemove2 = zeros(size(ParticlesCh2, 1),1);
                        for l = 1:size(matches, 1)
                            ToRemove1(matches(l,1),1) = 1;
                            ToRemove2(matches(l,2),1) = 1;
                        end
                        ParticlesCh1 = ParticlesCh1(ToRemove1 == 0);
                        ParticlesCh2 = ParticlesCh2(ToRemove2 == 0);
                    end

                    %%% Pass particles from Ch1 to Ch2
                    ParticlesCh1(cellfun(@isempty, ParticlesCh1)) = [];
                    NewParticlesCh2 = [obj.particles{2, 1}.List{1,i}, cell(size(obj.particles{2, 1}.List{1,i}, 1), 1)];
                    [data] = obj.getFrame(i, 2); % get frame from ch2
                    sig = [obj.info.sigma_px obj.info.sigma_px];
                    for l = 1:size(ParticlesCh1, 1)
                        CurrentParticle = ParticlesCh1{l, 1};
                        Coord = [CurrentParticle.col, CurrentParticle.row];
                        %Coord = Transformation.Coords1toCoords2.b*Coord*Transformation.Coords1toCoords2.T + Transformation.Coords1toCoords2.c;
                        Coord = transformPointsInverse(tform, Coord);
                        NewParticle = CurrentParticle;
                        NewParticle.row = Coord(:,2);
                        NewParticle.col = Coord(:,1);
                        for m = 1:size(NewParticle, 1)
                            if ~isnan(NewParticle.plane(m))
                                planeData = data(:,:, NewParticle.plane(m));
                                %AdjCoord = Transformation.Coords1toCoords2NotCorr.b*[NewParticle.rowNotCorr(m), NewParticle.colNotCorr(m)]*Transformation.Coords1toCoords2NotCorr.T + Transformation.Coords1toCoords2NotCorr.c;
                                AdjCoord = transformPointsInverse(tform, [NewParticle.colNotCorr(m), NewParticle.rowNotCorr(m)]);
                                [NewParticle.roiLims{m}] = EmitterSim.getROI(AdjCoord(1), AdjCoord(2),...
                                                    obj.info.detectParam.delta, size(planeData,2), size(planeData,1));
                                ROI = planeData(NewParticle.roiLims{m}(3):NewParticle.roiLims{m}(4), NewParticle.roiLims{m}(1):NewParticle.roiLims{m}(2));
                                [NewParticle.intensity(m), NewParticle.SNR(m), ~] = obj.getIntensityGauss(ROI,sig);
                            end
                        end
                        NewParticlesCh2{end+1,1 } = NewParticle;
                        NewParticlesCh2{end,2} = 'passed';
                    end

                    %%% Pass particles from Ch2 to Ch1
                    ParticlesCh2(cellfun(@isempty, ParticlesCh2)) = [];
                    NewParticlesCh1 = [obj.particles{1, 1}.List{1,i}, cell(size(obj.particles{1, 1}.List{1,i}, 1), 1)];
                    [data] = obj.getFrame(i, 1); % get frame from ch1
                    sig = [obj.info.sigma_px obj.info.sigma_px];
                    for l = 1:size(ParticlesCh2, 1)
                        CurrentParticle = ParticlesCh2{l, 1};
                        Coord = [CurrentParticle.col, CurrentParticle.row];
                        %Coord = Transformation.Coords2toCoords1.b*Coord*Transformation.Coords2toCoords1.T + Transformation.Coords2toCoords1.c;
                        Coord = transformPointsForward(tform, Coord);
                        NewParticle = CurrentParticle;
                        NewParticle.row = Coord(:,2);
                        NewParticle.col = Coord(:,1);
                        for m = 1:size(NewParticle, 1)
                            if ~isnan(NewParticle.plane(m))
                                planeData = data(:,:, NewParticle.plane(m));
                                AdjCoord = transformPointsForward(tform, [NewParticle.colNotCorr(m), NewParticle.rowNotCorr(m)]);
                                [NewParticle.roiLims{m}] = EmitterSim.getROI(AdjCoord(1), AdjCoord(2),...
                                                    obj.info.detectParam.delta, size(planeData,2), size(planeData,1));
                                ROI = planeData(NewParticle.roiLims{m}(3):NewParticle.roiLims{m}(4), NewParticle.roiLims{m}(1):NewParticle.roiLims{m}(2));
                                [NewParticle.intensity(m), NewParticle.SNR(m), ~] = obj.getIntensityGauss(ROI,sig);
                            end
                        end
                        NewParticlesCh1{end+1, 1} = NewParticle;
                        NewParticlesCh1{end, 2} = 'Passed';
                    end

                    obj.particles{1,1}.List{1,i} = NewParticlesCh1;
                    obj.particles{2,1}.List{1,i} = NewParticlesCh2;

                    obj.particles{1,1}.nParticles(1,i) = size(NewParticlesCh1, 1);
                    obj.particles{2,1}.nParticles(1,i) = size(NewParticlesCh2, 1);
                end
                close(f)
            end
            
            function [particle] = getParticles(obj,frames)
                for q = 1:obj.info.multiModal+1
                    %GetParticles
                    [idx] = obj.checkFrame(frames,obj.calibrated{1,q}.nFrames(1));
                    particle = obj.particles{q,1}.List{idx};
                    
                    if isempty(particle)
                        
                        warning('There was no particle found in this frames, check than you ran superResConsolidate method on this frame beforehand');
                        
                    end
                end
            end

            
            function showCandidate(obj,idx)
                %Display Candidate
                for q = 1:(obj.info.multiModal + 1)
                    % if obj.info.rotationalCalib ~= 1
                    assert(length(idx)==1, 'Only one frame can be displayed at once');
                    [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                    assert(~isempty(obj.candidatePos{q,1}{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
                    
                    [frame] = getFrame(obj,idx,q);
                    
                    nImages = size(frame,3);
                   
                    nsFig = ceil(nImages/4);
                    
                    if isempty(obj.corrLocPos)
                        candidate = obj.getCandidatePos(idx,q);
                    else
                        candidate = obj.corrLocPos{q, 1}{idx, 1};
                    end
                                               
                    rowPos    = candidate.row;
                    colPos    = candidate.col;
                    planeIdx  = candidate.plane;
                    if ismember('ParticlePassed', candidate.Properties.VariableNames)
                        passed    = candidate.ParticlePassed;
                    else
                        passed = zeros(size(candidate, 1) , 1);
                    end
                    h = figure();
                    h.Name = sprintf('Frame %d',idx);
                    for i = 1:nImages
                        passedIdx = passed(planeIdx==i);
                        colcoord = colPos(planeIdx==i);
                        rowcoord = rowPos(planeIdx==i);
                        
                        subplot(2,nImages/nsFig,i)
                        hold on
                        imagesc(frame(:,:,i))
                        colormap('hot')
                        hold on
                        for j = 1:size(colcoord, 1)
                            if passedIdx(j,1) == 0
                                color = 'g+';
                            elseif passedIdx(j,1) == 1
                                color = 'b+';
                            end
                            plot(colcoord(j),rowcoord(j),color,'MarkerSize',10)
                        end
                        axis image;
                        title({['Plane ' num2str(i)],sprintf(' Zpos = %0.3f',obj.calibrated{1,q}.oRelZPos(i))});
                        hold on
                    end
                    sgtitle(append('Channel ', num2str(q), ' - SR cal applied'))
                end

                for q = 1:(obj.info.multiModal + 1)
                    % if obj.info.rotationalCalib ~= 1
                    assert(length(idx)==1, 'Only one frame can be displayed at once');
                    [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                    assert(~isempty(obj.candidatePos{q,1}{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
                    
                    [frame] = getFrame(obj,idx,q);
                    
                    nImages = size(frame,3);
                   
                    nsFig = ceil(nImages/4);
                    
                    if isempty(obj.corrLocPos)
                        candidate = obj.getCandidatePos(idx,q);
                    else
                        candidate = obj.corrLocPos{q, 1}{idx, 1};
                    end
                                               
                    rowPos    = candidate.rowNotCorr;
                    colPos    = candidate.colNotCorr;
                    planeIdx  = candidate.plane;
                    if ismember('ParticlePassed', candidate.Properties.VariableNames)
                        passed    = candidate.ParticlePassed;
                    else
                        passed = zeros(size(candidate, 1) , 1);
                    end
                    h = figure();
                    h.Name = sprintf('Frame %d',idx);
                    for i = 1:nImages
                        passedIdx = passed(planeIdx==i);
                        colcoord = colPos(planeIdx==i);
                        rowcoord = rowPos(planeIdx==i);
                        
                        subplot(2,nImages/nsFig,i)
                        hold on
                        imagesc(frame(:,:,i))
                        colormap('hot')
                        hold on
                        for j = 1:size(colcoord, 1)
                            if passedIdx(j,1) == 0
                                color = 'g+';
                            elseif passedIdx(j,1) == 1
                                color = 'b+';
                            end
                            plot(colcoord(j),rowcoord(j),color,'MarkerSize',10)
                        end
                        axis image;
                        title({['Plane ' num2str(i)],sprintf(' Zpos = %0.3f',obj.calibrated{1,q}.oRelZPos(i))});
                        hold on
                    end
                    sgtitle(append('Channel ', num2str(q), ' - SR cal not applied'))
                end
            end
            
            function showParticles(obj,idx)
                %display particles (after consolidation), On top of the
                %localization, consolidated particles are circled.
                assert(length(idx)==1, 'Only one frame can be displayed at once');
                [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                % Show Candidate
                obj.showCandidate(idx);
                
                if isempty(obj.particles)
                    
                    warning('You did not consolidate the candidate yet, please use the consolidate method before showing the particles');
                    
                else
                    
                    if isempty(obj.particles.List(idx))
                        
                        warning('The candidates of the requested frame were not consolidated yet, only showing the candidate');
                        
                    else
                        
                        roiSize = obj.particles.roiSize;
                        nParticles = obj.particles.nParticles(idx);
                        h = gcf;
                        nPlanes = obj.calibrated.nPlanes;
                        colors = rand(nParticles,3);
                        %Display circled
                        for i = 1 : nPlanes
                            subplot(2,nPlanes/2,i)
                            hold on
                            for j = 1 : nParticles
                                currPart = obj.particles.List{idx}{j};
                                if(~isempty(currPart(currPart.plane == i,:)))
                                    part2Plot = currPart(currPart.plane == i,:);
                                    plot(part2Plot.col,part2Plot.row,'o',...
                                        'LineWidth',2, 'MarkerSize',10, 'MarkerEdgeColor',colors(j,:));
                                end
                            end
                            hold off
                        end
                        %Here we display a zoom onto the particle visible on
                        %the specific frame onto the consolidated planes
                        [frame] = getFrame(obj,idx);
                        assert(isstruct(frame),'Error unknown data format, data should be a struct');
                        for i = 1:nParticles
                            
                            currPart = obj.particles.List{idx}{i};
                            %Remove rows containing NaNs
                            idx2NaN = isnan(currPart.row);
                            currPart(idx2NaN,:) = [];
                            planes = currPart.plane;
                            figure(20+i)
                            hold on
                            for j = 1 : length(planes)
                                jdx = planes(j);
                                currFrame = frame.(sprintf('plane%d',jdx));
                                ROI = EmitterSim.getROI(currPart.col(j), currPart.row(j),...
                                    roiSize, size(currFrame,2), size(currFrame,1));
                                subplot(1,length(planes),j)
                                imagesc(currFrame(ROI(3):ROI(4),ROI(1):ROI(2)));
                                title({['Particle ' num2str(i)],[ ' Plane ' num2str(jdx)]});
                                axis image
                                colormap('jet')
                                
                            end
                            hold off
                        end
                    end
                end
            end
            
            function [candidateList] = planeConsolidation(obj,candMet,focusMetric,consThresh,q)
                %Loop through all candidate of a given frame and match them
                %between frame until none can be match or all are matched.
                nPlanes = obj.calibrated{1,q}.nPlanes;
                counter = 1;
                nPart = 0;
                maxIt = size(candMet,1);
                zMethod = obj.info.zMethod;
                candidateList = cell(max(size(find(~isnan(focusMetric)))),1);
                %continue until the list is not empty
                while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                    
                    if counter> maxIt
                        
                        error('While loop ran for an unexpectedly long time, something might be wrong');
                        
                    end
                    
                    %Find candidate in best focus
                    [~,idx] = max(focusMetric);
                    currentPlane = candMet.plane(idx);
                    
                    switch nPlanes
                        case 1
                            planes2Check = [];
                            
                        otherwise
                            
                            %Check which planes are to be checked (currently 2 planes
                            %above and 2 planes below the given plane
                            planes2Check = currentPlane-nPlanes:currentPlane-1;
                            planes2Check = planes2Check(planes2Check>0);
                            planes2Check = [planes2Check currentPlane+1:currentPlane+nPlanes];
                            planes2Check = planes2Check(planes2Check<nPlanes+1);
                            
                    end
                    currentCand = candMet(idx,:);
                    direction = -1;%Start by checking above

                    nCols = size(currentCand, 2);
                    
                    numericBlock1 = array2table(nan(nPlanes, nCols - 4), 'VariableNames', strcat("N1_", string(1:(nCols - 4))));
                    cellBlock = cell2table(cell(nPlanes, 2), 'VariableNames', strcat("C_", string(1:2)));
                    numericBlock2 = array2table(nan(nPlanes, 2), 'VariableNames', strcat("N2_", string(1:2)));

                    particle = [numericBlock1, cellBlock, numericBlock2];
                    particle.Properties.VariableNames = currentCand.Properties.VariableNames;
       
                    particle(currentCand.plane,:) = currentCand;
                    nCheck = length(planes2Check);
                    camConfig = obj.calibrated{1,q}.camConfig;
                    for i = 1:nCheck
                        
                        cand = candMet(candMet.plane == planes2Check(i),:);
                        if(planes2Check(i) > currentPlane)
                            direction = +1;%check below (Plane 1 is the uppest plane 8 is lowest)
                        end

                        [isPart] = Core.MPParticleMovie.isPartPlane(currentCand,cand,direction,consThresh,zMethod);
                        if ~all(isPart ==0)
                            id = cand.plane(isPart);
                            particle(id,:) = cand(isPart,:);
                        end
                        
                    end
                   
                    %We remove the particle(s) from the list
                    focusMetric(ismember(candMet(:,1), particle(:,1))) = [];
                    candMet(ismember(candMet(:,1), particle(:,1)),:) = [];
                    %format particle to be the same as before:
                    [particle] = obj.makeParticle(particle);
                    
                     %Check if the resulting configuration of the plane make
                    %sense e.g. no hole in the configuration
                     planeConfig = particle.plane;
                    [checkRes] = Core.MPParticleMovie.checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod);
                    
                    %Store
                    if checkRes
                        
                        nPart = nPart +1;
                        %store particle in a new list
                        candidateList{nPart,1} = particle;
                        
                    else
                        %Otherwise we remove it from the best focus search list
                        %by putting focus metric to NaN
                        %focusMetric(idx) = NaN;
                        
                    end
                    
                    counter = counter+1;
                    
                end
            end     
    end

     methods (Static)
       
         function [isPart]   = isPartPlane(current, next, direction,consThresh,zMethod)
            %This function aim at determining whether a candidate from one
            %plane and the another are actually the same candidate on
            %different plane or different candidate. The decision is based
            %on threshold on localization distance, ellipticity and focus
            %metric.
            
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');

            thresh = consThresh;
            [checkRes1] = Core.MPParticleMovie.checkEuDistRot([current.row, current.col],...
                [next.row, next.col],thresh);
            
            if strcmp(zMethod,'PSFE')
             % Test ellipticity
                [checkRes2] = Core.MPParticleMovie.checkEllipticity(current.ellip,...
                next.ellip,direction);
            
            elseif or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                %we do not test ellipticity here
                checkRes2 = checkRes1;
                
            else
                
                error('Unknown Z method for consolidation');
                
            end
            
            % Test focus Metric
            maxExpFM = current.fMetric+0.1*current.fMetric;
            checkRes3 = next.fMetric < maxExpFM;
            
            %isPart will only be true for particle that passes the 3 tests
            %isPart = checkRes1.*checkRes2.*checkRes3;
            isPart = checkRes1.*checkRes2.*checkRes3;
            
            PossCand = find(isPart);
            if(length(find(isPart))>1)
                if size(PossCand, 1) == 2
                    D = sqrt((next.row(PossCand(1)) - next.row(PossCand(2))).^2 + (next.col(PossCand(1)) - next.col(PossCand(2))).^2);
                    if D < consThresh
                        isPart(PossCand(2)) = 0;
                    else
                        warning('Could not choose which particle was the partner of the requested particle, killed them both');
                        isPart(isPart==1) = 0;
                    end
                else
                    warning('Could not choose which particle was the partner of the requested particle, killed them both');
                    isPart(isPart==1) = 0;
                end
            end
            
            if isempty(isPart)
                isPart = false;
            end
            
            isPart = logical(isPart);
            
        end
        
        function [checkRes] = checkEuDist(current,next,Thresh)
            %Use to check if the Euclidian distance is within reasonable
            %range
            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;
            % if rotational == 1
            %     [~,Idx] = min(EuDist);
            %     checkRes = zeros(size(EuDist));
            %     checkRes(Idx) = 1;
            % end
            if isempty(checkRes)
                checkRes = false;
            end
            
        end

        function [checkRes, EuDist] = checkEuDistRot(current,next,Thresh)
            %Use to check if the Euclidian distance is within reasonable
            %range
            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;

            if isempty(checkRes)
                checkRes = false;
            end
            
        end
        
        function [checkRes] = checkEllipticity(current, next, direction)
            %Use to check if the Ellipticity make sense with what we
            %expect from the behavior of the PSFEngineering plate
            switch direction
                
                case 1
                    
                    ellip = current < next+0.1*next;
                    
                case -1
                    
                    ellip = current +0.1 *current > next;
                    
            end
            
            checkRes = ellip;

            if isempty(checkRes)
                checkRes = false;
            end
            
          
            
        end
        
        function [checkRes] = checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            if or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                if nPlanes ==1
                    testPlanes = true;
                else
                    testPlanes = sum(~isnan(planeConfig))>=2;
                end
            elseif strcmp(zMethod,'PSFE')
                switch nPlanes
                    case 1
                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;
                    case 2

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;

                    case 4

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 2;

                    case 8

                        nPlanesEdgeFrange = 1;
                        nPlanesEdgeInterleaved = 2;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 3;
                    otherwise
                        error('Unknown number of planes, only expect 1,2,4,8')
                end
                %Let us test that we have consolidate the particle in at least
                %3 Planes
                isEdgePlane = or(~isempty(find(planeConfig==1,1)),~isempty(find(planeConfig==8,1)));


                switch camConfig
                    case 'fullRange'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

                    case 'interleaved'
                        if isEdgePlane
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeInterleaved;
                        else

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesInterleaved;
                        end
                    case 'equal'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

                    otherwise
                        error('unknown camera configuration');
                end
            
            end
            if testPlanes
                %We check that there is no "Gap" in the plane configuration
                % %as it would not make sense.
                % testConsec = diff(planeConfig(~isnan(planeConfig)));
                % checkRes = length(testConsec==1)>=2;
                checkRes = true;
            else
                
                checkRes = false;
                
            end
            
      

        end
        
        function [int,SNR]  = getIntensity(ROI,sig)
            %extract central position
            center = [ceil(size(ROI,1)/2),ceil(size(ROI,2)/2)];
            rowPos = center(1);
            colPos = center(2);
            %we integrate at 3*sigma (take ROI
            roiSignal = ceil(sig);
            if roiSignal(1) > center-1
                roiSignal = [center-1 center-1];
            end
            %get the idx for the ROI to integrate
            rowIdx = rowPos-roiSignal(1):rowPos+roiSignal(1);
            colIdx = colPos-roiSignal(2):colPos+roiSignal(2);
            
            %Pixel to integrate for signal
            px2SumInt = ROI(rowIdx,colIdx);
            
             %get background pixels
            bkg = ROI;
            bkg(rowIdx,colIdx) = 0;
            Padsize = (size(ROI, 1) - size(px2SumInt, 1))./2;
            if Padsize > 2
                CutIdx = Padsize - 2;
            else
                CutIdx = Padsize;
            end
            bkg(1:CutIdx, :) = 0;
            bkg(:, 1:CutIdx) = 0;
            bkg(end-(CutIdx-1):end, :) = 0;
            bkg(:, end-(CutIdx-1):end) = 0;
            px2SumBkg = bkg(bkg~=0);
            bkg = mean(px2SumBkg);
            bkgVar = std(px2SumBkg);
            %calculate signal
            int = px2SumInt - bkg;
            int = mean(int,'all');
            %SNR = max(max(px2SumInt))/bkgVar;
            SNR = sqrt(int);
            % if or(SNR<0,int<0)
            %     int = sum(sum(px2SumInt));
            %     SNR = sqrt(int);
            % end
           
        end

        function [int,SNR, gaussian_particle_mask]  = getIntensityGauss(ROI,sig)
            % --- Step 1: Create a Gaussian-weighted mask ---
            sz = 2*ceil(max(sig))+1;
            % [x, y] = meshgrid(-(sz(1)-1)/2:(sz(1)-1)/2, -(sz(1)-1)/2:(sz(1)-1)/2);
            % 
            % gauss_mask = exp(-(x.^2 / (2*sig(1)^2) + y.^2 / (2*sig(2)^2)));
            % gauss_mask = gauss_mask / sum(gauss_mask(:)); 
            gauss_mask = ones(sz, sz)./(sz*sz);
            % if size(ROI, 1) == size(ROI, 2)
            %     startIdx = (size(ROI, 1) - sz)/2+1;
            %     gaussian_particle_mask = zeros(size(ROI));
            %     gaussian_particle_mask(startIdx: startIdx+sz-1, startIdx:startIdx+sz-1) = gauss_mask;
            % else
                response = conv2(ROI, gauss_mask, 'same');
                [rows, cols] = size(response);
                if rows <= 5 && cols <= 5
                    M = response; % Keep as is
                else
                    M = zeros(size(response));
                    row_center = floor(rows/2) + 1;
                    col_center = floor(cols/2) + 1;
                    row_start = max(1, row_center - 2);
                    row_end   = min(rows, row_center + 2);
                    col_start = max(1, col_center - 2);
                    col_end   = min(cols, col_center + 2);
                    M(row_start:row_end, col_start:col_end) = response(row_start:row_end, col_start:col_end);
                end
                M = response;
                [maxVal, linearIdx] = max(M(:));
                [row, col] = ind2sub(size(M), linearIdx);

                %--- Step 3: Create shifted Gaussian mask at detected position ---
                half_sz = floor(sz/2);
                [x_full, y_full] = meshgrid(1:size(ROI,2), 1:size(ROI,1));
                gaussian_particle_mask = zeros(size(ROI));
    
                % Calculate bounds for placing the mask
                r_min = max(row - half_sz, 1);
                r_max = min(row + half_sz, size(ROI,1));
                c_min = max(col - half_sz, 1);
                c_max = min(col + half_sz, size(ROI,2));
    
                % Corresponding indices in the Gaussian mask
                gm_r_min = 1 + (r_min - (row - half_sz));
                gm_r_max = sz - ((row + half_sz) - r_max);
                gm_c_min = 1 + (c_min - (col - half_sz));
                gm_c_max = sz - ((col + half_sz) - c_max);
    
                % Ensure all indices are integers and scalars
                r_min = round(r_min); r_max = round(r_max);
                c_min = round(c_min); c_max = round(c_max);
                gm_r_min = round(gm_r_min); gm_r_max = round(gm_r_max);
                gm_c_min = round(gm_c_min); gm_c_max = round(gm_c_max);
    
                %Place the Gaussian mask inside the ROI
                gaussian_particle_mask(r_min:r_max, c_min:c_max) = ...
                    gauss_mask(gm_r_min:gm_r_max, gm_c_min:gm_c_max);
            % end

             % --- Step 5: Compute particle intensity ---
            int = sum(ROI .* gaussian_particle_mask, 'all');

            % --- Step 4: Background subtraction ---
            background_mask = gaussian_particle_mask == 0;
            background_mean = mean(ROI(background_mask));
            int = int - background_mean;

            SNR = int./background_mean;        
           
        end
       
     end
     
     methods (Access = protected)
        %method linked to candidate
        function [run,candidate] = existCandidate(obj,Path,ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'candidatePos')==true)
                        candidate = load([file2Analyze(1).folder filesep 'candidatePos.mat']);
                        candidate = candidate.candidate;
                        
                        if size(candidate,1)== obj.calibrated.nFrames
                            run = false;
                        else
                           disp('Detection missing in some frames, rerunning detection');
                           candidate = [];
                           run = true;
                        end
                            
                    else
                
                        run = true;
                        candidate =[];
                
                    end
                case 'run'
                    run = true;
                    candidate =[];
            end
        end
        
        %method linked to fitting
        function [run,SRLocPos] = existLocPos(obj,Path,ext) 
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'SRLocPos')==true)
                        SRLocPos = load([file2Analyze(1).folder filesep 'SRLocPos.mat']);
                        name = fieldnames(SRLocPos);
                        SRLocPos = SRLocPos.(name{1});
                        run = false;
                    else
                
                         run = true;
                        SRLocPos =[];
                
                    end
                case 'run'
                     run = true;
                     SRLocPos =[];
            end    
        end
        
        %method Linked to particles/planeConsolidation
        function [run, particle] = existParticles(obj,Path, ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        if isfield(particle, 'Particle')
                            particle = particle.Particle;
                        else
                            particle = particle.particle;
                        end
                        run = false;
                    else
                
                        run = true;
                        particle = [];
                
                    end
                    
                case 'run'
                    
                     run = true;
                     particle = [];
                     
            end
        end
        
        %Methods linked to Candidate
        function [candidate] = detectCandidate(obj,detectParam,frames,q)
            %Do the actual localization
            assert(~isempty(obj.calibrated),'Data should be calibrated to detect candidate');
            assert(isstruct(detectParam),'Detection parameter should be a struct with two fields');
            nFrames = length(frames);
            if q == 1
                if size(obj.candidatePos, 1) < 2
                    currentCandidate = obj.candidatePos;
                else
                    currentCandidate = obj.candidatePos{q,1};
                end
            elseif q == 2
                if size(obj.candidatePos, 1) == 1
                    currentCandidate = [];
                else
                    currentCandidate = obj.candidatePos{q,1};
                end
            end
            
            detectionMethod = obj.info.detectionMethod;

            if(isempty(currentCandidate))
                
                candidate = cell(obj.calibrated{1,q}.nFrames,1);
                
            else
                
                candidate = currentCandidate;
                
            end
           
            
            %parameter for localization
            h = waitbar(0,'detection of candidates...');

            for i = 1 : 1:nFrames
                
                position = table(zeros(10,1),zeros(10,1),zeros(10,1), zeros(10,1), zeros(10,1), zeros(10,1),...
                    zeros(10,1), zeros(10,1),'VariableNames',{'row', 'col', 'meanFAR','plane', 'AreaPx', 'Intensity', 'Eccentricity', 'Orientation'});
                [volIm] = obj.getFrame(frames(i),q);
                nPlanes = size(volIm,3);
                
                for j = 1:nPlanes

                        currentIM = volIm(:,:,j);
                        %localization occurs here
                        switch detectionMethod 
                            case 'Intensity'
                             
                                 bwImage = imbinarize(currentIM);

                                 % if strcmp(obj.info.file.MovToLoad, "Mask")
                                 %     bwImage = imcomplement(bwImage);
                                 % end

                                 % bwImage = im2gray(0 + bwImage);

                                 
                                 % SE = strel('disk',5);
                                 % bwImage = imopen(bwImage,SE);
                                 % bwImage = bwareaopen(bwImage,300);
        
                                 [ctr] = regionprops(bwImage, currentIM, 'Area','Centroid', 'MeanIntensity', 'Eccentricity', 'Orientation');
    
                                 %%% filter out too small and too big areas
                                    
                                 for k = 1:size(ctr, 1)
                                     if or(ctr(k).Area < detectParam.size(1), ctr(k).Area > detectParam.size(2))
                                         ctr(k).Area = [];
                                         ctr(k).Centroid = [];
                                         ctr(k).MeanIntensity = [];
                                         ctr(k).Eccentricity = [];
                                         ctr(k).Orientation = [];
                                     end
                                 end

                                 fun = @(ctr) all(structfun(@isempty,ctr));
                                 idx = arrayfun(fun,ctr);
                                 ctr(idx)=[];

                                 pos = cat(1,ctr.Centroid);
                                 if ~isempty(pos)
                                    startIdx = find(position.row==0,1,'First');
                                    if isempty(startIdx)
                                        startIdx = length(position.row)+1;
                                    end
                                     pos = flip(pos,2);
                                     pos(:,3) = NaN;
                                     pos(:,4) = j;
                                     pos(:,5) = [ctr.Area]';
                                     pos(:,6) = [ctr.MeanIntensity]';
                                     pos(:,7) = [ctr.Eccentricity]';
                                     pos(:,8) = [ctr.Orientation]';
                                     position(startIdx:startIdx+size(pos,1)-1,:) = array2table(pos);
                                 end
                        end    
                end
                
               

                idx = find(position.row==0,1,'First');
                if isempty(idx)
                    candidate{frames(i)} = position;                  
                else
                    candidate{frames(i)} = position(1:idx-1,:);
                    
                end

                waitbar(i/nFrames,h,...
                    sprintf('detection of candidates in Frame %d/%d done',i,nFrames));
            end 
            close(h);
        end
                
        function [candMet] = superResLocFit(obj,data,frameCandidate,roiSize)
            %Candidate metric are determined here (x,y,e,+focusmetric)
            delta = roiSize;
            
            %initialize table
            varNames = {'row','col','z','ellip','magX','magY','meanFAR','fMetric','gFitMet','plane'};
                candMet = table(zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    'VariableNames',varNames);
                sigSetup = [obj.info.sigma_px obj.info.sigma_px];
                
            for i = 1:size(frameCandidate,1)
                
                plane = frameCandidate.plane(i);
                planeData = double(data(:,:,plane));
                %Get the ROI
                [roi_lims] = EmitterSim.getROI(frameCandidate.col(i), frameCandidate.row(i),...
                    delta, size(planeData,2), size(planeData,1));
                roiLims{i,1} = roi_lims;
                ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                roiSize = size(ROI);
                if and(roiSize(1)==roiSize(2),mod(roiSize(1),2)==1)
                    if strcmpi(obj.info.fitMethod,'Phasor')
                        %Phasor fitting to get x,y,e
                        
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;

                    elseif strcmpi(obj.info.fitMethod,'Gauss')
                        [X,Y] = meshgrid(frameCandidate.col(i)-delta:frameCandidate.col(i)+...
                            delta,frameCandidate.row(i)-delta:frameCandidate.row(i)+delta);
                        domain(:,:,1) = X;
                        domain(:,:,2) = Y;

                        %Gauss (slower)
                        
                        x = [1:size(ROI,1)];
                        y = [1:size(ROI,1)];
                        meanROIhoriz = mean(ROI,1);
                        f = fit(x.',meanROIhoriz.','gauss2');
                        g = coeffvalues(f);
                        [~,idx] = min(abs([g(2), g(5)] - max(x)./2));
                        Width(1) = 2*sqrt(2*log(2))*(g(idx+2));
                        meanROIvert = mean(ROI,2);
                        f = fit(y.',meanROIvert,'gauss2');
                        g = coeffvalues(f);
                        [~,idx] = min(abs([g(2), g(5)] - max(x)./2));
                        Width(2) = 2*sqrt(2*log(2))*(g(idx+2));
                        Width = mean(Width);



                        [gPar] = Localization.Gauss.MultipleFitting(ROI,frameCandidate.col(i),...
                            frameCandidate.row(i),domain,1,Width);%data,x0,y0,domain,nbOfFit
                        colPos = gPar(5); %Should be directly the position of the particle as we
                            %gave above the domain of the ROI in the space of the image
                        rowPos = gPar(6);

                        row = rowPos - round(frameCandidate.row(i));
                        col = colPos - round(frameCandidate.col(i));

                        e = gPar(3)/gPar(2);

                        magX = gPar(1);
                        magY = gPar(1);
                    else
                        %Phasor fitting to get x,y,e
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;
                    end

                     %LRT focus metric
                    [fMetric,~] = Localization.likelihoodRatioTest(ROI,sigSetup,[row col]);

                    if magX>=magY
                        sig(1) = sigSetup(1) * magX/magY;
                        sig(2) = sigSetup(2);
                    else
                        sig(1) = sigSetup(1);
                        sig(2) = sigSetup(2) * magY/magX;
                    end

                    [int,SNR, gauss{i,1}] = obj.getIntensityGauss(ROI,sigSetup);
                    %LRT focus metric
                    [gFitMet,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);

                    %storing info
                    candMet.row(i) = rowPos;
                    candMet.col(i) = colPos;
                    candMet.z(i) = 0;
                    candMet.ellip(i) = e;
                    candMet.magX(i) = magX;
                    candMet.magY(i) = magY;
                    candMet.intensity(i) = int;
                    candMet.SNR(i) = SNR;
                    candMet.meanFAR(i) = frameCandidate.meanFAR(i);
                    candMet.fMetric(i) = fMetric;
                    candMet.gFitMet(i) = gFitMet;
                    candMet.plane(i) = plane;
                else
                    
                    %storing info
                    candMet.row(i)          = NaN;
                    candMet.col(i)          = NaN;
                    candMet.z(i)            = NaN;
                    candMet.ellip(i)        = NaN;
                    candMet.magX(i)         = NaN;
                    candMet.magY(i)         = NaN;
                    candMet.intensity(i)    = NaN;
                    candMet.SNR(i)          = NaN;
                    candMet.meanFAR(i)      = NaN;
                    candMet.fMetric(i)      = NaN;
                    candMet.gFitMet(i)      = NaN;
                    candMet.plane(i)        = NaN;
                end
            end
            %remove NaN's  
            idx = isnan(candMet.row);
            candMet(idx,:) = [];

            candMet.gauss = gauss;
            candMet.roiLims = roiLims;
        end
        
        function [doAvg]  = checkDoAverage(obj,ellip)
            camConfig = obj.calibrated{1,1}.camConfig;
            switch camConfig
                case 'fullRange'

                    if and(ellip>0.8,ellip<1.25)

                        doAvg = false;
                    else
                        doAvg = true;
                    end

                case 'interleaved'

                        doAvg = true;

                case 'equal'

                    if and(ellip>0.8,ellip<1.25)

                        doAvg = false;
                    
                    else
                        
                        doAvg = true;
                    
                    end

                otherwise
                    error('unknown camera config');
            end
        end
        
        function [newPart] = makeParticle(~,particleData)

            nCols = size(particleData, 2);
            nPlanes = 5;
            numericBlock1 = array2table(nan(nPlanes, nCols - 4), 'VariableNames', strcat("N1_", string(1:(nCols - 4))));
            cellBlock = cell2table(cell(nPlanes, 2), 'VariableNames', strcat("C_", string(1:2)));
            numericBlock2 = array2table(nan(nPlanes, 2), 'VariableNames', strcat("N2_", string(1:2)));

            newPart = [numericBlock1, cellBlock, numericBlock2];
            newPart.Properties.VariableNames = particleData.Properties.VariableNames;
     
            %store best focus in center
            [~,idx] = nanmax(particleData.fMetric);
            if idx-2 > 0
                newPart(1,:) = particleData(idx-2,:);
            end
            
            if idx-1 > 0
                newPart(2,:) = particleData(idx-1,:);
            end
            
            newPart(3,:) = particleData(idx,:);
            
            if idx+1 <= 8 && idx+1<= height(particleData)
                newPart(4,:) = particleData(idx+1,:);
            end
            
            if idx+2 <= 8 && idx+2 <= height(particleData)
                newPart(5,:) = particleData(idx+2,:);
            end
        
        end
       
     end
end

