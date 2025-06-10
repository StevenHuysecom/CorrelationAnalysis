classdef MPLocMovie < Core.MPParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'protected')
        
        SRCal
        ZCal
        corrected
    end
    
    methods
        
        function obj = MPLocMovie(raw, MPCal,info, SRCal, zCal)
            
            obj  = obj@Core.MPParticleMovie(raw,MPCal,info);
            
            switch nargin
                
                
                case 4
                    obj.SRCal
                    obj.ZCal = [];
                case 5
                    
                    obj.SRCal = SRCal;
                    obj.ZCal = zCal;
               
                otherwise
                    error('Too many input arguments');
            end
        end
        
        function set.SRCal(obj,SRCal)
            
            if ~isempty(SRCal)
                assert(isfolder(SRCal), 'The given path is not a folder');

                %Check Given path
                for q = 1:obj.info.multiModal + 1
                    [file2Analyze] = Core.Movie.getFileInPath(SRCal, append('SRCalibration', num2str(q), '.mat'));
    
                    if isempty(file2Analyze)
                        error('No SR calibration file found in the given folder');
                    else
                        fileName = [file2Analyze.folder filesep file2Analyze.name];
                        cal = load(fileName);
                        field = fieldnames(cal);
                        cal = cal.(field{1});
                        assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                            'SR calibration is supposed to be a struct with 2 fields');
    
                        obj.SRCal{q,1} = cal; 
                    end
                end
            else
                obj.SRCal = [];
            end
        end
        
        function set.ZCal(obj,zCal)
            if ~isempty(zCal)
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

                    obj.ZCal = cal; 
                end
            else
                obj.ZCal = [];
            end
        end
        
        function applyCorr(obj,rot,refPlane)
            
            %apply SRCal
            obj.applySRCal(rot,refPlane);
            
            %transform ellipticity into Z
            obj.applyZCal;
            
        end
        
        function applySRCal(obj, rot, refPlane)
            for q = 1: obj.info.multiModal + 1
                assert(~isempty(obj.unCorrLocPos{q,1}),'You need to find candidate and SR Localized them before applying corrections');

                if isempty(obj.corrLocPos{q,1})

                        obj.corrLocPos = obj.unCorrLocPos{q,1};

                end

                if(~isempty(obj.SRCal))
                    if size(obj.corrLocPos{q,1},2) == 1
                        if nargin <2
                            refPlane = 5;
                        end
                        for z = 1:size(size(obj.corrLocPos{q,1},2))
                            data = obj.unCorrLocPos{q,1};
                            nPlanesCal = size(obj.SRCal{q,1}.trans,1)+1;
                            nPlanesFile = obj.calibrated{1,q}.nPlanes;
                            assert(nPlanesCal == nPlanesFile,'Mismatch between number of planes in SR calibration and file');             

                            disp(['Applying SR calibration...']);
                            for i = 1 : length(data)
                                currData = data{i,1};
                                if ~isempty(currData)
                                    currPlanes = unique(currData.plane);
                                    for j = 1 : length(currPlanes)
                                        currentPlane = currPlanes(j);
                                        data2Corr = currData(currData.plane==currentPlane,{'row','col','plane'});

                                        if rot
                                            corrMat = obj.SRCal{q,1}.rot; %% if ch2 is corrected, we only need the SRCal of ch1
                                            [corrData] = Core.MPSRCalMovie.applyRot(data2Corr, corrMat,refPlane);

                                        else
                                            corrMat = obj.SRCal{q,1}.trans; %% if ch2 is corrected, we only need the SRCal of ch1
                                             [corrData] = Core.MPSRCalMovie.applyTrans(data2Corr,corrMat,refPlane);                    
                                        end

                                        obj.corrLocPos{q,1}{i}(currData.plane==currentPlane,{'row','col','plane', 'rowNotCorr', 'colNotCorr'}) = corrData;        
                                    end
                                end

                            end
                            obj.corrected.XY = true;
                            disp('========> DONE ! <=========');
                        end
                    end            
                else
                    obj.corrected.XY = false;
                    disp('========> DONE ! <=========');
                    warning('SR Calibration not found, no correction was applied');
                end
            end
        end

        function CalcChannelTransition(obj, tform, frame)
            %%% Calculate transformation 2 times: Once on the coordinates
            %%% that are not corrected => to get intensity, once on the
            %%% coordinates that are SR corrected => for tracking

            % matchedCoords1 = [];
            % matchedCoords2 = [];
            % matchedCoords1NotCorr = [];
            % matchedCoords2NotCorr = [];
            % 
            % 
            % for frameIdx = 1:size(obj.unCorrLocPos{1,1}, 1)
            %     data1 = obj.corrLocPos{1,1}{frameIdx};  
            %     data2 = obj.corrLocPos{2,1}{frameIdx};
            %     if ~or(isempty(data1), isempty(data2)) 
            %         nPlanes = max(max(data1.plane), max(data2.plane));
            %         for i = 4
            %             PartPlane1 = table2array(data1(data1.plane == i,{'row', 'col'}));
            %             PartPlane2 = table2array(data2(data2.plane == i,{'row', 'col'}));
            %             D = pdist2(PartPlane1, PartPlane2);
            %             D(D > threshold) = Inf;
            %             [matches, costs] = matchpairs(D, 15);
            %             matchedCoords1 = [matchedCoords1; PartPlane1(matches(:,1), :)];
            %             matchedCoords2 = [matchedCoords2; PartPlane2(matches(:,2), :)];
            % 
            %             PartPlane1NotCorr = table2array(data1(data1.plane == i,{'rowNotCorr', 'colNotCorr'}));
            %             PartPlane2NotCorr = table2array(data2(data2.plane == i,{'rowNotCorr', 'colNotCorr'}));
            %             DNotCorr = pdist2(PartPlane1NotCorr, PartPlane2NotCorr);
            %             DNotCorr(DNotCorr > threshold) = Inf;
            %             [matchesNotCorr, costsNotCorr] = matchpairs(DNotCorr, 15);
            %             matchedCoords1NotCorr = [matchedCoords1NotCorr; PartPlane1(matchesNotCorr(:,1), :)];
            %             matchedCoords2NotCorr = [matchedCoords2NotCorr; PartPlane2(matchesNotCorr(:,2), :)];
            %         end
            %     end
            % end
            % 
            % [~, ~, transform] = procrustes(matchedCoords1, matchedCoords2);
            % transform.c = transform.c(1,:);
            % transform2.T = transform.T';
            % transform2.b = 1/(transform.b);
            % transform2.c = -(1 /transform.b) * transform.c * transform.T';
            % 
            % [~, ~, transformNotCorr] = procrustes(matchedCoords1NotCorr, matchedCoords2NotCorr);
            % transformNotCorr.c = transformNotCorr.c(1,:);
            % transformNotCorr2.T = transformNotCorr.T';
            % transformNotCorr2.b = 1/(transformNotCorr.b);
            % transformNotCorr2.c = -(1 /transformNotCorr.b) * transformNotCorr.c * transformNotCorr.T';

            for frameIdx = 1:frame
                data1 = obj.corrLocPos{1,1}{frameIdx};  
                data2 = obj.corrLocPos{2,1}{frameIdx};
                if ~or(isempty(data1), isempty(data2)) 

                    data1.OriginChannel(:,1) = 1;
                    data2.OriginChannel(:,1) = 2;
                    combined = [data1; data2];
    
                    for q = 1:obj.info.multiModal + 1
                        CombinedLoc = combined;
                        CombinedLoc.OriginChannel(CombinedLoc.OriginChannel == q) = 0;
                        CombinedLoc.OriginChannel(CombinedLoc.OriginChannel ~= 0) = 1;
    
                        PassedPart = CombinedLoc(CombinedLoc.OriginChannel == 1, :);
    
                        for i = 1:size(PassedPart,1)
                            Coords = [PassedPart.col(i) PassedPart.row(i)];
                            CoordsNotCorr = [PassedPart.colNotCorr(i) PassedPart.rowNotCorr(i)];

                            if q == 1
                                Coordsnew = transformPointsForward(tform, Coords);
                                CoordsnewNotCorr = transformPointsForward(tform, CoordsNotCorr);
   
                            elseif q == 2
                                Coordsnew = transformPointsInverse(tform, Coords);
                                CoordsnewNotCorr = transformPointsInverse(tform, CoordsNotCorr);

                            end
                            PassedPart.col(i) = Coordsnew(:,1);
                            PassedPart.row(i) = Coordsnew(:,2);

                            PassedPart.colNotCorr(i) = CoordsnewNotCorr(:,1);
                            PassedPart.rowNotCorr(i) = CoordsnewNotCorr(:,2);
                        end
    
                        CombinedLoc(CombinedLoc.OriginChannel == 1, :) = PassedPart;
    
                        obj.corrLocPos{q,1}{frameIdx,1} = CombinedLoc;%(~toRemove, :);
                        obj.corrLocPos{q,1}{frameIdx,1}.Properties.VariableNames{'OriginChannel'} = 'ParticlePassed';
                    end
                end
            end

            transformation = tform;

            Filename = append(obj.raw.movInfo.Path, filesep, 'ChannelTransformations.mat');
            save(Filename, "transformation");
        end
        
        function applyZCal(obj)
            for q = 1:obj.info.multiModal+1
                disp('Applying Z Calibration... ');
                assert(~isempty(obj.unCorrLocPos{q,1}),'Need to fit before applying the calibration');
                if isempty(obj.ZCal)
                    
                    warning('Z Calibration needed to correct the data, using Intensity instead');
                    if strcmp(obj.info.zMethod,'PSFE')
                        error('zMethod is selected is PSFE while no z calibration was provided')
                    end
                 
                    obj.corrected.Z = false;
                    disp('========> DONE ! <=========');
                end
                
                if isempty(obj.corrLocPos{q,1})
                    obj.corrLocPos{q,1} = obj.unCorrLocPos{q,1};
                    warning('Z calibration is currently being applied on non-SRCorrected (X-Y) data');
                end
                
                if size(obj.corrLocPos{q,1}, 2) == 1
                    for z = 1:size(obj.corrLocPos{q,1}, 1)
                        data = obj.corrLocPos{q,1}{z,1}; 
                        zCal = obj.ZCal;
                        zMethod = obj.info.zMethod;
                        
                        if or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                            obj.corrected.Z = false;
                           
                        elseif strcmp(zMethod,'PSFE')
                        
                            %we check which method is best:
                            [method] = obj.pickZFitMethod;
                            
                            %Here we translate ellipticity into z position based on
                            %calibration
                            nPlanesCal = size(zCal.calib,1);
                            nPlanesFile = obj.calibrated{1,q}.nPlanes;
                            assert(nPlanesCal == nPlanesFile,'Mismatch between number of planes in Z calibration and file');
            
                            disp('Applying Z Calibration using PSFE and ZCal');
                            for i = 1 : length(data)
                                currData = data{i};
                                nPos = size(currData,1);
            
                                for j = 1 : nPos
            
                                    currentEllip = currData.ellip(j);
                                    currentPlane = currData.plane(j);
                                    [zPos] = obj.getZPosition(currentEllip,zCal,currentPlane,method);
            
                                    obj.corrLocPos{q,1}{z,1}{i}.z(j) = zPos;
                                end
                                
                            end
                                     %Here we translate the ellipticity range into zRange for each
                            %plane
            
                            ellipRange = zCal.fitZParam.ellipRange;
                            nPlanes = obj.calibrated{1,q}.nPlanes;
                            zRange = cell(nPlanes,1);
                            
                            for i = 1 : nPlanes
                                zRange{i} = obj.getZRange(ellipRange,zCal,i,method);
                            end
                            
                            obj.corrected.Z = true;
                            obj.calibrated{1,q}.zRange{z} = zRange;
                            
                        else
                            error('Unknown Z method');
                        end
                        
                        disp('=======> DONE ! <========');

                    end
                else
                    data = obj.corrLocPos{q,1}; 
                    zCal = obj.ZCal;
                    zMethod = obj.info.zMethod;
                    
                    if or(strcmp(zMethod,'Intensity'),strcmp(zMethod,'3DFit'))
                        obj.corrected.Z = false;
                       
                    elseif strcmp(zMethod,'PSFE')
                    
                        %we check which method is best:
                        [method] = obj.pickZFitMethod;
                        
                        %Here we translate ellipticity into z position based on
                        %calibration
                        nPlanesCal = size(zCal.calib,1);
                        nPlanesFile = obj.calibrated{1,q}.nPlanes;
                        assert(nPlanesCal == nPlanesFile,'Mismatch between number of planes in Z calibration and file');
        
                        disp('Applying Z Calibration using PSFE and ZCal');
                        for i = 1 : length(data)
                            currData = data{i};
                            nPos = size(currData,1);
        
                            for j = 1 : nPos
        
                                currentEllip = currData.ellip(j);
                                currentPlane = currData.plane(j);
                                [zPos] = obj.getZPosition(currentEllip,zCal,currentPlane,method);
        
                                obj.corrLocPos{q,1}{i}.z(j) = zPos;
                            end
                            
                        end
                                 %Here we translate the ellipticity range into zRange for each
                        %plane
        
                        ellipRange = zCal.fitZParam.ellipRange;
                        nPlanes = obj.calibrated{1,q}.nPlanes;
                        zRange = cell(nPlanes,1);
                        
                        for i = 1 : nPlanes
                            zRange{i} = obj.getZRange(ellipRange,zCal,i,method);
                        end
                        
                        obj.corrected.Z = true;
                        obj.calibrated{1,q}.zRange = zRange;
                        
                    else
                        error('Unknown Z method');
                    end
                    
                    disp('=======> DONE ! <========');
                end
            end
        end
        
        function [locPos] = getLocPos(obj,frames,q,z)
             %Extract the position of the candidate of a given frame
             if isnan(z)
                [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                locPos = obj.corrLocPos{q,1}{idx};
             elseif isnumeric(z)
                 [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                locPos = obj.corrLocPos{q,1}{z,1}{idx};
             else
                 [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                locPos = obj.corrLocPos{q,1}{idx};
             end
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
                
        function superResolve(obj)
            for q = 1:obj.info.multiModal+1
                disp('super resolving positions ... ');
               
                %Check if some particle were super resolved already:
                if q == 1
                    [run,SRList] = obj.existZResParticles(obj.info.runMethod,append(obj.raw.movInfo.Path, filesep, 'calibrated1'),'.mat',q);
                elseif q == 2
                    [run,SRList] = obj.existZResParticles(obj.info.runMethod,append(obj.raw.movInfo.Path, filesep, 'calibrated2'),'.mat',q);
                else
                    [run,SRList] = obj.existZResParticles(obj.info.runMethod,obj.raw.movInfo.Path,'.mat');
                end
                % run = 1;

                if run
                    
                    data2Resolve = obj.particles{q,1}.List;
                    nPlanes = obj.calibrated{1,q}.nPlanes;
                    nParticles = sum(obj.particles{q,1}.nParticles);
                    pxSize = obj.info.pxSize;
                    SRList = table(zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1),zeros(nParticles,1),...
                            'VariableNames',...
                            {'row','col','z','rowM','colM','zM','adjR','intensity','SNR','t'});
                    nFrames = length(data2Resolve);
                    h = waitbar(0,'SuperResolving position...');
    
                    for i = 1:nFrames
                      
                        frameData = data2Resolve{i};
                        frameData2Store = table(zeros(size(frameData, 1),1),...
                            zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),...
                            zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),...
                            zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),zeros(size(frameData, 1),1),...
                            'VariableNames',...
                            {'row','col','z','rowM','colM','zM','adjR','intensity','SNR','t'});
                        
                        if strcmp(obj.info.zMethod,'PSFE')
                            for j = 1:length(frameData)
                                partData = frameData{j};
                                [data] = obj.resolveXYZ(partData(:,{'row','col','z','ellip','plane'}));
                                frameData2Store(j,{'row','col','z','rowM','colM','zM'}) = data;
                                frameData2Store.intensity(j) = partData.intensity(3);
                                frameData2Store.SNR(j) = partData.SNR(3);
                                frameData2Store.t(j) = i;
                            end
                            
                        else
                            
                            ROIRad = ceil(obj.info.FWHM_px/2+1);
                            fDataAll{1} = obj.getFrame(i, 1);
                            fDataAll{2} = obj.getFrame(i, 2);
 
                            for j = 1:length(frameData)
                                partData = frameData{j,1};
                                if obj.info.rotational == 1
                                    if isempty(frameData{j, 2})
                                        r = q;
                                    elseif strcmp(frameData{j,2}, 'Passed')
                                        if q == 1
                                            r = 2;
                                        elseif q == 2
                                            r = 1;
                                        end
                                    end
                                else
                                    r = q;
                                end
                                fData = fDataAll{r};
                                
                                
                                
                                switch obj.info.zMethod
                                    case 'Intensity'
                                        if nPlanes ==1
                                            row  = partData.row(3)*pxSize;
                                            col  = partData.col(3)*pxSize;
                                            z    = partData.z(3);
                                            rowM = partData.row(3)*pxSize;
                                            colM = partData.col(3)*pxSize;
                                            zM   = partData.z(3);
                                            adjR = 0; 
                                            data = table(row,col,z,rowM,colM,zM,adjR,...
                               'VariableNames',{'row','col','z','rowM','colM','zM','adjR'});
    
                                        else
                                            PartVolIm = obj.getPartVolIm(partData,ROIRad,fData);
                                            [data] = obj.resolveXYZInt(partData,PartVolIm, q);                            
                                        end
    
                                    case '3DFit'
                                        if nPlanes ==1
                                            row  = partData.row(3)*pxSize;
                                            col  = partData.col(3)*pxSize;
                                            z    = partData.z(3);
                                            rowM = partData.row(3)*pxSize;
                                            colM = partData.col(3)*pxSize;
                                            zM   = partData.z(3);
                                            data = table(row,col,z,rowM,colM,zM,...
                               'VariableNames',{'row','col','z','rowM','colM','zM'});
    
                                        else
                                            [data] = obj.resolveXYZ3DFit(partData(:,{'row','col','z','ellip','plane'}),fData,q);
                                        end
                                end
    
                                frameData2Store(j,{'row','col','z','rowM','colM','zM','adjR', 'intensity', 'SNR'}) = data;
                                frameData2Store.t(j) = i;
    
                            end
                        end
                    startIdx = find(SRList.row==0,1);   
                    SRList(startIdx:startIdx+height(frameData2Store)-1,:) = frameData2Store;   
                    waitbar(i/nFrames,h,['SuperResolving positions: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                    end
                    close(h);
                    %clean up the list
                    SRList(isnan(SRList.row),:) = [];
                    obj.particles{q,1}.SRList = SRList;
                    particle = obj.particles{q,1};
                    
                else
                    if ~isfield(obj.particles{q,1}, 'SRList')
                        for z = 1:size(obj.particles{q,1}, 1)
                            obj.particles{q,1}{z,1} = SRList;
                        end
                    else
                        obj.particles{q,1}.SRList = SRList;
                    end                
                    particle = obj.particles;
                end
                             
                %Save the data
                folder = append('calibrated', num2str(q));

                fileName = sprintf('%s%s%s%sparticle.mat',obj.raw.movInfo.Path,'\', folder, '\');
                profile('off')
                save(fileName,'particle');
                disp('========> DONE ! <=========');
            end
        end
                   
        function showCorrLoc(obj,frames)
             part = obj.particles.SRList;
            switch nargin
                case 1
                    frames = min(part.t):max(part.t);
                case 2 
                    [frames] = obj.checkFrame(frames,obj.raw.maxFrame(1));
            end
           
           
            figure()
            hold on
            
            sizeMarker =5;
            scatter3(part.col,part.row,part.z,sizeMarker,part.z,'filled')
            axis ij;

            title('all Localization plotted');
            xlabel('x position in nm');
            ylabel('y position in nm');
            zlabel('z position in nm');
            
            
            hold off
            
        end
        
    end
    
    methods (Static)
        
                    
    end
    
    
    methods (Access = protected)
        
        function [corrData] = applyTrans(obj, data2Corr, transMat, refPlane, currentPlane)
            %act depending on whether the current plane is smaller or
            %bigger than the user-selected reference plane
            
            if currentPlane < refPlane

                idx2Corr = currentPlane:refPlane-1;
                sign = -1;

            elseif currentPlane > refPlane

                idx2Corr = refPlane:currentPlane-1;
                idx2Corr = fliplr(idx2Corr);
                sign = +1;

            else

                idx2Corr = [];

            end
            %1 Translation
            
            row = data2Corr.row;
            col = data2Corr.col;
            if ~isempty(idx2Corr)
                for j = 1:length(idx2Corr)
                                       
                    row = row + sign* transMat.rowTrans(idx2Corr(j));
                    col = col + sign* transMat.colTrans(idx2Corr(j));
                    
                end
                data2Corr.row = row;
                data2Corr.col = col;
                
            end
            
            corrData = data2Corr;
        
        end
        
        function [corrData] = applyRot(obj, data2Corr, corrMat, refPlane, currentPlane)
            
            %act depending on whether the current plane is smaller or
                %bigger than the user-selected reference plane
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = false;
                    
                elseif currentPlane > refPlane
                    
                    idx2Corr = refPlane:currentPlane-1;
                    idx2Corr = fliplr(idx2Corr);
                    sign = true;
                    
                else
                    
                    idx2Corr = [];
                    
                end
                
                %Rotation
                if ~isempty(idx2Corr)
                    %Pad Z coordinate
                    data2C(:,1) = data2Corr.row;
                    data2C(:,2) = data2Corr.col;
                    data2C(:,3) = 0;
                    %remove center of mass (CM)
                    CM = mean(data2C);
                    data2C = data2C - CM;
                    
                    %Change the orientation of the data (should be [x;y;z]
                    %not [x y z]
                    data2C =  data2C';
                    
                    %Correction occur here
                    for j = 1:length(idx2Corr)
                        
                        if sign
                            rot = corrMat.rot{idx2Corr(j)}';
                        else
                            rot = corrMat.rot{idx2Corr(j)};
                        end
                           
                            data2C  = (rot*data2C);
                            data2Store = data2C';
                            corrData   = data2Store(:,1:2)+CM(1:2);
                    end
                else
                    corrData = data2Corr(:,1:2);
            
                end
  
        end
    
        function [zPos,inRange] = getZPosition(obj,val2Z,zCal,currentPlane,method)
            
            relZ = obj.calibrated.oRelZPos;

            zRange = zCal.fitZParam.zRange;
            zRange = zRange{currentPlane};
            zVec = zRange(1):1:zRange(2); %Here we assume accuracy >= 1nm

            switch method
                case 'poly'

                    fit = polyval(zCal.calib{currentPlane,1},zVec);

                case 'spline'
                    fit = ppval(zCal.calib{currentPlane,2},zVec);
            end

            %find the index of the value the closest to the particle
            %ellipticity
             [~,idx] = min(abs(fit-val2Z));

             zPos = zVec(idx)+ relZ(currentPlane)*1000;          
             inRange = and(val2Z>=zCal.fitZParam(1).ellipRange(1),...
                 val2Z<=zCal.fitZParam(1).ellipRange(2));
            
             if isempty(zPos)
                 disp('ouuups zpos is empty');
             end

        end
        
        function [zRange] = getZRange(obj,ellipRange,zCal,currentPlane,method)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
            
            switch method
                case 'poly'
                    
                    fit = polyval(zCal.calib{currentPlane,1},zVec);
                
                case 'spline'
                    fit = ppval(zCal.calib{currentPlane,2},zVec);
            end
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx1] = min(abs(fit-ellipRange(1)));
             [~,idx2] = min(abs(fit-ellipRange(2)));
             
             zPos1 = zVec(idx1)+ relZ(currentPlane)*1000;    
             zPos2 = zVec(idx2)+ relZ(currentPlane)*1000;
             
             zRange = [zPos1, zPos2];
             
        end
        
        function [data]  = resolveXYZ(obj,partData)
         
            pxSize = obj.info.pxSize;
            ellipRange = obj.ZCal.fitZParam.ellipRange;  

            idx2Keep = and(partData.ellip > ellipRange(1), partData.ellip < ellipRange(2));
            partData(~idx2Keep,:) = table(nan);

            row  = partData.row(3)*pxSize;
            col  = partData.col(3)*pxSize;
            z    = partData.z(3);
            data = table(row,col,z,'VariableNames',{'row','col','z'});
            %check how to perform averaging depending on the camera config
            [doAvg]  = obj.checkDoAverage(partData.ellip(3));

            if doAvg
                elliptRange = ellipRange(1):0.001:ellipRange(2);
                %we weigh the average later base on how much out of focus the
                %plane was.
                wRange1 = length(elliptRange(elliptRange<=1));
                wRange2 = length(elliptRange(elliptRange>=1));
                weight1 = linspace(1,5,wRange1);
                weight2 = linspace(5,1,wRange2);
                finalWeight = [weight1 weight2];
                ellipKept = partData.ellip(idx2Keep);
                idx = ellipKept;
                for k = 1 :length(ellipKept)

                    [~,idx(k)] = min(abs(elliptRange-ellipKept(k)));

                end

                weight = finalWeight(idx);
                %Weighed average
                row = sum(diag(partData.row(idx2Keep)* weight))/sum(weight) * pxSize;
                col = sum(diag(partData.col(idx2Keep)* weight))/sum(weight) * pxSize;
                z   = sum(diag(partData.z(idx2Keep)* weight))/sum(weight) * pxSize;
            end

            data.rowM = row;
            data.colM = col;
            data.zM = z;
            
        end
         
        function [data] = resolveXYZInt(obj,partData,partVolIm, q)
          
            pxSize = obj.info.pxSize;
          
            [~, bestFocusIdx] = nanmax(partData.fMetric);
            bf = partData.plane(bestFocusIdx);
            planePos = obj.calibrated{1,q}.oRelZPos;
            %partVolIm(partVolIm == 0) = NaN;
            %Get ROI XZ, YZ scaled to same pixel size
            [Mag] = Core.MPLocMovie.getZPhasorMag(partVolIm);

            domain = planePos;
            data   = [Mag.x]+[Mag.y];
            guess.sig = 2*obj.info.FWHM_px*pxSize/1000;
            guess.mu  = planePos(bf);
            guess.A   = max(data)-min(data);
           
%             [Res,fitData] = SimpleFitting.gauss1D(data,domain,guess);
%             RMS = sqrt(sum((data(:)-fitData(:)).^2)/length(data));
%             adjR = 1 - RMS.^2/var(data);
%           
            params = [guess.sig guess.mu guess.A min(data)];
         
            fun = @(x) SimpleFitting.minGauss1D(domain,data,x);
            opt = optimset('Display','off');
            % then we can do:
            [out, RMSD] = fminsearch(fun,params,opt);
            %normalize RMSD with mean data
            adjR = 1 - RMSD.^2/var(data);

            %[~, gaussFit] = SimpleFitting.minGauss1D(domain,data,out);
      
            z = out(2);
            %if the z position is out of bound we do not consider the data
            if or(z<min(domain),z>max(domain))
                z   = NaN;                           
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
            else
                z = z*1000;
                row = partData.row(3)*pxSize;
                col = partData.col(3)*pxSize;
                zM = z;                      
                rowM = partData.row(3)*pxSize;
                colM = partData.col(3)*pxSize;

            end
           
            Int = partData.intensity(bestFocusIdx);
            SNR = partData.SNR(bestFocusIdx);
            %store the data
            data = table(row,col,z,rowM,colM,zM,adjR,Int, SNR,...
                   'VariableNames',{'row','col','z','rowM','colM','zM','adjR', 'Int', 'SNR'});
            
   
        end
        
        function [Int] = getXYZIntRot(obj, partData, partVolIm,q)
            Mag = [];
            Nplanes = 0;
            for i = 1:size(partVolIm,3)
                Mag(end+1) = nanmean(partVolIm(:,:,i),'all');
            end
            Mag(Mag < 0) = 0;
            Int = nansum(Mag);
        end

        function [data] = resolveXYZ3DFit(obj,partData,ROI,q)
             
            pxSize = obj.info.pxSize;
            bf = partData.plane(3);
            planePos = obj.calibrated{1,q}.oRelZPos;

            x0 = mean([ROIs(3) ROIs(4)])*pxSize;
            y0 = mean([ROIs(1) ROIs(2)])*pxSize;
      
            z0 = planePos(bf)*1000;
            width.xy = 200;
            width.z  = 600;
            
            x = (ROIs(3):ROIs(4))*pxSize;
            y = (ROIs(1):ROIs(2))*pxSize;
            z = planePos*1000;
         
            [domX,domY,domZ] = meshgrid(x,y,z);
           
            dom(:,:,:,1) = domX;
            dom(:,:,:,2) = domY;
            dom(:,:,:,3) = domZ;

            %Multiple gaussian fitting occurs here
            [gPar,resnorm,res] = Localization.Gauss.MultipleGFit3D(ROI,x0,y0,z0,dom,1,width);
            z = gPar(:,7);

            if or(z<min(planePos*1000),z>max(planePos*1000))
                z   = NaN;                     
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
            else
                row = partData.row(3)*pxSize;
                col = partData.col(3)*pxSize;
                zM = z;                      
                rowM = partData.row(3)*pxSize;
                colM = partData.col(3)*pxSize;
            end
            %store the data
            data = table(row,col,z,rowM,colM,zM,...
                   'VariableNames',{'row','col','z','rowM','colM','zM'});
            
            
        end
        
        function [method] = pickZFitMethod(obj)
            
            names = fieldnames(obj.ZCal.zAccuracy);
            nMethod = numel(names);
            
            if nMethod == 1
                method = names{1};
            else
                for i = 1: nMethod
                  currentMethod = names{i};
                  currentAccuracy = obj.ZCal.zAccuracy.(names{i}).BestFocus;
                  
                  %Here accuracy should be small (high accuracy mean small
                  %number)
                  if i==1
                    finalMethod = currentMethod;
                  elseif and(i>1, or(currentAccuracy > finalMethod,currentAccuracy==0))
                      
                      %finalMethod stay the same
                  
                  elseif and(i>1, and(currentAccuracy <  finalMethod,currentAccuracy>0))
                  
                      finalMethod = currentMethod;              
                       
                  end
                                                    
                end
                method = finalMethod;
            end
            
            
        end
        
    end
    
    methods (Static)
        function [run, SRList] = existZResParticles(runMethod,Path, ext,q)
            SRList = [];
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        Particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        if isfield(Particle, 'Particle')
                            Particle = Particle.Particle;
                        else
                            Particle = Particle.particle;
                        end

                        if size(Particle, 2) == 1
                            for z = 1:size(Particle,1)
                                %Particle = particle{z,1};
                                try
                                    if isfield(Particle,'SRList')
                                        if ~isempty(Particle.SRList)
                                            run = false;
                                            SRList = Particle.SRList;
                                        else
                                            run = true;
                                        end
                                    elseif isfield(Particle{q,1},'SRList')
                                        if ~isempty(Particle{q,1}.SRList)
                                            run = false;
                                            SRList = Particle{q,1}.SRList;
                                        else
                                            run = true;
                                        end
                                    else
                                        run = true;
                                    end
                                catch
                                    run = true;
                                end

                            end
                        else
                            if isfield(Particle,'SRList')
                                if ~isempty(Particle.SRList)
                                    run = false;
                                    SRList = Particle.SRList;
                                else
                                    run = true;
                                end
                            elseif isfield(Particle{q,1},'SRList')
                                if ~isempty(Particle{q,1}.SRList)
                                    run = false;
                                    SRList = Particle{q,1}.SRList;
                                else
                                    run = true;
                                end
                            else
                                run = true;
                            end
                        end
                        
                    else
                
                        run = true;
                        
                
                    end
                    
                case 'run'
                    
                     run = true;
                     
                     
            end
        end
        function [Mag] = getZPhasorMag(ROI)

            %Possible improvement : Translate the coordinate of the best
            %focus into the otherplanes to extract the exact value where
            %the PSF should be taken    
           
            Mag = struct('x',zeros(1,size(ROI,3)),'y',zeros(1,size(ROI,3)));
            for i =1:size(ROI,3)
                [~,~,~,Mag(i).x,Mag(i).y] = Localization.phasor(ROI(:,:,i));
            end

        end
        
        function [partVolIm] = getPartVolIm(partData,ROIRad,volIm)
            %Extract data from particle in the 8 planes
            imSize = size(volIm);
            
            pos = [round(nanmean(partData.row)),round(nanmean(partData.col))];

            ROIs = Misc.getROIs(pos,ROIRad,imSize(1:2));

            partVolIm = volIm(ROIs(1):ROIs(2),ROIs(3):ROIs(4),:);
            
        end

        function [partVolIm] = getPartVolImRot(partData,ROIRad,volIm)
            imSize = size(volIm);
            
            non_empty_indices = [];
            for i = 1:size(partData,1)
                pos{1,i} = [round(partData.row(i)), round(partData.col(i))];
                if or(isnan(partData.row(i)), isnan(partData.col(i)))
                    succes(i) = 0;
                else
                    %non_empty_indices(end+1) = i;
                    succes(i) = 1;
                end
            end
            partVolIm = [];
            for i = 1:size(pos,2)
                if succes(i) == 0
                    % distances = abs(non_empty_indices - i);
                    % [~, closest_idx] = min(distances);
                    % closest_non_empty_index = non_empty_indices(closest_idx);
                    % pos{i} = pos{closest_non_empty_index};
                else
                    ROI = Misc.getROIs(pos{1,i},ROIRad,imSize(1:2));
                    partVolIm(:,:,end+1) = volIm(ROI(1):ROI(2),ROI(3):ROI(4),i); 
                end                    
            end               

            partVolIm(partVolIm == 0) = NaN; 
        end
    end
end