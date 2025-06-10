clc
clear all
close all

%% Calibration info
path2ZCal = [];
path2SRCal = [];

%% Pathinfo
MainFolder = 'S:\DDM_TestData';
TimeFolders = {'Dancing cells'};
ProteinFolders = {'1_Vinculin'}; 
DiseaseFolders = {'CCM'};
SampleFolders = {'Set1'};

%% Storing info about the file
file.MovToLoad = ['DataMasked']; %always DataMasked
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = true; %true to recalibrate;
file.ext   = '.tif';
path2Cal =  [];
dimension = '2D';
correctDrift = false;
info.TotTime = [];
info.PxSize = 81;
info.FWHM = 3;
info.multiModal = 0;
info.detectionMethod = 'Intensity';

%% Detection parameters
detectParam.size = [10 100]; % min and mix size of FAs in pixels
trackParam.radius  = 3500;%nm
trackParam.memory  = 3;

%% Input For Correlation thingies
FrameTimes = struct('Vinculin_CT', 6.949, 'Vinculin_CCM', 5.489, 'Paxillin_CT', 6.073, 'Paxillin_CCM', 5.364, ...
    'VECadherin_CT', 5.673, 'VECadherin_CCM', 5.592, 'ActinSPY555_CT', 5.721, 'ActinSPY555_CCM', 5.390, ...
    'ActinVinculin_CT', 6.949, 'ActinVinculin_CCM', 5.489, 'ActinPaxillin_CT', 6.073, 'ActinPaxillin_CCM', 5.364, ...
    'ActinVECadherin_CT', 5.673, 'ActinVECadherin_CCM', 5.592);


for c = 1:numel(TimeFolders)

    for o = 1:numel(ProteinFolders)
        key1 = ProteinFolders{o}(3:end);

        Time = cell(numel(SampleFolders), 1);
        Corr = cell(numel(SampleFolders), 1);
        
        for r = 1:numel(DiseaseFolders)
            
            key2 = sprintf('%s_%s', key1, DiseaseFolders{r});
            info.ExpTime = FrameTimes.(key2);
            if ~isempty(info.TotTime)
                CorrelationInfo.nFrames = ceil((TotTime *60)/info.ExpTime);
            end

            MeanCorr.(DiseaseFolders{r}) = table(Time, Corr);
            for e = 1:numel(SampleFolders)
                file.path = append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}, filesep,...
                        DiseaseFolders{r}, filesep, SampleFolders{e});
                info.file = file;

                %% get TrackingData
                trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);
                trackingExp.retrieveMovies;
                trackingExp.retrieveTrackData(detectParam,trackParam);

            end
        end
    end
end
