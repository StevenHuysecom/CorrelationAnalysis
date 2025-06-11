clc
clear all
close all

%% Calibration info
path2ZCal = [];
path2SRCal = [];

%% Pathinfo
MainFolder = 'E:\DDM_TestData';
TimeFolders = {'Dancing cells'};
ProteinFolders = {'1_Vinculin'}; 
DiseaseFolders = {'CT'};
SampleFolders = {'Set2', 'Set3', 'Set4'};

%% Storing info about the file
file.MovToLoad = ['DataMasked']; %always DataMasked
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = true; %true to recalibrate;
file.ext   = '.tif'; %extenstion of video
path2Cal =  [];
dimension = '2D';
correctDrift = false;
info.TotTime = [];
info.PxSize = 81; %in nm
info.FWHM = 3; %Half with of gauss
info.multiModal = 0; 
info.detectionMethod = 'Intensity'; %always take Intensity
MakeMovie = 1; %Make movie from traces (takes long)

%% Detection parameters
detectParam.size = [5 100]; % min and mix size of FAs in pixels
trackParam.radius  = 300;%nm
trackParam.memory  = 3; %If trace is lost, how many frames to keep

%% Diffusion parameters:
diffInfo.fitRDiff = 4; %number of datapoints to plot MSD
diffInfo.T = 296.15; %temperature in Kelvin
diffInfo.minSize = 20; %min number of frames to be a trace


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
            diffInfo.ExpTime = FrameTimes.(key2);
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
                traces = trackingExp.getTraces3D;

                if MakeMovie == 1
                    %% get an amazing movie with the traces
                    trackingExp.MakeMovie(10, diffInfo.minSize, 500, 1000);
                end
                
                %% get the MSD, diffusion coefficient etc
                diffExp = Core.DiffusionMovie(diffInfo, file, traces);
                diffExp.CalcMSD;
                diffExp.CalcDiff;
                diffExp.CalcVelocity;
                diffExp.SaveData;
            end
        end
    end
end
