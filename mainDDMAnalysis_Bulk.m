clc
clear all
close all

%% Pathinfo
MainFolder = 'S:\DDM_TestData';
TimeFolders = {'Dancing cells'};
ProteinFolders = {'1_Vinculin'}; 
DiseaseFolders = {'CCM'};
SampleFolders = {'Set1\data'};

%% Storing info about the file
file.MovToLoad = ['Mask']; %For dancing cells: specify which mask to use, otherwise leave empty
                     %'RawData', 'DataMasked', 'Mask'
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = true; %true to recalibrate;
file.ext   = '.tif';
path2Cal =  [];
dimension = '2D';
correctDrift = false;
DDMInfo.TotTime = [];
DDMInfo.useROI = 'off'; %on
DDMInfo.PixelSize = 0.095; %µm
DDMInfo.Wavelength = 0.488; %µm
DDMInfo.NA = 1.20; %0.85 for 40x %1.2 for 63x
DDMInfo.Temp = 296.15; %Temp in Kelvin
DDMInfo.MinTimeLag = 1; %adjust timelag over which ISF is plotted, in number of frames % to skip frames
DDMInfo.Dimensions = [];
DDMInfo.MultiPlane = 1;

%% Parameters that vary sample to sample
time_frames_map = struct('Vinculin_CT', 6.949, 'Vinculin_CCM', 5.489, 'Paxillin_CT', 6.073, 'Paxillin_CCM', 5.364, ...
    'VECadherin_CT', 5.673, 'VECadherin_CCM', 5.592, 'ActinSPY555_CT', 5.721, 'ActinSPY555_CCM', 5.390, ...
    'ActinVinculin_CT', 6.949, 'ActinVinculin_CCM', 5.489, 'ActinPaxillin_CT', 6.073, 'ActinPaxillin_CCM', 5.364, ...
    'ActinVECadherin_CT', 5.673, 'ActinVECadherin_CCM', 5.592);
% time_frames_map = [];
QInfo_ParticleSize_map = struct('Vinculin', 1, 'Paxillin', 1, 'VECadherin', 10, ...
    'ActinSPY555', 5, 'ActinVinculin', 5, 'ActinPaxillin', 5, 'ActinVECadherin', 5); 
% QInfo_ParticleSize_map = [];
QInfo_Qmin_map = struct('Vinculin', 0, 'Paxillin', 0, 'VECadherin', 0, ...
    'ActinSPY555', 0, 'ActinVinculin', 0, 'ActinPaxillin', 0, 'ActinVECadherin', 0); 
% QInfo_Qmin_map = [];
QInfo_Qmax_map = struct('Vinculin', 10, 'Paxillin', 10, 'VECadherin', 10, ...
    'ActinSPY555', 10, 'ActinVinculin', 10, 'ActinPaxillin', 10, 'ActinVECadherin', 10);

%% Input For Correlation thingies
FrameTimes = struct('Vinculin_CT', 6.949, 'Vinculin_CCM', 5.489, 'Paxillin_CT', 6.073, 'Paxillin_CCM', 5.364, ...
    'VECadherin_CT', 5.673, 'VECadherin_CCM', 5.592, 'ActinSPY555_CT', 5.721, 'ActinSPY555_CCM', 5.390, ...
    'ActinVinculin_CT', 6.949, 'ActinVinculin_CCM', 5.489, 'ActinPaxillin_CT', 6.073, 'ActinPaxillin_CCM', 5.364, ...
    'ActinVECadherin_CT', 5.673, 'ActinVECadherin_CCM', 5.592);

for c = 1:numel(TimeFolders)

    for o = 1:numel(ProteinFolders)

        if isempty(QInfo_ParticleSize_map)
            DDMInfo.ParticleSize = 0.020;
        else
            key1 = ProteinFolders{o}(3:end); %-> go from 3rd character (skip '1_' prefix)
            DDMInfo.ParticleSize = QInfo_ParticleSize_map.(key1);
        end
        if ~isempty(QInfo_Qmin_map)
            DDMInfo.Qmin = QInfo_Qmin_map.(key1);
        else
            DDMInfo.Qmin = [];
        end

        if ~isempty(QInfo_Qmax_map)
            DDMInfo.Qmax = QInfo_Qmax_map.(key1);
        else
            DDMInfo.Qmax = []; 
        end

        key1 = ProteinFolders{o}(3:end);

        Time = cell(numel(SampleFolders), 1);
        Corr = cell(numel(SampleFolders), 1);
        
        for r = 1:numel(DiseaseFolders)
            
            if ~isempty(time_frames_map)
                key2 = sprintf('%s_%s', key1, DiseaseFolders{r});
                DDMInfo.ExpTime = time_frames_map.(key2);
            else 
                DDMInfo.ExpTime = 1;
            end
            
            if isempty(DDMInfo.TotTime)
                DDMInfo.FramesToAnalyze = 1;
            else
                DDMInfo.FramesToAnalyze = ceil((TotTime *60)/DDMInfo.ExpTime);
            end
            DDMInfo.nFrames = DDMInfo.FramesToAnalyze;

            MeanCorr.(DiseaseFolders{r}) = table(Time, Corr);
            for e = 1:numel(SampleFolders)
                file.path = append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}, filesep,...
                        DiseaseFolders{r}, filesep, SampleFolders{e});
                DDMInfo.file = file;

                %% Loading data
                DDMMovie = Core.DDMMovie(file,path2Cal,info,DDMInfo);
                DDMMovie.calibrate();
                DDMMovie.LoadAllFrames2(correctDrift);
                
                %% Extract correlation over time 
                [CorrelationOutput] = DDMMovie.main(DDMInfo, file,'NumBins',30);
                MeanCorr.(DiseaseFolders{r}).Time{e} = CorrelationOutput.Time;
                MeanCorr.(DiseaseFolders{r}).Corr{e} = CorrelationOutput.Correlation;
            end
        end
        CorrelationMovie.PlotCorrelation(MeanCorr, ProteinFolders{o}, append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}));
    end
end