clc
clear all
close all

%% Pathinfo
MainFolder = 'E:\DDM_TestData';
TimeFolders = {'Dancing cells'};
ProteinFolders = {'1_Vinculin'}; %'1_Vinculin', , '4_ActinSPY555', '6_ActinPaxillin', '7_ActinVECadherin'};'20250321', '20250328'
DiseaseFolders = {'CCM'};%, 'CCM'}; %};%
SampleFolders = {'Set2'};%,'Set3','Set4'}; %};%'data'

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = false; %true to recalibrate;
file.ext   = '.tif';
path2Cal =  [];
dimension = '2D';
correctDrift = false;

%% Input For DDM 
CorrelationInfo.nFrames = 400;% number of frames to load into memory
CorrelationInfo.FramesToAnalyze = 350 ; %Number of frames to analyze
CorrelationInfo.PixelSize = 0.081; %µm
% CorrelationInfo.ParticleSize = 0.050; %radius in µm
% CorrelationInfo.ExpTime = 0.030; %s
CorrelationInfo.Temp = 296.15; %Temp in Kelvin


FrameTimes = struct('Vinculin_CT', 6.949, 'Vinculin_CCM', 5.489, 'Paxillin_CT', 6.073, 'Paxillin_CCM', 5.364, ...
    'VECadherin_CT', 5.673, 'VECadherin_CCM', 5.592, 'ActinSPY555_CT', 5.721, 'ActinSPY555_CCM', 5.390, ...
    'ActinVinculin_CT', 6.949, 'ActinVinculin_CCM', 5.489, 'ActinPaxillin_CT', 6.073, 'ActinPaxillin_CCM', 5.364, ...
    'ActinVECadherin_CT', 5.673, 'ActinVECadherin_CCM', 5.592);
ParticleSizes = struct('Vinculin', 1, 'Paxillin', 1, 'VECadherin', 10, ...
    'ActinSPY555', 5, 'ActinVinculin', 5, 'ActinPaxillin', 5, 'ActinVECadherin', 5); 


for c = 1:numel(TimeFolders)

    for o = 1:numel(ProteinFolders)
        key1 = ProteinFolders{a}(3:end); %-> go from 3rd character (skip '1_' prefix)
        CorrelationInfo.ParticleSize = QInfo_ParticleSize_map.(key1);    

        for r = 1:numel(DiseaseFolders)
            key2 = sprintf('%s_%s', key1, DiseaseFolders{n});
            CorrelationInfo.ExpTime = time_frames_map.(key2);

            if isempty(TotTime)
                DDMInfo.FramesToAnalyze = 1;
            else
                DDMInfo.FramesToAnalyze = ceil((TotTime *60)/DDMInfo.ExpTime);
            end
            DDMInfo.nFrames = DDMInfo.FramesToAnalyze;

            for e = 1:numel(SampleFolders)
                file.path = append(MainFolder, filesep, TimeFolders{d}, filesep, ProteinFolders{a}, filesep,...
                        DiseaseFolders{n}, filesep, SampleFolders{c});
                CorrelationInfo.file = file;


            end
        end
    end
end




%% Loading data
CorrelationMovie = Core.CorrelationMovie(file,path2Cal,info,CorrelationInfo);
CorrelationMovie.calibrate();
CorrelationMovie.LoadAllFrames(correctDrift);

%% Extract DDM signal 
CorrelationMovie.main(CorrelationInfo, file);