clc
clear all
close all

%% Pathinfo
file.path  = ['E:\DDM_TestData\Dancing cells\1_Vinculin\CCM\Set1\CorrelationAnalysis' ];

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
CorrelationInfo.ParticleSize = 0.050; %radius in µm
CorrelationInfo.ExpTime = 0.030; %s
CorrelationInfo.Temp = 296.15; %Temp in Kelvin

%% Loading data
CorrelationMovie = Core.CorrelationMovie(file,path2Cal,info,CorrelationInfo);
CorrelationMovie.calibrate();
CorrelationMovie.LoadAllFrames(correctDrift);

%% Extract DDM signal 
CorrelationMovie.main(CorrelationInfo, file);

