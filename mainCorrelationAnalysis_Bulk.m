clc;
clear all;
close all;

%% Pathinfo
MainFolder = 'F:';
TimeFolders = {'Sarah'};
ProteinFolders = {'1_Vinculin'};%, '5_ActinVinculin', '6_ActinPaxillin', '7_ActinVECadherin'};
DiseaseFolders = {'CCM'};
SampleFolders = {'Set1'};

%% Storing info about the file
file.MovToLoad = ['DataSummedMask']; %For dancing cells: specify which mask to use, otherwise leave empty
                     %'RawData', 'DataMasked', 'Mask', 'DataSummedMask' % DO ALL THREE
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.calibrate = true; %true to recalibrate;
file.ext   = '.tif';
path2Cal =  [];
dimension = '2D';
correctDrift = false;
CorrelationInfo.TotTime = 40; % in minutes (if empty: [])

%% Input For Correlation thingies
FrameTimes = struct('Vinculin_CT', 6.949, 'Vinculin_CCM', 5.489, 'Paxillin_CT', 6.073, 'Paxillin_CCM', 5.364, ...
    'VECadherin_CT', 5.673, 'VECadherin_CCM', 5.592, 'ActinSPY555_CT', 5.721, 'ActinSPY555_CCM', 5.390, ...
    'ActinVinculin_CT', 6.949, 'ActinVinculin_CCM', 5.489, 'ActinPaxillin_CT', 6.073, 'ActinPaxillin_CCM', 5.364, ...
    'ActinVECadherin_CT', 5.673, 'ActinVECadherin_CCM', 5.592);

for c = 1:numel(TimeFolders)

    for o = 1:numel(ProteinFolders)
        key1 = ProteinFolders{o}(3:end); % key becomes 'Vinculin' etc.

        Time = cell(numel(SampleFolders), 1); %
        Corr = cell(numel(SampleFolders), 1); %
        
        for r = 1:numel(DiseaseFolders)
            
            key2 = sprintf('%s_%s', key1, DiseaseFolders{r});
            CorrelationInfo.ExpTime = 1;% FrameTimes.(key2);
            if ~isempty(CorrelationInfo.TotTime)
                CorrelationInfo.nFrames = ceil((CorrelationInfo.TotTime*60)/CorrelationInfo.ExpTime);
            end

            MeanCorr.(DiseaseFolders{r}) = table(Time, Corr);
            for e = 1:numel(SampleFolders)
                file.path = append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}, filesep,...
                        DiseaseFolders{r}, filesep, SampleFolders{e});
                CorrelationInfo.file = file;

                %% Loading data
                CorrelationMovie = Core.CorrelationMovie(file,path2Cal,info,CorrelationInfo);
                CorrelationMovie.calibrate();
                CorrelationMovie.LoadAllFrames(correctDrift);
                
                %% Extract correlation over time 
                [CorrelationOutput] = CorrelationMovie.main(CorrelationInfo, file, r, e);
                MeanCorr.(DiseaseFolders{r}).Time{e} = CorrelationOutput.Time;
                MeanCorr.(DiseaseFolders{r}).Corr{e} = CorrelationOutput.Correlation;

                %% Save info
                filenameRunInfo = append(file.path, filesep, sprintf('CorrelationInfo_%s.mat', file.MovToLoad));
                save(filenameRunInfo, 'CorrelationInfo');
                filenameRawCorr = append(file.path, filesep, sprintf('RawCorrelation_%s.mat', file.MovToLoad));
                save(filenameRawCorr, 'CorrelationOutput');
            end
        end
        filenameMeanCor = append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}, sprintf('MeanCorrelation_%s.mat', file.MovToLoad));
        save(filenameMeanCor, 'MeanCorr');
        CorrelationMovie.PlotCorrelation(MeanCorr, ProteinFolders{o}, append(MainFolder, filesep, TimeFolders{c}, filesep, ProteinFolders{o}));
    end
end




