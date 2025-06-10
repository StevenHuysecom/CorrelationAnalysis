clc ;
clear ;
close all;

%% USER INPUT
expTime = 0.030; %in sec
T = 296.15; %temperature in Kelvin
R = 0.050; %Radius of particle in um;
fitRDiff = 4; %in number of data
minSize = 20; %frames
ext = '.mat';
MultiModalChannels = 0; %0 if Channels 1-8, 1 if Channels 9-16

MainFolder = 'S:\Dual Color\20250220\Test_stickyness_2D\PS100 green';
SizeFolder = {'sample1', 'sample2', 'sample3'};% 
SampleFolder = {'3_min', '4_min', '5_min', '6_min', '7_min', '8_min', '9_min', '10_min',....
                '11_min',  '12_min', '13_min', '14_min', '15_min', '16_min','17_min', '18_min'}; % '5_min', '6_min', '8_min', '10_min', '12_min', '14_min', '16_min', 

for t = 1:numel(SizeFolder)
    for r = 1:numel(SampleFolder)
        try
            %% Loading
            path = append(MainFolder, filesep, SizeFolder{t}, filesep, SampleFolder{r});
            folder = dir(path);
            name = append('trackResults', num2str(MultiModalChannels+1), '.mat');
            idx = contains({folder.name},name);
            folder(~idx) = [];
            
            f2Load = [folder(1).folder filesep folder(1).name];
            
            tmpData = load(f2Load);
            name = fieldnames(tmpData);
            data = tmpData.(name{1});
            
            %% Processing
            currMov =  data.traces{1,1};
            allHeight = cellfun(@height,currMov(1,:));
            idx = allHeight>minSize;
            currMov = currMov(1, idx);
            allRes = struct('msdX',0,'msdY',0,'msdZ',0,'msdR',0,'tau',0,'DX',0,'DY',0,'DZ',0,'DR',0,...
                'nX',0,'nY',0,'nZ',0,'nR',0,'aX',0,'aY',0,'aZ',0,'aR',0,'vX',0,'vY',0,...
                'vZ',0,'vR',0);
            allRes(length(currMov)).msdX = [];
            maxLength = max(allHeight);
            allMSDX = zeros(length(currMov),maxLength-1);
            allMSDY = allMSDX;
            %allMSDZ = allMSDY;
            allMSDR = allMSDY;
            for i = 1:length(currMov)
                currPart = currMov{i};
            
                coord = [currPart.col, currPart.row, currPart.z];
                CM = mean(coord,1);
                coord = coord-CM;
            
                %in X
                msdx = MSD.calc(coord(:,1)/10^3);%convert to um;
                tau = (1:length(msdx))'*expTime;
                allMSDX(i,1:length(msdx)) = msdx;
                DX   = MSD.getDiffCoeff(msdx,tau,fitRDiff,'1D');
                nX   = MSD.getViscosity(DX,R,T);
                aX   = MSD.getDiffTypeAlpha(msdx,expTime);
                vX   = coord(1,1) - coord(end,1)/10^3/(length(coord)*expTime); %um/s
            
                %inY
                msdy = MSD.calc(coord(:,2)/10^3);%convert to um;
                allMSDY(i,1:length(msdy)) = msdy;
                DY   = MSD.getDiffCoeff(msdy,tau,fitRDiff,'1D');
                nY   = MSD.getViscosity(DY,R,T);
                aY   = MSD.getDiffTypeAlpha(msdy,expTime);
                vY   = coord(1,2) - coord(end,2)/10^3/(length(coord)*expTime); %um/s
            
                %inZ
                % msdz = MSD.calc(coord(:,3)/10^3);%convert to um;
                % allMSDZ(i,1:length(msdz)) = msdz;
                % DZ   = MSD.getDiffCoeff(msdz,tau,fitRDiff,'1D');
                % nZ   = MSD.getViscosity(DZ,R,T);
                % aZ   = MSD.getDiffTypeAlpha(msdz,expTime);
                % vZ   = coord(1,3) - coord(end,3)/10^3/(length(coord)*expTime); %um/s
            
            
                %inR
                msdr = MSD.calc(coord(:, 1:2)/10^3);%convert to um;
                allMSDR(i,1:length(msdr)) = msdr;
                DR   = MSD.getDiffCoeff(msdr,tau,fitRDiff,'2D');
                nR   = MSD.getViscosity(DR,R,T);
                aR   = MSD.getDiffTypeAlpha(msdr,expTime);
                dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2 +...
                    (coord(1,3)-coord(end,3))^2);
                vR = dR/10^3/(length(coord)*expTime); %um/s
            
            
                allRes(i).msdX = msdx;% in um^2
                allRes(i).msdY = msdy;
                % allRes(i).msdZ = msdz;
                allRes(i).msdR = msdr;
                allRes(i).tau = tau; % in sec
            
            
                allRes(i).DX   = DX;% in um^2 /sec
                allRes(i).DY   = DY;% in um^2 /sec
                % allRes(i).DZ   = DZ;% in um^2 /sec
                allRes(i).DR   = DR;% in um^2 /sec
            
                allRes(i).nX   = nX;
                allRes(i).nY   = nY;
                % allRes(i).nZ   = nZ;
                allRes(i).nR   = nR;
                
                allRes(i).aX   = aX;
                allRes(i).aY   = aY;
                % allRes(i).aZ   = aZ;
                allRes(i).aR   = aR;
                
                allRes(i).vX   = vX;
                allRes(i).vY   = vY;
                % allRes(i).vZ   = vZ;
                allRes(i).vR   = vR;
                
                
                allRes(i).num  = length(msdx);
            end
            
            %%
            
            meanMSDR = nanmean(allMSDR,1);
            tau = (1:length(meanMSDR))'*expTime;
            DR   = MSD.getDiffCoeff(meanMSDR,tau,fitRDiff,'3D');
            nR   = MSD.getViscosity(DR,R,T);
            
            disp(['Planes ' num2str(MultiModalChannels*8+1) '-' num2str(MultiModalChannels*8+8) ': ' 'The diffusion coefficient is ' num2str(DR) ' \mum^2/s and the viscosity is ' num2str(nR) ' cp']);
            %%
            name = append('msdRes', num2str(MultiModalChannels+1), '.mat');
            filename = [path filesep 'msdRes1.mat'];
            save(filename,'allRes');
            h = msgbox('Data succesfully saved');
            
            DiffList = [];
            for i = 1:size(allRes, 2)
                DiffList(i,1) = allRes(i).DR;
            end
            DiffConstant{t,r} = DiffList;
        catch
        end
    end
end

meanMSDR = nanmean(allMSDR,1);
tau = (1:length(meanMSDR))'*expTime;
DR   = MSD.getDiffCoeff(meanMSDR,tau,fitRDiff,'2D');
nR   = MSD.getViscosity(DR,R,T);

disp(['Planes ' num2str(MultiModalChannels*8+1) '-' num2str(MultiModalChannels*8+8) ': ' 'The diffusion coefficient is ' num2str(DR) ' \mum^2/s and the viscosity is ' num2str(nR) ' cp']);
%%
name = append('msdRes', num2str(MultiModalChannels+1), '.mat');
filename = [path filesep 'msdRes1.mat'];
save(filename,'allRes');
h = msgbox('Data succesfully saved');


