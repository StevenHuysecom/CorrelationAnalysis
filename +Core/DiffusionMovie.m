classdef DiffusionMovie < handle
    %DIFFUSIONMOVIE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        info;
        file;
        Alltraces;
        Traces
        MSD;
        Tau;
        Diff;
        Velocity;
        Results;
        Alpha;

    end
    
    methods

        function obj = DiffusionMovie(info, file, traces)
            obj.info = info;
            obj.file = file;
            obj.Alltraces = traces; 
        end

        function CalcMSD(obj)

            f = waitbar(0,'Initializing');
            currMov = obj.Alltraces{1, 1};
            allHeight = cellfun(@height,currMov(:, 1));
            idx = allHeight>obj.info.minSize;
            currMov = currMov(idx, 1);
            for i = 1:size(currMov, 1)
                waitbar(i./size(currMov, 1),f,'Calculating MSD');
                Trace = currMov{i, 1};
                Coord = [Trace.row, Trace.col]*10^-3; %in µm
                [MSD,D] = obj.calc(Coord);
                MSDList{i,1} = MSD;
                TauList{i,1} = [1:1:size(MSD, 1)]*obj.info.ExpTime;
            end
            obj.Traces = currMov;
            obj.MSD = MSDList;
            obj.Tau = TauList;
            close(f)
        end
        
        function CalcDiff(obj)
            f = waitbar(0,'Initializing');
            for i = 1:size(obj.MSD, 1)
                waitbar(i./size(obj.MSD, 1),f,'Calculating diffusion coefficient');
                MSD = obj.MSD{i,1};
                Tau = obj.Tau{i,1};
                DList(i,1) = obj.getDiffCoeff(MSD,Tau,obj.info.fitRDiff, '2D');
            end
            obj.Diff = DList;
            close(f)
        end

        function CalcAlpha(obj)
            f = waitbar(0,'Initializing');
            for i = 1:size(obj.MSD, 1)
                waitbar(i./size(obj.MSD, 1),f,'Calculating alpha diffusion exponent');
                MSD = obj.MSD{i,1};
                t = obj.Tau{i,1};
                toFit = log(MSD);
                fitPar = fit(log(t(:)),toFit(:),'a*x+b');
                AlphaList(i,1) = fitPar.a;
            end
            obj.Alpha = AlphaList;
            close(f)
        end

        function CalcVelocity(obj)
            f = waitbar(0,'Initializing');
            for i = 1:size(obj.Traces, 1)
                waitbar(i./size(obj.MSD, 1),f,'Calculating velocity');
                coord = [obj.Traces{i, 1}.row obj.Traces{i, 1}.row];
                dR   = sqrt((coord(1,1)-coord(end,1))^2 + (coord(1,2)-coord(end,2))^2);
                VelocityList(i,1) = dR/10^3/(size(coord,1)*obj.info.ExpTime);
            end
            obj.Velocity = VelocityList;
            close(f)
        end

        function SaveData(obj)
            for i = 1:size(obj.Traces, 1)
                AverageInt(i, 1) = mean(obj.Traces{i, 1}.Intensity);
                AverageArea(i, 1) = mean(obj.Traces{i, 1}.Area);
                AverageEccentricity(i, 1) = mean(obj.Traces{i, 1}.Eccentricity);
                AverageOrientation(i, 1) = mean(obj.Traces{i, 1}.Orientation);
            end
            Results = table(obj.Traces, obj.MSD, obj.Tau, obj.Diff, obj.Velocity, AverageInt, AverageArea, AverageEccentricity, AverageOrientation,...
                'VariableNames', {'Traces', 'MSD', 'TimeLag', 'D', 'v', 'avInt', 'avArea', 'avEcc', 'avOrient'});
            obj.Results = Results;

            Filename = append(obj.file.path, filesep, 'MSDResults.mat');
            save(Filename, "Results");

            disp('MSD and diffusion data is saved')
        end

        function [MSD, D] = calc(obj, Coord)
            dim = size(Coord,2);

            switch dim
                case 1
                    Coord(:,2:3) = 0;
                case 2
                    Coord(:,3) = 0;
                case 3
                otherwise
                    error('unexpected dimension for the vector')
            end
            
            %%%Fist calculate distance between each frame
            
            DX = diff(Coord(:,1).');
            DY = diff(Coord(:,2).');
            DZ = diff(Coord(:,3).');
            D = sqrt(DX.^2 + DY.^2 + DZ.^2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            MSD = zeros(size(Coord,1)-1,1);
            %Calculate mean-squere-displacement
            for i = 1:size(Coord,1)-1
                
                stp = i;
                cnt =  1;
                D1  = [];
                while cnt<=stp && cnt+stp<=size(Coord,1)
                    
                    idx = cnt:stp:size(Coord,1);
                    DX  = diff(Coord(idx,1).');
                    DY  = diff(Coord(idx,2).');
                    DZ  = diff(Coord(idx,3).');
                    D1  = [D1 sqrt(DX.^2 + DY.^2 + DZ.^2)];
                    cnt = cnt+1;
                    
                    if ~isempty(D1)
                        D2=D1(~isnan(D1));
                        if ~isempty(D2)
                            MSD(i) = mean(D2.^2);
                        else
                            MSD(i) = NaN;
                        end
                    end
                end %while
            end
            
            D   = D(:);
            MSD = MSD(:);

        end

        function D = getDiffCoeff(obj, msd,tau,fitRange,dim)
            switch dim
                case '1D'
                    div = 2;
                case '2D'
                    div = 4;
                case '3D'
                    div = 6;
                otherwise
                    error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
            end
        
            assert(min(size(msd))==1,'MSD needs to be provided as a vector')
            
            tofit = msd(1:fitRange);
            tau   = tau(1:fitRange);
            f     = fit(tau(:),tofit(:),'a*x');
            g = coeffvalues(f);
            D = g(1)/div;
        end
    end
end

