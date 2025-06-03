classdef MPCalibration < Core.Movie
    %This class will hold information about the 2D Cal movie and will be
    %able to calculate the calibration as well.
    
    properties (SetAccess = 'private')
        cal
    end
    
    methods
        function obj = MPCalibration(raw,info)
            %In this case raw path is used for both the loading of the
            %movie and the calibration calculation
            obj = obj@Core.Movie(raw,info);
          
        end
        
        function set.cal(obj,cal)
            assert(isstruct(cal),'2DCal expected to be a struct');
            obj.cal = cal;
                
        end
        
        function calc(obj, nChan)
            
            switch nargin
                case 1
                    nChan = 4;
                case 2
                otherwise
                    error('too many input arguments')
            end
            
            assert(~isempty(obj.raw),'Please provide a path');
            %Calculate from the raw path stored in the movie
            path = obj.raw.movInfo.Path;
            [file2Analyze] = obj.getFileInPath( path, 'calibration.mat');
            
            if (~isempty(file2Analyze))
                %If calibration already exist we load it
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                calibration = tmp.calibration;

            else
                %Otherwise we Calculate it
                path = obj.raw.fullPath;
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(path);
                assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
                assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');

                [calib, inform] = mpSetup.cali.calculate(path,nChan);

                calibration.info = inform;
                calibration.file = calib;
                [folder,~] = fileparts(path);
                filename = [folder filesep 'calibration.mat'];
                calibration.fullPath = filename;
                save(filename,'calibration');
                
            end
            
            obj.cal = calibration;
            disp('Done');
            
        end
        
        function cal = getCal(obj)
            
            cal = obj.cal;
            
        end
        
        function showCal(obj)
            focusMet = obj.cal.file.focusMet;
            fit = obj.cal.file.fit;
            color = rand(8,3);
            ZPos = obj.cal.file.Zpos;
            FocusZ = {obj.cal.file.inFocus.zpos};
            figure()
            hold on
            leg = cell(1,size(focusMet,2));
            height =  max(max(focusMet));
            y = 1:height;
            for i = 1 : size(focusMet,2)
                [~,idx] = max(fit(:,i));
                scatter(ZPos(:),focusMet(:,i),[],color(i,:),'filled')
                plot(ZPos(:),fit(:,i),'Color', color(i,:),'LineWidth',2.5,'HandleVisibility','off')
               
                x = ones(1,length(y))*FocusZ{i};
                plot(x(:),y(:),'k--','HandleVisibility','off');
                
                leg{i} = ['Cam' num2str(obj.cal.file.inFocus(i).cam) ' - Plane' num2str(obj.cal.file.inFocus(i).ch)];
                
            end
            ylim([min(min(focusMet)), max(max(focusMet))]);
            xlim([min(ZPos), max(ZPos)]);
            title('Setup Plane Calibration');
            legend(leg)
            hold off
            
        end
        
    end
end

