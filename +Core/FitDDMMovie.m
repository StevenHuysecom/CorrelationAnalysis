classdef FitDDMMovie < handle
%% Class to fit DDM output data
    
    properties (SetAccess ='public')
    FitResults
    DDMDataFit
    Model
    end
    
methods   
    
        function obj = FitDDMMovie(DataInput,Model)
            for c = 1:size(DataInput, 1)
                obj.DDMDataFit{c,1} = DataInput{c,1};
                obj.FitResults{c,1} = [];
                obj.Model{c,1} = Model;
            end
        end
        
        function [FitResults] = FitDDMMovieDancingCells(obj, DDMResults, DDMInfo, Path)
            DDMResults = DDMResults{1,1};
            f = waitbar(0, 'Intializing');
            Starts = find(DDMResults.QVector > DDMInfo.Qmin, 1, 'first');
            End = find(DDMResults.QVector < DDMInfo.Qmax & DDMResults.QVector ~= 0, 1, 'last');
            mkdir(append(Path, filesep, 'FiguresISFwithFit'));
            for i = Starts:End
                % try
                    waitbar((i-Starts+1)./(End - Starts+1),f,'Fitting');
                    q = DDMResults.QVector(i);
                    model = append('c*(1-((a+1)/(a*', num2str(q), '*b*x))*((sin(a*atan((', num2str(q), '*b*x)/(a+1))))/(1+((', num2str(q), '*b*x)/(a+1))^(a/2)))) + d');
                    Lower = [0, 0 , DDMResults.DDMSignalValue{i}(end) - DDMResults.DDMSignalValue{i}(1), 0];
                    Upper = [10^99, 10^99 , 10^99, DDMResults.DDMSignalValue{i}(1)];
                    Start = [25, 0.1, DDMResults.DDMSignalValue{i}(end) - DDMResults.DDMSignalValue{i}(1), DDMResults.DDMSignalValue{i}(1)];
                    [ff, gov] = fit(DDMResults.Time{i}', DDMResults.DDMSignalValue{i}', model, 'StartPoint', Start, 'Lower', Lower, 'Upper', Upper);
                    
                    if gov.adjrsquare > 0.99
                        coef = coeffvalues(ff);
                        v(i) = coef(2);
                        s(i) = (v(i))*(coef(1) + 1)^(-1/2);
                        A(i) = coef(3);
                        B(i) = coef(4);
                    else                            
                        v(i) = NaN;
                        s(i) = NaN;
                        A(i) = NaN;
                        B(i) = NaN;
                    end

                    Fig = figure();
                    plot(ff, DDMResults.Time{i}, DDMResults.DDMSignalValue{i})
                    title(append('adj r^2 = ', num2str(gov.adjrsquare)))

                    filename = append(Path, filesep, 'FiguresISFwithFit', filesep, 'ISFfit_q', num2str(q), '.png');
                    saveas(Fig, filename);
                    close(Fig)

                % catch
                %     v(i) = NaN;
                %     s(i) = NaN;
                %     A(i) = NaN;
                %     B(i) = NaN;
                % end
            end
            close(f)

            FitResults = [v', s'];
        end
        
        function FitDDMtoMSD(obj, DDMInfo, file, info)
            MSDfile = append(file.path, filesep, 'MSD.mat');
            if ~exist(MSDfile)
                info.runMethod  = 'run';
            end
            
                %if strcmp(info.runMethod, 'run')
                DDMFitOutput = table({},{},'VariableNames',{'Fit','CoeffVals'});
                StartIdx = find(obj.DDMDataFit.QVector >= DDMInfo.Qmin, 1); 
                EndIdx = find(obj.DDMDataFit.QVector <= DDMInfo.Qmax, 1, "last");
                MSD = [];
                QMatrix = [];
    
                for i = StartIdx:EndIdx
                    Idx = i;
                    time = obj.DDMDataFit.Time(Idx);
                    time = time{1,1};
                    DDMSignal = obj.DDMDataFit.DDMSignalValue(Idx);
                    DDMSignal = DDMSignal{1,1};
                    %%% Initial guess A
                    [pks, ~, w, prom] = findpeaks(DDMSignal);
                    MaxProm = max(prom);
                    IsProm = prom > MaxProm./2;
                    pks(IsProm == 0) = [];
                    A = mean(pks);
                    %%% Initial guess B
                    B = 0;
                    Q(i) = obj.DDMDataFit.QVector(Idx);


                    DDMSignal = (DDMSignal-B)./A;
                    
                    MSDRow = (-4./(Q(i)^2))*log(1-DDMSignal(1,:));
                    MSDRow = real(MSDRow);

                    if isempty(MSD)
                        MSD = MSDRow;
                    else
                        if size(MSD, 2) < size(MSDRow, 2)
                            ToAdd = NaN(size(MSD, 1), (size(MSDRow, 2) - size(MSD, 2)));
                            MSD = [MSD, ToAdd];
                            MSD = [MSD; MSDRow];
                        elseif size(MSD, 2) > size(MSDRow, 2)
                            ToAdd = NaN(1, (size(MSD, 2) - size(MSDRow, 2)));
                            MSDRow = [MSDRow, ToAdd];
                            MSD = [MSD; MSDRow];
                        else
                            MSD = [MSD; MSDRow];
                        end
                    end   
                end
                
                MSDMatrix = MSD;
                filename = append(file.path, filesep, 'MSDMatrix.mat');
                save(filename, 'MSDMatrix');
    
                MSDAv = nanmean(MSD, 1);
                time = [0, time];
                MSDAv(2,:) = time(1:size(MSDAv,2));
    
                filename = append(file.path, filesep, 'MSD.mat');
                save(filename, 'MSD') 
            % else
            %     disp('Found MSD file - loading it');
            %     MSD = load(MSDfile);
            %     MSD = MSD.MSD;
            %     disp('Found MSD file - Done');
            % end

            obj.FitResults.MSDResults = MSD;
        end

        function MSD_optimization(obj, DDMInfo, file, info, maxIter)
            for c = 1:size(obj.DDMDataFit, 1)
                for i = 1:size(obj.DDMDataFit{c,1},1)-1
                    D(i,:) = obj.DDMDataFit{c,1}.DDMSignalValue{i,1};
                    q(i) = obj.DDMDataFit{c,1}.QVector(i);
                end
                t = obj.DDMDataFit{c,1}.Time{1,1} ;
                q1 = DDMInfo.Qmin(c);
                q2 = DDMInfo.Qmax(c);
    
                q_idx = (q >= q1) & (q <= q2);
                q_sel = q(q_idx);
                D_sel = D(q_idx, :);
                
                A0 = max(D_sel, [], 2); % Initial guess for A(q)
                B0 = zeros(size(A0)); %min(D_sel, [], 2); % Initial guess for B(q)
                
                params0 = [A0; B0];
                
                options = optimset('MaxIter', maxIter, 'Display', 'iter');
                [params_opt, fval, exitflag] = fminsearch(@cost_function, params0, options);
    
                % Extract optimal parameters
                n_q = length(q_sel);
                A_opt = params_opt(1:n_q);
                B_opt = params_opt(n_q+1:end);
                
                obj.FitResults{c,1}.ParamA = A_opt;
                obj.FitResults{c,1}.ParamB = B_opt;
    
                [~, MSD_final] = cost_function(params_opt);
    
                obj.FitResults{c,1}.MSD = MSD_final;
                obj.FitResults{c,1}.MSDstddev = fval;
                obj.FitResults{c,1}.MSDExit = exitflag;
    
                Filename = append(file.path, filesep, 'MSD_', num2str(c), '.mat');
                save(Filename, 'MSD_final');
            end


            function [sigma2, MSD_avg] = cost_function(params)
                % COST_FUNCTION: Computes the dispersion \sigma^2 and MSD for given parameters.
                
                n_q = length(q_sel);
                A = params(1:n_q); % Amplitude A(q)
                B = params(n_q+1:end); % Noise baseline B(q)
                
                % Step 3: Calculate MSD(t|q)
                MSD_tq = zeros(n_q, length(t));
                for i = 1:n_q
                    q_i = q_sel(i);
                    % Ensure no invalid log values by checking (D(q,t) - B(q)) / A(q) < 1
                    valid_idx = (D_sel(i, :) - B(i)) < A(i);
                    MSD_tq(i, valid_idx) = (-4 ./ (q_i^2)) .* log(1 - (D_sel(i, valid_idx) - B(i)) ./ A(i));
                    % For invalid values, set to NaN to exclude from averages
                    MSD_tq(i, ~valid_idx) = NaN;
                end
                
                % Step 4: Determine subset J(t) and calculate N(t)
                J_t = MSD_tq < (4 ./ q_sel.^2)'; %< (4 ./ q_sel.^2)';
                N_t = sum(J_t, 1);
                
                % Step 5: Calculate MSD(t)
                MSD_avg = zeros(1, length(t));
                for j = 1:length(t)
                    valid_idx = J_t(:, j);
                    if N_t(j) > 0
                        MSD_avg(1,j) = mean(MSD_tq(valid_idx, j), 'omitnan');
                    else
                        MSD_avg(j) = NaN; % No valid q values for this time
                    end
                end
                
                % Step 6: Calculate dispersion \sigma^2_t and \sigma^2
                sigma2_t = zeros(1, length(t));
                for j = 1:length(t)
                    valid_idx = J_t(:, j);
                    if N_t(j) > 1
                        log_ratios = log(MSD_tq(valid_idx, j) ./ MSD_avg(j));
                        sigma2_t(j) = sum(log_ratios.^2) / (N_t(j) - 1);
                    else
                        sigma2_t(j) = NaN;
                    end
                end
                
                sigma2 = nansum(sigma2_t); % Total dispersion
            end
                    
        end      
              
        function [D, a, Visc] = CalcDFromMSDAnom(obj, DDMInfo, file, model, StartValues)
            MSD = obj.FitResults.MSDResults;
            ft = fittype(model);
            options = fitoptions(ft);
            options.Lower = [0 0];
            options.Upper = [StartValues(1)*100 1];
            try
                f = fit(MSD(2,1:round(size(MSD,2)./2)).', MSD(1,1:round(size(MSD,2)./2)).', model, 'StartPoint',StartValues);
                
                MSDFig = figure()
                plot(f, MSD(2,1:round(size(MSD,2)./2)).', MSD(1,1:round(size(MSD,2)./2)).');
                xlabel('Time (s)')
                ylabel('MSD (µm^2)')
    
                filename = append(file.path, filesep, 'MSDfigAnom.png');
                saveas(MSDFig, filename)
                
                coef = coeffvalues(f);
                D = coef(1);
                a = coef(2);
                Visc = (((1.380649*10^(-23))*DDMInfo.Temp)./(6*pi*DDMInfo.ParticleSize*10^(-6)*D*10^(-12)))*10^3;

                disp(sprintf('Linear diffusion coefficient is %f µm^2s^-^1', D));
                %disp(sprintf('Anomalous exponent is %f', a));
                disp(sprintf('Viscosity is %f cP', Visc));
            catch
                disp('Fitting raised error')
                D = NaN;
                a = NaN;
                Visc = NaN;
            end
        end

        function [D, a, Visc] = CalcDFromMSDAnomLog(MSD)
            MSD = obj.FitResults.MSDResults;
            ft = fittype(model);
            options = fitoptions(ft);
            options.Lower = [0 0];
            options.Upper = [StartValues(1)*100 1];
            try
                f = fit(MSD(2,1:round(size(MSD,2)./2)).', MSD(1,1:round(size(MSD,2)./2)).', model, 'StartPoint',StartValues);
                
                MSDFig = figure()
                plot(f, MSD(2,1:round(size(MSD,2)./2)).', MSD(1,1:round(size(MSD,2)./2)).');
                xlabel('Time (s)')
                ylabel('MSD (µm^2)')
    
                filename = append(file.path, filesep, 'MSDfigAnom.png');
                saveas(MSDFig, filename)
                
                coef = coeffvalues(f);
                D = coef(1);
                a = coef(2);
                Visc = (((1.380649*10^(-23))*DDMInfo.Temp)./(6*pi*DDMInfo.ParticleSize*10^(-6)*D*10^(-12)))*10^3;

                disp(sprintf('Linear diffusion coefficient is %f µm^2s^-^1', D));
                %disp(sprintf('Anomalous exponent is %f', a));
                disp(sprintf('Viscosity is %f cP', Visc));
            catch
                disp('Fitting raised error')
                D = NaN;
                a = NaN;
                Visc = NaN;
            end
        end

        function CalcDFromMSDLin(obj, DDMInfo, file, model)
            MSDMatrix = load(append(file.path, filesep, 'MSDMatrix.mat'));
            MSDMatrix = MSDMatrix.MSDMatrix;
            Time = obj.DDMDataFit.Time{1,1};

            ResultMatrix = [];
            
            for i = 1:size(MSDMatrix, 1)
                MSD = MSDMatrix(i,:);
                ft = fittype(model);
                try
                    f = fit(Time(1,1:4).', MSD(1,1:4).', ft);
                    
                    % MSDFig = figure()
                    % plot(f, MSD(2,1:round(size(MSD,2)./2)).', MSD(1,1:round(size(MSD,2)./2)).');
                    % xlabel('Time (s)')
                    % ylabel('MSD (µm^2)')
                    % 
                    % filename = append(file.path, filesep, 'MSDfigLin.png');
                    % saveas(MSDFig, filename)
                    
                    coef = coeffvalues(f);
                    D = coef(1);
                    Visc = (((1.380649*10^(-23))*DDMInfo.Temp)./(6*pi*DDMInfo.ParticleSize*10^(-6)*D*10^(-12)))*10^3;
                    %A = coef(1);
    
                    disp(sprintf('Linear diffusion coefficient is %f µm^2s^-^1 - calculated from first 4 datapoints', D));
                    %disp(sprintf('Anomalous exponent is %f', a));
                    disp(sprintf('Viscosity is %f cP - calculated from first 4 datapoints', Visc));
                catch
                    disp('Fitting raised error')
                    D = NaN;
                    Visc = NaN;
                end

                ResultMatrix(i,1) = obj.DDMDataFit.QVector(i);
                ResultMatrix(i,2) = D;
                ResultMatrix(i,3) = 1;
                ResultMatrix(i,4) = Visc;
            end

            ResultMatrix = table(ResultMatrix(:,1), ResultMatrix(:,2), ResultMatrix(:,3), ResultMatrix(:,4),...
                'VariableNames', {'QVector', 'Diffusion (µm^2/s)', 'Anom coeff', 'Visc (cP)'});

            Filename = append(file.path, filesep, 'ResultMatrix.mat');
            save(Filename, 'ResultMatrix');
        end

        function ShowFitResult(obj)
            colors = jet(size(obj.DDMDataFit,1));
            AllParametersInFunctionOFQ = [];
             for i=1:size(obj.DDMDataFit,1)
                try
                time =obj.FitResults.Time(i);
                DDMSignal = obj.FitResults.DDMSignalValue(i);
                f = obj.FitResults.Fit(i);
                CFVals = obj.FitResults.CoeffVals(i);
                QVector = obj.FitResults.QVector(i);
                
                time = time{1,1};
                DDMSignal = DDMSignal{1,1};
                f = f{1,1};
                CFVals = CFVals{1,1};
                if ~ischar(CFVals)
                     AllParametersInFunctionOFQ = [ AllParametersInFunctionOFQ ; [QVector , CFVals]];
                     NumPlots =numel(CFVals);
                end
                subplot(2,NumPlots,[1 NumPlots])
                xAxis = time(1):(time(end)/(5*length(time))):time(end);
                
                scatter(time,DDMSignal,5,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
                hold on
                plot(xAxis,f(xAxis).*max(DDMSignal),'Color',colors(i,:));
                catch
                     
                end

             end 
             
            set(gca,'LineWidth',2)
            set(gca,'FontSize',24)
            box on
            set(gca,'YScale','log')
            for i=1:NumPlots
                subplot(2,NumPlots, NumPlots+i)
                hold on
                box on
                scatter(AllParametersInFunctionOFQ(:,1),AllParametersInFunctionOFQ(:,i+1),50,'MarkerFaceColor','k','MarkerEdgeColor','k')
                set(gca,'LineWidth',2)
                set(gca,'FontSize',24)
                xlabel('Q-Vector')
                ylabel(['Parameter ' 'c_{' num2str(i) '}'])
            end
        end
        
        function ShowParameterResult(obj)

        end
        
        function [D] = CalcD(obj, DDMInfo)
            K = [];
            Q = [];
            for i = 1:height(obj.FitResults)
                if strcmp(obj.FitResults.CoeffVals(i),'Could Not Compute')
                    % K(i,1) = 0;
                    % Q(i,1) = obj.FitResults.QVector(i);
                else
                    Coeff = cell2mat(obj.FitResults.CoeffVals(i));
                    K(end+1,1) = Coeff(2);
                    Q(end+1,1) = obj.FitResults.QVector(i);
                end
            end
            QMinIdx = find(Q(:,1)>DDMInfo.Qmin,1,'first');
            QMaxIdx = find(Q(:,1)>DDMInfo.Qmax,1,'first') - 1;
            if isempty(QMaxIdx)
                QMaxIdx = length(K)
            end
            Q = log(Q(QMinIdx:QMaxIdx));
            K = log(K(QMinIdx:QMaxIdx));
            %only choose points that are close to slope -2
            dx = diff(Q);
            dy = diff(K);
            slope = dy./dx;
            tolerance = 1.00;
            indices = find(abs(slope + 2) < tolerance);

            start_idx = indices(1);
            contiguous_indices = start_idx;
            for i = 2:length(indices)
                if indices(i) == indices(i-1) + 1
                    contiguous_indices = [contiguous_indices, indices(i)];
                else
                    break;  % Stop at the first gap
                end
            end

            Q_filtered = Q(contiguous_indices);
            K_filtered = K(contiguous_indices);

            f = fit(Q_filtered, K_filtered, '(-2)*x+a')
            %f = fit(log(Q(QMinIdx:QMaxIdx,1)), log(K(QMinIdx:QMaxIdx,1)), '(-2)*x+a');
            %f = fit(log(Q), log(K), '(-2)*x+a');

            close all
            figure()
            plot(Q, K)
            hold on
            plot(f)

            g = coeffvalues(f);
            D = exp(-g(1));

        end
    
end  
end