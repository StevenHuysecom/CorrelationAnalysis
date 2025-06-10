function [MSD] = calc(coord, tau, expTime)

    dim = size(coord,2);
    
    switch dim
        case 1
            coord(:,2) = 0;
        case 2        
        otherwise
            error('unexpected dimension for the vector')
    end
    
    % MSAD = zeros(size(coord,1)-1,1);
    %Calculate mean-square-displacement
    maxIdx = round((max(tau)-min(tau))/expTime);
    for dt = 1:size(coord,1)-1
        
        cnt =  1;
        AngDisp = [];
        while cnt<=dt && cnt+dt<=size(coord,1)
            
            idx = cnt:dt:size(coord,1);
            
            SelectTheta = coord(idx, 1);
            SelectPhi = coord(idx, 2);
            for i = 1:size(SelectTheta, 1) - 1
                Theta1 = SelectTheta(i, 1);
                Theta2 = SelectTheta(i+1, 1);
                Phi1 = SelectPhi(i, 1);
                Phi2 = SelectPhi(i+1, 1);

                r1 = [sin(Theta1)*cos(Phi1), sin(Theta1)*sin(Phi1), cos(Theta1)];
                r2 = [sin(Theta2)*cos(Phi2), sin(Theta2)*sin(Phi2), cos(Theta2)];
                AngDisp = [AngDisp; acos(dot(r1, r2))];
            end

            cnt = cnt + 1;
            AngDisp(AngDisp == 0) = [];
        end

        if ~isempty(AngDisp)
            Dr = AngDisp(~isnan(AngDisp));
            
            if ~isempty(Dr)
                MSAD(dt) = mean(Dr.^2);
            else
                MSAD(dt) = NaN;
            end
        else
            MSAD(dt) = NaN;
        end
        TimeLag(dt) = dt*expTime;
    end

    MSD = [MSAD; TimeLag];
end