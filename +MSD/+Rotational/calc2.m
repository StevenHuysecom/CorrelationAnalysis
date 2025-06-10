function [MSAD, TimeLag] = calc2(coord, tau, expTime)

    dim = size(coord,2);
    
    switch dim
        case 1
            coord(:,2) = 0;
        case 2        
        otherwise
            error('unexpected dimension for the vector')
    end
    
    MSAD = zeros(size(coord,1)-1,1);
    %Calculate mean-square-displacement
    maxIdx = round((max(tau)-min(tau))/expTime);
    AngDisp = [];

    Theta = coord(:,1).';
    Phi = coord(:,2).';

    AngDisp = cell(round((tau(end) - tau(1))./expTime), 1);

    for dt = 1:size(coord,1)-1
        for z = 1:size(Theta, 2)-dt
            TimeLagCurr = round((tau(z+dt) - tau(z))./expTime);

            theta1 = Theta(z);
            phi1 = Phi(z);
            r1 = [cos(theta1)*cos(phi1) cos(theta1)*sin(phi1) sin(theta1)];

            theta2 = Theta(z+dt);
            phi2 = Phi(z+dt);
            r2 = [cos(theta2)*cos(phi2) cos(theta2)*sin(phi2) sin(theta2)];

            AngDisp{TimeLagCurr, 1} = [AngDisp{TimeLagCurr, 1}; real(acos(dot(r1, r2)))];
        end
    end

    MSAD = [];
    TimeLag = [];
    for dt = 1:size(AngDisp,1)
        Dr = AngDisp{dt, 1};
        Dr(isnan(Dr)) = [];
        Dr(Dr == 0) = [];
        if ~isempty(Dr)
            MSAD = [MSAD, mean(Dr.^2)];
            TimeLag = [TimeLag, dt*expTime];
        end
    end
end
