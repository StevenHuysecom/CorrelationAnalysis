function Dr = getDiffCoeff(msad,tau,fitRange,dim)
    
    switch dim
        case '1D'
            error('Cannot do rotational tracking in 1D')
        case '2D'
            div = 1;
        case '3D'
            div = 2;
        otherwise
            error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
    end

    try
        % [f, gov] = fit(msad(2,:)',msad(1,:)','a*(1-(1-b^2)*exp((-1.6*c*x)^0.95))');
        [f, gov] = fit(tau(1,1:fitRange)',msad(1,1:fitRange)','a*x');
        % figure()
        % plot(f, tau(1,:)', msad(1,:)');

        if gov.rsquare > 0.03
            g = coeffvalues(f);
            Dr = g(1)/div;
        else
            Dr = NaN;
        end
    catch
        Dr = NaN;
    end
end