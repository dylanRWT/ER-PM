%% HUFP_osscilation_fitter_auto
%
% Created by Dylan Greening for use in RWTurner Lab UCALGARY
%
% A function nested inside of HUFP_runone - > HUFP_ladder_analysis
%
% Fits data with a damped osscilatory model using hands-free (double
%   smoothed period, and linear decay) inital values 
% 
% Input Arg:
% xCull = pixel positions which correspond to distance in the
%   centroid/feature images
% y_raw = input y_value data which corresponds to values a xCull position
%
% Output Arg:
% best_coefficients = the coefficients which generated the best fit after 4
%   iterations of lsqnonlin fitting of the damped model
% best_jacobian = not used at this time
% resnorm = normalized residuals 
% e_flag = the exit flag from lsqnonlin, defines why the solver stopped 
% solout = the y_values of the best fit solution for the damped model
%   fitting
%%
function [best_coefficients,best_jacobian,resnorm,e_flag,solout] = HUFP_osscilation_fitter_auto(xCull,y_raw,initial_period)
    %% Normalize the data
    y1 = y_raw/max(y_raw);
    coefficients = zeros(5,6); % coefficients of the damped oscillatory model
    solout = zeros(4,width(y1));
    resnorm = zeros(4,1);
    residual = zeros(4,width(y1));
    e_flag = zeros(4,1);
    %
    ossfit1 = @(x)  x(1)*exp(-x(2) * xCull) .* cos((2*pi/initial_period) * (xCull + x(3))) + x(4) *exp(xCull * -1 * x(5));
    coefficients = polyfit(xCull,(y1),1); % this determines a good guess for baseline
    %                 amp, amp decay, shift, baseline, baseline decay
    initial1 = [0.3,0.019,-2,coefficients(2),0.0045]  ;
    initial_out1 = ossfit1(initial1);
    % Lsqnonlin in two steps
    ub = [1.4,0.03,10,1,0.01] ;
    lb = [0.1,0.0001,-10,0,0.00001] ; 
    options =  optimoptions('lsqnonlin','MaxIter',1000,'TolFun', 1e-9, 'MaxFunEvals', 1000 ,'TolX', 1e-9); 
    [solval1,resnorm(1,1),residual(1,:),e_flag(1,1),~,~,jacobian1] = lsqnonlin(@(x) (ossfit1(x) - y1).^2 ,initial1,lb,ub,options) ;
    solout(1,:) = ossfit1(solval1);
    %% Do a second solve where only period and phase are manipulated
    ossfit2 = @(x)  solval1(1)*exp(-1 * solval1(2) * xCull) .* cos((2*pi/x(1)) * (xCull + x(2))) + solval1(4) *exp(xCull * -1 * solval1(5));
    initial2 = [initial_period, solval1(3)] ;
    ub = [40, 10] ;
    lb = [13, -10] ;
    [solval2,resnorm(2,1),residual(2,:),e_flag(2),~,~,jacobian2] = lsqnonlin(@(x) (ossfit2(x)-y1).^2, initial2, lb, ub, options) ;
    solout(2,:) = ossfit2(solval2) ; 
    %% ensure that you are at a global min by solving with final coefficients
    ossfit3 = @(x)  x(1)*exp(-1 * x(2) * xCull) .* cos((2*pi/x(3)) * (xCull + x(4))) + x(5) *exp(xCull * -1 * x(6));
    initial3 = [solval1(1),solval1(2),solval2(1),solval2(2),solval1(4),solval1(5)];
    ub = [0.8,0.03,24,5,1.2,0.02] ;
    lb = [0.15,0.0001,13,-5,0.2,0.0001] ; 
    [solval3,resnorm(3,1),residual(3,:),e_flag(3),~,~,jacobian3] = lsqnonlin(@(x) (ossfit3(x)-y1).^2, initial3, lb, ub, options) ;
    solout(3,:) = ossfit3(solval3);
    initial4 = [solval3(1),solval3(2),solval3(3),solval3(4),solval3(5),solval3(6)];
    [solval4,resnorm(4,1),residual(4,:),e_flag(4),~,~,jacobian4] = lsqnonlin(@(x) (ossfit3(x)-y1).^2, initial4, lb, ub, options) ;
    solout(4,:) = ossfit3(solval4);
    %% only return the coefficients of the least squares fit
    if (resnorm(1) < resnorm(2) && resnorm(1) < resnorm(3)) && resnorm(1) < resnorm(4)
        best_jacobian = jacobian1;
        best_coefficients = solval1;
    elseif (resnorm(2) < resnorm(1) && resnorm(2) < resnorm(3)) && resnorm(2) < resnorm(4)
        best_jacobian = jacobian2;
        best_coefficients = [solval1(1),solval1(2),solval2(1),solval2(2),solval1(4),solval1(5)];
    elseif (resnorm(3) < resnorm(1) && resnorm(3) < resnorm(2)) && resnorm(3) < resnorm(4)
        best_jacobian = jacobian3;
        best_coefficients = solval3;
    else 
        best_jacobian = jacobian4;
        best_coefficients = solval4;
    end
    %%
end