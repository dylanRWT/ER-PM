function [best_coefficients,best_jacobian,resnorm,e_flag,solout,user_input] = HUFP_osscilation_fitter_manual(xCull,y_raw)
    % Normalize the data
    y1 = y_raw/max(y_raw);
    coefficients = zeros(5,6);
    solout = zeros(4,width(y1));
    resnorm = zeros(4,1);
    residual = zeros(4,width(y1));
    e_flag = zeros(4,1);
    %
    satisfied = 'n' ; % initialize a loop to take user input for initial model values
    plot(xCull,y1,'r') % plot the graph
    set(gcf, "units", 'inches', 'position', [2 2 7 7])
    drawnow()
    while satisfied == 'n'
        prompt = "What values do you want to use for amplitude, amp decay, period, shift, baseline and baseline decay    " ;
        user_input = input(prompt) ;
        if isa(user_input,'double')
            ossfit1 = @(x)  x(1)*exp(-x(2) * xCull) .* cos((2*pi/user_input(3)) * (xCull + x(3))) + x(4) *exp(xCull * -1 * x(5));
            %                 amp, amp decay, shift, baseline, baseline decay
            initial1 = [user_input(1),user_input(2),user_input(4),user_input(5),user_input(6)]  ;
            initial_out1 = ossfit1(initial1);
            close all
            plot(xCull,y1,'r',xCull,initial_out1,'b')
            set(gcf, "units", 'inches', 'position', [2 2 7 7])
            drawnow
        end
        prompt = "Are you satisfied with the initial fit?   y or n " ; 
        satisfied = input(prompt, 's') ; 
    end
    coefficients(1,:) = user_input;
    close all
    % Lsqnonlin in two steps
    %     amp +- 10%           amp decay +- 100%
    ub = [(user_input(1)*1.1),(user_input(2)*2),(user_input(4)+3),(user_input(5)+2),(user_input(6)*2)] ;
    lb = [(user_input(1)*0.9),(user_input(2)*0.5),(user_input(4)-3),(user_input(5)-2),(user_input(6)*0.5)] ; 
    options =  optimoptions('lsqnonlin','MaxIter',1000,'TolFun', 1e-9, 'MaxFunEvals', 1000 ,'TolX', 1e-9); 
    [solval1,resnorm(1,1),residual(1,:),e_flag(1,1),~,~,jacobian1] = lsqnonlin(@(x) (ossfit1(x) - y1).^2 ,initial1,lb,ub,options) ;
    solout(1,:) = ossfit1(solval1);
        % exit flags 
        % 1 == Function converged to a solution x
        % 2 == Change in x is less than the specified tolerance, or Jacobian at x is undefined.
        % 3 == Change in the residual is less than the specified tolerance.
        % 4 == Relative magnitude of search direction is smaller than the step tolerance.
        % 0 == max iterations or evaluations exceeded change options.MaxFunctionEvaluations or options.MaxIterations
    %% Do a second solve where only period and phase are manipulated
    ossfit2 = @(x)  solval1(1)*exp(-1 * solval1(2) * xCull) .* cos((2*pi/x(1)) * (xCull + x(2))) + solval1(4) *exp(xCull * -1 * solval1(5));
    initial2 = [user_input(3), solval1(3)] ;
    if solval1(3) > 0
        ub = [(user_input(3)*1.1), (solval1(3)*2)]  ;
        lb = [(user_input(3)*0.9), (solval1(3)*0.5)] ;
    else 
        ub = [(user_input(3)*1.1), (solval1(3)*0.5)]  ;
        lb = [(user_input(3)*0.9), (solval1(3)*2)] ;
    end
    [solval2,resnorm(2,1),residual(2,:),e_flag(2),~,~,jacobian2] = lsqnonlin(@(x) (ossfit2(x)-y1).^2, initial2, lb, ub, options) ;
    solout(2,:) = ossfit2(solval2) ; 
    % ensure that you are at a global min by solving with final
    % coefficients twice
    ossfit3 = @(x)  x(1)*exp(-1 * x(2) * xCull) .* cos((2*pi/x(3)) * (xCull + x(4))) + x(5) *exp(xCull * -1 * x(6));
    initial3 = [solval1(1),solval1(2),solval2(1),solval2(2),solval1(4),solval1(5)];
    ub = [0.8,0.03,24,5,1.2,0.02] ;
    lb = [0.15,0.0001,13,-5,0.2,0.0001] ; 
    [solval3,resnorm(3,1),residual(3,:),e_flag(3),~,~,jacobian3] = lsqnonlin(@(x) (ossfit3(x)-y1).^2, initial3, lb, ub, options) ;
    solout(3,:) = ossfit3(solval3);
    initial4 = [solval3(1),solval3(2),solval3(3),solval3(4),solval3(5),solval3(6)];
    [solval4,resnorm(4,1),residual(4,:),e_flag(4),~,~,jacobian4] = lsqnonlin(@(x) (ossfit3(x)-y1).^2, initial4, lb, ub, options) ;
    solout(4,:) = ossfit3(solval4);
    % only return the coefficients of the least squares fit
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
    %
end