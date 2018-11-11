function conc = cfunction(time,selector)
%cfunction defines the concentration input function used in analysis
%   INPUTS: time - the time period to be considered in seconds.
%   selector - an integer between 1 - 7 to specify the type of function.
%   OUTPUT: the concetration as a function of time, excpet if selector ==
%   3, the output is the intergral of the function
switch selector
    case 1
        %Biexponential Pk model
        a1 = 0.833;
        a2 = 0.336;
        t1 = 0.171;
        t2 = 0.364;
        sigma1 = 0.055;
        sigma2 = 0.134;
        alpha = 1.064;
        beta = 0.166;
        s = 37.772;
        tau = 0.482;
        tMin = time/60;
        r1 = (a1 / (sigma1*sqrt(2*pi)))*exp(-(tMin-t1).^2/(2*sigma1^2));
        r2 = (a2 / (sigma2*sqrt(2*pi)))*exp(-(tMin-t2).^2/(2*sigma2^2));
        wo = alpha*exp(-beta*tMin)/(1+exp(-s*(tMin-tau)));
        conc = r1 + r2 + wo;
    case 2
        %constant function, with an integral == 1
        conc = (1/numel(time));
        
    case 3
        %integral of case 1
        fun = @(tp) (2252549397762555*exp(-(144115188075855872*(tp/60 - 91/250).^2)/5175464634180137))/2251799813685248 + (1700717457231141*exp(-(20000*(tp/60 - 171/1000).^2)/121))/281474976710656 + (133*exp(-(83*tp)/30000))/(125*(exp(2275763/125000 - (9443*tp)/15000) + 1));
        conc = integral(fun,0,1000);
        
    case 4
        %heaviside function 
        conc = heaviside(time);
    case 5
        %heaviside function combined with a sign function. small offset
        %needed for analysis 
        conc = heaviside(time).*sin(time)+ 1e-10;
    case 6
        %exponential decay Pk model
        conc = (0.81*exp(-0.081*t)+0.34)./(1+4*exp(0.99*t));
    case 7
       %rectangular pulse with width 2. small offset needed for analysis
        conc = rectpuls(time-1.5,2) + 1e-10;
end

