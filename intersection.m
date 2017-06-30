%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           
%                               AEM 413 
%                             Problem 4.9
%                          Kellen Schroeter
%
%              From Modern Compressible Flow: Third Edition
%                          John D. Aanderson
%
%                 Ref for analytical soln. for Beta
%            <http://arc.aiaa.org/doi/abs/10.2514/2.2349>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------- Input Values -----------------

mach1 = 3.0;             % Initial Mach number (greater than 1)
pres1 = 1.0;             % Initial Pressure
theta2 = 20;             % Lower Deflection Angle
theta3 = 15;             % Upper Deflection Angle
phi = .1;                % Initial Guess for Phi
phiIncrement = 0.00001;  % Increment for each iteration
presConvCrit = .0001;    % Convergence Criterion for Pressure
n = 0;                   % 0 for Weak Shock, 1 for Strong Shock  
gamma = 1.4;             % Ratio of Specific Heats


%---------------- Loop Setup -----------------

theta2 = theta2 * pi / 180;       %Convert Angles to Radians
theta3 = theta3 * pi / 180; 
phi = phi * pi / 180;

presConverge = 1;                 %Initialize Pressure Convergence

while presConverge > presConvCrit
    
    theta4 = phi + theta3;        %Flow directions across refracted shocks
    theta4Prime = theta2 - phi; 
    
    %---Solve From Bottom
         %Analytical Solution for Beta (L. Rudd & M.J. Lewis 1998)
    mu=asin(1/mach1);             % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta2);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta2);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    beta2=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)));
    
    Mn1 = mach1*sin(beta2);       % Shock Relations 
    presRatio1 = 1 + (2*gamma)/(gamma+1)*(Mn1^2 - 1);
    Mn2 = sqrt((1 + Mn1^2*(gamma-1)/2)/(gamma*Mn1^2 - (gamma-1)/2));
    mach2 = Mn2/sin(beta2 - theta2);
          
          %Analytical Solution for Beta (L. Rudd & M.J. Lewis 1998)
    mu=asin(1/mach2);             % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta4Prime);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta4Prime);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    beta4Prime=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)));
    
    Mn2p = mach2*sin(beta4Prime); % Shock Relations 
    presRatio2 = 1 + (2*gamma)/(gamma+1)*(Mn2p^2 - 1);
    Mn4Prime = sqrt((1 + Mn2p^2*(gamma-1)/2)/(gamma*Mn2p^2 - (gamma-1)/2));
    mach4Prime = Mn4Prime/sin(beta4Prime - theta4Prime);
    
    pres4Prime = pres1 * presRatio1 * presRatio2;
    
    %---Solve From Top
         %Analytical Solution for Beta (L. Rudd & M.J. Lewis 1998)
    mu=asin(1/mach1);             % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta3);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta3);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    beta3=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)));
    
    Mn1 = mach1*sin(beta3);       % Shock Relations 
    presRatio1 = 1 + (2*gamma)/(gamma+1)*(Mn1^2 - 1);
    Mn3 = sqrt((1 + Mn1^2*(gamma-1)/2)/(gamma*Mn1^2 - (gamma-1)/2));
    mach3 = Mn3/sin(beta3 - theta3);
          
          %Analytical Solution for Beta (L. Rudd & M.J. Lewis 1998)
    mu=asin(1/mach3);             % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta4);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta4);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    beta4=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)));
    
    Mn3p = mach3*sin(beta4);      % Shock Relations 
    presRatio3 = 1 + (2*gamma)/(gamma+1)*(Mn3p^2 - 1);
    Mn4 = sqrt((1 + Mn3p^2*(gamma-1)/2)/(gamma*Mn3p^2 - (gamma-1)/2));
    mach4 = Mn4/sin(beta4 - theta4);
    
    pres4 = pres1 * presRatio1 * presRatio3;
    
    %---Compare
    presConverge = abs(pres4 - pres4Prime); 
    phi = phi + phiIncrement;
    
end

fprintf('Pressure4 = %.5f atm\n\n' , pres4)
fprintf('Pressure4prime = %.5f atm\n\n' , pres4Prime)
fprintf('Mach4 = %.5f\n\n' , mach4)
fprintf('Mach4prime = %.5f\n\n' , mach4Prime)
fprintf('Phi = %.5f deg\n\n' , phi*180/pi)

