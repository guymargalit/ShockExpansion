function [c_l,c_d] = ShockExpansion(M,alpha,epsilon)
clearvars -except M alpha epsilon

% Oblique shock relations calculations adapted from:
% https://www.mathworks.com/matlabcentral/fileexchange/28242-oblique-shock-relations-solver

%% Define constants %%
gamma = 1.4;
if(epsilon == 0)
    %% Calculate nu %%
    nu_inf = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M^2 - 1))) - atand(sqrt(M^2 - 1));
    nu = nu_inf + alpha;

    %% Solve for mach, T1/T_inf, and p1/p_inf %%
    syms x
    M1 = abs(vpasolve(nu == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T1 = (1+((1.4-1)/2)*(M^2))/(1+((1.4-1)/2)*M1^2);
    p1 = double((T1)^((1.4)/(1.4-1)));

    %% Solve for beta with oblique shock relations %%
    theta = alpha;
    LHS = tan(theta*pi/180);
    beta = 1;
    RHS = inline('2*cot(beta*pi/180)*((M*sin(beta*pi/180))^2-1)/(M^2*(1.4+cos(2*beta*pi/180))+2)','beta','M');
    counter = 0;
    while (RHS(beta,M) < LHS*0.98 | RHS(beta,M) > LHS*1.02)
        if (beta > 90 | beta < 0 | counter >200)
            break;
        end
        if (RHS(beta,M) < LHS*0.70)
            beta=beta+1;
        elseif(RHS(beta,M) < LHS*0.90)
            beta=beta+0.5;
        elseif(RHS(beta,M) < LHS*0.95)
            beta=beta+0.3;
        elseif(RHS(beta,M) < LHS*0.98)
            beta=beta+0.1;
        end
        
        if(RHS(beta,M) > LHS*1.03)
            beta=beta-0.1;
        elseif(RHS(beta,M) > LHS*1.10)
            beta = beta-0.5;
        end
        counter = counter+1;   
    end
    
    counter = 0;
    
    while (RHS(beta,M) < LHS*0.99 | RHS(beta,M) > LHS*1.01)
        if (beta > 90 | beta < 0 | counter > 100)
            break;
        end
        if (RHS(beta,M) < LHS*0.98)
            beta=beta+0.05;
        elseif(RHS(beta,M) < LHS*0.985)
            beta=beta+0.03;
        elseif(RHS(beta,M) < LHS*0.99)
            beta=beta+0.01;
        end
        
        if(RHS(beta,M) > LHS*1.01)
            beta=beta-0.01;
        elseif(RHS(beta,M) > LHS*1.02)
            beta = beta-0.1;
        end
        counter = counter+1;
    end
    
    counter = 0;
    
    while (RHS(beta,M) ~= LHS)
        if (beta > 90 | beta < 0 | counter > 40)
            break;
        end
        
        if (RHS(beta,M) < LHS)
            beta=beta+0.005;
        end
        if(RHS(beta,M) > LHS)
            beta=beta-0.005;
        end
        counter = counter+1;
    end
    
    if (beta < 0)
        fprintf('Beta not in bounds');
    end

    %% Solve for mach and p2/p_inf %%
    M2 = M*sind(beta);
    p2 = 1 + ((2*gamma)/(gamma+1))*(M2^2 - 1);

    %% Calculate c_l and c_d %%
    c_l = (2/(gamma*(M^2)))*(p2 - p1)*cosd(alpha);
    c_d = (2/(gamma*(M^2)))*(p2 - p1)*sind(alpha);
    
elseif(alpha >= 0 && alpha < epsilon)
    %% Find p1/p_inf with oblique shock relations %%
    theta = epsilon-alpha;
    LHS = tan(theta*pi/180);
    beta = 1;
    RHS = inline('2*cot(beta*pi/180)*((M*sin(beta*pi/180))^2-1)/(M^2*(1.4+cos(2*beta*pi/180))+2)','beta','M');
    counter = 0;
    while (RHS(beta,M) < LHS*0.98 | RHS(beta,M) > LHS*1.02)
        if (beta > 90 | beta < 0 | counter >200)
            break;
        end
        if (RHS(beta,M) < LHS*0.70)
            beta=beta+1;
        elseif(RHS(beta,M) < LHS*0.90)
            beta=beta+0.5;
        elseif(RHS(beta,M) < LHS*0.95)
            beta=beta+0.3;
        elseif(RHS(beta,M) < LHS*0.98)
            beta=beta+0.1;
        end
        
        if(RHS(beta,M) > LHS*1.03)
            beta=beta-0.1;
        elseif(RHS(beta,M) > LHS*1.10)
            beta = beta-0.5;
        end
        counter = counter+1;   
    end
    
    counter = 0;
    
    while (RHS(beta,M) < LHS*0.99 | RHS(beta,M) > LHS*1.01)
        if (beta > 90 | beta < 0 | counter > 100)
            break;
        end
        if (RHS(beta,M) < LHS*0.98)
            beta=beta+0.05;
        elseif(RHS(beta,M) < LHS*0.985)
            beta=beta+0.03;
        elseif(RHS(beta,M) < LHS*0.99)
            beta=beta+0.01;
        end
        
        if(RHS(beta,M) > LHS*1.01)
            beta=beta-0.01;
        elseif(RHS(beta,M) > LHS*1.02)
            beta = beta-0.1;
        end
        counter = counter+1;
    end
    
    counter = 0;
    
    while (RHS(beta,M) ~= LHS)
        if (beta > 90 | beta < 0 | counter > 40)
            break;
        end
        
        if (RHS(beta,M) < LHS)
            beta=beta+0.005;
        end
        if(RHS(beta,M) > LHS)
            beta=beta-0.005;
        end
        counter = counter+1;
    end
    
    if (beta < 0)
        fprintf('Beta not in bounds');
    end
    Mn_inf = M*sind(beta);
    p1 = double(1 + ((2*gamma)/(gamma+1))*(Mn_inf^2 - 1));
    Mn1 = sqrt((1+((gamma-1)/2)*Mn_inf^2)/(gamma*Mn_inf^2 - (gamma-1)/2));
    M1 = Mn1/(sind(beta-(theta)));
    nu1 = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M1^2 - 1))) - atand(sqrt(M1^2 - 1));
    
    %% Find p2/p_inf with expansion fan %%
    nu2 = nu1 + (2*epsilon);
    syms x
    M2 = abs(vpasolve(nu2 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T2 = (1+((1.4-1)/2)*(M1^2))/(1+((1.4-1)/2)*M2^2);
    p2p1 = double((T2)^((1.4)/(1.4-1)));
    p2 = (p2p1)*p1;
    
    %% Find p3/p_inf with oblique shock relations %%
    theta = epsilon + alpha;
    LHS = tan(theta*pi/180);
    beta = 1;
    RHS = inline('2*cot(beta*pi/180)*((M*sin(beta*pi/180))^2-1)/(M^2*(1.4+cos(2*beta*pi/180))+2)','beta','M');
    counter = 0;
    while (RHS(beta,M) < LHS*0.98 | RHS(beta,M) > LHS*1.02)
        if (beta > 90 | beta < 0 | counter >200)
            break;
        end
        if (RHS(beta,M) < LHS*0.70)
            beta=beta+1;
        elseif(RHS(beta,M) < LHS*0.90)
            beta=beta+0.5;
        elseif(RHS(beta,M) < LHS*0.95)
            beta=beta+0.3;
        elseif(RHS(beta,M) < LHS*0.98)
            beta=beta+0.1;
        end
        
        if(RHS(beta,M) > LHS*1.03)
            beta=beta-0.1;
        elseif(RHS(beta,M) > LHS*1.10)
            beta = beta-0.5;
        end
        counter = counter+1;   
    end
    
    counter = 0;
    
    while (RHS(beta,M) < LHS*0.99 | RHS(beta,M) > LHS*1.01)
        if (beta > 90 | beta < 0 | counter > 100)
            break;
        end
        if (RHS(beta,M) < LHS*0.98)
            beta=beta+0.05;
        elseif(RHS(beta,M) < LHS*0.985)
            beta=beta+0.03;
        elseif(RHS(beta,M) < LHS*0.99)
            beta=beta+0.01;
        end
        
        if(RHS(beta,M) > LHS*1.01)
            beta=beta-0.01;
        elseif(RHS(beta,M) > LHS*1.02)
            beta = beta-0.1;
        end
        counter = counter+1;
    end
    
    counter = 0;
    
    while (RHS(beta,M) ~= LHS)
        if (beta > 90 | beta < 0 | counter > 40)
            break;
        end
        
        if (RHS(beta,M) < LHS)
            beta=beta+0.005;
        end
        if(RHS(beta,M) > LHS)
            beta=beta-0.005;
        end
        counter = counter+1;
    end
    
    if (beta < 0)
        fprintf('Beta not in bounds');
    end
    Mn_inf = M*sind(beta);
    p3 = double(1 + ((2*gamma)/(gamma+1))*(Mn_inf^2 - 1));
    Mn3 = sqrt((1+((gamma-1)/2)*Mn_inf^2)/(gamma*Mn_inf^2 - (gamma-1)/2));
    M3 = Mn3/(sind(beta-(theta)));
    
    %% Find p4/p_inf with expansion fan %%
    nu3 = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M3^2 - 1))) - atand(sqrt(M3^2 - 1));
    nu4 = nu3 + (2*epsilon);
    syms x
    M4 = abs(vpasolve(nu4 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T4 = (1+((1.4-1)/2)*(M3^2))/(1+((1.4-1)/2)*M4^2);
    p4p3 = double((T4)^((1.4)/(1.4-1)));
    p4 = (p4p3)*p3;
    
    %% Compute c_n and c_a %%
    c_n = (1/(gamma*M^2))*(-p1 - p2 + p3 + p4);
    c_a = (1/(gamma*M^2))*(p1 - p2 + p3 - p4)*tand(epsilon);
    
    %% Compute c_l and c_d %%
    c_l = c_n*cosd(alpha) - c_a*sind(alpha);
    c_d = c_n*sind(alpha) + c_a*cosd(alpha);

elseif(alpha == epsilon)
    %% Find p1/p_inf %%
    nu_inf = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M^2 - 1))) - atand(sqrt(M^2 - 1));
    nu1 = nu_inf;
    syms x
    M1 = abs(vpasolve(nu1 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    p1 = 1;
    %% Find p2/p_inf with expansion fan %%
    nu2 = nu1 + (2*epsilon);
    syms x
    M2 = abs(vpasolve(nu2 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T2 = (1+((1.4-1)/2)*(M1^2))/(1+((1.4-1)/2)*M2^2);
    p2p1 = double((T2)^((1.4)/(1.4-1)));
    p2 = (p2p1)*p1;
    
    %% Find p3/p_inf with oblique shock relations %%
    theta = (2*epsilon);
    LHS = tan(theta*pi/180);
    beta = 1;
    RHS = inline('2*cot(beta*pi/180)*((M*sin(beta*pi/180))^2-1)/(M^2*(1.4+cos(2*beta*pi/180))+2)','beta','M');
    counter = 0;
    while (RHS(beta,M) < LHS*0.98 | RHS(beta,M) > LHS*1.02)
        if (beta > 90 | beta < 0 | counter >200)
            break;
        end
        if (RHS(beta,M) < LHS*0.70)
            beta=beta+1;
        elseif(RHS(beta,M) < LHS*0.90)
            beta=beta+0.5;
        elseif(RHS(beta,M) < LHS*0.95)
            beta=beta+0.3;
        elseif(RHS(beta,M) < LHS*0.98)
            beta=beta+0.1;
        end
        
        if(RHS(beta,M) > LHS*1.03)
            beta=beta-0.1;
        elseif(RHS(beta,M) > LHS*1.10)
            beta = beta-0.5;
        end
        counter = counter+1;   
    end
    
    counter = 0;
    
    while (RHS(beta,M) < LHS*0.99 | RHS(beta,M) > LHS*1.01)
        if (beta > 90 | beta < 0 | counter > 100)
            break;
        end
        if (RHS(beta,M) < LHS*0.98)
            beta=beta+0.05;
        elseif(RHS(beta,M) < LHS*0.985)
            beta=beta+0.03;
        elseif(RHS(beta,M) < LHS*0.99)
            beta=beta+0.01;
        end
        
        if(RHS(beta,M) > LHS*1.01)
            beta=beta-0.01;
        elseif(RHS(beta,M) > LHS*1.02)
            beta = beta-0.1;
        end
        counter = counter+1;
    end
    
    counter = 0;
    
    while (RHS(beta,M) ~= LHS)
        if (beta > 90 | beta < 0 | counter > 40)
            break;
        end
        
        if (RHS(beta,M) < LHS)
            beta=beta+0.005;
        end
        if(RHS(beta,M) > LHS)
            beta=beta-0.005;
        end
        counter = counter+1;
    end
    
    if (beta < 0)
        fprintf('Beta not in bounds');
    end
    Mn_inf = M*sind(beta);
    p3 = double(1 + ((2*gamma)/(gamma+1))*(Mn_inf^2 - 1));
    Mn3 = sqrt((1+((gamma-1)/2)*Mn_inf^2)/(gamma*Mn_inf^2 - (gamma-1)/2));
    M3 = Mn3/(sind(beta-(theta)));
    
    %% Find p4/p_inf with expansion fan %%
    nu3 = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M3^2 - 1))) - atand(sqrt(M3^2 - 1));
    nu4 = nu3 + (2*epsilon);
    syms x
    M4 = abs(vpasolve(nu4 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T4 = (1+((1.4-1)/2)*(M3^2))/(1+((1.4-1)/2)*M4^2);
    p4p3 = double((T4)^((1.4)/(1.4-1)));
    p4 = (p4p3)*p3;
    
    %% Compute c_n and c_a %%
    c_n = (1/(gamma*M^2))*(-p1 - p2 + p3 + p4);
    c_a = (1/(gamma*M^2))*(p1 - p2 + p3 - p4)*tand(epsilon);
    
    %% Compute c_l and c_d %%
    c_l = c_n*cosd(alpha) - c_a*sind(alpha);
    c_d = c_n*sind(alpha) + c_a*cosd(alpha);
    
elseif(alpha > epsilon)
    %% Find p1/p_inf with expansion fan %%
    nu_inf = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M^2 - 1))) - atand(sqrt(M^2 - 1));
    nu1 = nu_inf + (alpha-epsilon);
    syms x
    M1 = abs(vpasolve(nu1 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T1 = (1+((1.4-1)/2)*(M^2))/(1+((1.4-1)/2)*M1^2);
    p1 = double((T1)^((1.4)/(1.4-1)));
    
    %% Find p2/p_inf with expansion fan %%
    nu2 = nu1 + (2*epsilon);
    syms x
    M2 = abs(vpasolve(nu2 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T2 = (1+((1.4-1)/2)*(M1^2))/(1+((1.4-1)/2)*M2^2);
    p2p1 = double((T2)^((1.4)/(1.4-1)));
    p2 = (p2p1)*p1;
    
    %% Find p3/p_inf with oblique shock relations %%
    theta = (alpha+epsilon);
    LHS = tan(theta*pi/180);
    beta = 1;
    RHS = inline('2*cot(beta*pi/180)*((M*sin(beta*pi/180))^2-1)/(M^2*(1.4+cos(2*beta*pi/180))+2)','beta','M');
    counter = 0;
    while (RHS(beta,M) < LHS*0.98 | RHS(beta,M) > LHS*1.02)
        if (beta > 90 | beta < 0 | counter >200)
            break;
        end
        if (RHS(beta,M) < LHS*0.70)
            beta=beta+1;
        elseif(RHS(beta,M) < LHS*0.90)
            beta=beta+0.5;
        elseif(RHS(beta,M) < LHS*0.95)
            beta=beta+0.3;
        elseif(RHS(beta,M) < LHS*0.98)
            beta=beta+0.1;
        end
        
        if(RHS(beta,M) > LHS*1.03)
            beta=beta-0.1;
        elseif(RHS(beta,M) > LHS*1.10)
            beta = beta-0.5;
        end
        counter = counter+1;   
    end
    
    counter = 0;
    
    while (RHS(beta,M) < LHS*0.99 | RHS(beta,M) > LHS*1.01)
        if (beta > 90 | beta < 0 | counter > 100)
            break;
        end
        if (RHS(beta,M) < LHS*0.98)
            beta=beta+0.05;
        elseif(RHS(beta,M) < LHS*0.985)
            beta=beta+0.03;
        elseif(RHS(beta,M) < LHS*0.99)
            beta=beta+0.01;
        end
        
        if(RHS(beta,M) > LHS*1.01)
            beta=beta-0.01;
        elseif(RHS(beta,M) > LHS*1.02)
            beta = beta-0.1;
        end
        counter = counter+1;
    end
    
    counter = 0;
    
    while (RHS(beta,M) ~= LHS)
        if (beta > 90 | beta < 0 | counter > 40)
            break;
        end
        
        if (RHS(beta,M) < LHS)
            beta=beta+0.005;
        end
        if(RHS(beta,M) > LHS)
            beta=beta-0.005;
        end
        counter = counter+1;
    end
    
    if (beta < 0)
        fprintf('Beta not in bounds');
    end
    
    Mn_inf = M*sind(beta);
    p3 = double(1 + ((2*gamma)/(gamma+1))*(Mn_inf^2 - 1));
    Mn3 = sqrt((1+((gamma-1)/2)*Mn_inf^2)/(gamma*Mn_inf^2 - (gamma-1)/2));
    M3 = Mn3/(sind(beta-(theta)));
    
    %% Find p4/p_inf with expansion fan %%
    nu3 = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(M3^2 - 1))) - atand(sqrt(M3^2 - 1));
    nu4 = nu3 + (2*epsilon);
    syms x
    M4 = abs(vpasolve(nu4 == sqrt(6)*(atan(sqrt((1/6)*(x^2 - 1)))*(180/pi)) - (atan(sqrt(x^2 - 1))*(180/pi)),x));
    T4 = (1+((1.4-1)/2)*(M3^2))/(1+((1.4-1)/2)*M4^2);
    p4p3 = double((T4)^((1.4)/(1.4-1)));
    p4 = (p4p3)*p3;
    
    %% Compute c_n and c_a %%
    c_n = (1/(gamma*M^2))*(-p1 - p2 + p3 + p4);
    c_a = (1/(gamma*M^2))*(p1 - p2 + p3 - p4)*tand(epsilon);
    
    %% Compute c_l and c_d %%
    c_l = c_n*cosd(alpha) - c_a*sind(alpha);
    c_d = c_n*sind(alpha) + c_a*cosd(alpha);
end
end