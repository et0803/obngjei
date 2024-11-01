function [optimalPsiInitial,optimalQ5Initial] = localOptimizationForJointErrorInitials(K,psiSeed,q5Seed)
% maximize function f0 = [sin(psi), cos(psi),1]*K*[sin(q5), cos(q5),1]'
% minimize function f = -[sin(psi), cos(psi),1]*K*[sin(q5), cos(q5),1]' =  [sin(psi), cos(psi),1]*K_N*[sin(q5), cos(q5),1]', where K_N = -K

% minimize function f = psiFeature'*K_N*q5Feature
K_N = -K;

x0 = [psiSeed,q5Seed];

% d(psiFeature')/d(psi) = [cos(psi), -sin(psi), 0 ] =[sin(psi), cos(psi),1] * G_psi = psiFeature'* G_psi;
G_psi = [0,-1,0;1,0,0;0,0,0];
df_dPsi_matrix = G_psi*K_N;             % df_dPsi = psiFeature' * df_dPsi_matrix * q5Feature

% d(q5Feature)/d(q5) = [cos(q5); -sin(q5); 0 ] = G_q5 * [sin(q5); cos(q5);1]= G_q5 * q5Feature ;
G_q5 = [0,1,0;-1,0,0;0,0,0];
df_dq5_matrix =K_N*G_q5;           % df_dq5 = psiFeature' * df_dq5_matrix * q5Feature


% https://www.mathworks.com/matlabcentral/answers/787539-steepest-descent-algorithm-in-matlab
stepSize = 0.1;
stepSizeLowBound =0.0001;
tol = 1e-6;

while 1
    x_prev = x0;
    grad_prev = [func(x_prev,df_dPsi_matrix), func(x_prev,df_dq5_matrix)];
    temp = x_prev -  stepSize*grad_prev;
    
    if(func(x_prev,K_N)<func(temp,K_N))
        stepSize = stepSize/2;
        if stepSize <stepSizeLowBound
            break;
        end
    else
        x0=temp;
        if norm(x0-x_prev)<tol
            break;
        end
    end
    %fprintf("iteration psi:%f, q5:%f, stepSize:%f,loss:%f\n",x0(1),x0(2),stepSize,func(x0,K_N))
end

% map angle to [0,2*pi)
x0 = x0 -2*pi*fix(x0/(2*pi));
x0 = x0+2*pi;
x0 = x0 -2*pi*fix(x0/(2*pi));

optimalPsiInitial = x0(1);
optimalQ5Initial = x0(2);
end


