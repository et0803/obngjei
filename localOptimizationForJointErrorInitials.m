function [psiInitialCandidate,phiInitialCandidate] = localOptimizationForJointErrorInitials(K,psiSeed,phiSeed)
% maximize function f0 = [sin(psi), cos(psi),1]*K*[sin(phi), cos(phi),1]'
% minimize function f = -[sin(psi), cos(psi),1]*K*[sin(phi), cos(phi),1]' =  [sin(psi), cos(psi),1]*K_N*[sin(phi), cos(phi),1]', where K_N = -K

% minimize function f = psiFeature'*K_N*phiFeature
K_N = -K;

x0 = [psiSeed,phiSeed];

% d(psiFeature')/d(psi) = [cos(psi), -sin(psi), 0 ] =[sin(psi), cos(psi),1] * G_psi = psiFeature'* G_psi;
G_psi = [0,-1,0;1,0,0;0,0,0];
df_dPsi_matrix = G_psi*K_N;             % df_dPsi = psiFeature' * df_dPsi_matrix * phiFeature

% d(phiFeature)/d(phi) = [cos(phi); -sin(phi); 0 ] = G_phi * [sin(phi); cos(phi);1]= G_phi * phiFeature ;
G_phi = [0,1,0;-1,0,0;0,0,0];
df_dPhi_matrix =K_N*G_phi;           % df_dPhi = psiFeature' * df_dPhi_matrix * phiFeature


% https://www.mathworks.com/matlabcentral/answers/787539-steepest-descent-algorithm-in-matlab
stepSize = 0.1;
stepSizeLowBound =0.0001;
tol = 1e-6;

while 1
    x_prev = x0;
    grad_prev = [func(x_prev,df_dPsi_matrix), func(x_prev,df_dPhi_matrix)];
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
    %fprintf("iteration psi:%f, phi:%f, stepSize:%f,loss:%f\n",x0(1),x0(2),stepSize,func(x0,K_N))
end

% map angle to [0,2*pi)
x0 = x0 -2*pi*fix(x0/(2*pi));
x0 = x0+2*pi;
x0 = x0 -2*pi*fix(x0/(2*pi));

psiInitialCandidate = x0(1);
phiInitialCandidate = x0(2);
end


