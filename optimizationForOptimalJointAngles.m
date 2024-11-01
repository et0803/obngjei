function [finalOptimizedJointAngles] = optimizationForOptimalJointAngles(roboticArm, jointConfigStruct,sensedFrameCToBaseLink, optimalJointInitials)
% modified from:
% https://github.com/et0803/memikrr/blob/main/src/hybridIK.cpp
% https://github.com/orocos/orocos_kinematics_dynamics/blob/master/orocos_kdl/src/chainiksolverpos_lma.cpp

jointPosition = optimalJointInitials;

W = diag([1,1,1,10,10,10]);
stopCriterion =5e-4;


tau=10;
lambda = tau;
v=2;
eps_joint = 1e-15;
maxiter=10000;

for i=1:size(jointConfigStruct,2)
    if i==1        
        jointConfigStruct(i).JointPosition = 0;  % simplified the jacobian matrix computation.
    else
        if i<7
            jointConfigStruct(i).JointPosition =jointPosition(i-1);  % 5 joint for 7-DoF arm
        else
            jointConfigStruct(i).JointPosition =0;  % joint for fingers
        end
    end
end

currentPose = getTransform(roboticArm,jointConfigStruct,'link6','world');
deltaPose_twist = twist_from_pose_diff(currentPose,sensedFrameCToBaseLink); %[wx, wy, wz, vx, vy, vz]
weightedDeltaPose_twist = W* deltaPose_twist;
weightedDeltaPose_twist_norm = norm(weightedDeltaPose_twist);

currentGeoJacob  = geometricJacobian(roboticArm, jointConfigStruct,'link6'); %V_ee = J * \dot(q)  V_ee=[wx, wy, wz, vx, vy, vz]
weightedGeoJacob =  W*currentGeoJacob(:,2:6);

if weightedDeltaPose_twist_norm<stopCriterion
    fprintf('init value is accepted! position error(x,y,z) /m: (%f,%f,%f); orientation error /rad:%f \n',deltaPose_twist(4),deltaPose_twist(5),deltaPose_twist(6),norm(deltaPose_twist(1:3)));
else
    for i=1:maxiter
        [U,SS,V] = svd(weightedGeoJacob);
        S=zeros(size(SS,2),size(SS,1));
        for k=1:5
            S(k,k) = SS(k,k)/(SS(k,k)*SS(k,k)+lambda); %LM算法,更新Hessian
        end
        diffq =V*S*U'* weightedDeltaPose_twist;
        grad = weightedGeoJacob'*weightedDeltaPose_twist;
        diffqNormLInfinity = norm(diffq, 'inf');  %关节增量最大的分量
        
        if diffqNormLInfinity<eps_joint
            fprintf('%dth LM iteration. The joint position increments are to small, lambda:%.3f; v:%.3f\n',i,lambda,v);
            break;
        end
        if grad'*grad < eps_joint*eps_joint
            fprintf('%dth LM iteration. The gradient of E towards the joints is to small, lambda:%.3f; v:%.3f\n',i,lambda,v);
            break;
        end
        
        jointPosition_new = jointPosition+diffq';
        for j=2:6
            jointConfigStruct(j).JointPosition =jointPosition_new(j-1);  % 5 joint for 7-DoF arm
        end
        currentPose = getTransform(roboticArm,jointConfigStruct,'link6','world');
        deltaPose_twist_new = twist_from_pose_diff(currentPose,sensedFrameCToBaseLink); %[wx, wy, wz, vx, vy, vz]
        weightedDeltaPose_twist_new = W* deltaPose_twist_new;
        weightedDeltaPose_twist_new_norm = norm(weightedDeltaPose_twist_new);
        rho = weightedDeltaPose_twist_norm^2 - weightedDeltaPose_twist_new_norm^2;
        rho = rho / (diffq'*(lambda*diffq+grad));
        
        if rho>0
            jointPosition = jointPosition_new;
            weightedDeltaPose_twist = weightedDeltaPose_twist_new;
            weightedDeltaPose_twist_norm = weightedDeltaPose_twist_new_norm;
            
            if weightedDeltaPose_twist_norm < stopCriterion
                fprintf('%dth LM iteration. Return with no error, lambda:%.3f; v:%.3f\n',i,lambda,v);
                break;
            end
            
            currentGeoJacob  = geometricJacobian(roboticArm, jointConfigStruct,'link6'); %V_ee = J * \dot(q)  V_ee=[wx, wy, wz, vx, vy, vz]
            weightedGeoJacob =  W*currentGeoJacob(:,2:6);
            temp = 2*rho -1;
            lambda = lambda*max(1/3.0, 1-temp^3);
            v=2;
        else
            lambda = lambda*v;
            v=v*2;
        end
        
        if i == maxiter
            fprintf('%dth LM iteration. Max iteration exceeded, lambda:%.3f; v:%.3f\n',i,lambda,v);
        end
    end
end

finalOptimizedJointAngles = jointPosition;
end

