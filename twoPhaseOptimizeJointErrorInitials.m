function [jointErrorInitials] = twoPhaseOptimizeJointErrorInitials(robotBasePoseToWorld, sensedFrameCToBaseLink, upperArmLength, forearmLength, commandedJointPosition)
% 1, assume that there is no position error, and compute the  BC_R_0 (geometrical relation between the frame {C} and the frame {B} when ψ(psi) and φ(phi) are both zeros)
gravityDirection = [1,0,0]';
wristPosition = sensedFrameCToBaseLink(1:3,4);
armSwivelReferencePlaneNormalDirection =  cross(wristPosition,gravityDirection);

distanceFromShoulderToWrist = norm(wristPosition);
upperArmWristDirectionAngle = acos((upperArmLength^2 + distanceFromShoulderToWrist^2 - forearmLength^2)/(2*upperArmLength*distanceFromShoulderToWrist));
elbowPositionAtReferencePlane =axang2rotm([armSwivelReferencePlaneNormalDirection', upperArmWristDirectionAngle])*(wristPosition/norm(wristPosition) * upperArmLength) ;

upperArmForeamrAngle = pi - acos((upperArmLength^2 + forearmLength^2 - distanceFromShoulderToWrist^2)/(2*upperArmLength*forearmLength)); %in [0, pi]

BC_R_0 = zeros([3,3]);
BC_R_0(1:3,1) = ((wristPosition-elbowPositionAtReferencePlane)/norm(wristPosition-elbowPositionAtReferencePlane));
BC_R_0(1:3,2) = (armSwivelReferencePlaneNormalDirection/norm(armSwivelReferencePlaneNormalDirection));
BC_R_0(1:3,3) = cross(BC_R_0(1:3,1),BC_R_0(1:3,2))/norm(cross(BC_R_0(1:3,1),BC_R_0(1:3,2)));

% show BC_R_0
axisLength = 0.15;
wristPositionInWorld = robotBasePoseToWorld*[wristPosition;1];
elbowPositionAtReferencePlaneInWorld = robotBasePoseToWorld*[elbowPositionAtReferencePlane;1];
plot3([0;elbowPositionAtReferencePlaneInWorld(1);wristPositionInWorld(1);0],[0;elbowPositionAtReferencePlaneInWorld(2);wristPositionInWorld(2);0],[0;elbowPositionAtReferencePlaneInWorld(3);wristPositionInWorld(3);0]);

BC_R_0_inWorld = robotBasePoseToWorld(1:3,1:3)*BC_R_0;
quiver3(wristPositionInWorld(1),wristPositionInWorld(2),wristPositionInWorld(3),BC_R_0_inWorld(1,1)*axisLength,BC_R_0_inWorld(2,1)*axisLength,BC_R_0_inWorld(3,1)*axisLength,'r-','LineWidth',2); %x
quiver3(wristPositionInWorld(1),wristPositionInWorld(2),wristPositionInWorld(3),BC_R_0_inWorld(1,2)*axisLength,BC_R_0_inWorld(2,2)*axisLength,BC_R_0_inWorld(3,2)*axisLength,'g-','LineWidth',2); %y
quiver3(wristPositionInWorld(1),wristPositionInWorld(2),wristPositionInWorld(3),BC_R_0_inWorld(1,3)*axisLength,BC_R_0_inWorld(2,3)*axisLength,BC_R_0_inWorld(3,3)*axisLength,'b-','LineWidth',2); %z

% 2, calculate the constant K for [sin(psi), cos(psi),1]*K*[sin(phi), cos(phi),1]'
% For equation 4,  arg max Tr(inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 * e^([n_psi_In_BC_R_0]*psi) * e^([xi5]*phi))
% inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 = L

% e^([n_psi_In_BC_R_0]*psi) = (I+[n_psi_In_BC_R_0]^2) + sin(psi)*[n_psi_In_BC_R_0] + cos(psi)*(-[n_psi_In_BC_R_0]^2) = M1 + sin(psi)*M2 + cos(psi)*M3
% e^([xi5]*phi) = (I+[xi5]^2) + sin(phi)*[xi5] + cos(phi)*(-[xi5]^2) = N1 + sin(phi)*N2 + cos(phi)*N3 (https://thenumb.at/Exponential-Rotations/)

% let: inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 * e^([n_psi_In_BC_R_0]*psi) * e^([xi5]*phi) =
% C0 + C1*cos(psi) + C2*sin(psi) + C3*cos(phi) + C4*sin(phi) + C5*cos(psi)*cos(phi) + C6*cos(psi)*sin(phi) + C7*sin(psi)*cos(phi) + C8*sin(psi)*sin(phi)

n_psi_In_BC_R_0 = BC_R_0 \ (-wristPosition/norm(wristPosition));
M2=[0, -n_psi_In_BC_R_0(3),n_psi_In_BC_R_0(2);  n_psi_In_BC_R_0(3), 0, -n_psi_In_BC_R_0(1); -n_psi_In_BC_R_0(2),n_psi_In_BC_R_0(1),0];  %sin(psi) coefficient
M1=eye([3,3])+M2*M2;
M3=-M2*M2;       %cos(psi) coefficient

xi5=[-1;0;0];
N2=[0, -xi5(3),xi5(2);  xi5(3), 0, -xi5(1); -xi5(2),xi5(1),0]; %sin(phi) coefficient
N1=eye([3,3])+N2*N2;
N3=-N2*N2; %cos(phi) coefficient

L=sensedFrameCToBaseLink(1:3,1:3)\ BC_R_0;

C0 = L*M1*N1;
C1 = L*M3*N1;
C2 = L*M2*N1;
C3 = L*M1*N3;
C4 = L*M1*N2;
C5 = L*M3*N3;
C6 = L*M3*N2;
C7 = L*M2*N3;
C8 = L*M2*N2;

% [sin(psi), cos(psi),1]*K*[sin(phi), cos(phi),1]'
tr_constant = trace(C0);
tr_c_psi = trace(C1);
tr_s_psi = trace(C2);
tr_c_phi = trace(C3);
tr_s_phi = trace(C4);
tr_c_psi_c_phi = trace(C5);
tr_c_psi_s_phi = trace(C6);
tr_s_psi_c_phi = trace(C7);
tr_s_psi_s_phi = trace(C8);

K = [tr_s_psi_s_phi,tr_s_psi_c_phi,tr_s_psi; tr_c_psi_s_phi, tr_c_psi_c_phi,tr_c_psi;tr_s_phi,tr_c_phi,tr_constant];

% 3,two-phase global maximization algorithm to compute psi and phi
angleInterval=5;
stopCriteron = 0.01;

while 1
    trialsNum= 360/angleInterval;
    foundLocalMaximums = -1*ones([3,trialsNum*trialsNum]); %-1 indicates not assigned
    foundLocalMaximumsNum=0;
    for i=1:trialsNum
        for j=1:trialsNum
            psiSeed = i*angleInterval/180*pi;
            phiSeed = j*angleInterval/180*pi;
            [psiInitialCandidate,phiInitialCandidate] = localOptimizationForJointErrorInitials(K,psiSeed,phiSeed);
            
            if foundLocalMaximums(3,1) <0
                foundLocalMaximums(3,1) = 1;
                foundLocalMaximums(1:2,1) = [psiInitialCandidate,phiInitialCandidate]';
                foundLocalMaximumsNum = 1;
                fprintf("when angle interval is %f degrees, find %dth local minimum, psi:%f rad, phi:%f rad\n",angleInterval,1, foundLocalMaximums(1,1),foundLocalMaximums(2,1));
            else
                k=1;
                existTag = 0;
                while foundLocalMaximums(3,k)>0
                    if(norm( foundLocalMaximums(1:2,k) - [psiInitialCandidate,phiInitialCandidate]')<1e-5)
                        existTag =1;
                        break;
                    else
                        k=k+1;
                    end
                end
                if existTag < 1
                    foundLocalMaximumsNum = foundLocalMaximumsNum+1;
                    foundLocalMaximums(1:2,foundLocalMaximumsNum) = [psiInitialCandidate,phiInitialCandidate]';
                    foundLocalMaximums(3,foundLocalMaximumsNum) = 1;
                    fprintf("when angle interval is %f degrees, find %dth local minimum, psi:%f rad, phi:%f rad\n",angleInterval,foundLocalMaximumsNum, foundLocalMaximums(1,foundLocalMaximumsNum),foundLocalMaximums(2,foundLocalMaximumsNum));
                end
            end
        end
    end
    
    if(((foundLocalMaximumsNum*(foundLocalMaximumsNum+1))/(trialsNum*(trialsNum-1)))<stopCriteron)
        break;
    else
        angleInterval = angleInterval/2;
    end
end


finalPsiInintial = foundLocalMaximums(1,1);
finalPhiInintial = foundLocalMaximums(2,1);  % [0,2pi)
if foundLocalMaximumsNum>1
    objectValue = func([finalPsiInintial,finalPhiInintial],K);
    for i=2:foundLocalMaximumsNum
        objectValueCandidate = func([foundLocalMaximums(1,i),foundLocalMaximums(2,i)],K);
        if objectValueCandidate>objectValue
            finalPsiInintial = foundLocalMaximums(1,i);
            finalPhiInintial = foundLocalMaximums(2,i);
            objectValue = objectValueCandidate;
        end
    end
end

% 4,map psi and phi to 5 joint position theta, and 5 joint position error epsilon
elbowPosition = axang2rotm([(-wristPosition/norm(wristPosition))', finalPsiInintial]) * elbowPositionAtReferencePlane;

% compute the parent frame of the elbow considering the colinear case of shoulder center, elbow center and wrist center
elbowJointParentFrame = eye([3,3]);
elbowJointParentFrame(:,1) = elbowPosition/norm(elbowPosition);
elbowJointParentFrame(:,2) = cross(wristPosition,elbowJointParentFrame(:,1));
if norm(elbowJointParentFrame(:,2)) < 1e-5 % shoulder center, elbow center, wrist center are colinear.
    elbowJointParentFrame(:,2) = cross(sensedFrameCToBaseLink(1:3,3), elbowJointParentFrame(:,1));
    elbowJointParentFrame(:,2) = elbowJointParentFrame(:,2)/norm(elbowJointParentFrame(:,2));
    elbowJointParentFrame(:,3) = cross(elbowJointParentFrame(:,1),elbowJointParentFrame(:,2));
    elbowJointParentFrame(:,3) = elbowJointParentFrame(:,3)/norm(elbowJointParentFrame(:,3));
else
    elbowJointParentFrame(:,2) = elbowJointParentFrame(:,2)/norm(elbowJointParentFrame(:,2));
    elbowJointParentFrame(:,3) = cross(elbowJointParentFrame(:,1),elbowJointParentFrame(:,2));
    elbowJointParentFrame(:,3) = elbowJointParentFrame(:,3)/norm(elbowJointParentFrame(:,3));
end

% compute the euler angle of order of Y-Z-X (-jointAngle1,-jointAngle2,-jointAngle3) into [-pi, pi)
eulerAngles = zeros([3,2]);
eulerAngles(2,1) = asin(elbowJointParentFrame(2,1));
if eulerAngles(2,1) >0
    eulerAngles(2,2) = pi-eulerAngles(2,1);
else
    eulerAngles(2,2) = -pi-eulerAngles(2,1);
end

if norm(cos(eulerAngles(2,1)))<1e-10   % gimbal lock, joint angle 2 is 90 degrees or -90 degrees
    if sin(eulerAngles(2,1)) >0  %90 degree, each angle have half an rotation angle with same sign
        eulerAngles(1,1) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(3,2))/2;
        eulerAngles(1,2) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(3,2))/2;
        
        eulerAngles(3,1) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(3,2))/2;
        eulerAngles(3,2) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(3,2))/2;
    else  % -90 degree,   each angle have half an rotation angle with opposite sign
        eulerAngles(1,1) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(1,3))/2;
        eulerAngles(1,2) = aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(1,3))/2;
        
        eulerAngles(3,1) = -aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(1,3))/2;
        eulerAngles(3,2) = -aCosSin(elbowJointParentFrame(3,3),elbowJointParentFrame(1,3))/2;
    end
else
    eulerAngles(1,1) = aCosSin(elbowJointParentFrame(1,1)/cos(eulerAngles(2,1)),-elbowJointParentFrame(3,1)/cos(eulerAngles(2,1)));
    eulerAngles(1,2) = aCosSin(elbowJointParentFrame(1,1)/cos(eulerAngles(2,2)),-elbowJointParentFrame(3,1)/cos(eulerAngles(2,2)));
    
    eulerAngles(3,1) = aCosSin(elbowJointParentFrame(2,2)/cos(eulerAngles(2,1)),-elbowJointParentFrame(2,3)/cos(eulerAngles(2,1)));
    eulerAngles(3,2) = aCosSin(elbowJointParentFrame(2,2)/cos(eulerAngles(2,2)),-elbowJointParentFrame(2,3)/cos(eulerAngles(2,2)));
end

% get the joint angle considering difference in reference directions
jointAngle =diag([-1,-1,-1])* eulerAngles;
errorInLimitsTag = [1,1];
errorLimit = 20/180*pi;
for i=1:3
    for j=1:2
        if abs(jointAngle(i,j) - commandedJointPosition(i)) >errorLimit
            errorInLimitsTag(j)=0;
        end
    end
end

if sum(errorInLimitsTag)>0
    validIKSolutionIndex = find(errorInLimitsTag,1);
else
    normError1 = norm(jointAngle(:,1)-commandedJointPosition(1:3));
    normError2 = norm(jointAngle(:,2)-commandedJointPosition(1:3));
    if normError1<normError2
        validIKSolutionIndex = 1;
    else
        validIKSolutionIndex = 2;
    end
end

jointErrorInitials(1) = commandedJointPosition(1) - jointAngle(1,validIKSolutionIndex);  % joint 1, 2, 3 is in [-pi, pi)
jointErrorInitials(2) = commandedJointPosition(2) - jointAngle(2,validIKSolutionIndex);
jointErrorInitials(3) = commandedJointPosition(3) - jointAngle(3,validIKSolutionIndex);
jointErrorInitials(4) = commandedJointPosition(4) - upperArmForeamrAngle; % upperArmForeamrAngle is in [0,pi]
jointErrorInitials(5) = commandedJointPosition(5) - (finalPhiInintial -2*pi*fix((finalPhiInintial+pi)/(2*pi))); %map finalPhiInintial from [0,2pi) to [-pi, pi)

end

