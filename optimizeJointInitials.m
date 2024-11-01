function [optimizedJointInitials] = optimizeJointInitials(robotBasePoseToWorld, simulatedElbowFrameToBaseLink, sensedFrameCToBaseLink, upperArmLength, forearmLength, nominalJointInitial)
% 1, assume that there is no position error, and compute the  BC_R_0 (geometrical relation between the frame {C} and the frame {B} when ψ(psi) and q5 are both zeros)
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

% 2, calculate the constant K for [sin(psi), cos(psi),1]*K*[sin(q5), cos(q5),1]'
% For equation 18,  arg max Tr(inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 * e^([n_psi_In_BC_R_0]*psi) * e^([xi5]*q5))
% inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 = L

% e^([n_psi_In_BC_R_0]*psi) = (I+[n_psi_In_BC_R_0]^2) + sin(psi)*[n_psi_In_BC_R_0] + cos(psi)*(-[n_psi_In_BC_R_0]^2) = M1 + sin(psi)*M2 + cos(psi)*M3
% e^([xi5]*q5) = (I+[xi5]^2) + sin(q5)*[xi5] + cos(q5)*(-[xi5]^2) = N1 + sin(q5)*N2 + cos(q5)*N3 (https://thenumb.at/Exponential-Rotations/)

% let: inv(frameCToBaseLink_sensorMeasurement(1:3,1:3) * BC_R_0 * e^([n_psi_In_BC_R_0]*psi) * e^([xi5]*q5) =
% C0 + C1*cos(psi) + C2*sin(psi) + C3*cos(q5) + C4*sin(q5) + C5*cos(psi)*cos(q5) + C6*cos(psi)*sin(q5) + C7*sin(psi)*cos(q5) + C8*sin(psi)*sin(q5)

n_psi_In_BC_R_0 = BC_R_0 \ (-wristPosition/norm(wristPosition));
M2=[0, -n_psi_In_BC_R_0(3),n_psi_In_BC_R_0(2);  n_psi_In_BC_R_0(3), 0, -n_psi_In_BC_R_0(1); -n_psi_In_BC_R_0(2),n_psi_In_BC_R_0(1),0];  %sin(psi) coefficient
M1=eye([3,3])+M2*M2;
M3=-M2*M2;       %cos(psi) coefficient

xi5=[-1;0;0];
N2=[0, -xi5(3),xi5(2);  xi5(3), 0, -xi5(1); -xi5(2),xi5(1),0]; %sin(q5) coefficient
N1=eye([3,3])+N2*N2;
N3=-N2*N2; %cos(q5) coefficient

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

% [sin(psi), cos(psi),1]*K*[sin(q5), cos(q5),1]'
tr_constant = trace(C0);
tr_c_psi = trace(C1);
tr_s_psi = trace(C2);
tr_c_q5 = trace(C3);
tr_s_q5 = trace(C4);
tr_c_psi_c_q5 = trace(C5);
tr_c_psi_s_q5 = trace(C6);
tr_s_psi_c_q5 = trace(C7);
tr_s_psi_s_q5 = trace(C8);

K = [tr_s_psi_s_q5,tr_s_psi_c_q5,tr_s_psi; tr_c_psi_s_q5, tr_c_psi_c_q5,tr_c_psi;tr_s_q5,tr_c_q5,tr_constant];



% 3 compute nominal initial for q5 and q5
q5_initial = nominalJointInitial(5);
elbowPosition = simulatedElbowFrameToBaseLink(1:3,4);
armSwivalAnglePlaneNormalDirection = cross(wristPosition,elbowPosition);
psi_initial = acos(dot(armSwivelReferencePlaneNormalDirection/norm(armSwivelReferencePlaneNormalDirection), armSwivalAnglePlaneNormalDirection/norm(armSwivalAnglePlaneNormalDirection)));

if dot(gravityDirection,armSwivalAnglePlaneNormalDirection)>0
    psi_initial = -psi_initial;
end

[optimalPsiInitial,optimalQ5Initial] = gradientDecentForJointInitials(K,psi_initial,q5_initial);


% 4,map psi and q5 to 5 joint angles
elbowPosition = axang2rotm([(-wristPosition/norm(wristPosition))', optimalPsiInitial]) * elbowPositionAtReferencePlane;

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

% choose the solution closest to the nominal one
normError1 = norm(jointAngle(:,1)-nominalJointInitial(1:3));
normError2 = norm(jointAngle(:,2)-nominalJointInitial(1:3));
if normError1<normError2
    validIKSolutionIndex = 1;
else
    validIKSolutionIndex = 2;
end

optimizedJointInitials(1) = jointAngle(1,validIKSolutionIndex);  % joint 1, 2, 3 is in [-pi, pi)
optimizedJointInitials(2) = jointAngle(2,validIKSolutionIndex);
optimizedJointInitials(3) = jointAngle(3,validIKSolutionIndex);
optimizedJointInitials(4) = upperArmForeamrAngle; % upperArmForeamrAngle is in [0,pi]
optimizedJointInitials(5) = (optimalQ5Initial -2*pi*fix((optimalQ5Initial+pi)/(2*pi))); %map finalq5Inintial from [0,2pi) to [-pi, pi)

end

