%% show the model of the 7-DoF anthropomorphic robotic arm
clc;
clear;
roboticArmModelFileName='./x_arm_aligned.urdf'; 

roboticArm = importrobot(roboticArmModelFileName);
jointConfigStruct =randomConfiguration(roboticArm);

jointErrorLevelForEpsilonMeas = 8/180*pi;
jointErrorLevelForEpsilon3 = 8/180*pi;

commandedJointPosition =[0.71, 0.58, 1.04, 1.91, -0.99];
epsilonMeas = jointErrorLevelForEpsilonMeas*randn(1,5);
epsilon3 = jointErrorLevelForEpsilon3*randn(1,5);
acturalJointPosition_GroudTruth = commandedJointPosition + epsilonMeas + epsilon3;
nominalJointInitial = commandedJointPosition + epsilonMeas;

for i=1:size(jointConfigStruct,2)
    if i==1         % robot base to world, for a better plot view
        jointConfigStruct(i).JointPosition = pi/2;  
    else
        if i<7
            jointConfigStruct(i).JointPosition = acturalJointPosition_GroudTruth(i-1);  % 5 joint for 7-DoF arm
        else
            jointConfigStruct(i).JointPosition =0;  % joint for fingers
        end
    end
end

robotBasePoseToWorld  = getTransform(roboticArm,jointConfigStruct,'base_link','world');
actualFrameCToRobotBase  = getTransform(roboticArm,jointConfigStruct,'link6','base_link');
simulatedElbowFrameToBaseLink = getTransform(roboticArm,jointConfigStruct,'link4','base_link');

upperArmLength = 0.33;
forearmLength = 0.257;

show(roboticArm,jointConfigStruct,'PreservePlot',0 ,'Frames', 'off', 'visuals','on','Collisions','off');
axis equal
axis([-0.8,0.8,-0.8,0.8,-0.8,0.8])
axis on
grid off
hold on
view([68,12]);

%% error generation for 6D pose sensor meaesurement 
noiseLevel_SO3 = 0.02/180*pi; % 0.02 degree (vive tracker accuracy reported in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5439658/)
noiseLevel_position = 0.0002;  % 0.2 millimeters
noiseType = 'G';

sensedFrameCToBaseLink = eye([4,4]);   
sensedFrameCToBaseLink(1:3,1:3) = addNoise(actualFrameCToRobotBase(1:3,1:3) ,noiseLevel_SO3,'right',noiseType);
sensedFrameCToBaseLink(1:3,4) = actualFrameCToRobotBase(1:3,4) + noiseLevel_position*randn(3,1);

%% compute optimal joint intials
optimalJointInitials = optimizeJointInitials(robotBasePoseToWorld, simulatedElbowFrameToBaseLink, sensedFrameCToBaseLink, upperArmLength, forearmLength, nominalJointInitial);
fprintf("\nthe actual joint angles are:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", acturalJointPosition_GroudTruth(1),acturalJointPosition_GroudTruth(2),acturalJointPosition_GroudTruth(3),acturalJointPosition_GroudTruth(4),acturalJointPosition_GroudTruth(5));

fprintf("\noptmized initial for joint angles:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", optimalJointInitials(1),optimalJointInitials(2),optimalJointInitials(3),optimalJointInitials(4),optimalJointInitials(5));

fprintf("\nnominal initial for joint angles:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", nominalJointInitial(1),nominalJointInitial(2),nominalJointInitial(3),nominalJointInitial(4),nominalJointInitial(5));


%% compute optimized joint
finalOptimizedJointAngles = optimizationForOptimalJointAngles(roboticArm, jointConfigStruct,sensedFrameCToBaseLink, optimalJointInitials);
fprintf("\noptimized joint angles:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", finalOptimizedJointAngles(1),finalOptimizedJointAngles(2),finalOptimizedJointAngles(3),finalOptimizedJointAngles(4),finalOptimizedJointAngles(5));


finalError =finalOptimizedJointAngles -  acturalJointPosition_GroudTruth;
fprintf("\njoint angle errors of epsilon 3 are (in deg):\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", epsilon3(1)/pi*180,epsilon3(2)/pi*180,epsilon3(3)/pi*180,epsilon3(4)/pi*180,epsilon3(5)/pi*180);
fprintf("\nfinal joint angle errors are (in deg):\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", finalError(1)/pi*180,finalError(2)/pi*180,finalError(3)/pi*180,finalError(4)/pi*180,finalError(5)/pi*180);


%% show robot joint axis and frames of two 6D pose sensor 
link6ToLink5AtHomeConfig = getTransform(roboticArm,jointConfigStruct,'link6','link5');

linkFrameInWorld = zeros(4,4,8);
linkFrameInWorld(:,:,1) = getTransform(roboticArm,jointConfigStruct,'base_link','world');

for i=1:7
    linkFrameInWorld(:,:,i+1) = getTransform(roboticArm,jointConfigStruct,['link' num2str(i)],'world');
end

robotJointAxisInLastParentLink = zeros(7,3);
for i=1:7
    robotJointAxisInLastParentLink(i,:) = roboticArm.Bodies{1, i+1}.Joint.JointAxis;    
end

% draw joint axis frame
axisLine = zeros([7,6]);
axisLength = 0.08;
for i=1:7
    startPoint = linkFrameInWorld(1:3,4,i+1);
    endPoint = linkFrameInWorld(:,:,i+1)*[robotJointAxisInLastParentLink(i,:),1]';
    delta = endPoint(1:3) - startPoint;
    %quiver3(startPoint(1),startPoint(2),startPoint(3),delta(1),delta(2),delta(3),axisLength,'m-','LineWidth',2);
    axisLine(i,:)=[startPoint',delta'];
end

% draw link line
plot3(axisLine(:,1),axisLine(:,2),axisLine(:,3),'Color','k','LineWidth',1);

% draw cylinder for joint
cylinderHeight = ones([1,7])*0.08;
cylinderAxialOffset = ones([1,7])*0;

cylinderRadius=0.01;
cylinderBin=1000;
endPlaneColor='m';
withEndPlane=1;
withCylinderBin=0;
cylinderFaceAlpha=1;
endPlaneFaceAlpha=1;
cylinderColor=['b','g','b','g','b','g','b'];

for i=1:7
    X1= axisLine(i,1:3 ) - axisLine(i,4:6)*cylinderHeight(i)/2 + cylinderAxialOffset(i); 
    X2= axisLine(i,1:3 ) + axisLine(i,4:6)*cylinderHeight(i)/2 + cylinderAxialOffset(i); 
    Cylinder(X1,X2,cylinderRadius,cylinderBin,cylinderColor(i),endPlaneColor,withEndPlane,withCylinderBin,cylinderFaceAlpha,endPlaneFaceAlpha);
end

baseLinkToTracker1 = [0.0, 0.0, -1.0, -0.039; -1.0, 0.0, 0.0, -0.113;   0.0, 1.0, 0.0, -0.0905; 0,0,0,1]; % in x_arm_gaze_control/data/transformationConfig/baseLinkToBaseTracker.yaml
tracker1ToWorld = linkFrameInWorld(:,:,1) /baseLinkToTracker1;  %linkFrameInWorld(:,:,1) * inv(baseLinkToTracker1);

link6ToTracker2AtHomeConfig = [0.0,-1.0, 0.0, 0.0; -1.0, 0.0, 0.0, -0.099;  0.0, 0.0, -1.0, 0.052; 0,0,0,1]; % zero joint position. In  x_arm_gaze_control/data/transformationConfig/calibLink8ToWristTracker.yaml
tracker2ToWorld = linkFrameInWorld(:,:,6) * link6ToLink5AtHomeConfig / link6ToTracker2AtHomeConfig;
compensationFrameCenteredAtWristInWorld= linkFrameInWorld(:,:,6) *link6ToLink5AtHomeConfig; % compensation frame is aligned with link6 at zero position and fixed with the link5(forearm)

% draw frame for tracker 1# and tracker 2#
tracker1AxisLength = 0.08;
quiver3(tracker1ToWorld(1,4),tracker1ToWorld(2,4),tracker1ToWorld(3,4),tracker1ToWorld(1,1)*tracker1AxisLength,tracker1ToWorld(2,1)*tracker1AxisLength,tracker1ToWorld(3,1)*tracker1AxisLength,'r-','LineWidth',2); %x
quiver3(tracker1ToWorld(1,4),tracker1ToWorld(2,4),tracker1ToWorld(3,4),tracker1ToWorld(1,2)*tracker1AxisLength,tracker1ToWorld(2,2)*tracker1AxisLength,tracker1ToWorld(3,2)*tracker1AxisLength,'g-','LineWidth',2); %y
quiver3(tracker1ToWorld(1,4),tracker1ToWorld(2,4),tracker1ToWorld(3,4),tracker1ToWorld(1,3)*tracker1AxisLength,tracker1ToWorld(2,3)*tracker1AxisLength,tracker1ToWorld(3,3)*tracker1AxisLength,'b-','LineWidth',2); %z

tracker2AxisLength = 0.08;
quiver3(tracker2ToWorld(1,4),tracker2ToWorld(2,4),tracker2ToWorld(3,4),tracker2ToWorld(1,1)*tracker2AxisLength,tracker2ToWorld(2,1)*tracker2AxisLength,tracker2ToWorld(3,1)*tracker2AxisLength,'r-','LineWidth',2); %x
quiver3(tracker2ToWorld(1,4),tracker2ToWorld(2,4),tracker2ToWorld(3,4),tracker2ToWorld(1,2)*tracker2AxisLength,tracker2ToWorld(2,2)*tracker2AxisLength,tracker2ToWorld(3,2)*tracker2AxisLength,'g-','LineWidth',2); %y
quiver3(tracker2ToWorld(1,4),tracker2ToWorld(2,4),tracker2ToWorld(3,4),tracker2ToWorld(1,3)*tracker2AxisLength,tracker2ToWorld(2,3)*tracker2AxisLength,tracker2ToWorld(3,3)*tracker2AxisLength,'b-','LineWidth',2); %z