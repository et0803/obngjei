%% show the collision model of the ARA
clc;
clear;
roboticArmModelFileName='./x_arm_aligned.urdf'; 

roboticArm = importrobot(roboticArmModelFileName);
jointConfigStruct =randomConfiguration(roboticArm);

commandedJointPosition =[0.71, 0.58, 1.04, 1.91, -0.99];
jointErrorLevelForSimulation = 10/180*pi;
jointErrorsForSimulation = jointErrorLevelForSimulation*randn(1,5);
acturalJointPosition = commandedJointPosition - jointErrorsForSimulation;

for i=1:size(jointConfigStruct,2)
    if i==1         % robot base to world, for a better plot view
        jointConfigStruct(i).JointPosition = pi/2;  
    else
        if i<7
            jointConfigStruct(i).JointPosition = acturalJointPosition(i-1);  % 5 joint for 7-DoF arm
        else
            jointConfigStruct(i).JointPosition =0;  % joint for fingers
        end
    end
end

robotBasePoseToWorld  = getTransform(roboticArm,jointConfigStruct,'base_link','world');
actualFrameCToRobotBase  = getTransform(roboticArm,jointConfigStruct,'link6','base_link');
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
noiseLevel_position = 0.00002;  % 0.02 millimeters
noiseType = 'G';

sensedFrameCToBaseLink = eye([4,4]);   
sensedFrameCToBaseLink(1:3,1:3) = addNoise(actualFrameCToRobotBase(1:3,1:3) ,noiseLevel_SO3,'right',noiseType);
sensedFrameCToBaseLink(1:3,4) = actualFrameCToRobotBase(1:3,4) + noiseLevel_position*randn(3,1);

%% compute initials for joint error
%jointErrorInitials = twoPhaseOptimizeJointErrorInitials(robotBasePoseToWorld, actualFrameCToRobotBase, upperArmLength, forearmLength, commandedJointPosition);
jointErrorInitials = twoPhaseOptimizeJointErrorInitials(robotBasePoseToWorld, sensedFrameCToBaseLink, upperArmLength, forearmLength, commandedJointPosition);

%% compute optimized joint error
mu = 4; % 5mm and 1 degree contributes to equivalent error.
optimizedJointError = localOptimizationForOptimalJointError(roboticArm, jointConfigStruct,sensedFrameCToBaseLink, commandedJointPosition, jointErrorInitials, mu);

%% print joint errors
fprintf("\nsimulated joint errors:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", jointErrorsForSimulation(1),jointErrorsForSimulation(2),jointErrorsForSimulation(3),jointErrorsForSimulation(4),jointErrorsForSimulation(5));

fprintf("optimized initials for joint errors:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", jointErrorInitials(1),jointErrorInitials(2),jointErrorInitials(3),jointErrorInitials(4),jointErrorInitials(5));

fprintf("optimized joint errors:\n");
fprintf("%.6f, %.6f, %.6f, %.6f, %.6f\n", optimizedJointError(1),optimizedJointError(2),optimizedJointError(3),optimizedJointError(4),optimizedJointError(5));


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