function deltaPose_twist = twist_from_pose_diff(currentPose,targetPose)
    positionDiff = targetPose(1:3,4) - currentPose(1:3,4); %旋转矩阵对于一个要变换的向量，都是premultify
    currentOrientation = currentPose(1:3,1:3);
    targetOrientation = targetPose(1:3,1:3);
    targetToCurrentOrientationInAxang = rotm2axang(currentOrientation'*targetOrientation)';
    orientationDiff = currentOrientation * targetToCurrentOrientationInAxang(1:3)* targetToCurrentOrientationInAxang(4);
    deltaPose_twist=zeros(6,1);
    %twist 为6个元素的向量，依照matlab的jacobian矩阵顺序，前面是角速度，后面的距离。
    deltaPose_twist(1) = orientationDiff(1);
    deltaPose_twist(2) = orientationDiff(2);
    deltaPose_twist(3) = orientationDiff(3);   
    deltaPose_twist(4) = positionDiff(1);
    deltaPose_twist(5) = positionDiff(2);
    deltaPose_twist(6) = positionDiff(3);
end

