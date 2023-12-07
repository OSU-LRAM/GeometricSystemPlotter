function plotCurrentSalpGeometry_GeneralCurvature(ax,base_anim_geom,salp,pose,forcedef,t)
    
    %Get transform that moves all salp elements to current pose in world frame
    poseTransform = vec_to_mat_SE2(pose);

    %Plot each link
    B = base_anim_geom.links.points;
    inWorldSalp = poseTransform*B;
    fill(ax,inWorldSalp(1,:),inWorldSalp(2,:),base_anim_geom.links.color);

    %Plot each thruster
    thrusterAngles = salp.geometry.thrusters.angles;
    for thruster_ind = 1:salp.geometry.thrusters.nThrusters

        rotateThrusterTransform = vec_to_mat_SE2([0,0,thrusterAngles(thruster_ind)]);
        toLinkTransform = base_anim_geom.thrusters.transforms{thruster_ind};

        currentThruster = base_anim_geom.thrusters.points{thruster_ind};
        rotatedThruster = rotateThrusterTransform*currentThruster;
        onLinkThruster = toLinkTransform*rotatedThruster;
        inWorldThruster = poseTransform*onLinkThruster;
        fill(ax,inWorldThruster(1,:),inWorldThruster(2,:),base_anim_geom.thrusters.color);

    end

    %Plot each force jet
    maxThrustAmp = forcedef.maxThrust;
    thrustAmps = forcedef.amplitudes(t);
    for thruster_ind = 1:salp.geometry.thrusters.nThrusters

        rotateThrusterTransform = vec_to_mat_SE2([0,0,thrusterAngles(thruster_ind)]);
        toLinkTransform = base_anim_geom.thrusters.transforms{thruster_ind};

        currentJetOut = base_anim_geom.out_force_jets.points{thruster_ind};
        currentJetIn = base_anim_geom.in_force_jets.points{thruster_ind};

        jetScaling = thrustAmps(thruster_ind)/maxThrustAmp;
        scaledJetOut = scaleThrusterForcesAnimation(currentJetOut,jetScaling);
        scaledJetIn = scaleThrusterForcesAnimation(currentJetIn,jetScaling);

        rotatedJetOut = rotateThrusterTransform*scaledJetOut;
        rotatedJetIn = rotateThrusterTransform*scaledJetIn;

        onLinkJetOut = toLinkTransform*rotatedJetOut;
        onLinkJetIn = toLinkTransform*rotatedJetIn;

        inWorldJetOut = poseTransform*onLinkJetOut;
        inWorldJetIn = poseTransform*onLinkJetIn;

        fill(ax,inWorldJetOut(1,:),inWorldJetOut(2,:),base_anim_geom.out_force_jets.color);
        fill(ax,inWorldJetIn(1,:),inWorldJetIn(2,:),base_anim_geom.in_force_jets.color);

    end

end