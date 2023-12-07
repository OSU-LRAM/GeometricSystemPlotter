%Scales down jet shape to match thrust output
function scaledShape = scaleThrusterForcesAnimation(inShape,scaleFactor)

    scaledShape = inShape*scaleFactor;
    scaledShape(3,:) = ones(1,size(scaledShape,2));

end