%Calls the force definition to get force vectors acting at each thruster
%location in the salp body frame at that point
function f_vec = getForceVectors(salp,forcedef,t)

    nThrusters = salp.geometry.thrusters.nThrusters;
    f_vec = zeros(nThrusters*3,1);
    amplitudes = forcedef.amplitudes(t);
    thrusterAngles = salp.geometry.thrusters.angles;

    for i = 1:nThrusters
        thrusterAngle = thrusterAngles(i);
        amp = amplitudes(i);
        f_vec(3*i-2:3*i) = [amp*cos(thrusterAngle);amp*sin(thrusterAngle);0];
    end

end