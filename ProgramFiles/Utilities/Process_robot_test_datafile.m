function Process_robot_test_datafile

load SerpenoidGranularLocomotion2016

FAx1 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,3));
FAy1 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,4));
FAtheta1 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,5));
FAx2 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,6));
FAy2 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,7));
FAtheta2 = scatteredInterpolant(data_log(:,1),data_log(:,2),data_log(:,8));

[alpha1,alpha2] = ndgrid(linspace(-3.1,3.1,51));

Ax1 = FAx1(alpha1,alpha2);
Ax2 = FAx2(alpha1,alpha2);
Ay1 = FAy1(alpha1,alpha2);
Ay2 = FAy2(alpha1,alpha2);
Atheta1 = FAtheta1(alpha1,alpha2);
Atheta2 = FAtheta2(alpha1,alpha2);


save('SerpenoidGranularConnection2016','alpha1','alpha2','Ax1','Ax2','Ay1','Ay2','Atheta1','Atheta2')