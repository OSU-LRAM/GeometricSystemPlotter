function efficiencyplot(vp)


syslist = {'links','cc','serp'};

f = figure(1337);
clf(f,'reset')
cost = 'pathlength';
for idx = 1:numel(syslist)
    line([vp.(syslist{idx}).(cost).normalizedPeriod{:}],[vp.(syslist{idx}).(cost).netDisp{:}])
end
lims = axis;
lims(1) = 0;
lims(3) = 0;
axis(lims);

f = figure(1338);
clf(f,'reset')
cost = 'covaccel';
for idx = 1:numel(syslist)
    line([vp.(syslist{idx}).(cost).normalizedPeriod{:}],[vp.(syslist{idx}).(cost).netDisp{:}])
end
lims = axis;
lims(1) = 0;
lims(3) = 0;
axis(lims);