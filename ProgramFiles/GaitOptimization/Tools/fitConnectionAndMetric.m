%Stores connection, metric, and metric derivatives as interpolated
%functions that can be called for any 2D joint pair
function s = fitConnectionAndMetric(s);

%Get metric and grid data
metric_cell = s.metricfield.metric_eval.content.metric;
grid = s.grid.metric_eval;
r1 = grid{1};
r2 = grid{2};

%Fit metric and coriolis data to interpolant functions
metric_fun = matrixFit(r1,r2,metric_cell);
dmdr1_fun = matrixFit(r1,r2,s.coriolisfield.coriolis_eval.content.dM{1});
dmdr2_fun = matrixFit(r1,r2,s.coriolisfield.coriolis_eval.content.dM{2});

%Fit connection data to interpolant functions
A_fun = connectionFit(r1,r2,s.vecfield.eval.content.A_num);

%Put all functions in the same structure
funs.metric_fun = metric_fun;
funs.dmdr1_fun = dmdr1_fun;
funs.dmdr2_fun = dmdr2_fun;
funs.A_fun = A_fun;

s.funs = funs;

end


%Fits connection data to a 3x2 interpolated function
function myfun = connectionFit(r1,r2,A_cell)

    %Get interpolation data for each entry in the matrix
    fs = {};
    for i = 1:6
        ai = A_cell{i};
        fi = fit([r1(:),r2(:)],ai(:),'linearinterp');
        fs{i} = fi;
    end
   
    %Compose matrix of individual functions
    myfun = @(r1,r2) [fs{1}(r1,r2),fs{4}(r1,r2);...
        fs{2}(r1,r2),fs{5}(r1,r2);...
        fs{3}(r1,r2),fs{6}(r1,r2)];
end

%Fits metric/coriolis data to a 2x2 interpolated function
function myfun = matrixFit(r1,r2,m_cell)

    %Get interpolation data for each entry in the matrix
    fs = {};
    for i = 1:4
        mii = m_cell{i};
        fi = fit([r1(:),r2(:)],mii(:),'linearinterp');
        fs{i} = fi;
    end

    %Compose matrix of individual functions
    myfun = @(r1,r2) [fs{1}(r1,r2),fs{2}(r1,r2);...
        fs{3}(r1,r2),fs{4}(r1,r2)];
end
    