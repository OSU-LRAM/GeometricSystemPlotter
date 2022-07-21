function hessianfourier = evaluate_hessian_fourier(f,s,n,dimension,direction,invert)

y = path_from_fourier(f,n,dimension);
y = y(1:end-1,:); % Remove the last points because path_from_fourier returns self-connected gait

hessian=evaluate_hessian(y,s,n,dimension,direction,invert);

chy=chy_generator(f,n,dimension);
dy = [];

for i = 1:dimension
    dy=blkdiag(dy,chy{i});
end

hessian.disp=cell2mat(hessian.disp);
hessian.stroke=cell2mat(hessian.stroke);
% hessianstroke=hessianstrokecalculator(y,chy,s,n,dimension,totalstroke);
% hessiandisp=numelhessdisp(y,s,n,dimension,direction);

hessianfourier.disp=dy*hessian.disp*dy.';
hessianfourier.stroke=dy*hessian.stroke*dy.';

% hessianfourier = struct();
% num_f = numel(f)-dimension;
% hessianfourier.stroke = zeros(num_f,num_f);
% hessianfourier.disp = zeros(num_f,num_f);
% jacobian=evaluate_jacobian_fourier(f,s,n,dimension,direction);
% df = zeros(size(f));
% for i = 1:size(f,1) - 1
%     for j = 1:size(f,2)
%         df(i,j) = 1e-6;
%         djacobian = evaluate_jacobian_fourier(f+df,s,n,dimension,direction);
%         temphessstroke = reshape((djacobian.stroke-jacobian.stroke)/1e-6,[num_f 1]);
%         temphessdisp = reshape((djacobian.disp-jacobian.disp)/1e-6,[num_f 1]);
%         hessianfourier.stroke((j-1)*(size(f,1) - 1)+i,:) = temphessstroke.';
%         hessianfourier.disp((j-1)*(size(f,1) - 1)+i,:) = temphessdisp.';
%     end
% end
% 
% djacobian = [];

end