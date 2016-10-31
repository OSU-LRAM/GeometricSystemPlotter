function [p,t]=distmeshnd_hypercube(fdist,fh,h,box,fix,display_fig_handle,varargin)
%DISTMESHND N-D Mesh Generator in a hypercubic box. More specifically, this
%code assumes that the distance function can return a second output giving
%the component-wise distances from the target region to the point
%   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%      P:           Node positions (NxNDIM)
%      T:           Triangle indices (NTx(NDIM+1))
%      FDIST:       Distance function
%      FH:          Edge length function
%      H:           Smallest edge length
%      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
%      FIX:         Fixed node positions (NFIXxNDIM)
%	   display_fig_handle: handle of figure to display into
%      FDISTPARAMS: Additional parameters passed to FDIST
%
%   Example: Unit ball
%      dim=3;
%      d=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
%
%   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
%   Modified 2011 by Ross L. Hatton





if ~isempty(display_fig_handle)
	
	% Ensure figure existence
	f = figure(display_fig_handle);
	clf(f)
	
	ax = axes('parent',f);
	
else
	
	f = NaN;
	
end


dim=size(box,2);
ptol=.001; ttol=.1; Fscale=1+.4/2^(dim-1); deltat=.2; geps=0.001*h; deps=sqrt(eps)*h;

% 1. Create initial distribution in bounding box (note, in the original
% code, this was a grid and could jam if placed in a square domain

% Build the sequence of "heights" of an n-d regular simplex
hsequence = zeros(dim,1);
hsequence(1) = h;
for i = 2:dim
	hsequence(i) = sqrt(i^2 - 1)/i * hsequence(i-1);
end

if dim==1
  p=(box(1):h:box(2))';
else
  cbox=cell(1,dim);
  cgrid = cell(1,dim); % Hold a grid to determine odd and even points along each dimension
  for ii=1:dim
    cbox{ii}=box(1,ii):hsequence(ii):box(2,ii);
	cgrid{ii} = 1:length(cbox{ii});
  end
  pp=cell(1,dim);
  pgrid = cell(1,dim);
  [pp{:}]=ndgrid(cbox{:});
  [pgrid{:}] = ndgrid(cgrid{:});
  p=zeros(numel(pp{1}),dim);
  for ii=1:dim
    p(:,ii)=pp{ii}(:);
  end
  for ii = 1:(dim-1)
	p(:,ii) = p(:,ii) + (hsequence(ii)/2 * ~mod(pgrid{ii+1}(:),2)); 
	% Add an offset to arrange the points in a face-centered lattice
  end
  
end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p,varargin{:})<geps,:);
r0=feval(fh,p);
p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
N=size(p,1);

wb = waitbar2a(0,'Mesh completion');

count=0;
pold=inf;
while 1
  % 3. Retriangulation by Delaunay
  if max(sqrt(sum((p-pold).^2,2)))>ttol*h
    pold=p;
    warning('off','MATLAB:delaunayn:DuplicateDataPoints') % Suppress duplicate warnings
    t=delaunayn(p);
    pmid=zeros(size(t,1),dim);
    for ii=1:dim+1
      pmid=pmid+p(t(:,ii),:)/(dim+1);
    end
    t=t(feval(fdist,pmid,varargin{:})<-geps,:);

	
	% 4. Describe each edge by a unique pair of nodes
	pair=zeros(0,2);
    localpairs=nchoosek(1:dim+1,2);
    for ii=1:size(localpairs,1)
      pair=[pair;t(:,localpairs(ii,:))];
    end
    pair=unique(sort(pair,2),'rows');

	% 5. Graphical output of the current mesh
	if ~isnan(f)
		if dim==2
		  trimesh(t,p(:,1),p(:,2),zeros(N,1),'parent',ax,'EdgeColor',[235 14 30]/255)
		  view(ax,2),axis(ax,'equal'),axis(ax,'off'),drawnow
		elseif dim==3
		  if mod(count,5)==0
			simpplot(p,t,'p(:,2)>0');
			title(['Retriangulation #',int2str(count)])
			drawnow
		  end
		else
		  disp(sprintf('Retriangulation #%d',count))
		end
	end
    count=count+1;
  end

  % 6. Move mesh points based on edge lengths L and forces F
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));
  hbars=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
  L0=hbars*Fscale*(sum(L.^dim)/sum(hbars.^dim))^(1/dim); % note that ^(1/dim) corresponds to a square root in the 2d version
  F=max(L0-L,0);
  Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
  dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...
                 ones(size(pair,1),1)*[1:dim,1:dim], ...
                 Fbar,N,dim));
  dp(1:size(fix,1),:)=0;
  p=p+deltat*dp;

  % 7. Bring outside points back to the boundary
  [d, d_component, closest_face] = feval(fdist,p,varargin{:}); % Get the full and component-wise distances of points to the hypercube
  ix=d>0;
  n_outside = sum(ix);
  kickfactor = 1;
  while any(ix)

	  p(ix,:)=p(ix,:)-( d_component(ix,:).*(d_component(ix,:)>0).*closest_face(ix,:));
	  d=feval(fdist,p,varargin{:});
	  ix=d>0;
	  
	  if sum(ix) >= n_outside
		  keyboard
		  %kickfactor = kickfactor*(1+1e-1); % Exponentially growing kickfactor to bring points inside
	  end
	  n_outside = sum(ix);
	  
  end

  % 8. Termination criterion
  %maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)))<ptol*h;
  maxdp = max(sqrt(sum(deltat*dp(d<-geps,:).^2,2)));
  if (maxdp/h)<ptol
	  break; 
  end
  
  waitbar2a(ptol/(maxdp/h),wb);
  
end

close(wb)

end
