function flatprint(f,destination,printoptions)
%Print out the figure as a flattened pdf or eps

%make sure all my scripts are available
setenv('PATH','/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/usr/texbin:/opt/local/bin:/opt/local/sbin:/Users/rlhatton/scripts')

%make sure that the tempdir for building the composite image exists and is
%empty

recycle('on'); %move files to trash rather than rm-ing them

if ~exist('~/Temp/flatprint/','dir')
	mkdir('~/Temp/flatprint/');
else
	delete('~/Temp/flatprint/*')
end

%Set the figure to print with no background
ihc = get(f,'InvertHardcopy');
c = get(f,'Color');

%Print the figure to an eps file
print(f,'-depsc2','-painters','-loose','~/Temp/flatprint/flatprint',printoptions{:});

%Flatten the output, and send the desired copy to the destination
unix(['cd ~/Temp/flatprint/; epsflatten flatprint.eps; cp flatprint_flat.' destination(end-2:end) ' ' destination]);

%restore the figure properties
set(f,'InvertHardCopy',ihc,'Color',c)