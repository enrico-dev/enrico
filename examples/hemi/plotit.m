
load circ.dat

%%% aviobj = avifile('mv1.avi','compression','None');
%%% fig=figure;

x = [-1 -1]; y=[-.25 .25]; z=[.13 .13];


for ii=001:161
    iname=sprintf('part%05d.3D',ii);
    part = load(iname);

    h=plot3(x,y,z,'k-',...
            circ(:,1),circ(:,2),circ(:,3),'r-',...
            circ(:,2),circ(:,1),circ(:,3),'r-',...
            part(:,1),part(:,2),part(:,3),'k.','markersize',3);
    axis equal;axis ([-1 1.5 -1.1 1.1 0 .7]);view([90,00]);
    drawnow

%   F = getframe(fig);
%   aviobj = addframe(aviobj,F);
end
%%% aviobj = close(aviobj);
