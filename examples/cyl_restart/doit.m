format compact

load 'c.b'; b=c; load 'c.a'; a=c; 
skp=1; sa=5000+skp; sb=skp;
for k=2:3;
   n=length(a); a1=a(sa:n,k);
   m=length(b); b1=b(sb:m,k);
   m=min(length(b1),length(a1)); a2=a1(1:m); i=1:m; b2=b1(1:m);
   d=a2-b2; max(abs(d))
   d(1:7)'
%  if k==2; pause(3); else; plot(i,d,'b-'); end;
end;

dc = d;

load 'p.b'; b=p; load 'p.a'; a=p; 
skp=1; sa=5000+skp; sb=skp;
for k=2:3;
   n=length(a); a1=a(sa:n,k);
   m=length(b); b1=b(sb:m,k);
   m=min(length(b1),length(a1)); a2=a1(1:m); i=1:m; b2=b1(1:m);
   d=a2-b2; max(abs(d))
   d(1:7)'
%  if k==2; pause(3); else; plot(i,d,'b-'); end;
end;

dp = d;


