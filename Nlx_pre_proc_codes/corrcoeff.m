
function c = corrcoeff(x,y)
%syntax :
% function c = corrcoeff(x,y)
%
% -- David B. Omer 2012*
%tic
[strn, coln]= size(y);
[strnx, colnx]= size(x);

if( strnx~=strn | colnx~=coln | max(strn,coln) == 1)
  error('matrix dimensions must agree');
end;

meanx=mean(x);
meany=mean(y);
x=x-meanx(ones(size(x,1),1),:);
y=y-meany(ones(size(y,1),1),:);
stdx = std(x);
stdy = std(y);

c=zeros(size(stdy));

k = find(stdx == 0 | stdy == 0 );

stdx(k) = k;
stdy(k) = k;
  
if(strn == 1 ),
  strn = coln;
end;  
c = sum(x.*y)./(stdx.*stdy.*(strn-1));
c(k)=k*0-2;

%t = toc;
%disp(['corrcoeff() took: ',num2str(t),' sec']); 

