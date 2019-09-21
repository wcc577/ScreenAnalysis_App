function limits=autoScale(Lower,Upper,data)
%Returns upper and lower limits for a data set. 


[m n]=size(data);

lineData=reshape(data,1,m*n);
lineData=sort(lineData);

[y,x]=hist(lineData,5001);

nPix=0;
count1=0;
while nPix<(Lower*m*n)
    count1=count1+1;
    nPix=nPix+y(count1);
end
limits(1)=x(count1);

nPix=0;
count2=0;
while nPix<(Upper*m*n)
    count2=count2+1;
    nPix=nPix+y(count2);
end

limits(2)=x(count2);

if limits(1)==limits(2)
    error('Lower Limit and Upper Limit results from autoScale are the same.  Consider Revising the Lower and Upper input values.')
    
end