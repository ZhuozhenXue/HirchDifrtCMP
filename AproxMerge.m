function f = AproxMerge(X,e)
%This function can merge similar values from the original set of data into their original positions
%e represents the allowable error
S=size(X);
X=round(X,4);
    
for i=1:S(2)-1
    %Each column of elements only discards similar elements after itself
    for j=1:S(2)-i
        %The similar standards are:
        if norm(X(:,i)-X(:,i+j))<=e && norm(X(:,i))~=0
        X(:,i+j)=0;
        end
    end
end

f=X;
end
