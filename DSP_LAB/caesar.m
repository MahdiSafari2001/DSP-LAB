function coded = caesar(v,n)
coded=zeros([length(v) 1]);
for i=1:length(v)
    if double(v(i))+n>126
        coded(i)=(mod(double(v(i))+n-127,126-31)+32);
    elseif  double(v(i))+n<=126 & double(v(i))+n>=32
        coded(i)=(double(v(i))+n);
    else
        coded(i)=126-mod(31-double(v(i))+n,95);
    end
end
coded=(char(coded))';
end