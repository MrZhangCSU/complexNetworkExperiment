clear;clc;
filename = 'degree.csv';
M= csvread(filename,1,0,[1 0 1 999]);
a=[];b=[];count=0;
for i=1:1000
    a(i)=i;
end
for i=1:1000
    for j=1:1000
        if (M(j) == i)
            count=count+1;
        end
    end
    b(i)=count/1000;
    count = 0;
end

loglog(a,b,'-s');
grid on;