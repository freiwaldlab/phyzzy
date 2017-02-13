function [ v, lambda ] = dpss( len, bw )
%DPSS copies the eponymous method from the signal processing toolbox
%   in other words, rtns eigenvectors of M(t,t') = modified sinc of t-t'
rtnNum = round(2*bw);
bw = bw/len;
M = zeros(len);
for i = 1:len
  for j = 1:len
    M(i,j) = sin(2*pi*bw*(i-j))/(pi*(i-j));
  end
end
for i = 1:len
  M(i,i) = 2*bw;
end
[v,lambda] = eig(M,'vector');
v = fliplr(v);
v = v(:,1:rtnNum);
lambda = flipud(lambda);
lambda = lambda(1:rtnNum);
end

