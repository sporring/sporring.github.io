function J = quick_linear_diffusion(I, time, order)

v = 2.^[0:order];
invA = inv((ones(order+1,1)*v).^([0:order]'*ones(1,order+1)));

J = zeros(size(I));
for i = 1:size(invA,1)
  j = v(size(v,2)-i+1);
  J = J+invA(i,1)*linear_diffusion(I,[0,linspace(time/j,time,j)]);
end
