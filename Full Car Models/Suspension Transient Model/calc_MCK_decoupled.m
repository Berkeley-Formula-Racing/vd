function [m,c,k] = calc_MCK(car,tirePos,tireVel)

for i=1:2 %calculate motion ratios by interp of data
    tireX = tirePos(i)*39.3701; %get current tire pos (m to in)
    MRArr(i) = lininterp1(car.MR_F(:,1),car.MR_F(:,2),tireX); %MR data within damping mat
end
for i=3:4 %calculate motion ratios by interp of data
    tireX = tirePos(i)*39.3701; %get current tire pos (m to in)
    MRArr(i) = lininterp1(car.MR_R(:,1),car.MR_R(:,2),tireX); %MR data within damping mat
end

c_c = car.c_compression;
c_r = car.c_rebound;
for i =1:4 %calculate damping coeff by interp of dampingcurves.mat data
    v = tireVel(i)*MRArr(i)*39.3701; %get current tire vel (m/s to in/s)    
    if v == 0
        cArr(i) = (c_c(2,2)-c_c(1,2))/(c_c(2,1)-c_c(1,1))*175.126835;
    elseif v > 0
        cArr(i) = lininterp1(c_c(:,1),c_c(:,2),v)/v*175.126835;
    else
        cArr(i) = lininterp1(c_r(:,1),-c_r(:,2),-v)/v*-175.126835;
    end
end

k_tf = car.k_tf;
k_tr = car.k_tr;
m = car.M;
m_f = m/10;
m_r = m_f;
k_rf = car.k_rf;
k_rr = car.k_rr;
Ix = car.Ixx;
Iy = car.Iyy;
a_1 = car.l_f;
a_2 = car.l_r;
b_1 = car.t_f/2;
b_2 = car.t_f/2;
w = car.t_f;

k_1 = car.k*MRArr(1)^2; %tire k = spring k * MR^2
k_2 = car.k*MRArr(2)^2;
k_3 = car.k*MRArr(3)^2;
k_4 = car.k*MRArr(4)^2;
c_1 = cArr(1)*MRArr(1)^2;
c_2 = cArr(2)*MRArr(2)^2;
c_3 = cArr(3)*MRArr(3)^2;
c_4 = cArr(4)*MRArr(4)^2;

m = diag([m Ix Iy m_f m_f m_r m_r]);

c11 = c_1+c_2+c_3+c_4;
c21 = b_1*c_1-b_2*c_2+b_1*c_3-b_2*c_4;
c12 = c21;
c31 = a_2*(c_3+c_4)-a_1*(c_1+c_2);
c13 = c31;
c22 = b_1^2*c_1+b_2^2*c_2+b_1^2*c_3+b_2^2*c_4;
c32 = a_1*b_2*c_2-a_1*b_1*c_1+a_2*b_1*c_3-a_2*b_2*c_4;
c23 = c32;
c33 = (c_1+c_2)*a_1^2+(c_4+c_3)*a_2^2;

c = [c11 c12 c13 -c_1 -c_2 -c_3 -c_4;
    c21 c22 c23 -b_1*c_1 b_2*c_2 -b_1*c_3 b_2*c_4;
    c31 c32 c33 a_1*c_1 a_1*c_2 -a_2*c_3 -a_2*c_4;
    -c_1 -b_1*c_1 a_1*c_1 c_1 0 0 0;
    -c_2 b_2*c_2 a_1*c_2 0 c_2 0 0;
    -c_3 -b_1*c_3 -a_2*c_3 0 0 c_3 0;
    -c_4 b_2*c_4 -a_2*c_4 0 0 0 c_4];

k11 = k_1+k_2+k_3+k_4;
k21 = b_1*k_1-b_2*k_2+b_1*k_3-b_2*k_4;
k12 = k21;
k31 = a_2*(k_3+k_4)-a_1*(k_1+k_2);
k13 = k31;
k22 = k_rf+k_rr+b_1^2*k_1+b_2^2*k_2+b_1^2*k_3+b_2^2*k_4;
k32 = a_1*b_2*k_1-a_1*b_2*k_2+a_2*b_1*k_3-a_2*b_2*k_4;
k23 = k32;
k42 = -b_1*k_1-k_rf/w;
k24 = k42;
k52 = b_2*k_2+k_rf/w;
k25 = k52;
k62 = -b_1*k_3-k_rr/w;
k26 = k62;
k72 = b_2*k_4+k_rr/w;
k27 = k72;
k33 = (k_1+k_2)*a_1^2+(k_1+k_2)*a_2^2;
k44 = k_1+k_tf+k_rf/w^2;
k55 = k_2+k_tf+k_rf/w^2;
k66 = k_3+k_tr+k_rr/w^2;
k77 = k_4+k_tr+k_rr/w^2;

k = [k11 k12 k13 -k_1 -k_2 -k_3 -k_4;
    k21 k22 k23 k24 k25 k26 k27;
    k31 k32 k33 a_1*k_1 a_1*k_2 -a_2*k_3 -a_2*k_4;
    -k_1 k42 a_1*k_1 k44 -k_rf/w^2 0 0;
    -k_2 k52 a_1*k_2 -k_rf/w^2 k55 0 0;
    -k_3 k62 -a_2*k_3 0 0 k66 -k_rr/w^2;
    -k_4 k72 -a_2*k_4 0 0 -k_rr/w^2 k77];