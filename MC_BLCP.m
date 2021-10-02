clc; clear all;


lambda = 0.1;
nB = 10;
iterations = 2000;
xt = 150;
t_vec = linspace(0,100,20);
R = 100;
L = 1000;

for j = 1:length(t_vec)
    F(j) = 0;
    t = t_vec(j);
    
    for i = 1:iterations
        N_points = nB;
        r_vec = R * (rand(1,nB));
        theta_vec = 2*pi *rand(1,nB);
        d_vec = abs(xt * cos(theta_vec) - r_vec);
        
        for k = 1:nB
            N_points(k) = poissrnd(2*L*lambda);
            near(k) = min(L*rand(1,N_points(k)));
        end
        
        if min(sqrt(d_vec.^2 + near.^2)) < t
            F(j) = F(j) + 1/iterations;
        end
    end
  j       
end

plot(t_vec, F,'*')
hold on