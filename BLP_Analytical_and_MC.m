clc; clear all;


lambda = 0.1;
nB = 20;
iterations = 2000;
xt = 500;
t_vec = linspace(0,50,20);
R = 100;

for j = 1:length(t_vec)
    F(j) = 0;
    t = t_vec(j);

    %%
    % Monte Carlo Simulations
    
    for i = 1:iterations
        N_points = nB;
        r_vec = R * (rand(1,nB));
        theta_vec = 2*pi *rand(1,nB);
        d_vec = abs(xt * cos(theta_vec) - r_vec);
        if min(d_vec) < t
            F(j) = F(j) + 1/iterations;
        end
    end
    
  %%
    %Analytical Expression

    if xt + t < R
        A_a(j) = 2*pi*t;
    elseif (xt - t < R) && (xt + t >= R)
        A_a(j) = 2* (pi *t - xt * sqrt(1 - ((R - t)/xt)^2) + (R - t) * acos((R - t)/xt));
    elseif xt - t >= R
        A_a(j) = 2* (pi *t - xt * (sqrt(1 - ((R - t)/xt)^2) -  ...
            sqrt(1 - ((R + t)/xt)^2)) + (R - t) * acos((R - t)/xt) - ...
            (R + t) * acos((R + t)/xt));
    end
    
    F_a(j) = 1 - ((2 * pi * R - A_a(j))/(2 * pi * R))^nB;
         
end


plot(t_vec, F,'*')
hold on
plot(t_vec, F_a)