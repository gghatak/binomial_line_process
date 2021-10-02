clc; clear all;

x_t = 150;
R = 100;
lambda = 0.1;
nB = 10;

r_vec = linspace(-R,R, 500);
theta_vec = linspace(0, pi, 500);
t_vec = linspace(eps, 100, 50);

for q = 1:length(t_vec)
    r = x_t;
    t = t_vec(q);
    x = linspace(eps, pi ,100);
    y1 = (sin(x)).^2.*(r * cot(x) .* csc(x) + t * (csc(x)).^2);
    y2 = (sin(x)).^2.*(r * cot(x) .* csc(x) - t * (csc(x)).^2);
    
    y1_bound = min(R, max(-R , y1));
    y2_bound = min(R, max(-R , y2));
    
    Area(q) = trapz(x,y1_bound) - trapz(x,y2_bound);
    
    
    
for i = 1:length(theta_vec)
    theta = theta_vec(i);
    for j = 1:length(r_vec)
        
        r = r_vec(j);
        
        x_1 = (1/(1 + (cot(theta))^2)) * (r * cot(theta) * csc(theta) + x_t + ...
            sqrt((-r^2 * (csc(theta))^2 + 2 * r * x_t * csc(theta) * cot(theta)+...
            (t^2 - x_t^2) * (cot(theta))^2 + t^2))); 

        x_2 = (1/(1 + (cot(theta))^2)) * (r * cot(theta) * csc(theta) + x_t - ...
            sqrt((-r^2 * (csc(theta))^2 + 2 * r * x_t * csc(theta) * cot(theta)+...
            (t^2 - x_t^2) * (cot(theta))^2 + t^2)));

        y_1 = (r - x_1 * cos(theta))/sin(theta);

        y_2 = (r - x_2 * cos(theta))/sin(theta);
        
        
        
        if isnan(x_1) || isnan(x_2) || ~isreal(x_1) || ~isreal(x_2)
            C(j) = 0;
        else
            C(j) = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2 );
        end
        
        
    end
    Temp(i) = trapz(r_vec, exp( - lambda * C));
end
Final1(q) = ((1/(2* pi * R))*(trapz(theta_vec, Temp)))^nB;

%Final(q) = (Final1(q)^nB - (2 * pi * R - Area(q))^nB )/(Final1(q) - (2 * pi * R - Area(q)));

q
end
plot(t_vec, 1 - Final1)
hold on