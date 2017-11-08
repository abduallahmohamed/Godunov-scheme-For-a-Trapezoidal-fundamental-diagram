%We want to solve the classical LWR PDE on a grid, using the Godunov scheme.
%Consider a Trapezoidal fundamental diagram, defined by ?(k)=k v for k<kc1 ,
%?(k)=kc1 v for k>kc1 and k<kc2  and ?(k)=-w(k-kmax) for k>kc2 ,
%with parameters kc1=0.02/m, kmax=0.2/m, v=30 m/s and w=5 m/s

kc1 = 0.02;
kmax = 0.2;
v = 30;
w = 5;
kc2 = ((kc1*v)/(-1*w))+ kmax;
default_k = 0;
%Draw TPFD
figure;
density_axis = 0:0.0001:kmax;
tpfd =0:0.0001:kmax;
for j = 1:length(density_axis)
    tpfd(j) =  trapezoidal_fundamental_diagram(density_axis(j),kc1,kc2,v,w,kmax);
end
plot(density_axis,tpfd);
title('Trapezoidal fundamental diagram');
xlabel('Density k');
ylabel('Flow ?(k)');
saveas(gcf,'TPFD.png')

%Draw Demand and Supply
figure;
demnd =0:0.0001:kmax;
for j = 1:length(density_axis)
    demnd(j) =  demand(density_axis(j),kc1,kc2,v,w,kmax);
end
plot(density_axis,demnd);
title('Demand');
xlabel('Density k');
ylabel('Demand D(k)');
saveas(gcf,'Demand.png')

figure;
supp =0:0.0001:kmax;
for j = 1:length(density_axis)
    supp(j) =  supply(density_axis(j),kc1,kc2,v,w,kmax);
end
plot(density_axis,supp);
title('Supply');
xlabel('Density k');
ylabel('Supply S(k)');
saveas(gcf,'Supply.png')

%Godunov scheme
cell_width = 100; %Width of each cell in x, space step
delta_x = cell_width;
cell_x = 10; %Number of cells
delta_t = 3; %Minimal time step that satisifies CFL condition
time_scale = 60; %Width of time
cell_t = time_scale / delta_t; %Number of cells in time

k = zeros(cell_x,cell_t); %Grid
for x = 1:cell_x
    k(x,1) = density_mapper(x*cell_width);
end

for t = 2:cell_t
    for x = 1:cell_x
        %??.?(?+??,?)=??.?(?,?)+??.(min(D(k(t, ?- ??)),S(k(t,?))?min(D(k(t," " ?)),S(k(t,?+??)))
        
        if x-1 >= 1
            k_t_1_x_1 = k(x-1,t-1);
        else
            k_t_1_x_1 = default_k;
            
        k_t_1_x = k(x,t-1);
        
        end
        if x+1 <= cell_x
            k_t_1_x_p_1 = k(x+1,t-1);
        else
            k_t_1_x_p_1 = default_k;
        end
        
        d0 = demand(k_t_1_x_1,kc1,kc2,v,w,kmax);
        s1 = supply(k_t_1_x,kc1,kc2,v,w,kmax);
        d1 = demand(k_t_1_x,kc1,kc2,v,w,kmax);
        s2 = supply(k_t_1_x_p_1,kc1,kc2,v,w,kmax);

        
    k(x,t) =  k_t_1_x + (delta_t/delta_x)*(min(d0,s1)-min(d1,s2));
    end
end

%a.	Initial condition: k(0,x)=0.04 if 0<x<300, k(0,x)=0.1/m 
%for 300<x<500, and k(0,x)=0.02/m for x>500
figure;
%surf(k);
plot(k);
colormap('hot');
title('Godunov Scheme with initial condition only');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G1.png')

%b.	Upstream boundary demand: d(t)=0.3/s if 0s<t<40s, and d(t)=0.5 for 40s<t<60s
k = zeros(cell_x,cell_t); %Grid
for x = 1:cell_x
    k(x,1) = density_mapper(x*cell_width);
end

for t = 2:cell_t
    for x = 1:cell_x
        %??.?(?+??,?)=??.?(?,?)+??.(min(D(k(t, ?- ??)),S(k(t,?))?min(D(k(t," " ?)),S(k(t,?+??)))
        
        if x-1 >= 1
            k_t_1_x_1 = k(x-1,t-1);
        else
            k_t_1_x_1 = default_k;
            
        k_t_1_x = k(x,t-1);
        
        end
        if x+1 <= cell_x
            k_t_1_x_p_1 = k(x+1,t-1);
        else
            k_t_1_x_p_1 = default_k;
        end
        d0 = upstream_boundary(t*delta_t);
        d0_tmp = demand(k_t_1_x_1,kc1,kc2,v,w,kmax);
        d0 = min(d0,d0_tmp);
        
        s1 = supply(k_t_1_x,kc1,kc2,v,w,kmax);
        
        d1 = upstream_boundary(t*delta_t);
        d1_tmp = demand(k_t_1_x,kc1,kc2,v,w,kmax);
        d1 = min(d1,d1_tmp);
        
        s2 = supply(k_t_1_x_p_1,kc1,kc2,v,w,kmax);

        
    k(x,t) =  k_t_1_x + (delta_t/delta_x)*(min(d0,s1)-min(d1,s2));
    end
end

figure;
%surf(k);
plot(k);
colormap('hot');
title('Godunov Scheme with initial condition and upstream boundary demand');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G2.png')


%c.	Donwstream boundary supply: s(t)=0.6/s if 0s<t<20s, and s(t)=0.05/s for 20s<t<60s
k = zeros(cell_x,cell_t); %Grid
for x = 1:cell_x
    k(x,1) = density_mapper(x*cell_width);
end

for t = 2:cell_t
    for x = 1:cell_x
        %??.?(?+??,?)=??.?(?,?)+??.(min(D(k(t, ?- ??)),S(k(t,?))?min(D(k(t," " ?)),S(k(t,?+??)))
        
        if x-1 >= 1
            k_t_1_x_1 = k(x-1,t-1);
        else
            k_t_1_x_1 = default_k;
            
        k_t_1_x = k(x,t-1);
        
        end
        if x+1 <= cell_x
            k_t_1_x_p_1 = k(x+1,t-1);
        else
            k_t_1_x_p_1 = default_k;
        end
        d0 = upstream_boundary(t*delta_t);
        d0_tmp = demand(k_t_1_x_1,kc1,kc2,v,w,kmax);
        d0 = min(d0,d0_tmp);
        
        s1 = downstream_boundary(t*delta_t);
        s1_tmp = supply(k_t_1_x,kc1,kc2,v,w,kmax);
        s1 = min(s1,s1_tmp);
        
        d1 = upstream_boundary(t*delta_t);
        d1_tmp = demand(k_t_1_x,kc1,kc2,v,w,kmax);
        d1 = min(d1,d1_tmp);
        
        s2 = downstream_boundary(t*delta_t);
        s2_tmp = supply(k_t_1_x_p_1,kc1,kc2,v,w,kmax);
        s2 = min(s2,s2_tmp);

        
    k(x,t) =  k_t_1_x + (delta_t/delta_x)*(min(d0,s1)-min(d1,s2));
    end
end

figure;
%surf(k);
plot(k);
colormap('hot');
title('Godunov Scheme with initial condition , upstream boundary demand and downstream boundary supply');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G3.png')

function dpd = downstream_boundary(t)
    if (t>=0)&&(t<20)
        dpd = 0.6;
    elseif (t >= 20) && (t<=60)
        dpd = 0.05;
    end
end
function upd = upstream_boundary(t)
    if (t>=0)&&(t<40)
        upd = 0.3;
    elseif (t >= 40) && (t<=60)
        upd = 0.5;
    end
end
function k_0_x = density_mapper(x)
%This function is a helper function to map different constant densities
%correspending to different space values 
%x is a space value, k_0_x = intial  density value at point x
    if (x>=0)&&(x<300)
        k_0_x = 0.04;
    elseif (x>=300)&&(x<500)
        k_0_x = 0.1;
    elseif (x>=500)
        k_0_x = 0.02;
    end

end

function tfd = trapezoidal_fundamental_diagram(k,kc1,kc2,v,w,kmax)
%This function is the  Trapezoidal fundamental diagram, defined by ?(k)=k v for k<=kc1 ,
%?(k)=kc1 v for k>kc1 and k<=kc2  and ?(k)=-w(k-kmax) for k>kc2 
    if k <= kc1
       tfd = k*v;
    elseif (k>kc1)&&(k<=kc2)
        tfd = kc1*v;
    elseif k > kc2
        tfd = -1*w*(k-kmax);
    else
        tfd = inf;
    end
end

function dmn = demand(k,kc1,kc2,v,w,kmax)
%This function is the  demand function, defined by 
%D(k)= (?(k),if k?k_c1, ?(k_c1 )=vk_c1,if k>k_c1  and k?k_c2, ?(k_c2 )=vk_c1,if k>k_c2 
    if k <= kc1
       dmn = trapezoidal_fundamental_diagram(k,kc1,kc2,v,w,kmax);
    elseif (k>kc1)&&(k<=kc2)
        dmn = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        dmn = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);

    end
end

function spy = supply(k,kc1,kc2,v,w,kmax)
%This function is the  supply function, defined by 
%S(k)= ?(k_c1 )=vk_c1,if k?k_c1,?(k_c2 )=vk_c1,if k>k_c1  and k?k_c2 , ?(k_ ),if k>k_c2 
    if k <= kc1
       spy = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    elseif (k>kc1)&&(k<=kc2)
        spy = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        spy = trapezoidal_fundamental_diagram(k,kc1,kc2,v,w,kmax);

    end
end



