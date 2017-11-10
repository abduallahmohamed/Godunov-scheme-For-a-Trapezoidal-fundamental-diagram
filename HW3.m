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

%Godunov scheme - Solution #4
cell_width = 100; %Width of each cell in x, space step
delta_x = cell_width;
cell_x = 10; %Number of cells
delta_t = 3; %Minimal time step that satisifies CFL condition
time_scale = 60; %Width of time
cell_t = time_scale / delta_t; %Number of cells in time
fac_t_x = delta_t/delta_x;
k = zeros(cell_t,cell_x); %Grid
k_inv = zeros(cell_t,cell_x);

%Apply initial condition 
for x = 1:cell_x
    k(1,x) = density_mapper(x*delta_x);
end


for t = 2:cell_t
    
    for x = 1:cell_x
      
        k_t_nv_one_x = k(t-1,x);
        s_k_t_nv_one_x = supply(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        d_k_t_nv_one_x = demand(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        
        
        if x-1 < 1
            d_k_t_nv_one_x_nv_one = upstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_nv_one = k(t-1,x-1);
            d_k_t_nv_one_x_nv_one = demand(k_t_nv_one_x_nv_one,kc1,kc2,v,w,kmax);

        end
        if x+1 > cell_x
            s_t_nv_one_x_pl_one = downstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_pl_one = k(t-1,x+1);
            s_t_nv_one_x_pl_one = supply(k_t_nv_one_x_pl_one,kc1,kc2,v,w,kmax);

        end
        

        k(t,x) = k_t_nv_one_x +fac_t_x * (min(d_k_t_nv_one_x_nv_one,s_k_t_nv_one_x) - min(d_k_t_nv_one_x,s_t_nv_one_x_pl_one));
             
    end  
end


figure;
imagesc(k);
ax = gca;
ax.XAxisLocation = 'top';
ax.YAxisLocation = 'left';
title('Godunov Scheme 20*10 Grid');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G1.png')

%Godunov scheme - Solution #5
cell_width = 10; %Width of each cell in x, space step
delta_x = cell_width;
cell_x = 100; %Number of cells
delta_t = 0.3; %Minimal time step that satisifies CFL condition
time_scale = 60; %Width of time
cell_t = time_scale / delta_t; %Number of cells in time
fac_t_x = delta_t/delta_x;
k = zeros(cell_t,cell_x); %Grid
k_inv = zeros(cell_t,cell_x);

%Apply initial condition 
for x = 1:cell_x
    k(1,x) = density_mapper(x*delta_x);
end


for t = 2:cell_t
    
    for x = 1:cell_x
      
        k_t_nv_one_x = k(t-1,x);
        s_k_t_nv_one_x = supply(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        d_k_t_nv_one_x = demand(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        
        
        if x-1 < 1
            d_k_t_nv_one_x_nv_one = upstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_nv_one = k(t-1,x-1);
            d_k_t_nv_one_x_nv_one = demand(k_t_nv_one_x_nv_one,kc1,kc2,v,w,kmax);

        end
        if x+1 > cell_x
            s_t_nv_one_x_pl_one = downstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_pl_one = k(t-1,x+1);
            s_t_nv_one_x_pl_one = supply(k_t_nv_one_x_pl_one,kc1,kc2,v,w,kmax);

        end
        

        k(t,x) = k_t_nv_one_x +fac_t_x * (min(d_k_t_nv_one_x_nv_one,s_k_t_nv_one_x) - min(d_k_t_nv_one_x,s_t_nv_one_x_pl_one));
             
    end  
end


figure;
imagesc(k);
ax = gca;
ax.XAxisLocation = 'top';
ax.YAxisLocation = 'left';
title('Godunov Scheme 200*100 Grid');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G2.png')


%Godunov scheme - Solution #5
cell_width = 10; %Width of each cell in x, space step
delta_x = cell_width;
cell_x = 100; %Number of cells
delta_t = 0.3; %Minimal time step that satisifies CFL condition
time_scale = 60; %Width of time
cell_t = time_scale / delta_t; %Number of cells in time
fac_t_x = delta_t/delta_x;
k = zeros(cell_t,cell_x); %Grid
k_inv = zeros(cell_t,cell_x);

%Apply initial condition 
for x = 1:cell_x
    k(1,x) = density_mapper(x*delta_x);
end


for t = 2:cell_t
    
    for x = 1:cell_x
      
        k_t_nv_one_x = k(t-1,x);
        s_k_t_nv_one_x = supply(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        d_k_t_nv_one_x = demand(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        
        
        if x-1 < 1
            d_k_t_nv_one_x_nv_one = upstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_nv_one = k(t-1,x-1);
            d_k_t_nv_one_x_nv_one = demand(k_t_nv_one_x_nv_one,kc1,kc2,v,w,kmax);

        end
        if x+1 > cell_x
            s_t_nv_one_x_pl_one = downstream_boundary_zero(t*delta_t);

        else
            k_t_nv_one_x_pl_one = k(t-1,x+1);
            s_t_nv_one_x_pl_one = supply(k_t_nv_one_x_pl_one,kc1,kc2,v,w,kmax);

        end
        

        k(t,x) = k_t_nv_one_x +fac_t_x * (min(d_k_t_nv_one_x_nv_one,s_k_t_nv_one_x) - min(d_k_t_nv_one_x,s_t_nv_one_x_pl_one));
             
    end  
end


figure;
imagesc(k);
ax = gca;
ax.XAxisLocation = 'top';
ax.YAxisLocation = 'left';
title('Godunov Scheme 200*100 Grid with ds boundary supply of 0');
xlabel('Space (x)');
ylabel('Time (t)');
saveas(gcf,'G3.png')


%FUNCTIONS

function dpd = downstream_boundary_zero(t)
   dpd = 0;
end
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



