%----------------------2-D Transient Heat Conduction---------------------- 
%----------------------No Heat Generation----------------------
%Dr.M.binazadeh
%Korosh Agha Mohammad Ghasemi
%Chemical Engineering at Shiraz University
clc;
clear all;
%% Variable Declaration
n = 21;                 %number of nodes
L = 1;                  %length of domain
W = 1;                  %width of domain
alpha = 1e-4;           %thermal diffusivity (m^2/s)

m =(n-2)*(n-2);         %variable used to construct penta diagonal matrix
sm = sqrt(m);

dx = L/(n-1);           %delta x domain 
dy = W/(n-1);           %delta y domain
x = linspace(0,L,n);    %linearly spaced vectors x direction
y = linspace(0,W,n);    %linearly spaced vectors y direction
[X,Y]=meshgrid(x,y);

Tin = 200;              %internal temperature
T  = ones(m,1)* Tin;    %initilizing Space
Ta = zeros(n,n);
A  = zeros(m,m);         
Ax = zeros(m,1);
B  = zeros(sm,sm);

dt = 0.1;                %time step
tmax = 500;              %total Time steps (s)
t = 0 : dt : tmax;
r = alpha * dt /(dx^2);  %for stability, must be 0.5 or less

Tt = 100;                %Top Wall
Tb = 100;                %Bottom Wall
Tl = 800;                %Left Wall
Tr = 100;                %Right Wall

%% Boundry Conditions
Ta(1,1:n) = Tt;         %Top Wall
Ta(n,1:n) = Tb;         %Bottom Wall
Ta(1:n,1) = Tl;         %Left Wall
Ta(1:n,n) = Tr;         %Right Wall

%% Setup Matrix
for ix = 1 : 1 : m
    
    for jx =1 : 1 : m
        
        if (ix == jx)
                A(ix,jx) = (1 + 4*r); 
                
        elseif ( (ix == jx + 1) && ( (ix - 1) ~= sm * round( (ix-1)/sm) ) ) %RHS
                A(ix,jx) = -r;
                
        elseif ( (ix == jx - 1) && ( ix ~= sm*round(ix/sm) ) )
                A(ix,jx) = -r;
                
        elseif (ix == jx + sm)
                A(ix,jx) = -r;
                
        elseif (jx == ix + sm)
                A(ix,jx) = -r;
        else
                A(ix,jx) = 0;
        end
    end
end

 
for iy = 1 : 1 : sm
    
    for jy = 1 : 1 : sm
        
        if (iy == 1) && (jy == 1)
            B(iy,jy) = r * (Tl + Tt);
            
        elseif (iy == 1) && (jy == sm)
            B(iy,jy) = r * (Tt + Tr);                           %LHS
            
        elseif (iy == sm) && (jy == sm)
            B(iy,jy) = r * (Tb + Tr);
            
        elseif (iy == sm) && (jy == 1)
            B(iy,jy) = r * (Tb + Tl);
            
        elseif (iy == 1) && (jy == sm)
            B(iy,jy) = r * (Tt + Tr);
            
        elseif (iy == 1)&&(jy > 1 || jy < sm)
            B(iy,jy) = r * Tt;
            
        elseif (jy == sm) && (iy > 1 || iy < sm)
            B(iy,jy) = r * Tr;
            
        elseif (iy == sm) && ( jy > 1 || jy < sm)
            B(iy,jy) = r * Tb;
            
        elseif (jy == 1) && ( iy > 1 || jy < sm)
            B(iy,jy) = r * Tl;
        else
            B(iy,jy) = 0;
        end
                    
    end
end

 Bx = reshape(B,[],1);  %Convert matrix to vector
 
%% Solution
for l = 2 : length(t)   %time steps
    
     Xx = ( T + Bx );
     
           Ax = A \ Xx;
           
                 T( 1 : m ) = Ax( 1 : m );

                       fprintf('Time t=%d\n',l-1);
    
end

Tx = reshape( Ax , sm , sm); %convert vector to matrix

for i = 2 : 1 : n-1
    
    for j = 2 : 1 :n-1
        
        Ta(i,j) = Tx ( i-1 , j-1 );
        
    end
    
end

%% Plot

   contourf(X,Y,Ta,50,'edgecolor','none');
        h = colorbar;
        ylabel(h, 'Temperature °C')
        colormap jet
        axis equal
            
        title(['Top (Tt)= ',num2str(Tt),'°C']);
            xlabel(['Bottom (Tb)= ',num2str(Tb),'°C'])
            
            yyaxis left
            ylabel(['Left (Tl)= ',num2str(Tl),'°C'])
            
            yyaxis right
            ylabel(['Right (Tr)= ',num2str(Tr),'°C'])
        
        