% Loop for results and comparisons

clear all; 
clc; 

Main_path= ['/Users/johanatch/Desktop/Quant Econ/HW3/Final'];
cd(Main_path)

NN = [31 50 100 300]; % Number of initial grid points

sizeNN = size(NN,2);

for comparison=1:sizeNN

    N = NN(comparison);

    [kgrid, T, VALUEFUNCTION, K_TOMORR, Kstar] = ngm(N);
    TT(comparison)=T;
    KKGRID(comparison).kgrid = kgrid;
    VALFUNCTION(comparison).val = VALUEFUNCTION;
    KTOM(comparison).kdeci = K_TOMORR;

    BIG_GRID = linspace(0.75*Kstar, 1.25*Kstar,1000)';
  
    POLICY(comparison).policy = interp1(KKGRID(comparison).kgrid, KTOM(comparison).kdeci, BIG_GRID, 'cubic');
    VALUE(comparison).value = interp1(KKGRID(comparison).kgrid, VALFUNCTION(comparison).val, BIG_GRID, 'cubic');

end

 v3=figure; % Policy Function
    plot(BIG_GRID,VALUE(1).value); hold;
    plot(BIG_GRID,VALUE(2).value);
    plot(BIG_GRID,VALUE(3).value);
    plot(BIG_GRID,VALUE(4).value);
    legend('N=31','N=50','N=100','N=300');
    title('Value Function');
    xlabel('kt'); 
    ylabel('v(kt,z=1)');
    
    saveas(v3,'VAL_COMPARISON.png')

 f3=figure; % Policy Function
    plot(BIG_GRID,POLICY(1).policy); hold;
    plot(BIG_GRID,POLICY(2).policy); 
    plot(BIG_GRID,POLICY(3).policy);
    plot(BIG_GRID,POLICY(4).policy);
    legend('N=31','N=50','N=100','N=300');
    title('Policy Function');
    xlabel('kt'); 
    ylabel('kt+1(kt,z=1)');
    
    saveas(f3,'POL_COMPARISON.png');

dif=max(abs(POLICY(4).policy-POLICY(1).policy));
