    clear;
    close all;
    clc;
    %% system parameters
    system.chip.N = 5;
    system.layer.N = 5; %two packages, PoP.
     
    %% for the interposer dimension
    system.hsp.Xsize = 40e-3;
    system.hsp.Ysize = 35e-3;
    system.hsp.Xgrid = 70e-6;
    system.hsp.Ygrid = 70e-6;
%     system.hsp.Xgrid = 100e-6;
%     system.hsp.Ygrid = 100e-6;
    
    system.pack.Xsize = 8.22e-3;
    system.pack.Ysize = 8.22e-3;
    system.pack.Xgrid = 7e-6;
    system.pack.Ygrid = 7e-6;
%     system.pack.Xgrid = 50e-6;
%     system.pack.Ygrid = 50e-6;

    %% for the chip grid information; size is assigned by each chip
    system.chip.Xgrid = 7e-6; %x grid size of chip
    system.chip.Ygrid = 7e-6; %y grid size of chip
%     system.chip.Xgrid = 50e-6; %x grid size of chip
%     system.chip.Ygrid = 50e-6; %y grid size of chip
    
    %% thickness, material identifier
%     system.Nchip = system.chip.N;
%     system.pack.metal_portion = 70*5/1600 * 0.5; % 5 routing layers, 60um per routing layer, total 1600, 50% metal per layer
    system.layer.hsp = [1000e-6, 4];
    system.layer.tim2 = [30e-6, 10];
%     system.layer.mold_a = [1e-6, 1];
%     system.layer.ILD = [2e-6, 5];
%     system.layer.ubumps = [2e-6, 3];
    system.layer.metal_portion = 0.5*0.3; %copper portion of the metal layer
    system.num.hsp = 1; %the layer of hsp modeled
    system.num.tim2 = 1; %the layer of TIM2 modeled
    
    %% transient analysis parameter
    system.tran = 0;
    system.T = 10;
    system.dt = 0.025;
    system.su = 1;
    system.su_limit = 10;
    system.Ntime = system.T/system.dt + 1;
        
    Material.K = [3 %0.9 %Molding 1:
                  149 %silicon 2:
                  1.6 %0.9 %underfill 3:
                  400 %copper 4:
                  1.38 %sio2 5:
                  60 %microbump 6:
                  149 %interposer 7:
                  0.024 %air 8:
                  0.1 %ceramic package property 9:
                  3 %TIM2 10
                  ]; 
    
    Material.C = [1000 % TIM1 1 : Thermal Management Products & Custom Solutions Catalog
                  705 % CHIP 2 : wikipedia
                  1000 % UNDER 3 : Thermal Management Products & Custom Solutions Catalog
                  385 % COPPER 4 : http://hyperphysics.phy-astr.gsu.edu/hbase/tables/sphtt.html
                  705 % SIO2 5 : http://www.azom.com/properties.aspx?ArticleID=1114
                  227 % UBUMP 6 : PAGE 11: http://www.almit.com/dloads/Agents/SAC%20Alloy%20Comparison.pdf
                  705 % SILICON INTERPOSER 7 :
                  1003.5 % air 8 :
                  600 % FR4 package 9 : http://www.electronics-cooling.com/2002/08/simplified-transient-model-for-ic-packages/
                  1000 %TIM2 10
                  ]; 
    
    Material.D = [2900 %TIM1 1
                  2329 %CHIP 2: http://en.wikipedia.org/wiki/Silicon
                  2100 %UNDER 3
                  8690 %COPPER 4: http://en.wikipedia.org/wiki/Copper
                  2648 %SIO2 5: http://en.wikipedia.org/wiki/Silicon_dioxide
                  12000 %UBUMP 6: http://en.wikipedia.org/wiki/Solder
                  2329 %SILICON INTERPOSER 7
                  1.225e-3 % air @ 15 centigrade 8
                  1850 % FR4 package 9
                  2900 %TIM2 10
                  ];     
    
    %% figure controller
    draw.granularity = 100;
    % thermal map, number of color used
    draw.write = 1;
    %write chip temperature files
    draw.T = 0;
    % whether to draw the thermal map; 1 yes; 0 no
    draw.P = 0;
    draw.pdrange = [0 250];
    draw.pdclamp = 1;
    %whether to draw the conductivity distribution
    draw.C = 0;
    % whether to draw the power maps; 1 yes, 0 no
    % only draw the die with blk_num > 1 (non-uniform cases)
    draw.displayT = 1;
    % print the temperature information
    draw.absolutely = 1;
    %print Tjunc - Ta or Tjunc; 1 for Tjunc
    draw.gif = 0;
    draw.range = [55 106];
    draw.clamp = 0;

    %% for die1 parameter
    chip(1).layer.id = 1;
    chip(1).model = 4; % for each die, how many layers we model
    
    % flip chip package; order:
    %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
    %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
    %chip thick format: thickness, (material identifier #)
    %% die layer information
    
    chip1_tim_thickness = 20e-6;
    chip1_thickness = 150e-6;
    
    chip(1).layer.mold_b = [chip1_tim_thickness, 1];
    chip(1).layer.die = [chip1_thickness, 2]; %die thickness
    chip(1).layer.ild = [3e-6, 5]; %metal layer thickness
    chip(1).layer.under = [2e-6, 5]; %underfill bonding thickness; between two dies
    chip(1).metal_portion = 0.5; %copper portion of the metal layer
    
    %% die size
    chip(1).Xsize = 8.22e-3; %x dimension of chip 1
    chip(1).Ysize = 8.22e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(1).bump.d = 2.5e-6;
    chip(1).bump.px = 5e-6;
    chip(1).bump.py = 5e-6;
    chip(1).bump.map = [0 0 chip(1).Xsize chip(1).Ysize];
%     chip(1).bump.map = [];  
    chip(1).bump.material = 4;

    %% power information
%     chip(1).power = [35.7552];  %power dissipation of each die
    chip(1).power = [5E+01];  %power dissipation of each die
    chip(1).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(1).map = [ ];
    chip(1).name = cellstr('Logic Tier #3');    
%     chip(1).blk_num = [12];
    chip(1).blk_num = [0];
   
    %%%%%%%%%%%%%%%%%%%%%%%%for die #2%%%%%%%%%%%%%%%%%%%%%%%%
    %% die layer information
    chip(2).layer.id = 2;
    chip(2).model = 3; % for each die, how many layers we model
    
%     chip1_chip2_middle_ILD_thickness = 0.1e-6;
    chip2_thickness = 25e-6;
    
%     chip(2).layer.mold_b = [chip1_chip2_middle_ILD_thickness, 5];
    chip(2).layer.die = [chip2_thickness, 2]; %die thickness
    chip(2).layer.ild = [3e-6, 5]; %metal layer thickness
    chip(2).layer.under = [2e-6, 5]; %underfill bonding thickness; between two dies
    chip(2).metal_portion = 0.5; %copper portion of the metal layer
    
    %% die size
    chip(2).Xsize = 8.22e-3; %x dimension of chip 1
    chip(2).Ysize = 8.22e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(2).bump.d = 2.5e-6;
    chip(2).bump.px = 5e-6;
    chip(2).bump.py = 5e-6;
    chip(2).bump.map = [0 0 chip(2).Xsize chip(2).Ysize];  
    chip(2).bump.material = 4;

    %% TSV geometry    
    chip(2).TSV.d = 2.5e-6; % tsv diameter including the liner thickness
    chip(2).TSV.liner = 0.2e-6; %liner thickness
    chip(2).TSV.px = 5e-6; %x direction pitch
    chip(2).TSV.py = 5e-6; %y direction pitch
    chip(2).TSV.material = [4, 5];
    x_TSV = 4.16e-3; %m
    xl_TSV = (chip(2).Xsize-x_TSV)/2;
    
    chip(2).TSV.map = [xl_TSV 0 x_TSV chip(2).Ysize];
%     chip(2).TSV.map = [0 0 chip(2).Xsize chip(2).Ysize];

    %% power information
    chip(2).power = [];  %power dissipation of each die
    chip(2).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(2).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(1).blk_name = cell(12, 1);
%     chip(1).blk_name(:) = cellstr('');
%     chip(1).blk_name(3) = cellstr('Cache');
%     chip(1).blk_name(5) = cellstr('Core');
%     chip(1).blk_name(8) = cellstr('Core');
    chip(2).name = cellstr('Memory Tier #2');    
    chip(2).blk_num = [0]; 

    %%%%%%%%%%%%%%%%%%%%%%%%for die #2%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%for die #3%%%%%%%%%%%%%%%%%%%%%%%%
    %% die layer information
    chip(3).layer.id = 3;
    chip(3).model = 3; % for each die, how many layers we model
    
%     chip1_chip2_middle_ILD_thickness = 0.1e-6;
    chip3_thickness = 25e-6;
    
%     chip(3).layer.mold_b = [chip1_chip2_middle_ILD_thickness, 5];
    chip(3).layer.die = [chip3_thickness, 2]; %die thickness
    chip(3).layer.ild = [3e-6, 5]; %metal layer thickness
    chip(3).layer.under = [2e-6, 5]; %underfill bonding thickness; between two dies
    chip(3).metal_portion = 0.5; %copper portion of the metal layer
    
    %% die size
    chip(3).Xsize = 8.22e-3; %x dimension of chip 1
    chip(3).Ysize = 8.22e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(3).bump.d = 2.5e-6;
    chip(3).bump.px = 5e-6;
    chip(3).bump.py = 5e-6;
    chip(3).bump.map = [0 0 chip(3).Xsize chip(3).Ysize];  
    chip(3).bump.material = 4;

    %% TSV geometry    
    chip(3).TSV.d = 2.5e-6; % tsv diameter including the liner thickness
    chip(3).TSV.liner = 0.2e-6; %liner thickness
    chip(3).TSV.px = 5e-6; %x direction pitch
    chip(3).TSV.py = 5e-6; %y direction pitch
    chip(3).TSV.material = [4, 5];
    x_TSV = 4.16e-3; %m
    xl_TSV = (chip(3).Xsize-x_TSV)/2;
    
    chip(3).TSV.map = [xl_TSV 0 x_TSV chip(3).Ysize];
%     chip(3).TSV.map = [0 0 chip(3).Xsize chip(3).Ysize];

    %% power information
    chip(3).power = [];  %power dissipation of each die
    chip(3).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(3).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(1).blk_name = cell(12, 1);
%     chip(1).blk_name(:) = cellstr('');
%     chip(1).blk_name(3) = cellstr('Cache');
%     chip(1).blk_name(5) = cellstr('Core');
%     chip(1).blk_name(8) = cellstr('Core');
    chip(3).name = cellstr('Logic Tier #2');    
    chip(3).blk_num = [0]; 

    %%%%%%%%%%%%%%%%%%%%%%%%for die #3%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%for die #4%%%%%%%%%%%%%%%%%%%%%%%%
    %% die layer information
    chip(4).layer.id = 4;
    chip(4).model = 3; % for each die, how many layers we model
    
%     chip1_chip2_middle_ILD_thickness = 0.1e-6;
    chip4_thickness = 25e-6;
    
%     chip(4).layer.mold_b = [chip1_chip2_middle_ILD_thickness, 5];
    chip(4).layer.die = [chip4_thickness, 2]; %die thickness
    chip(4).layer.ild = [3e-6, 5]; %metal layer thickness
    chip(4).layer.under = [2e-6, 5]; %underfill bonding thickness; between two dies
    chip(4).metal_portion = 0.5; %copper portion of the metal layer
    
    %% die size
    chip(4).Xsize = 8.22e-3; %x dimension of chip 1
    chip(4).Ysize = 8.22e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(4).bump.d = 2.5e-6;
    chip(4).bump.px = 5e-6;
    chip(4).bump.py = 5e-6;
    chip(4).bump.map = [0 0 chip(4).Xsize chip(4).Ysize];  
    chip(4).bump.material = 4;

    %% TSV geometry    
    chip(4).TSV.d = 2.5e-6; % tsv diameter including the liner thickness
    chip(4).TSV.liner = 0.2e-6; %liner thickness
    chip(4).TSV.px = 5e-6; %x direction pitch
    chip(4).TSV.py = 5e-6; %y direction pitch
    chip(4).TSV.material = [4, 5];
    x_TSV = 4.16e-3; %m
    xl_TSV = (chip(4).Xsize-x_TSV)/2;
    
    chip(4).TSV.map = [xl_TSV 0 x_TSV chip(4).Ysize];
%     chip(4).TSV.map = [0 0 chip(4).Xsize chip(4).Ysize];

    %% power information
    chip(4).power = [];  %power dissipation of each die
    chip(4).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(4).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(1).blk_name = cell(12, 1);
%     chip(1).blk_name(:) = cellstr('');
%     chip(1).blk_name(3) = cellstr('Cache');
%     chip(1).blk_name(5) = cellstr('Core');
%     chip(1).blk_name(8) = cellstr('Core');
    chip(4).name = cellstr('Memory Tier #1');    
    chip(4).blk_num = [0]; 

    %%%%%%%%%%%%%%%%%%%%%%%%for die #4%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%for die #5%%%%%%%%%%%%%%%%%%%%%%%%
    %% die layer information
    chip(5).layer.id = 5;
    chip(5).model = 3; % for each die, how many layers we model
    
%     chip1_chip2_middle_ILD_thickness = 0.1e-6;
    chip5_thickness = 25e-6;
    
%     chip(5).layer.mold_b = [chip1_chip2_middle_ILD_thickness, 5];
    chip(5).layer.die = [chip5_thickness, 2]; %die thickness
    chip(5).layer.ild = [3e-6, 5]; %metal layer thickness
    chip(5).layer.under = [2e-6, 5]; %underfill bonding thickness; between two dies
    chip(5).metal_portion = 0.5; %copper portion of the metal layer
    
    %% die size
    chip(5).Xsize = 8.22e-3; %x dimension of chip 1
    chip(5).Ysize = 8.22e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(5).bump.d = 2.5e-6;
    chip(5).bump.px = 5e-6;
    chip(5).bump.py = 5e-6;
    chip(5).bump.map = [0 0 chip(5).Xsize chip(5).Ysize];  
    chip(5).bump.material = 4;

    %% TSV geometry    
    chip(5).TSV.d = 2.5e-6; % tsv diameter including the liner thickness
    chip(5).TSV.liner = 0.2e-6; %liner thickness
    chip(5).TSV.px = 5e-6; %x direction pitch
    chip(5).TSV.py = 5e-6; %y direction pitch
    chip(5).TSV.material = [4, 5];
    x_TSV = 4.16e-3; %m
    xl_TSV = (chip(5).Xsize-x_TSV)/2;
    
    chip(5).TSV.map = [xl_TSV 0 x_TSV chip(5).Ysize];
%     chip(5).TSV.map = [0 0 chip(5).Xsize chip(5).Ysize];

    %% power information
    chip(5).power = [];  %power dissipation of each die
    chip(5).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(5).map = [ ];
    chip(5).name = cellstr('Logic Tier #1');    
    chip(5).blk_num = [0]; 

    %%%%%%%%%%%%%%%%%%%%%%%%for die #5%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%die, interposer, package information%%%%%%%%%%%
    %% adopt the method from PDN modeling
    %use the heat spreader as the reference plane
    system.pack.xl = (system.hsp.Xsize - system.pack.Xsize)/2;
    system.pack.yb = (system.hsp.Ysize - system.pack.Ysize)/2;
    
    chip(1).xl = (system.hsp.Xsize - chip(1).Xsize)/2;
    chip(1).yb = (system.hsp.Ysize - chip(1).Ysize)/2;
    
    chip(2).xl = (system.hsp.Xsize - chip(2).Xsize)/2;
    chip(2).yb = (system.hsp.Ysize - chip(2).Ysize)/2;
    
    chip(3).xl = (system.hsp.Xsize - chip(3).Xsize)/2;
    chip(3).yb = (system.hsp.Ysize - chip(3).Ysize)/2;
    
    chip(4).xl = (system.hsp.Xsize - chip(4).Xsize)/2;
    chip(4).yb = (system.hsp.Ysize - chip(4).Ysize)/2;
    
    chip(5).xl = (system.hsp.Xsize - chip(5).Xsize)/2;
    chip(5).yb = (system.hsp.Ysize - chip(5).Ysize)/2;
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the thermal resistivity of each boundary
     % r = 1/(hA); A is the size of top surface area
    % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
    % 0.5 W*cm^2/K
    h.up = 1/(2.25e-4);
%     h.up = 1/(0.03*1e-4);
    
    r_liquid = 1000*1e-4; %K-m2/W
%     r_liquid = 3.21*1e-4; %K-m2/W
    h.down = 1/(r_liquid); % the cooling of bottom surface; 
    %(only the area with the same size of chip;)
    %microfluidic is assumed to be as large as chip in the interposer
    
    h.side = 10;
    % side surface cooling, usually near adiabatic
    
    h.d = 1/(r_liquid);
    %the cooling of the bottom surface except for the MFHS area
    
    h.Ta = 300;
    %the ambient temperature; for CPU is 38C which is 311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%% main program to run

diary(['results_',datestr(now,'dd-mm-yy','local'),'_',datestr(now,'hh-MM-ss','local'),'.txt'])
diary on;
% format shortEng;

% cell_count = 4;
% sim_case = 3;
% T = ThermSim(system, chip, Material, draw, h, cell_count, sim_case);

%% simple case
cell_count = 4;
sim_case = [1];
substrate_thick = [150 50] * 1e-6; %units: m
tsv_diameter = [10] * 1e-6; %units: m
% hsp_size = [40e-3, 50e-3, 75e-3, 100e-3, 200e-3, 500e-3];


% RRAM
chip_size = [8.22] * 1e-3; %units: m
x_TSV = [0.28] * 1e-3; %m
pwr_logic = [27.42]; % W
pwr_mem = [99.99]; % W

% ADC variables
const = 1e-3;
x_margin_adc = 0.10;
y_margin_adc = 0.10;
nADC = 4500;
nyADC = 100;
nxADC = 45;
block_dimension_adc = 0.0632;
y_pitch_adc = 0.0171;
x_pitch_adc = 0.0143;
per_ADC_power = 3.02E-03;
x_dim_adc = 7.35;
pool_gb_xsize_adc = 0.5840;
y_dim_pool = 6.85;
y_dim_gb = 1.37;
tsv_xsize_adc = 0.2798;
pool_power = 1.40E-01;
gb_power = 2.00E-03;

% MEM variables
x_margin_MEM = 0.10;
y_margin_MEM = 0.10;
nMEM = 6325;
nyMEM = 115;
nxMEM = 55;
block_dimension_mem = 0.0616;
mem_y_dimension = 0.0167;
sar_y_dimension = 0.044877106;
y_pitch_mem = 0.0082;
x_pitch_mem = 0.0071;
per_MEM_power = 9.75E-04;
per_sar_power = 6.90E-03;
tsv_xsize_mem = 0.2798;
x_dim_mem = 7.94;


for i = 1:1:length(sim_case)
        
    system.pack.Xsize = chip_size(i);
    system.pack.Ysize = chip_size(i);
    
    chip(1).Xsize = chip_size(i);
    chip(1).Ysize = chip_size(i);
    chip(2).Xsize = chip_size(i);
    chip(2).Ysize = chip_size(i);
    chip(3).Xsize = chip_size(i);
    chip(3).Ysize = chip_size(i);
    chip(4).Xsize = chip_size(i);
    chip(4).Ysize = chip_size(i);
    chip(5).Xsize = chip_size(i);
    chip(5).Ysize = chip_size(i);
        
    %Bumps and TSVs
    %Tier 1
    chip(1).bump.d = tsv_diameter(i);
    chip(1).bump.px = chip(1).bump.d*2;
    chip(1).bump.py = chip(1).bump.d*2;
    chip(1).bump.map = [0 0 chip(1).Xsize chip(1).Ysize];
    
    %Tier 2
    % Bump geometry
    chip(2).bump.d = tsv_diameter(i);
    chip(2).bump.px = tsv_diameter(i)*2;
    chip(2).bump.py = tsv_diameter(i)*2;
    chip(2).bump.map = [0 0 chip(2).Xsize chip(2).Ysize];  
    
    % TSV geometry    
    chip(2).TSV.d = tsv_diameter(i); % tsv diameter including the liner thickness
    chip(2).TSV.px = tsv_diameter(i)*2; %x direction pitch
    chip(2).TSV.py = tsv_diameter(i)*2; %y direction pitch
%     chip(2).TSV.map = [0 0 chip(2).Xsize chip(2).Ysize];
    chip(2).TSV.map = [((chip(2).Xsize-x_TSV(i))/2) 0 x_TSV(i) chip(2).Ysize];
    
    %Tier 3
    % Bump geometry
    chip(3).bump.d = tsv_diameter(i);
    chip(3).bump.px = tsv_diameter(i)*2;
    chip(3).bump.py = tsv_diameter(i)*2;
    chip(3).bump.map = [0 0 chip(3).Xsize chip(3).Ysize];  
    
    % TSV geometry    
    chip(3).TSV.d = tsv_diameter(i); % tsv diameter including the liner thickness
    chip(3).TSV.px = tsv_diameter(i)*2; %x direction pitch
    chip(3).TSV.py = tsv_diameter(i)*2; %y direction pitch
    chip(3).TSV.map = [0 0 chip(3).Xsize chip(3).Ysize];
    chip(3).TSV.map = [((chip(3).Xsize-x_TSV(i))/2) 0 x_TSV(i) chip(3).Ysize];
    
    %Tier 4
    % Bump geometry
    chip(4).bump.d = tsv_diameter(i);
    chip(4).bump.px = tsv_diameter(i)*2;
    chip(4).bump.py = tsv_diameter(i)*2;
    chip(4).bump.map = [0 0 chip(4).Xsize chip(4).Ysize];  
    
    % TSV geometry    
    chip(4).TSV.d = tsv_diameter(i); % tsv diameter including the liner thickness
    chip(4).TSV.px = tsv_diameter(i)*2; %x direction pitch
    chip(4).TSV.py = tsv_diameter(i)*2; %y direction pitch
    chip(4).TSV.map = [0 0 chip(4).Xsize chip(4).Ysize];
    chip(4).TSV.map = [((chip(4).Xsize-x_TSV(i))/2) 0 x_TSV(i) chip(4).Ysize];
    
    %Tier 5
    % Bump geometry
    chip(5).bump.d = tsv_diameter(i);
    chip(5).bump.px = tsv_diameter(i)*2;
    chip(5).bump.py = tsv_diameter(i)*2;
    chip(5).bump.map = [0 0 chip(5).Xsize chip(5).Ysize];  
    
    % TSV geometry    
    chip(5).TSV.d = tsv_diameter(i); % tsv diameter including the liner thickness
    chip(5).TSV.px = tsv_diameter(i)*2; %x direction pitch
    chip(5).TSV.py = tsv_diameter(i)*2; %y direction pitch
    chip(5).TSV.map = [0 0 chip(5).Xsize chip(5).Ysize];
    chip(5).TSV.map = [((chip(5).Xsize-x_TSV(i))/2) 0 x_TSV(i) chip(5).Ysize];
    
        
    % Substrate thickness
    chip(1).layer.die(1) = substrate_thick(1); %die thickness
    chip(2).layer.die(1) = substrate_thick(2); %die thickness
    chip(3).layer.die(1) = substrate_thick(2); %die thickness
    chip(4).layer.die(1) = substrate_thick(2); %die thickness
    chip(5).layer.die(1) = substrate_thick(2); %die thickness
    
%     for m = 6%1:1:length(hsp_size)
%         system.hsp.Xsize = hsp_size(m);
%         system.hsp.Ysize = hsp_size(m);

        % adjust package and chip xl and yb
        system.pack.xl = (system.hsp.Xsize - system.pack.Xsize)/2;
        system.pack.yb = (system.hsp.Ysize - system.pack.Ysize)/2;

        chip(1).xl = (system.hsp.Xsize - chip(1).Xsize)/2;
        chip(1).yb = (system.hsp.Ysize - chip(1).Ysize)/2;
        chip(2).xl = (system.hsp.Xsize - chip(2).Xsize)/2;
        chip(2).yb = (system.hsp.Ysize - chip(2).Ysize)/2;
        chip(3).xl = (system.hsp.Xsize - chip(3).Xsize)/2;
        chip(3).yb = (system.hsp.Ysize - chip(3).Ysize)/2;
        chip(4).xl = (system.hsp.Xsize - chip(4).Xsize)/2;
        chip(4).yb = (system.hsp.Ysize - chip(4).Ysize)/2;
        chip(5).xl = (system.hsp.Xsize - chip(5).Xsize)/2;
        chip(5).yb = (system.hsp.Ysize - chip(5).Ysize)/2;

        chip(1).power = pwr_logic;
        chip(2).power = pwr_mem;
        chip(3).power = pwr_logic;
        chip(4).power = pwr_mem;
        chip(5).power = pwr_logic;

        chip(1).map = logic_power_map_gen(const, x_margin_adc, y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        A = logic_power_map_gen(const, (x_dim_adc/2 + pool_gb_xsize_adc + tsv_xsize_adc + x_margin_adc), y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        chip(1).map = [chip(1).map ; A];
        P1 = [(x_dim_adc/2)*const 0 (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB1 = [(x_dim_adc/2)*const (y_dim_pool)*const (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        P2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const (y_dim_gb)*const (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const 0 (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        chip(1).map = [chip(1).map ; P1 ; GB1 ; P2 ; GB2];
        chip(1).blk_name = cell((nADC*2)+4, 1);
        chip(1).blk_name(:) = cellstr('');
        chip(1).blk_num = [(nADC*2)+4];

        chip(2).map = mem_power_map_gen(const, x_margin_MEM, y_margin_MEM, nMEM, nyMEM, nxMEM, block_dimension_mem, mem_y_dimension, sar_y_dimension, y_pitch_mem, x_pitch_mem, per_MEM_power, per_sar_power);
        B = mem_power_map_gen(const, (x_dim_mem/2 + tsv_xsize_mem + x_margin_MEM), y_margin_MEM, nMEM, nyMEM, nxMEM, block_dimension_mem, mem_y_dimension, sar_y_dimension, y_pitch_mem, x_pitch_mem, per_MEM_power, per_sar_power);;
        chip(2).map = [chip(2).map ; B];
        chip(2).blk_name = cell(nMEM*4, 1);
        chip(2).blk_name(:) = cellstr('');
        chip(2).blk_num = [nMEM*4];    

        chip(3).map = logic_power_map_gen(const, x_margin_adc, y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        C = logic_power_map_gen(const, (x_dim_adc/2 + pool_gb_xsize_adc + tsv_xsize_adc + x_margin_adc), y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        chip(3).map = [chip(3).map ; C];
        P1 = [(x_dim_adc/2)*const 0 (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB1 = [(x_dim_adc/2)*const (y_dim_pool)*const (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        P2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const (y_dim_gb)*const (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const 0 (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        chip(3).map = [chip(3).map ; P1 ; GB1 ; P2 ; GB2];
        chip(3).blk_name = cell((nADC*2)+4, 1);
        chip(3).blk_name(:) = cellstr('');
        chip(3).blk_num = [(nADC*2)+4];

        chip(4).map = mem_power_map_gen(const, x_margin_MEM, y_margin_MEM, nMEM, nyMEM, nxMEM, block_dimension_mem, mem_y_dimension, sar_y_dimension, y_pitch_mem, x_pitch_mem, per_MEM_power, per_sar_power);
        D = mem_power_map_gen(const, (x_dim_mem/2 + tsv_xsize_mem + x_margin_MEM), y_margin_MEM, nMEM, nyMEM, nxMEM, block_dimension_mem, mem_y_dimension, sar_y_dimension, y_pitch_mem, x_pitch_mem, per_MEM_power, per_sar_power);
        chip(4).map = [chip(4).map ; D];
        chip(4).blk_name = cell(nMEM*4, 1);
        chip(4).blk_name(:) = cellstr('');
        chip(4).blk_num = [nMEM*4];

        chip(5).map = logic_power_map_gen(const, x_margin_adc, y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        E = logic_power_map_gen(const, (x_dim_adc/2 + pool_gb_xsize_adc + tsv_xsize_adc + x_margin_adc), y_margin_adc, nADC, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
        chip(5).map = [chip(5).map ; E];
        P1 = [(x_dim_adc/2)*const 0 (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB1 = [(x_dim_adc/2)*const (y_dim_pool)*const (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        P2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const (y_dim_gb)*const (pool_gb_xsize_adc/2)*const y_dim_pool*const pool_power 0];
        GB2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc)*const 0 (pool_gb_xsize_adc/2)*const y_dim_gb*const gb_power 0];
        chip(5).map = [chip(5).map ; P1 ; GB1 ; P2 ; GB2];
        chip(5).blk_name = cell((nADC*2)+4, 1);
        chip(5).blk_name(:) = cellstr('');
        chip(5).blk_num = [(nADC*2)+4];

        fprintf('\n\nDie size: %12.3e m\n', chip_size(i));
        fprintf('Die area: %12.3e m2\n', chip(1).Xsize*chip(1).Ysize);
        fprintf('Power for Logic: %12.3e W\n', chip(1).power);
        fprintf('Power for Memory: %12.3e W\n', chip(2).power);
        fprintf('Total stack power: %12.3e W\n', (chip(1).power*3 + chip(2).power*2));
        fprintf('Die bulk thickness: %12.3e m\n', substrate_thick(i));
        fprintf('TSV diameter: %12.3e m\n', chip(2).TSV.d);
        fprintf('TSV pitch (x and y): %12.3e m\n', chip(2).TSV.px);
        fprintf('ubump diameter: %12.3e m\n', chip(2).bump.d);
        fprintf('ubump pitch (x and y): %12.3e m\n', chip(2).bump.px);
%         fprintf('HSP size: %12.3e m\n', system.hsp.Xsize);

        T = ThermSim(system, chip, Material, draw, h, cell_count, sim_case(i));
%     end

end

%%
diary off;