    clear;
    close all;
    clc;
    %% system parameters
    system.chip.N = 8;
    system.layer.N = 3; %two packages, PoP.
    
    %% for the interposer dimension
    system.hsp.Xsize = 25e-3;
    system.hsp.Ysize = 35e-3;
    system.hsp.Xgrid = 250e-6;
    system.hsp.Ygrid = 350e-6;
%     system.hsp.Xsize = 12e-3;
%     system.hsp.Ysize = 12e-3;
%     system.hsp.Xgrid = 12e-3;
%     system.hsp.Ygrid = 12e-3;
    
    system.pack.Xsize = 10e-3;
    system.pack.Ysize = 10e-3;
    system.pack.Xgrid = 25e-6;
    system.pack.Ygrid = 25e-6;
%     system.pack.Xsize = 10e-3;
%     system.pack.Ysize = 10e-3;
%     system.pack.Xgrid = 10e-3;
%     system.pack.Ygrid = 10e-3;

    %% for the chip grid information; size is assigned by each chip
    system.chip.Xgrid = 50e-6; %x grid size of chip
    system.chip.Ygrid = 50e-6; %y grid size of chip
%     system.chip.Xgrid = 10e-3; %x grid size of chip
%     system.chip.Ygrid = 10e-3; %y grid size of chip
    
    %% thickness, material identifier
    %five layers of package
    %EMIB is located in the first 100 um.
    system.pack.metal_portion = 70*5/1600 * 0.5; % 5 routing layers, 60um per routing layer, total 1600, 50% metal per layer
    system.layer.hsp = [1000e-6, 4];
    system.layer.tim2 = [50e-6, 10];
%     system.layer.mold_a = [1e-6, 1];
    system.layer.ILD = [50e-6, 5];
%     system.layer.ubumps = [2e-6, 3];
    system.layer.metal_portion = 0.5*0.3; %copper portion of the metal layer
    system.num.hsp = 1; %the layer of hsp modeled
    system.num.tim2 = 1; %the layer of TIM2 modeled
    
    %% transient analysis parameter
    system.tran = 0;
    system.T = 10;
    system.dt = 1e-3;
    system.su = 1;
    system.su_limit = 10; 
        
    Material.K = [3 %Molding 1:
                  149 %silicon 2:
                  0.9 %underfill 3:
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
    draw.T = 1;
    % whether to draw the thermal map; 1 yes; 0 no
    draw.P = 1;
    %whether to draw the conductivity distribution
    draw.C = 0;
    % whether to draw the power maps; 1 yes, 0 no
    % only draw the die with blk_num > 1 (non-uniform cases)
    draw.displayT = 1;
    % print the temperature information
    draw.absolutely = 1;
    %print Tjunc - Ta or Tjunc; 1 for Tjunc
    draw.gif = 0;
    draw.range = [55 105];
    draw.clamp = 0;

    %% for die1 parameter
    chip(1).layer.id = 1;
    chip(1).model = 4; % for each die, how many layers we model
      
    %chip thick format: thickness, (material identifier #)
    %% die layer information
    
    chip1_tim_thickness = 50e-6;
    chip1_thickness = 150e-6;

    chip(1).layer.mold_b = [chip1_tim_thickness, 1];
    chip(1).layer.die = [chip1_thickness, 2]; %die thickness
    chip(1).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(1).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(1).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(1).Xsize = 10e-3; %x dimension of chip 1
    chip(1).Ysize = 10e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(1).bump.d = 40e-6;
    chip(1).bump.px = 200e-6;
    chip(1).bump.py = 200e-6;
    chip(1).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(1).bump.map = [ ];
    chip(1).bump.material = 6;

    %% power information
    chip(1).power = [100];  %power dissipation of each die
    chip(1).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(1).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(1).blk_name = cell(12, 1);
%     chip(1).blk_name(:) = cellstr('');
%     chip(1).blk_name(3) = cellstr('Cache');
%     chip(1).blk_name(5) = cellstr('Core');
%     chip(1).blk_name(8) = cellstr('Core');
    chip(1).name = cellstr('CHIP # 1');    
    chip(1).blk_num = [0];
   
    %%%%%%%%%%%%%%%%%%%%%%%%for die #2%%%%%%%%%%%%%%%%%%%%%%%%
    %% die layer information
    chip(2).layer.id = 2;
    chip(2).model = 4; % for each die, how many layers we model
    
    chip1_chip2_middle_ILD_thickness = 50e-6;
    chip2_thickness = 50e-6;

    chip(2).layer.mold_b = [chip1_chip2_middle_ILD_thickness, 5];
    chip(2).layer.die = [chip2_thickness, 2]; %die thickness
    chip(2).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(2).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(2).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(2).Xsize = 10e-3; %x dimension of chip 1
    chip(2).Ysize = 10e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(2).bump.d = 40e-6;
    chip(2).bump.px = 200e-6;
    chip(2).bump.py = 200e-6;
    chip(2).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(2).bump.map = [ ];
    chip(2).bump.material = 6;
    
    chip(2).TSV.map = [0     0  system.pack.Xsize system.pack.Ysize];
%     chip(2).TSV.map = [ ];
    chip(2).TSV.d = 20e-6;
    chip(2).TSV.liner = 0e-6;
    chip(2).TSV.material = [4, 5];
    chip(2).TSV.px = 200e-6;
    chip(2).TSV.py = 200e-6;    
    %% power information
    chip(2).power = [5];  %power dissipation of each die
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
    chip(2).name = cellstr('CHIP # 2');    
    chip(2).blk_num = [0]; 

    %%%%%%%%%%%%%%%%%%%%%%%%for die #2%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #3%%%%%%%%%%%%%%%%%%%%%%%%
    chip(3).layer.id = 3;
    chip(3).model = 4; % for each die, how many layers we model
    
    tier2_tier3_middle_ILD_thickness = 50e-6;
    chip3_thickness = 50e-6;

    chip(3).layer.mold_b = [tier2_tier3_middle_ILD_thickness, 5];
    chip(3).layer.die = [chip3_thickness, 2]; %die thickness
    chip(3).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(3).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(3).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(3).Xsize = 2.5e-3; %x dimension of chip 1
    chip(3).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
    chip(3).bump.d = 40e-6;
    chip(3).bump.px = 200e-6;
    chip(3).bump.py = 200e-6;
    chip(3).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(3).bump.map = [ ];  
    chip(3).bump.material = 6;

% chip(3).TSV.map = [0     0  chip(3).Xsize chip(3).Ysize];
%     chip(3).TSV.map = [];
    chip(3).TSV.d = 20e-6;
    chip(3).TSV.liner = 0e-6;
    chip(3).TSV.material = [4, 5];
    chip(3).TSV.px = 200e-6;
    chip(3).TSV.py = 200e-6;  
    %% power information
    chip(3).power = [1];  %power dissipation of each die
    chip(3).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(3).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(3).blk_name = cell(12, 1);
%     chip(3).blk_name(:) = cellstr('');
%     chip(3).blk_name(3) = cellstr('Cache');
%     chip(3).blk_name(5) = cellstr('Core');
%     chip(3).blk_name(8) = cellstr('Core');
    chip(3).name = cellstr('CHIP # 3');    
    chip(3).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #3%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #4%%%%%%%%%%%%%%%%%%%%%%%%
    chip(4).layer.id = 3;
    chip(4).model = 4; % for each die, how many layers we model
    
    chip4_mold_thickness = 50e-6;
    chip4_thickness = 50e-6;
    
    chip(4).layer.mold_b = [chip4_mold_thickness, 5];
    chip(4).layer.die = [chip4_thickness, 2]; %die thickness
    chip(4).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(4).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(4).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(4).Xsize = 2.5e-3; %x dimension of chip 1
    chip(4).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
%     chip(4).bump.d = 500e-9; %Modified this with current ASCENT theme 3 goal of 1um pitch, with which a diameter of 500nm is assumed
%     chip(4).bump.px = 1e-6;
%     chip(4).bump.py = 1e-6;
    chip(4).bump.d = 40e-6;
    chip(4).bump.px = 200e-6;
    chip(4).bump.py = 200e-6;
%     chip(4).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(4).bump.map = [ ];  
    chip(4).bump.material = 6;
    
%     chip(4).TSV.map = [0     0  chip(4).Xsize chip(4).Ysize];
    chip(4).TSV.d = 20e-6;
    chip(4).TSV.liner = 0e-6;
    chip(4).TSV.material = [4, 5];
    chip(4).TSV.px = 200e-6;
    chip(4).TSV.py = 200e-6;  
    %% power information
    chip(4).power = [3.755];  %power dissipation of each die
    chip(4).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(4).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(4).blk_name = cell(12, 1);
%     chip(4).blk_name(:) = cellstr('');
%     chip(4).blk_name(3) = cellstr('Cache');
%     chip(4).blk_name(5) = cellstr('Core');
%     chip(4).blk_name(8) = cellstr('Core');
    chip(4).name = cellstr('CHIP # 4');    
    chip(4).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #4%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #5%%%%%%%%%%%%%%%%%%%%%%%%
    chip(5).layer.id = 3;
    chip(5).model = 4; % for each die, how many layers we model
    
    chip5_mold_thickness = 50e-6;
    chip5_thickness = 50e-6;
    
    chip(5).layer.mold_b = [chip5_mold_thickness, 5];
    chip(5).layer.die = [chip5_thickness, 2]; %die thickness
    chip(5).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(5).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(5).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(5).Xsize = 2.5e-3; %x dimension of chip 1
    chip(5).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
%     chip(5).bump.px = 1e-6;
%     chip(5).bump.py = 1e-6;
    chip(5).bump.d = 40e-6;
    chip(5).bump.px = 200e-6;
    chip(5).bump.py = 200e-6;
%     chip(5).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(5).bump.map = [ ];  
    chip(5).bump.material = 6;
    
%     chip(5).TSV.map = [0     0  chip(5).Xsize chip(5).Ysize];
    chip(5).TSV.d = 20e-6;
    chip(5).TSV.liner = 0e-6;
    chip(5).TSV.material = [4, 5];
    chip(5).TSV.px = 200e-6;
    chip(5).TSV.py = 200e-6;  
    %% power information
    chip(5).power = [2];  %power dissipation of each die
    chip(5).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(5).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(5).blk_name = cell(12, 1);
%     chip(5).blk_name(:) = cellstr('');
%     chip(5).blk_name(3) = cellstr('Cache');
%     chip(5).blk_name(5) = cellstr('Core');
%     chip(5).blk_name(8) = cellstr('Core');
    chip(5).name = cellstr('CHIP # 5');    
    chip(5).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #5%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #6%%%%%%%%%%%%%%%%%%%%%%%%
    chip(6).layer.id = 3;
    chip(6).model = 4; % for each die, how many layers we model
    
    chip6_mold_thickness = 50e-6;
    chip6_thickness = 50e-6;
    
    chip(6).layer.mold_b = [chip6_mold_thickness, 5];
    chip(6).layer.die = [chip6_thickness, 2]; %die thickness
    chip(6).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(6).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(6).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(6).Xsize = 2.5e-3; %x dimension of chip 1
    chip(6).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
%     chip(6).bump.px = 1e-6;
%     chip(6).bump.py = 1e-6;
    chip(6).bump.d = 40e-6;
    chip(6).bump.px = 200e-6;
    chip(6).bump.py = 200e-6;
%     chip(6).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(6).bump.map = [ ];  
    chip(6).bump.material = 6;
    
%     chip(6).TSV.map = [0     0  chip(6).Xsize chip(6).Ysize];
    chip(6).TSV.d = 20e-6;
    chip(6).TSV.liner = 0e-6;
    chip(6).TSV.material = [4, 5];
    chip(6).TSV.px = 200e-6;
    chip(6).TSV.py = 200e-6;  
    %% power information
    chip(6).power = [3.75];  %power dissipation of each die
    chip(6).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(6).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(3).blk_name = cell(12, 1);
%     chip(3).blk_name(:) = cellstr('');
%     chip(3).blk_name(3) = cellstr('Cache');
%     chip(3).blk_name(5) = cellstr('Core');
%     chip(3).blk_name(8) = cellstr('Core');
    chip(6).name = cellstr('CHIP # 6');    
    chip(6).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #6%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #7%%%%%%%%%%%%%%%%%%%%%%%%
    chip(7).layer.id = 3;
    chip(7).model = 4; % for each die, how many layers we model
    
    chip7_mold_thickness = 50e-6;
    chip7_thickness = 50e-6;
    
    chip(7).layer.mold_b = [chip7_mold_thickness, 5];
    chip(7).layer.die = [chip7_thickness, 2]; %die thickness
    chip(7).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(7).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(7).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(7).Xsize = 2.5e-3; %x dimension of chip 1
    chip(7).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
%     chip(7).bump.px = 1e-6;
%     chip(7).bump.py = 1e-6;
    chip(7).bump.d = 40e-6;
    chip(7).bump.px = 200e-6;
    chip(7).bump.py = 200e-6;
%     chip(7).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(7).bump.map = [ ];  
    chip(7).bump.material = 6;
    
%     chip(7).TSV.map = [0     0  chip(7).Xsize chip(7).Ysize];
    chip(7).TSV.d = 20e-6;
    chip(7).TSV.liner = 0e-6;
    chip(7).TSV.material = [4, 5];
    chip(7).TSV.px = 200e-6;
    chip(7).TSV.py = 200e-6;  
    %% power information
    chip(7).power = [3.75];  %power dissipation of each die
    chip(7).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(7).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(3).blk_name = cell(12, 1);
%     chip(3).blk_name(:) = cellstr('');
%     chip(3).blk_name(3) = cellstr('Cache');
%     chip(3).blk_name(5) = cellstr('Core');
%     chip(3).blk_name(8) = cellstr('Core');
    chip(7).name = cellstr('CHIP # 7');    
    chip(7).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #7%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%for die #8%%%%%%%%%%%%%%%%%%%%%%%%
    chip(8).layer.id = 3;
    chip(8).model = 4; % for each die, how many layers we model
    
    chip8_mold_thickness = 50e-6;
    chip8_thickness = 50e-6;
    
    chip(8).layer.mold_b = [chip8_mold_thickness, 5];
    chip(8).layer.die = [chip8_thickness, 2]; %die thickness
    chip(8).layer.ild = [50e-6, 5]; %metal layer thickness
    chip(8).layer.under = [50e-6, 3]; %underfill bonding thickness; between two dies
    chip(8).metal_portion = 0.5*0.3; %copper portion of the metal layer
    
    %% die size
    chip(8).Xsize = 2.5e-3; %x dimension of chip 1
    chip(8).Ysize = 3e-3; %y dimension of chip 1
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %% Bump geometry
%     chip(8).bump.px = 1e-6;
%     chip(8).bump.py = 1e-6;
    chip(8).bump.d = 40e-6;
    chip(8).bump.px = 200e-6;
    chip(8).bump.py = 200e-6;
%     chip(8).bump.map = [0     0  system.pack.Xsize system.pack.Ysize];  
%     chip(8).bump.map = [ ];  
    chip(8).bump.material = 6;
    
%     chip(8).TSV.map = [0     0  chip(8).Xsize chip(8).Ysize];
    chip(8).TSV.d = 20e-6;
    chip(8).TSV.liner = 0e-6;
    chip(8).TSV.material = [4, 5];
    chip(8).TSV.px = 200e-6;
    chip(8).TSV.py = 200e-6;  
    %% power information
    chip(8).power = [2.75];  %power dissipation of each die
    chip(8).p_metal = [0.05];
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    chip(8).map = [ ];
    %blk_num is for splitting the power maps of each die
%     chip(3).blk_name = cell(12, 1);
%     chip(3).blk_name(:) = cellstr('');
%     chip(3).blk_name(3) = cellstr('Cache');
%     chip(3).blk_name(5) = cellstr('Core');
%     chip(3).blk_name(8) = cellstr('Core');
    chip(8).name = cellstr('CHIP # 8');    
    chip(8).blk_num = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%for die #8%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%

    %%%%%%%%%%%%%%%%%%die, interposer, package information%%%%%%%%%%%
    %% adopt the method from PDN modeling
    %use the heat spreader as the reference plane
    system.pack.xl = (system.hsp.Xsize - system.pack.Xsize)/2;
    system.pack.yb = (system.hsp.Ysize - system.pack.Ysize)/2;
    
    chip(1).xl = system.pack.xl;
    chip(1).yb = system.pack.yb;
    
    chip(2).xl = system.pack.xl;
    chip(2).yb = system.pack.yb;
    
    chip(3).xl = system.pack.xl + 0.5e-3;
    chip(3).yb = system.pack.yb + 1e-3;
    
    chip(4).xl = chip(3).xl + chip(3).Xsize + 0.5e-3;
    chip(4).yb = chip(3).yb;
    
    chip(5).xl = chip(4).xl + chip(4).Xsize + 0.5e-3;
    chip(5).yb = chip(4).yb;
    
    chip(6).xl = system.pack.xl + 0.5e-3;
    chip(6).yb = chip(3).yb + chip(3).Ysize + 1.5e-3;
    
    chip(7).xl = chip(6).xl + chip(6).Xsize + 0.5e-3;
    chip(7).yb = chip(6).yb;
    
    chip(8).xl = chip(7).xl + chip(7).Xsize + 0.5e-3;
    chip(8).yb = chip(6).yb;
    
    %% Update bump maps
    chip(1).bump.map = [chip(1).xl chip(1).yb chip(1).Xsize chip(1).Ysize];
    
    chip(2).bump.map = [chip(2).xl chip(2).yb chip(2).Xsize chip(2).Ysize];
    
    chip(3).bump.map = [chip(3).xl chip(3).yb chip(3).Xsize chip(3).Ysize;
                       chip(4).xl chip(4).yb chip(4).Xsize chip(4).Ysize;
                       chip(5).xl chip(5).yb chip(5).Xsize chip(5).Ysize;
                       chip(6).xl chip(6).yb chip(6).Xsize chip(6).Ysize;
                       chip(7).xl chip(7).yb chip(7).Xsize chip(7).Ysize;
                       chip(8).xl chip(8).yb chip(8).Xsize chip(8).Ysize];
    
    chip(3).TSV.map = [chip(3).xl chip(3).yb chip(3).Xsize chip(3).Ysize;
                       chip(4).xl chip(4).yb chip(4).Xsize chip(4).Ysize;
                       chip(5).xl chip(5).yb chip(5).Xsize chip(5).Ysize;
                       chip(6).xl chip(6).yb chip(6).Xsize chip(6).Ysize;
                       chip(7).xl chip(7).yb chip(7).Xsize chip(7).Ysize;
                       chip(8).xl chip(8).yb chip(8).Xsize chip(8).Ysize];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% the thermal resistivity of each boundary
     % r = 1/(hA); A is the size of top surface area
    % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
    % 0.5 W*cm^2/K
    h.up = 1/(2.25e-4);
%     h.up = 1/(100e-4);
%     h.up = 1/(0.63*1e-4);
%     h.up = 1/(0.03*1e-4);
    
    r_liquid = 1000*1e-4; %K-m2/W
%     r_liquid = 0.05*2000*1e-4; %K-m2/W
%     r_liquid = 0.25*1e-4; %K-m2/W
    h.down = 1/(r_liquid); % the cooling of bottom surface; 
    %(only the area with the same size of chip;)
    %microfluidic is assumed to be as large as chip in the interposer
    
    h.side = 0;
    % side surface cooling, usually near adiabatic
    
    h.d = 1/(r_liquid);
    %the cooling of the bottom surface except for the MFHS area
    
    h.Ta = 293.15;
    %the ambient temperature; for CPU is 38C which is 311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


%% assigning embedded power dies

% chip(3).power = [13.3];
% chip(4).power = [50];
% chip(5).power = [26.67];
% chip(6).power = [50];
% chip(7).power = [50];
% chip(8).power = [36.67];

chip(1).power = [100];
chip(2).power = [5];
chip(3).power = [3.75];
chip(4).power = [3.75];
chip(5).power = [3.75];
chip(6).power = [3.75];
chip(7).power = [37.5];
chip(8).power = [3.75];

%% main program to run
T = ThermSim(system, chip, Material, draw, h);
    

%     % main program to run
% %     for p_d7 = [0.375, 0.75, 1.125, 1.5, 1.875, 2.25, 3.75, 5.625, 7.5, 15, 22.5, 30, 37.5]
%     for p_d7 = [0.0375, 0.075, 0.15, 0.225, 0.3]
% %         disp('die 7 power = ',p_d7);
%         disp(p_d7);
%         chip(7).power = [p_d7]; %chip 7 power
%         T = ThermSim(system, chip, Material, draw, h);
% %         disp('Min Max Temperatures: ',T);
% %         disp(T);
%     end
    
%     for d_t3 = [1e-6, 2e-6, 5e-6, 10e-6, 25e-6, 50e-6, 100e-6, 150e-6, 200e-6, 250e-6, 300e-6]
% %         disp('die thickness = ',t)
%         disp(d_t3);
%         chip(3).layer.die(1) = d_t3; %die thickness
% %         chip(4).layer.die(1) = d_t3;
% %         chip(5).layer.die(1) = d_t3;
% %         chip(6).layer.die(1) = d_t3;
% %         chip(7).layer.die(1) = d_t3;
% %         chip(8).layer.die(1) = d_t3;
%         T = ThermSim(system, chip, Material, draw, h);
% %         disp('Min Max Temperatures: ',T);
% %         disp(T);
%     end

%% Mesh size variation

for m_hsp = [5.00e-06 1.00e-05 1.50e-05]
    system.hsp.Xgrid = m_hsp;
    system.hsp.Ygrid = m_hsp;
    for m_pack = [5.00e-06 1.00e-05 1.50e-05]
        %         disp('die 7 power = ',p_d7);
        disp(m_hsp);
        disp(m_pack);
        system.pack.Xgrid = m_pack;
        system.pack.Ygrid = m_pack;
        system.chip.Xgrid = m_pack; %x grid size of chip
        system.chip.Ygrid = m_pack;
        T = ThermSim(system, chip, Material, draw, h);
        %         disp('Min Max Temperatures: ',T);
        %         disp(T);
    end
end