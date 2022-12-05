    %clear;
    close all;
    clc;
    %% Input excel details
    system_sheet = 'TSV_stacked_3D_v4/air_cooling/inputs/system.csv';
    draw_sheet = 'TSV_stacked_3D_v4/air_cooling/inputs/draw.csv';
    bc_sheet = 'TSV_stacked_3D_v4/air_cooling/inputs/thermal_boundary_conditions.csv';
    block_level_specs_sheet = 'TSV_stacked_3D_v4/air_cooling/inputs/block_layout_specs';
    
    %% system parameters
    system.chip.N = readvars(system_sheet, 'Range', 'D10:D10');
    system.layer.N = readvars(system_sheet, 'Range', 'D14:D14');
    
    %% for the heat spreader, package dimensions
    system.hsp.Xsize = readvars(system_sheet, 'Range', 'D2:D2');
    system.hsp.Ysize = readvars(system_sheet, 'Range', 'D3:D3');
    system.hsp.Xgrid = readvars(system_sheet, 'Range', 'D4:D4');
    system.hsp.Ygrid = readvars(system_sheet, 'Range', 'D5:D5');
    
    system.pack.Xsize = readvars(system_sheet, 'Range', 'D6:D6');
    system.pack.Ysize = readvars(system_sheet, 'Range', 'D7:D7');
    system.pack.Xgrid = readvars(system_sheet, 'Range', 'D8:D8');
    system.pack.Ygrid = readvars(system_sheet, 'Range', 'D9:D9');

    %% for the chip grid information; size is assigned by each chip
    system.chip.Xgrid = readvars(system_sheet, 'Range', 'D11:D11'); %x grid size of chip
    system.chip.Ygrid = readvars(system_sheet, 'Range', 'D12:D12'); %y grid size of chip
    
    %% thickness, material identifier
    system.layer.hsp = readmatrix(system_sheet,'Range','D15:E15');
    system.layer.tim2 = readmatrix(system_sheet,'Range','D16:E16');
    system.layer.pack = readmatrix(system_sheet,'Range','D17:E17');
    system.layer.metal_portion = readvars(system_sheet, 'Range', 'D18:D18');% 5 routing layers, 60um per routing layer, total 1600, 50% metal per layer
    system.num.hsp = readvars(system_sheet, 'Range', 'D20:D20'); %the layer of hsp modeled
    system.num.tim2 = readvars(system_sheet, 'Range', 'D21:D21'); %the layer of TIM2 modeled
    system.num.pack = readvars(system_sheet, 'Range', 'D23:D23'); %the layer of TIM2 modeled
    
   %% transient analysis parameter
    system.tran = readvars(system_sheet, 'Range', 'D24:D24');
    system.T = readvars(system_sheet, 'Range', 'D25:D25');
    system.dt = readvars(system_sheet, 'Range', 'D26:D26');
    system.su = readvars(system_sheet, 'Range', 'D27:D27');
    system.su_limit = readvars(system_sheet, 'Range', 'D28:D28');
    system.Ntime = system.T/system.dt + 1;
    
    %% Material Parameters
    Material.K = readmatrix('material_properties.csv','Range','B2:B11');
    Material.C = readmatrix('material_properties.csv','Range','C2:C11');
    Material.D = readmatrix('material_properties.csv','Range','D2:D11');
    
    %% figure controller
    draw.granularity = readvars(draw_sheet, 'Range', 'C2:C2');
    % thermal map, number of color used
    draw.write = readvars(draw_sheet, 'Range', 'C3:C3');
    %write chip temperature files
    draw.T = readvars(draw_sheet, 'Range', 'C4:C4');
    % whether to draw the thermal map; 1 yes; 0 no
    draw.P = readvars(draw_sheet, 'Range', 'C5:C5');
    draw.pdrange = readmatrix(draw_sheet,'Range', 'C6:D6');
    draw.pdclamp = readvars(draw_sheet, 'Range', 'C7:C7');
    % whether to draw the power maps; 1 yes, 0 no
    % only draw the die with blk_num > 1 (non-uniform cases)
    draw.C = readvars(draw_sheet, 'Range', 'C8:C8');
    %whether to draw the conductivity distribution
    draw.displayT = readvars(draw_sheet, 'Range', 'C9:C9');
    % print the temperature information
    draw.absolutely = readvars(draw_sheet, 'Range', 'C10:C10');
    %print Tjunc - Ta or Tjunc; 1 for Tjunc
    draw.gif = readvars(draw_sheet, 'Range', 'C11:C11');
    draw.range = readmatrix(draw_sheet,'Range','C12:D12');
    draw.clamp = readvars(draw_sheet, 'Range', 'C13:C13');

    %% Die layer information
    % Loop through all chip files
    % Must follow naming convention: chip#.csv
    for i = 1:system.chip.N
        myfilename = sprintf('TSV_stacked_3D_v4/air_cooling/inputs/chips/chip%d.csv', i);
        %% for die1 parameter
        chip(i).model = readvars(myfilename,'Range', 'D4:D4'); % for each die, how many layers we model
        chip(i).layer.id = readvars(myfilename,'Range', 'D6:D6');
        
        % flip chip package; order:
        %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
        %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
        %chip thick format: thickness, (material identifier #)
        %% die layer information
        if i == 1
            chip(i).layer.tim = readmatrix(myfilename, 'Range', 'D7:E7');
        end
        
        chip(i).layer.die = readmatrix(myfilename, 'Range', 'D8:E8'); %die thickness
        chip(i).layer.ild = readmatrix(myfilename, 'Range', 'D9:E9'); %metal layer thickness
        chip(i).layer.under = readmatrix(myfilename, 'Range', 'D11:E11'); %underfill bonding thickness; between two dies
        chip(i).metal_portion = readvars(myfilename,'Range', 'D12:D12'); %copper portion of the metal layer
        
        %% die size
        chip(i).Xsize = readvars(myfilename,'Range', 'D14:D14'); %x dimension of chip 1
        chip(i).Ysize = readvars(myfilename,'Range', 'D15:D15'); %y dimension of chip 1
        %assumed TSV starting from top metal layer of a top die
        %            to the first metal layer of a bottom die
        %Thus the TSV passes through bonding layer & bulk of a botom die
        %% Bump geometry
        chip(i).bump.d = readvars(myfilename,'Range', 'D17:D17'); 
        chip(i).bump.px = readvars(myfilename,'Range', 'D18:D18'); 
        chip(i).bump.py = readvars(myfilename,'Range', 'D19:D19'); 
        chip(i).bump.map = readmatrix(myfilename,'Range', 'D20:G20');
        chip(i).bump.material = readvars(myfilename,'Range', 'E21:E21');
        
        %% TSV geometry
        if i > 1
            chip(i).TSV.d = readvars(myfilename,'Range', 'D23:D23'); % tsv diameter including the liner thickness
            chip(i).TSV.liner = readvars(myfilename,'Range', 'D24:D24'); %liner thickness
            chip(i).TSV.px = readvars(myfilename,'Range', 'D25:D25'); %x direction pitch
            chip(i).TSV.py = readvars(myfilename,'Range', 'D26:D26'); %y direction pitch
            chip(i).TSV.map = readmatrix(myfilename,'Range', 'D27:G27');
            chip(i).TSV.material = readmatrix(myfilename,'Range', 'D28:E28');
        end

        %% power information
        chip(i).power = readvars(myfilename,'Range', 'D30:D30');   %power dissipation of each die
        chip(i).p_metal = readvars(myfilename,'Range', 'D31:D31'); 
        % from top to bottom; unit: watt

        %power blocks by each die
        %format: bottom left point bl_x, bl_y, width, height, power
        %list blocks in die1 and then die2, die3 ....
        chip(i).name = cellstr(readmatrix(myfilename,'OutputType', 'string','Range','D32:D32'));
        chip(i).name2 = readmatrix(myfilename,'OutputType', 'string','Range','D32:D32');
        chip(i).map = [ readmatrix(myfilename,'Range','D34:I43') ];
        chip(i).blk_num = readvars(myfilename,'Range', 'D35:D35'); 
        
    end
    
        %% %%%%%%%%%%%%%%%%die, interposer, package information%%%%%%%%%%%
    %% adopt the method from PDN modeling
    %use the heat spreader as the reference plane
    system.pack.xl = (system.hsp.Xsize - system.pack.Xsize)/2;
    system.pack.yb = (system.hsp.Ysize - system.pack.Ysize)/2;
    
    for i = 1:system.chip.N
        chip(i).xl = (system.hsp.Xsize - chip(i).Xsize)/2;
        chip(i).yb = (system.hsp.Ysize - chip(i).Ysize)/2;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% the thermal resistivity of each boundary
    % r = 1/(hA); A is the size of top surface area
    % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
    % 0.5 W*cm^2/K
    h.up = readvars(bc_sheet, 'Range', 'C2:C2');
    
    h.down = readvars(bc_sheet, 'Range', 'C3:C3'); % the cooling of bottom surface; 
    %(only the area with the same size of chip;)
    %microfluidic is assumed to be as large as chip in the interposer
    
    h.side = readvars(bc_sheet, 'Range', 'C4:C4');
    % side surface cooling, usually near adiabatic
    
    h.d = readvars(bc_sheet, 'Range', 'C5:C5');
    %the cooling of the bottom surface except for the MFHS area
    
    h.Ta = readvars(bc_sheet, 'Range', 'C6:C6');
    %the ambient temperature; for CPU is 38C which is 311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simple case
cell_count = 4;
sim_case = [1];
substrate_thick = readvars(block_level_specs_sheet, 'Range', 'B3:B3'); %units: m
tsv_diameter = readvars(block_level_specs_sheet, 'Range', 'B4:B4'); %units: m

%scaled powers
% RRAM
chip_size = readvars(block_level_specs_sheet, 'Range', 'B7:B7'); %units: m
x_TSV = readvars(block_level_specs_sheet, 'Range', 'B8:B8'); %m
pwr_logic = readvars(block_level_specs_sheet, 'Range', 'B9:B9'); % W
pwr_mem = readvars(block_level_specs_sheet, 'Range', 'B10:B10'); % W

% ADC variables
x_margin_adc = readvars(block_level_specs_sheet, 'Range', 'B14:B14');
y_margin_adc = readvars(block_level_specs_sheet, 'Range', 'B15:B15');
nADC = readvars(block_level_specs_sheet, 'Range', 'B16:B16');
nyADC = readvars(block_level_specs_sheet, 'Range', 'B17:B17');
nxADC = readvars(block_level_specs_sheet, 'Range', 'B18:B18');
block_dimension_adc = readvars(block_level_specs_sheet, 'Range', 'B19:B19');
y_pitch_adc = readvars(block_level_specs_sheet, 'Range', 'B20:B20');
x_pitch_adc = readvars(block_level_specs_sheet, 'Range', 'B21:B21');
per_ADC_power = readvars(block_level_specs_sheet, 'Range', 'B22:B22');
x_dim_adc = readvars(block_level_specs_sheet, 'Range', 'B23:B23');
pool_gb_xsize_adc = readvars(block_level_specs_sheet, 'Range', 'B24:B24');
y_dim_pool = readvars(block_level_specs_sheet, 'Range', 'B25:B25');
y_dim_gb = readvars(block_level_specs_sheet, 'Range', 'B26:B26');
tsv_xsize_adc = readvars(block_level_specs_sheet, 'Range', 'B27:B27');
pool_power = readvars(block_level_specs_sheet, 'Range', 'B28:B28');
gb_power = readvars(block_level_specs_sheet, 'Range', 'B29:B29');


% MEM variables
x_margin_MEM = readvars(block_level_specs_sheet, 'Range', 'B32:B32');
y_margin_MEM = readvars(block_level_specs_sheet, 'Range', 'B33:B33');
nMEM = readvars(block_level_specs_sheet, 'Range', 'B34:B34');
nyMEM = readvars(block_level_specs_sheet, 'Range', 'B35:B35');
nxMEM = readvars(block_level_specs_sheet, 'Range', 'B36:B36');
block_dimension_mem = readvars(block_level_specs_sheet, 'Range', 'B37:B37');
y_pitch_mem = readvars(block_level_specs_sheet, 'Range', 'B38:B38');
x_pitch_mem = readvars(block_level_specs_sheet, 'Range', 'B39:B39');
per_MEM_power = readvars(block_level_specs_sheet, 'Range', 'B40:B40');
tsv_xsize_mem = readvars(block_level_specs_sheet, 'Range', 'B41:B41');
x_dim_mem = readvars(block_level_specs_sheet, 'Range', 'B42:B42');


for i = 1:1:length(sim_case)
        
    %underfill, Bumps and TSVs
    chip(1).bump.d = chip(1).layer.under(1);
    chip(1).bump.px = chip(1).bump.d*2;
    chip(1).bump.py = chip(1).bump.d*2;
    chip(1).bump.map = [0 0 chip(1).Xsize chip(1).Ysize];
    
    %Tier 2
    chip(2).bump.d = chip(2).layer.under(1);
    chip(2).bump.px = chip(2).bump.d*2;
    chip(2).bump.py = chip(2).bump.d*2;
    chip(2).bump.map = [0 0 chip(2).Xsize chip(2).Ysize];  
    
    % TSV geometry    
    chip(2).TSV.d = tsv_diameter(i); % tsv diameter including the liner thickness
    chip(2).TSV.px = tsv_diameter(i)*2; %x direction pitch
    chip(2).TSV.py = tsv_diameter(i)*2; %y direction pitch
    chip(2).TSV.map = [((chip(2).Xsize-x_TSV(i))/2) 0 x_TSV(i) chip(2).Ysize];
    
        
    % Substrate thickness
    chip(2).layer.die(1) = substrate_thick(i); %die thickness

        for j = 1:system.chip.N
            if chip(j).name2 == 'Memory Tier'
                %Memory (chip 1) power map definition
                chip(j).power = pwr_mem;
        %         chip(1).map = [];
                chip(j).map = logic_power_map_gen(x_margin_MEM, y_margin_MEM, nyMEM, nxMEM, block_dimension_mem, y_pitch_mem, x_pitch_mem, per_MEM_power);
                B = logic_power_map_gen((x_dim_mem/2 + tsv_xsize_mem + x_margin_MEM), y_margin_MEM, nyMEM, nxMEM, block_dimension_mem, y_pitch_mem, x_pitch_mem, per_MEM_power);
                chip(j).map = [chip(j).map ; B];
                chip(j).blk_name = cell(nMEM*2, 1);
                chip(j).blk_name(:) = cellstr('');
                chip(j).blk_num = [nMEM*2];
            else
                %Logic (chip 2) power map definition
                chip(j).power = pwr_logic;
        %         chip(2).map = [];
                chip(j).map = logic_power_map_gen(x_margin_adc, y_margin_adc, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
                A = logic_power_map_gen((x_dim_adc/2 + pool_gb_xsize_adc + tsv_xsize_adc + x_margin_adc), y_margin_adc, nyADC, nxADC, block_dimension_adc, y_pitch_adc, x_pitch_adc, per_ADC_power);
                chip(j).map = [chip(j).map ; A];
                P1 = [(x_dim_adc/2) 0 (pool_gb_xsize_adc/2) y_dim_pool pool_power 0];
                GB1 = [(x_dim_adc/2) (y_dim_pool) (pool_gb_xsize_adc/2) y_dim_gb gb_power 0];
                P2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc) (y_dim_gb) (pool_gb_xsize_adc/2) y_dim_pool pool_power 0];
                GB2 = [(x_dim_adc/2 + pool_gb_xsize_adc/2 + tsv_xsize_adc) 0 (pool_gb_xsize_adc/2) y_dim_gb gb_power 0];
                chip(j).map = [chip(2).map ; P1 ; GB1 ; P2 ; GB2];
                chip(j).blk_name = cell((nADC*2)+4, 1);
                chip(j).blk_name(:) = cellstr('');
                chip(j).blk_num = [(nADC*2)+4];
            end
        end

        fprintf('\n\nDie size: %12.3e m\n', chip_size(i));
        fprintf('Die area: %12.3e m2\n', chip(1).Xsize*chip(1).Ysize);
        fprintf('Die bulk thickness: %12.3e m\n', substrate_thick(i));
        fprintf('TSV diameter: %12.3e m\n', chip(2).TSV.d);
        fprintf('TSV pitch (x and y): %12.3e m\n', chip(2).TSV.px);
        fprintf('ubump diameter: %12.3e m\n', chip(2).bump.d);
        fprintf('ubump pitch (x and y): %12.3e m\n', chip(2).bump.px);
        
        total_power = 0;
        for j = 1:system.chip.N
            fprintf('Power for Chip %d: %12.3e W\n', j, chip(j).power);
            total_power = total_power + chip(j).power;
        end
        fprintf('Total stack power: %12.3e W\n', total_power);

        T = ThermSim(system, chip, Material, draw, h, cell_count, sim_case(i));

end

%%
% diary off;
