function vp = generate_locomotor_race_final(makevideo)

load sysplotter_config.mat

info_needed = struct('Framerate',30,...
                     'Duration',3,...
                     'Coordinates','minperturbation_coords',...
                     'sysplotterpath', sysplotterpath,...
                     'datapath', datapath,...
                     'UserFile_path', [syspath '/..'],...
                     'syspath', syspath,...
                     'current_system2', {{}},...
                     'current_shch2', {{}},...
                     'optcostfunction',{{}});
                 
info_needed.Movies = struct('system', 1,...
           'CCFx', 0,...
           'CCFy', 0,...
       'CCFtheta', 0,...
        'vfieldx', 0,...
        'vfieldy', 0,...
    'vfieldtheta', 0);


% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% seedlist = {'opt_181021dcd641302d8c00a64cae111ac7','opt_a5004c97c4eebdb9b2658e839901b216','opt_36f65157e963179771c57115cf0271bb'};
% optcostlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};
% 
% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% seedlist = {'opt_81e2fb7d09f4b17ab39587bc9fb83f59','opt_8045ea7a2eea0e4a3f6be1ed182ddf27','opt_07a0c5552b3d0852d5f87f6fcad6d30c'};
% optcostlist = {'pathlength metric','pathlength metric','pathlength metric'};

% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% seedlist = {'opt_81e2fb7d09f4b17ab39587bc9fb83f59','opt_8045ea7a2eea0e4a3f6be1ed182ddf27_2','opt_07a0c5552b3d0852d5f87f6fcad6d30c_2'};
% optcostlist = {'pathlength metric','pathlength metric','pathlength metric'};

% 3-, 4-, and 5-link swimmers, pathlength
vp.links.pathlength.systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe','six_link_HighRe'};
vp.links.pathlength.seedlist = {'circle_1p0_2serial','circle_1p0_3serial','circle_1p0_4serial','circle_1p0_5serial'};
vp.links.pathlength.optcostlist = {'pathlength metric','pathlength metric','pathlength metric','pathlength metric'};

% 3-, 4-, and 5-link swimmers, covariant acceleration
vp.links.covaccel.systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe','six_link_HighRe'};
vp.links.covaccel.seedlist = {'circle_1p0_2serial','circle_1p0_3serial','circle_1p0_4serial','circle_1p0_5serial'};
vp.links.covaccel.optcostlist = {'covariant acceleration','covariant acceleration','covariant acceleration','covariant acceleration'};

% Comparisons of link systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.links.comp1.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp1.seedlist = {vp.links.pathlength.seedlist{sys_idx},vp.links.covaccel.seedlist{sys_idx}};
vp.links.comp1.optcostlist = {'pathlength metric','covariant acceleration'};
vp.links.comp1.evalcostlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.links.comp2.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp2.seedlist = {vp.links.pathlength.seedlist{sys_idx},vp.links.covaccel.seedlist{sys_idx}};
vp.links.comp2.optcostlist = {'pathlength metric','covariant acceleration'};
vp.links.comp2.evalcostlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 3;
vp.links.comp3.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp3.seedlist = {vp.links.pathlength.seedlist{sys_idx},vp.links.covaccel.seedlist{sys_idx}};
vp.links.comp3.optcostlist = {'pathlength metric','covariant acceleration'};
vp.links.comp3.evalcostlist = {'covariant acceleration','covariant acceleration'};


% 2-, 3-, and 4-segment constant-curvature swimmers, pathlength
vp.cc.pathlength.systemlist = {'two_mode_CC','three_mode_CC','four_mode_CC','five_mode_CC'};
vp.cc.pathlength.seedlist = {'circle_4p0_2serial','circle_4p0_3serial','circle_4p0_4serial','circle_4p0_5serial'};
vp.cc.pathlength.optcostlist = {'pathlength metric','pathlength metric','pathlength metric','pathlength metric'};


% 2-, 3-, and 4-segment constant-curvature swimmers, covariant acceleration
vp.cc.covaccel.systemlist = {'two_mode_CC','three_mode_CC','four_mode_CC','five_mode_CC'};
vp.cc.covaccel.seedlist = {'circle_4p0_2serial','circle_4p0_3serial','circle_4p0_4serial','circle_4p0_5serial'};
vp.cc.covaccel.optcostlist = {'covariant acceleration','covariant acceleration','covariant acceleration','covariant acceleration'};



% Comparisons of piecewise-constant systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.cc.comp1.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp1.seedlist = {vp.cc.pathlength.seedlist{sys_idx},vp.cc.covaccel.seedlist{sys_idx}};
vp.cc.comp1.optcostlist = {'pathlength metric','covariant acceleration'};
vp.cc.comp1.evalcostlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.cc.comp2.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp2.seedlist = {vp.cc.pathlength.seedlist{sys_idx},vp.cc.covaccel.seedlist{sys_idx}};
vp.cc.comp2.optcostlist = {'pathlength metric','covariant acceleration'};
vp.cc.comp2.evalcostlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 3;
vp.cc.comp3.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp3.seedlist = {vp.cc.pathlength.seedlist{sys_idx},vp.cc.covaccel.seedlist{sys_idx}};
vp.cc.comp3.optcostlist = {'pathlength metric','covariant acceleration'};
vp.cc.comp3.evalcostlist = {'covariant acceleration','covariant acceleration'};



% 2- and 4-mode serpenoid swimmers, pathlength
vp.serp.pathlength.systemlist = {'two_mode_serpenoid','four_mode_serpenoid'};
vp.serp.pathlength.seedlist = {'circle_6p0_2oddeven','circle_6p0_4oddeven'};
vp.serp.pathlength.optcostlist = {'pathlength metric','pathlength metric'};


% 2- and 4-mode serpenoid swimmers, covariant acceleration
vp.serp.covaccel.systemlist = {'two_mode_serpenoid','four_mode_serpenoid'};
vp.serp.covaccel.seedlist = {'circle_6p0_2oddeven','circle_6p0_4oddeven'};
vp.serp.covaccel.optcostlist = {'covariant acceleration','covariant acceleration'};

% Comparisons of piecewise-constant systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.serp.comp1.systemlist = {vp.serp.pathlength.systemlist{sys_idx}, vp.serp.covaccel.systemlist{sys_idx}};
vp.serp.comp1.seedlist = {vp.serp.pathlength.seedlist{sys_idx},vp.serp.covaccel.seedlist{sys_idx}};
vp.serp.comp1.optcostlist = {'pathlength metric','covariant acceleration'};
vp.serp.comp1.evalcostlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.serp.comp2.systemlist = {vp.serp.pathlength.systemlist{sys_idx}, vp.serp.covaccel.systemlist{sys_idx}};
vp.serp.comp2.seedlist = {vp.serp.pathlength.seedlist{sys_idx},vp.serp.covaccel.seedlist{sys_idx}};
vp.serp.comp2.optcostlist = {'pathlength metric','covariant acceleration'};
vp.serp.comp2.evalcostlist = {'covariant acceleration','covariant acceleration'};

% vp.serp.deg.systemlist = {vp.serp.pathlength.systemlist{1}, vp.serp.covaccel.systemlist{2}};
% vp.serp.deg.seedlist = {vp.serp.covaccel.seedlist{1},vp.serp.covaccel.seedlist{1}};
% vp.serp.deg.optcostlist = {'covariant acceleration','covariant acceleration'};
% 
% Comparisons of two-mode performance, pathlength
vp.twomodes.pathlength.systemlist = {vp.links.pathlength.systemlist{1}, vp.cc.pathlength.systemlist{1}, vp.serp.pathlength.systemlist{1}};
vp.twomodes.pathlength.seedlist = {vp.links.pathlength.seedlist{1}, vp.cc.pathlength.seedlist{1}, vp.serp.pathlength.seedlist{1}};
vp.twomodes.pathlength.optcostlist = {'pathlength metric','pathlength metric','pathlength metric'};

% Comparisons of two-mode performance, covariant acceleration
vp.twomodes.covaccel.systemlist = {vp.links.covaccel.systemlist{1}, vp.cc.covaccel.systemlist{1}, vp.serp.covaccel.systemlist{1}};
vp.twomodes.covaccel.seedlist = {vp.links.covaccel.seedlist{1}, vp.cc.covaccel.seedlist{1}, vp.serp.covaccel.seedlist{1}};
vp.twomodes.covaccel.optcostlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};


% Comparisons of four-mode performance, pathlength
vp.fourmodes.pathlength.systemlist = {vp.links.pathlength.systemlist{3}, vp.cc.pathlength.systemlist{3}, vp.serp.pathlength.systemlist{2}};
vp.fourmodes.pathlength.seedlist = {vp.links.pathlength.seedlist{3}, vp.cc.pathlength.seedlist{3}, vp.serp.pathlength.seedlist{2}};
vp.fourmodes.pathlength.optcostlist = {'pathlength metric','pathlength metric','pathlength metric'};

% Comparisons of four-mode performance, covariant acceleration
vp.fourmodes.covaccel.systemlist = {vp.links.covaccel.systemlist{3}, vp.cc.covaccel.systemlist{3}, vp.serp.covaccel.systemlist{2}};
vp.fourmodes.covaccel.seedlist = {vp.links.covaccel.seedlist{3}, vp.cc.covaccel.seedlist{3}, vp.serp.covaccel.seedlist{2}};
vp.fourmodes.covaccel.optcostlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};


syslist = fieldnames(vp);

for idx = 1:numel(syslist)
    sys = syslist{idx};
    
    optcostlist = fieldnames(vp.(syslist{idx}));
    
    for idx2 = 1:numel(optcostlist)
        cost = optcostlist{idx2};

        info_needed.current_system2 = vp.(sys).(cost).systemlist;
        info_needed.current_shch2 = vp.(sys).(cost).seedlist;
        info_needed.optcostfunction = vp.(sys).(cost).optcostlist;
        if isfield(vp.(sys).(cost),'evalcostlist')
            info_needed.evalcostfunction = vp.(sys).(cost).evalcostlist;
        else
            info_needed.evalcostfunction = vp.(sys).(cost).optcostlist;
        end
        info_needed.moviename = [sys '_' cost];


        [vp.(sys).(cost).normalizedPeriod, vp.(sys).(cost).netDisp] = animate_locomotor_race(0,info_needed,makevideo);
    end
    
end

% vp.links6.pathlength.systemlist = {'six_link_HighRe'};
% vp.links6.pathlength.seedlist = {'opt_7aaf7309edd7d3778771a86c163869a5'};
% vp.links6.pathlength.optcostlist = {'pathlength metric'};
% 
% vp.links5.pathlength.systemlist = {'five_link_HighRe','five_link_HighRe'};
% vp.links5.pathlength.seedlist = {'circle_1p0_4serial','opt_175268cfa591e8fdb4ea2c9348a503b7'};
% vp.links5.pathlength.optcostlist = {'pathlength metric','pathlength metric'};
% 
% vp.serp4.pathlength.systemlist = {'four_mode_serpenoid','four_mode_serpenoid'};
% vp.serp4.pathlength.seedlist = {'circle_6p0_4oddeven','opt_3b51a14e53557544541f4431219ea451'};
% vp.serp4.pathlength.optcostlist = {'pathlength metric','pathlength metric'};
% 

% vp.test.pathlength.systemlist = {'three_mode_CC'};
% vp.test.pathlength.seedlist = {'circle_4p0_3serial'};
% vp.test.pathlength.optcostlist = {'pathlength metric'};
% 
% 
% sys = 'fourmodes';
% cost = 'pathlength';
% 
%         info_needed.current_system2 = vp.(sys).(cost).systemlist;
%         info_needed.current_shch2 = vp.(sys).(cost).seedlist;
%         info_needed.optcostfunction = vp.(sys).(cost).optcostlist;
%         if isfield(vp.(sys).(cost),'evalcostlist')
%             info_needed.evalcostfunction = vp.(sys).(cost).evalcostlist;
%         else
%             info_needed.evalcostfunction = vp.(sys).(cost).optcostlist;
%         end
%         info_needed.moviename = [sys '_' cost];
% 
% 
%         [vp.(sys).(cost).normalizedPeriod, vp.(sys).(cost).netDisp] = animate_locomotor_race(0,info_needed,makevideo);
