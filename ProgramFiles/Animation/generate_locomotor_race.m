function vp = generate_locomotor_race(makevideo)

load sysplotter_config.mat

info_needed = struct('Framerate',15,...
                     'Duration',3,...
                     'Coordinates','minperturbation_coords',...
                     'sysplotterpath', sysplotterpath,...
                     'datapath', datapath,...
                     'UserFile_path', [syspath '/..'],...
                     'syspath', syspath,...
                     'current_system2', {{'three_link_lowRe','three_link_lowRe','serpenoid_lowRe'}},...
                     'current_shch2', {{'LowRe_MaxDisplacement','LowRe_MaxEfficiency','circle_6'}},...
                     'costfunction',{{'pathlength metric','pathlength metric','pathlength metric'}});
                 
info_needed.Movies = struct('system', 1,...
           'CCFx', 0,...
           'CCFy', 0,...
       'CCFtheta', 0,...
        'vfieldx', 0,...
        'vfieldy', 0,...
    'vfieldtheta', 0);


% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% gaitlist = {'opt_181021dcd641302d8c00a64cae111ac7','opt_a5004c97c4eebdb9b2658e839901b216','opt_36f65157e963179771c57115cf0271bb'};
% costlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};
% 
% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% gaitlist = {'opt_81e2fb7d09f4b17ab39587bc9fb83f59','opt_8045ea7a2eea0e4a3f6be1ed182ddf27','opt_07a0c5552b3d0852d5f87f6fcad6d30c'};
% costlist = {'pathlength metric','pathlength metric','pathlength metric'};

% systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe'};
% gaitlist = {'opt_81e2fb7d09f4b17ab39587bc9fb83f59','opt_8045ea7a2eea0e4a3f6be1ed182ddf27_2','opt_07a0c5552b3d0852d5f87f6fcad6d30c_2'};
% costlist = {'pathlength metric','pathlength metric','pathlength metric'};

% 3-, 4-, and 5-link swimmers, pathlength
vp.links.pathlength.systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe','six_link_HighRe'};
vp.links.pathlength.gaitlist = {'opt_81dc07b81f8b5e5976cbdbcfbd5bc222','opt_78744a0b8b59ededa15fcc94a80c8050','opt_175268cfa591e8fdb4ea2c9348a503b7','opt_7aaf7309edd7d3778771a86c163869a5'};
vp.links.pathlength.costlist = {'pathlength metric','pathlength metric','pathlength metric','pathlength metric'};

% 3-, 4-, and 5-link swimmers, covariant acceleration
vp.links.covaccel.systemlist = {'three_link_HighRe','four_link_HighRe','five_link_HighRe','six_link_HighRe'};
vp.links.covaccel.gaitlist = {'opt_074cc8a8fe958909eda886934de4aa0d','opt_e367cce7e4018f9fd1a898ac68837563','opt_d50c438fe7b2469e5ffc67ecd329204a','opt_981d5a890c2f08705c8871d5546b37a0'};
vp.links.covaccel.costlist = {'covariant acceleration','covariant acceleration','covariant acceleration','covariant acceleration'};

% Comparisons of link systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.links.comp1.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp1.gaitlist = {vp.links.pathlength.gaitlist{sys_idx},vp.links.covaccel.gaitlist{sys_idx}};
vp.links.comp1.costlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.links.comp2.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp2.gaitlist = {vp.links.pathlength.gaitlist{sys_idx},vp.links.covaccel.gaitlist{sys_idx}};
vp.links.comp2.costlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 3;
vp.links.comp3.systemlist = {vp.links.pathlength.systemlist{sys_idx}, vp.links.covaccel.systemlist{sys_idx}};
vp.links.comp3.gaitlist = {vp.links.pathlength.gaitlist{sys_idx},vp.links.covaccel.gaitlist{sys_idx}};
vp.links.comp3.costlist = {'covariant acceleration','covariant acceleration'};


% 2-, 3-, and 4-segment constant-curvature swimmers, pathlength
vp.cc.pathlength.systemlist = {'two_mode_CC','three_mode_CC','four_mode_CC','five_mode_CC'};
vp.cc.pathlength.gaitlist = {'opt_14da671eded328d3635cc4b3b8a2e7c2','opt_08a9f975bb9ae77540d7a6ecdf728d3b','opt_74d3b3d293d9707eb8afbca804ee22ff','opt_50381ee845fff0a72bef8c64221e6927'};
vp.cc.pathlength.costlist = {'pathlength metric','pathlength metric','pathlength metric','pathlength metric'};


% 2-, 3-, and 4-segment constant-curvature swimmers, covariant acceleration
vp.cc.covaccel.systemlist = {'two_mode_CC','three_mode_CC','four_mode_CC','five_mode_CC'};
vp.cc.covaccel.gaitlist = {'opt_becaae5749ff7603a88791d0df178ef7','opt_91c1f2fc8d41ae5783b5486b6073aa54','opt_bc286eed659533fc39d891b0dbb19ddb','opt_cfc708f21b8915e131239b2a4c8f0f0b'};
vp.cc.covaccel.costlist = {'covariant acceleration','covariant acceleration','covariant acceleration','covariant acceleration'};



% Comparisons of piecewise-constant systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.cc.comp1.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp1.gaitlist = {vp.cc.pathlength.gaitlist{sys_idx},vp.cc.covaccel.gaitlist{sys_idx}};
vp.cc.comp1.costlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.cc.comp2.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp2.gaitlist = {vp.cc.pathlength.gaitlist{sys_idx},vp.cc.covaccel.gaitlist{sys_idx}};
vp.cc.comp2.costlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 3;
vp.cc.comp3.systemlist = {vp.cc.pathlength.systemlist{sys_idx}, vp.cc.covaccel.systemlist{sys_idx}};
vp.cc.comp3.gaitlist = {vp.cc.pathlength.gaitlist{sys_idx},vp.cc.covaccel.gaitlist{sys_idx}};
vp.cc.comp3.costlist = {'covariant acceleration','covariant acceleration'};



% 2- and 4-mode serpenoid swimmers, pathlength
vp.serp.pathlength.systemlist = {'two_mode_serpenoid','four_mode_serpenoid'};
vp.serp.pathlength.gaitlist = {'opt_33bf8d4335665da9b9362bc55f4a5c4b','opt_0f798d8a3cb62f89a869c7990e77220e'};
vp.serp.pathlength.costlist = {'pathlength metric','pathlength metric'};


% 2- and 4-mode serpenoid swimmers, covariant acceleration
vp.serp.covaccel.systemlist = {'two_mode_serpenoid','four_mode_serpenoid'};
vp.serp.covaccel.gaitlist = {'opt_c7573722a48f8767276a7cdee661e66b','opt_a033f8e1a5a5750d2b8f640a22ecf2dd'};
vp.serp.covaccel.costlist = {'covariant acceleration','covariant acceleration'};

% Comparisons of piecewise-constant systems between pathlength and covariant-acceleration
% gaits

sys_idx = 1;
vp.serp.comp1.systemlist = {vp.serp.pathlength.systemlist{sys_idx}, vp.serp.covaccel.systemlist{sys_idx}};
vp.serp.comp1.gaitlist = {vp.serp.pathlength.gaitlist{sys_idx},vp.serp.covaccel.gaitlist{sys_idx}};
vp.serp.comp1.costlist = {'covariant acceleration','covariant acceleration'};

sys_idx = 2;
vp.serp.comp2.systemlist = {vp.serp.pathlength.systemlist{sys_idx}, vp.serp.covaccel.systemlist{sys_idx}};
vp.serp.comp2.gaitlist = {vp.serp.pathlength.gaitlist{sys_idx},vp.serp.covaccel.gaitlist{sys_idx}};
vp.serp.comp2.costlist = {'covariant acceleration','covariant acceleration'};

% vp.serp.deg.systemlist = {vp.serp.pathlength.systemlist{1}, vp.serp.covaccel.systemlist{2}};
% vp.serp.deg.gaitlist = {vp.serp.covaccel.gaitlist{1},vp.serp.covaccel.gaitlist{1}};
% vp.serp.deg.costlist = {'covariant acceleration','covariant acceleration'};
% 
% Comparisons of two-mode performance, pathlength
vp.twomodes.pathlength.systemlist = {vp.links.pathlength.systemlist{1}, vp.cc.pathlength.systemlist{1}, vp.serp.pathlength.systemlist{1}};
vp.twomodes.pathlength.gaitlist = {vp.links.pathlength.gaitlist{1}, vp.cc.pathlength.gaitlist{1}, vp.serp.pathlength.gaitlist{1}};
vp.twomodes.pathlength.costlist = {'pathlength metric','pathlength metric','pathlength metric'};

% Comparisons of two-mode performance, covariant acceleration
vp.twomodes.covaccel.systemlist = {vp.links.covaccel.systemlist{1}, vp.cc.covaccel.systemlist{1}, vp.serp.covaccel.systemlist{1}};
vp.twomodes.covaccel.gaitlist = {vp.links.covaccel.gaitlist{1}, vp.cc.covaccel.gaitlist{1}, vp.serp.covaccel.gaitlist{1}};
vp.twomodes.covaccel.costlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};


% Comparisons of four-mode performance, pathlength
vp.fourmodes.pathlength.systemlist = {vp.links.pathlength.systemlist{3}, vp.cc.pathlength.systemlist{3}, vp.serp.pathlength.systemlist{2}};
vp.fourmodes.pathlength.gaitlist = {vp.links.pathlength.gaitlist{3}, vp.cc.pathlength.gaitlist{3}, vp.serp.pathlength.gaitlist{2}};
vp.fourmodes.pathlength.costlist = {'pathlength metric','pathlength metric','pathlength metric'};

% Comparisons of four-mode performance, covariant acceleration
vp.fourmodes.covaccel.systemlist = {vp.links.covaccel.systemlist{3}, vp.cc.covaccel.systemlist{3}, vp.serp.covaccel.systemlist{2}};
vp.fourmodes.covaccel.gaitlist = {vp.links.covaccel.gaitlist{3}, vp.cc.covaccel.gaitlist{3}, vp.serp.covaccel.gaitlist{2}};
vp.fourmodes.covaccel.costlist = {'covariant acceleration','covariant acceleration','covariant acceleration'};


% syslist = fieldnames(vp);
% 
% for idx = 1:numel(syslist)
%     sys = syslist{idx};
%     
%     costlist = fieldnames(vp.(syslist{idx}));
%     
%     for idx2 = 1:numel(costlist)
%         cost = costlist{idx2};
% 
%         info_needed.current_system2 = vp.(sys).(cost).systemlist;
%         info_needed.current_shch2 = vp.(sys).(cost).gaitlist;
%         info_needed.costfunction = vp.(sys).(cost).costlist;
%         info_needed.moviename = [sys '_' cost];
% 
% 
%         [vp.(sys).(cost).normalizedPeriod, vp.(sys).(cost).netDisp] = animate_locomotor_race(0,info_needed,makevideo);
%     end
%     
% end

% vp.links6.pathlength.systemlist = {'six_link_HighRe'};
% vp.links6.pathlength.gaitlist = {'opt_7aaf7309edd7d3778771a86c163869a5'};
% vp.links6.pathlength.costlist = {'pathlength metric'};
% 
% vp.links5.pathlength.systemlist = {'five_link_HighRe','five_link_HighRe'};
% vp.links5.pathlength.gaitlist = {'circle_1p0_4serial','opt_175268cfa591e8fdb4ea2c9348a503b7'};
% vp.links5.pathlength.costlist = {'pathlength metric','pathlength metric'};
% 
% vp.serp4.pathlength.systemlist = {'four_mode_serpenoid','four_mode_serpenoid'};
% vp.serp4.pathlength.gaitlist = {'circle_6p0_4oddeven','opt_3b51a14e53557544541f4431219ea451'};
% vp.serp4.pathlength.costlist = {'pathlength metric','pathlength metric'};
% 
sys = 'serp';
cost = 'covaccel';

        info_needed.current_system2 = vp.(sys).(cost).systemlist;
        info_needed.current_shch2 = vp.(sys).(cost).gaitlist;
        info_needed.costfunction = vp.(sys).(cost).costlist;
        info_needed.moviename = [sys '_' cost];


        [vp.(sys).(cost).normalizedPeriod, vp.(sys).(cost).netDisp] = animate_locomotor_race(0,info_needed,makevideo);
