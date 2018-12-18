function [delta_T, L, x0, u0, num_timesteps, Qtrue, Rtrue, yreal] = get_common_variables()
    delta_T = 0.1; L = 0.5;
    x0 = [10; 0; pi/2; -60; 0; -pi/2];
    u0 = [2; -pi/18; 12; pi/25];
    num_timesteps = 1000;
    given = load('cooplocalization_finalproj_KFdata.mat');
    Qtrue = given.Qtrue; Rtrue = given.Rtrue;
    yreal = given.ydata;
end