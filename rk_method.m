% Generates the Butcher tableau for several Runge-Kutta methods.
% INPUTS: 
% name: the name of the Runge-Kutta method
% OUTPUTS:
% BT_struct: A struct of the Butcher tableau for the specified method
function BT_struct = rk_method(name)
    BT_struct = struct();
    switch name
        case "midpoint"
            BT_struct.A = [0, 0; 0.5, 0];
            BT_struct.B = [0, 1];
            BT_struct.C = [0; 0.5];
        case "kutta3rd"
            BT_struct.A = [0, 0, 0;
                           0.5, 0, 0;
                           -1, 2, 0];
            BT_struct.B = [1/6, 2/3, 1/6];
            BT_struct.C = [0; 0.5; 1];
        case "nystrom5th"
            BT_struct.A = [0, 0, 0, 0, 0, 0; 
                           1/4, 0, 0, 0, 0, 0;
                           4/25, 6/25, 0, 0, 0, 0;
                           1/4, -3, 15/4, 0, 0, 0;
                           2/27, 10/9, -50/81, 8/81, 0, 0;
                           2/25, 12/25, 2/15, 8/75, 0, 0];
            BT_struct.B = [23/192, 0, 125/192, 0, -27/64, 125/192];
            BT_struct.C = [0; 1/3; 2/5; 1; 2/3; 4/5];
        case "dormandprince"
            BT_struct.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
            BT_struct.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
            5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
            BT_struct.A = [0,0,0,0,0,0,0;
            1/5, 0, 0, 0,0,0,0;...
            3/40, 9/40, 0, 0, 0, 0,0;...
            44/45, -56/15, 32/9, 0, 0, 0,0;...
            19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
            9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
            35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
        case "fehlberg"
            BT_struct.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
            BT_struct.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
            25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
            BT_struct.A = [0,0,0,0,0,0;...
            1/4, 0,0,0,0,0;...
            3/32, 9/32, 0,0,0,0;...
            1932/2197, -7200/2197, 7296/2197, 0,0,0;...
            439/216, -8, 3680/513, -845/4104, 0,0;...
            -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
        case "heuneuler"
            BT_struct.C = [0,1];
            BT_struct.B = [1/2,1/2;1,0];
            BT_struct.A = [0,0;1,0];
        case "fehlbergrk"
            BT_struct.C = [0,1/2,1];
            BT_struct.B = [1/512, 255/256, 1/512;...
            1/256, 255/256, 0];
            BT_struct.A = [0,0,0;1/2,0,0;1/256,255/256,0];
        case "bogacki"
            BT_struct.C = [0,1/2, 3/4, 1];
            BT_struct.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
            BT_struct.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];
    end
end