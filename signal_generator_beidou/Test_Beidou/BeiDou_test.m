% /*!
%  * \file BeiDou_test.m

%  create a file containing the desired B1 code sequence of the selected satellite 
%
%
%  * -------------------------------------------------------------------------

clc;
clear all;
close all;

%% Satellite B1 code

No = 1;  % select the PRN number of the satellite

B1code = generateB1code(No);
B1str = num2str(B1code);

settings = initSettings();

B1CodesTable = makeB1table(settings);


% B1codeT = B1code.';
% B1strT = num2str(B1codeT)

name = strcat('/home/giorgio/Desktop/GSOC_2015/Materiale/Test_BeiDou/data/BeiDou_B1_',num2str(No),'_ranging_code_', num2str(No),'.txt');

f = fopen(name, 'w');

fwrite(f, B1str, 'char');
fclose(f);


%fileID = fopen('nine.bin');
%A = fread(fileID)