%% testing shifting matrix
rst
c = (1:10);
c = c'*c;
J = fShiftingMatrix(max(size(c, 1)));

u = J*c; %shift once downwards
v = c*(J'); %shift once rightwards
x = J'*c; %shift once upwards
y = c*J; %shift once leftwards

%% triangular matrix
rst
N = 10;
T = triu(ones(N))

%% testing of matrix division
rst
N = 5;
I = eye(N);
J = fShiftingMatrix(N);
K = magic(N);
ANS = K/(I-J);

%% what comes out
rst
N = 5;
R = magic(N);
G = 1./R;

G.*R

%% is operation correct
rst
N = 5;
I = eye(N);
J = fShiftingMatrix(N);
K = magic(N);
A = K - J*K;
B = (I-J)* K;
sum(A(:)==B(:))/length(A(:))

%% test
rst
N = 5;
I = eye(N);
J = fShiftingMatrix(N);
A = I - J;
Ainv = A^-1
A = Ainv^-1;

%% replication test
rst
V = 1:10;
T = ones(length(V), 1);
T*V;

%% test of hardaman multiplicative distribution property
rst
N = 2;
G = magic(N);
U = magic(N);
Rc = rand(1, 1); %magic(N);
Io = magic(N);

% is (G(o)(UnRc))Io same as G(o)(Un(RcIo))
A = G.*(U*(Rc*Io));
% B = G.*(((U*Rc))*(Io));
B = G.*((U./G)*(Rc./G)*(Io./G));
sum(A(:)==B(:))/length(A(:))

% testing if (UnRc)Io is same as (Un(RcIo))
% A = (U*(Rc*Io));
% B = (U*Rc)*Io;
% sum(A(:)==B(:))/length(A(:))

%% test of Matrix distribution property
N = 4;
Io = magic(N);
J = magic(N);
A = Io - (J*Io);
B = (eye(N) - J)*Io;
sum(A(:)==B(:))/length(A(:))

%% haar wavelets
rst
fsamp = 1e3;
fsig = 500;
t = 0:2*pi/(31*fs):2*pi;
x = square(2*t);
figure;
plot(t, x)
title('actual waveform')

figure;
hft = haar_fft(x);
plot(t, hft)
title('haar transform')

%% random resistance value generator
R = 1e6 * randn(3, 3)

%% comparing X methods
N = 3;
Vc_1 = 0.5*ones(N);
Io = ones(N);
U = fUTriangularMatrix(N);
X_1 = Vc_1 + U*Io

X_2 = zeros(N);

%%
Vc_2 = [zeros(N-1, N);
    Vc_1(end,:)]
X_2 = (J'*X_2 + Vc_2) + Io

%% hardman property test
N = 10
X = magic(N); Y = magic(N);
A = X + (Y.*X);
B = ((ones(N) + Y).*X);
sum(A(:)==B(:))/length(A(:))

%% hardaman property test
N = 10
X = magic(N); Y = magic(N); Z = magic(N);
A = X*(Y.*Z);
B = (X.*Y)*(Z);
sum(A(:)==B(:))/length(A(:))

%% second model version 3
rst
N = 3;
%define array
R = 10e3; R = repmat(R, [N, N]); %R(3, 3) = 50;% R(3, 2) = 5e3; R(3, 3) = 30e3;
% R = 10e3*abs(randn(size(R)))
G = 1./R;
% Rc = R; Rr = R;
Rc = 1e3*ones(size(R)); Rr = 1e3*ones(size(R));
%define operation matrices
J = fShiftingMatrix(N);
U = fUTriangularMatrix(N);
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]);
Vs = [vs zeros(N, N-1)];
vc = zeros(N, 1);
Vc = ones(N, 1)*transpose(vc);
Vo = ones(N, 1)*transpose(vs); %Vo = [Vo(:, 1), zeros(N, N-1)]
Vi = Vo*J' + Vs;
X = Vc;
%currents
Is = zeros(N);
Io = zeros(N);
Ii = J*Io;
iter = 0
total_iter = 0;
iter_limit = 2000;
% %% Give answer:
% Io(N,:) = [1.03249e-3 1.03271e-3 2.7e-3];

% %% calculations for second model v2
iter = 0;
hasConverged = false;
while ((~hasConverged) && (iter<iter_limit))||(iter<5)
    Io_prev = Io;
    Vo_prev = Vo;
    Vo = ((Vo*J') + Vs) - Rr.*((Io - (J*Io) + (Is*J)) - Io + (J*Io));
    Io = G.*(Vi -(Vc + U*(Rc.*Io))) + (J*Io);
    Is = (Io - (J*Io) + (Is*J));
    iter = iter + 1;
    %recalculated inputs
    % X = (Vc + U*(Rc.*Io));        % with x method 2 (Vc has all non-zero elements) better!
    % % X = (J'*X) + Vc + Rc.*Io; % with x method 2 (Vc has zero elements)
    % Ii = (J*Io);
    % Vi = ((Vo*J') + Vs);
    % Is = (Io - Ii + (Is*J));
    %check for convergence
    hasConverged = fHasConverged([Io Vo], [Io_prev Vo_prev], 1e-3);
    [Io Vo]
end
total_iter = total_iter + iter
%display results
total_iter_ = fUnits(total_iter, 'iterations')
Vo_ = fUnits(Vo, 'V')
Io_prev_ = fUnits(Io_prev, 'A')
Io_ = fUnits(Io, 'A')


%% Remove whitespace from file in matlab
fileDir = '/Users/IOmotade/Desktop/School/4th_Year/FYP/Sims/';
fid_in  = fopen([fileDir 'out.txt']);
fid_out = fopen('text_out.txt','w+');
text_in = fgetl(fid_in);
while ischar(text_in)
    %     text_out = [regexprep(text_in,' ','')];
    %     text_out = [regexprep(text_in,'Cir(\w*)','')]
    
    %Check if line is empty
    valid = ~isempty(text_in)
    
    %Remove unneccessary lines
    %%Circuit*
    strtIdx = regexp(text_in,'Cir(\w*)')
    valid = valid & isempty(strtIdx);
    
    %%Date*
    strtIdx = regexp(text_in,'Date(\w*)')
    valid = valid & isempty(strtIdx);
    
    %%Date*
    strtIdx = regexp(text_in,'Oper(\w*)')
    valid = valid & isempty(strtIdx);
    
    %%Date*
    strtIdx = regexp(text_in,'Node(\w*)')
    valid = valid & isempty(strtIdx);
    
    text_out = [regexprep(text_in,'\r','a')]
    if  valid
        fprintf(fid_out,'%s\n',text_out);
        %         disp(text_in)
        %     disp(text_out)
    end
    text_in = fgetl(fid_in);
end
% disp(text_in)
disp(text_out)
disp('run')
fclose(fid_in);
fclose(fid_out);

%% Test line removal function
rst
fileDir = '/Users/IOmotade/Desktop/School/4th_Year/FYP/Sims/';
fileName = string([fileDir 'out.txt']);
tmpFileName = string(['../FYP/' 'tmp.txt']);
unwantedLineStarts = ["Cir";
    'Date';
    'Oper';
    'Node';
    '---';
    'Sou';
    'Resist';
    'mode';
    'rsh';
    'narr';
    'tc';
    'def';
    'resistan';
    'p';
    'acm';
    'dc';
    'CPU';
    'Tota';
    ];

copyfile(char(fileName), char(tmpFileName))
fRemoveLines(tmpFileName, unwantedLineStarts)
fRemoveExtraWhitespace(tmpFileName);
tmpFileID = fopen(tmpFileName);
a = textscan(tmpFileID, '%s %f')
fclose(tmpFileID);

%% Read in Simulation Results File

rst
N = 3;
MemR = 10e3*ones(N); LRowR = 1e3*ones(N); LColR = 1e3*ones(N);
%define signals
%voltages
vs = 5; vs = repmat(vs, [N, 1]); %Source Voltage Vector(Nx1)

%Open desired text file
fileDir = '/Users/IOmotade/Google Drive/FYP/';
filename = 'mat_array';
fileLoc = [fileDir, filename, '.cir'];

Circuit = fGenerateSpiceFile(N, vs, MemR, LRowR, LColR, fileLoc);

fileDir = '/Users/IOmotade/Desktop/School/4th_Year/FYP/Sims/';
fileLoc = string([fileDir 'out.txt']);

Circuit = fReadSpiceSimResults(fileLoc, Circuit)

%% Full-Simulation
%% Define Array Envirionment
rst
N = 3;

%resistors
MemR = 10e3*ones(N); LRowR = 1e3*ones(N); LColR = 1e3*ones(N);
%voltages
vs = 5; vs = repmat(vs, [N, 1]); %Source Voltage Vector(Nx1)

% %% fMacSpiceSim(vs, MemR, Rc, Rr)
%Open desired text file
% fileDir = '/Users/IOmotade/Google Drive/FYP/';
filename = 'mat_array';
fileLoc = [filename, '.cir']; %[fileDir, filename, '.cir'];

Circuit = fGenerateSpiceFile(N, vs, MemR, LRowR, LColR, fileLoc);

command = '/Applications/MacSpice.app/Contents/MacOS/MacSpice -b mat_array.cir > out.txt';
status = system(command);
if status ==0
    fileLoc = "out.txt";
end

Circuit = fReadSpiceSimResults(fileLoc, Circuit)

%% Full-Simulation with function
rst
N = 50;
MemR = 10e3*ones(N); LRowR = 1e3*ones(N); LColR = 1e3*ones(N);
vs = 5; vs = repmat(vs, [N, 1]); %Source Voltage Vector(Nx1)
Circuit = fMacSpiceSim(N, vs, MemR, LRowR, LColR)

%% Uniqueness study
for i = 1:50
    for j = 1:50
        val(i, j) = 10*i+j;
    end
end
val = val(:)
size(unique(val)) == size(val)