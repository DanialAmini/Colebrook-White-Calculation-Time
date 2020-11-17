clear;
clc;
close all;


name_='virtRome'; %Node: 150634   ,Pipe:  157044   ,10iter:    1.004   sec
%name_='exeter'; %Node: 1893   ,Pipe:  3032   ,10iter:    0.015   sec
%name_='nid1944'; %Node: 3977   ,Pipe:  5704   ,10iter:    0.026   sec
%name_='esf'; %Node: 39913   ,Pipe:  40064   ,10iter:    0.157   sec
%name_='synth1';  %Node: 21058   ,Pipe:  27038   ,10iter:    0.124   sec
%name_='synth3';  %Node: 19742   ,Pipe:  25434   ,10iter:    0.117  sec
name_='pescara'; %node=71 pipe=99 10iter=0.00144 sec
name_='Balerma_irrig'; %node=457   pipe=454  10iter=0.0043


file_pipe=strcat(name_,'_pipe.csv');
file_node=strcat(name_,'_node.csv');
file_reservoir=strcat(name_,'_reservoir.csv');



node_=csvread(file_node,1,0);
node_=int64(node_(:,1));

reservoir_=csvread(file_reservoir,1,0);
reservoir_=int64(reservoir_(:,1));

pipe_=csvread(file_pipe,1,0);
pipe_=int64(pipe_(:,1:3));

[N_node,~]=size(node_);
[N_reservoir,~]=size(reservoir_);
[N_pipe,~]=size(pipe_);


%get rid of unnecessary variables
clear name_ file_reservoir file_pipe file_node;


fprintf('node=%d, reservoir=%d, pipe=%d \n',N_node,N_reservoir,N_pipe);

% pipe_unique=unique(pipe_(:,1));
% [N_pipe,~]=size(pipe_unique);
% fprintf('unique pipe=%d \n',N_pipe);
% clear pipe_unique;

node_unique=unique(pipe_(:,2:3));
[N_node,~]=size(node_unique);
fprintf('unique node (from pipe list ends)=%d \n',N_node);
clear node_unique

node_=[node_;reservoir_]; %put reservoir# at the end of node#
%node_unique=unique(node_);
[N_node,~]=size(node_);
fprintf('node list+reservoir=%d \n',N_node);
clear reservoir_;


[~,idx] = ismember(pipe_(:,2:3),node_);
pipe_=[pipe_,idx];
clear idx;

t=-0.01*rand(N_pipe,1);


b=rand(N_node,1);

A=speye(N_node);
A=A+sparse(pipe_(:,4),pipe_(:,5),t,N_node,N_node);
A=A+sparse(pipe_(:,5),pipe_(:,4),t,N_node,N_node);

A=sparse(A);



fprintf("\n");

fprintf('time of symrcm sort & reorder\n');
tic
for i=1:10
    p= symrcm(A);
    AA=sparse(A(p,p));
    bb=b(p);
end
toc

[lower,upper] = bandwidth(AA);
fprintf('%d = bandwidth of symrcm\n',lower+upper);

fprintf("\n");
fprintf('time of symamd sort & reorder\n');
tic
for i=1:10
    p= symamd(A);
    AAA=sparse(A(p,p));
    bbb=b(p);
end
toc


fprintf("\n");

%warmup
for i=1:10
    x=A\b;
    x=AA\bb;
    x=AAA\bbb;
end

fprintf('time of normal solve\n');
tic
for i=1:10
    x=A\b;
end
toc

fprintf('time of symrcm solve\n');
tic
for i=1:10
    x=AA\bb;
end
toc

fprintf('time of symamd solve\n');
tic
for i=1:10
    x=AAA\bbb;
end
toc

figure();
spy(A);

figure();
spy(AA);

figure();
spy(chol(AA));

figure();
spy(AAA);

figure();
spy(chol(AAA));


