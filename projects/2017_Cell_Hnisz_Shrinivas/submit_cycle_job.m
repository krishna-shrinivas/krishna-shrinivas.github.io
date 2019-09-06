function submit_cycle_job()
tic;

%No of chains
params.Nc = 50;

%Number of modifiable sites
% per chain (valency)
params.f =6;

%Cross-link kinetics

%On-rate
params.kcr_on =0.5;
%Off-rate
params.kcr_off =0.5;

% Number of Modifier and Demodifier molecules
params.Mod =3;

params.Demod = 30;

%Define modification rates/off rates
params.k_mod = 0.05;

params.k_demod = 0.05;

%Define simulation_time
params.t_end =100;


%Flag for changing Modifier levels
%if flag=0, no change happens
% if flag=1, change in modifier levels
% is enforced at t_change to Mod_change

params.change_flag=0;
params.Mod_change =0.5;
params.t_change = 30;


%Define number of trajectories/replicates

Ntraj =1;


%% Define range of parameter
%If you wish to sweep across a parameter range
% for example, a range of valencies or equilibrium constants
% or systems size

%In this particular case, we change valencies
% without any other changes
Nc = [5:1:5];

for i =1:1:length(Nc)
    params.f =Nc(i);
    
    %Pass reaction conditions and get the following parameters returned
    %[Cell of reactions times, where each time{i} has Ntraj different cells
    % and i represents the counter of which parameter is passed
    %trajectory, cluster_max_size returns the cell of largest cluster size
    % Str, Val are returned for the last trajectory only: They are the
    % adjacency matrix of connections and distribution of modified sites
    
    [time{i},cluster_max_size{i},no_of_clusters{i},Str{i},Val{i}]= cyclic_clusters(params,Ntraj);
    
    %Saving the important parameters into a file of choice
    save('test_valency_sweep','Nc','params','time','cluster_max_size','no_of_clusters','Str','Val');
end

%Stores simulation time
sim_time =toc;
save('test_valency_sweep','Nc','params','time','cluster_max_size','no_of_clusters','Str','Val','sim_time');

end

function [t,cluster_max_size,no_of_clusters,Str,Val] = cyclic_clusters(params,Ntraj)

cluster_max_size =[];
t = [];

for Nruns =1:Ntraj
    
    %Define initial cluster structure
    % which has non-modified, non-cross-linked chains
    %No of chains per cluster, net valency per cluster, no of bonds per cluster
    
    Val = zeros(params.Nc,1);
    
    %Adjacency matrix
    Str = eye(params.Nc);
    Str(1:params.Nc+1:end) = 1:1:params.Nc;
    
    t_end =params.t_end;
    t{Nruns}(1) =0;
    
    
    k=1;
    cluster_max_size{Nruns}(1) =1;
    while(t{Nruns}(end)<t_end)
        
        if params.change_flag
            if t{Nruns}(k)>params.t_change
                params.Mod =params.Mod_change;
            end
        end
        
        [Str,Val,deltaT] = cyclic_reaction_cluster(Str,Val,params);
        k = k+1;
        t{Nruns}(k) = t{Nruns}(k-1)+deltaT;
        cluster_max_size{Nruns}(k) =1;
        no_of_clusters(Nruns) = max(diag(Str));
        for p=1:no_of_clusters(Nruns)
            if cluster_max_size{Nruns}(k)<length(find(diag(Str)==p))
                cluster_max_size{Nruns}(k) =length(find(diag(Str)==p));
            end
        end
        
    end
    
end
end

function [Str,Val,deltaT] = cyclic_reaction_cluster(Str,Val,params)

bonds=sum(Str,2)-diag(Str);

for i=1:1:length(Str)
    
    reac_prop{i} = [];
    reac_type{i} = [];
    reac_with{i} = [];
    reac_from{i} = [];
    reac_with_cluster{i} = [];
    reac_from_cluster{i} = [];
    
    %First intra-cluster phosphorylation reactions
    if  Val(i)< params.f
        k_mod_cl = params.k_mod*params.Mod*(params.f-Val(i));
        reac_prop{i} =horzcat(reac_prop{i},k_mod_cl);
        reac_type{i} = horzcat(reac_type{i},-2);
        reac_from{i} = horzcat(reac_from{i},i);
        reac_with{i} = horzcat(reac_with{i},i);
        reac_with_cluster{i} = horzcat(reac_with_cluster{i},Str(i,i));
        reac_from_cluster{i} = horzcat(reac_from_cluster{i},Str(i,i));
        
    end
    
    %Intra-cluster dephosphorylation reactions
    
    if Val(i) - bonds(i) > 0
        k_demod_cl = params.k_demod*params.Demod*(Val(i) - bonds(i));
        reac_prop{i} =horzcat (reac_prop{i},k_demod_cl);
        reac_type{i} = horzcat(reac_type{i},-1);
        reac_from{i} = horzcat(reac_from{i},i);
        reac_with{i} = horzcat(reac_with{i},i);
        
        reac_with_cluster{i} = horzcat(reac_with_cluster{i},Str(i,i));
        reac_from_cluster{i} = horzcat(reac_from_cluster{i},Str(i,i));
    end
    
    %Intra-cluster cross-link break --> Cluster breakdown
    %Remember to decide which cross-link to break
    if bonds(i)>0
        kcr_off_cl = params.kcr_off*(bonds(i));
        reac_prop{i} =horzcat(reac_prop{i},kcr_off_cl);
        reac_type{i} = horzcat(reac_type{i},2);
        reac_from{i} = horzcat(reac_from{i},i);
        reac_with{i} = horzcat(reac_with{i},i);
        reac_with_cluster{i} = horzcat(reac_with_cluster{i},Str(i,i));
        reac_from_cluster{i} = horzcat(reac_from_cluster{i},Str(i,i));
    end
    
    %Across cluster cross-link formation --> Cluster joining
    if bonds(i) < Val(i)
        for l=i+1:1:length(Str)
            if bonds(l) < Val(l)
                kacross_cr_on_cl = params.kcr_on*(Val(i)-bonds(i))*(Val(l)-bonds(l));
                reac_prop{i} =horzcat(reac_prop{i},kacross_cr_on_cl);
                reac_type{i} = horzcat(reac_type{i},3);
                reac_from{i} = horzcat(reac_from{i},i);
                reac_with{i} = horzcat(reac_with{i},l);
                reac_from_cluster{i} = horzcat(reac_from_cluster{i},Str(i,i));
                reac_with_cluster{i} = horzcat(reac_with_cluster{i},Str(l,l));
            end
        end
        
    end
    
    
    reac_prop{i} = cumsum(reac_prop{i});
    total_reac_prop(i) = reac_prop{i}(end);
    
end

total_reac_prop = cumsum(total_reac_prop);

p = rand(2,1);
reacting_chain = min(find(total_reac_prop>p(1)*total_reac_prop(end)));

%Time step
deltaT = -log(p(2))/total_reac_prop(end);

if reacting_chain~=1
    sub_reaction_index = min(find(reac_prop{reacting_chain}+total_reac_prop(reacting_chain-1)>p(1)*total_reac_prop(end)));
    sub_reaction_type = reac_type{reacting_chain}(sub_reaction_index);
    sub_reaction_with = reac_with{reacting_chain}(sub_reaction_index);
    sub_reaction_from = reac_from{reacting_chain}(sub_reaction_index);
    sub_reaction_with_cluster = reac_with_cluster{reacting_chain}(sub_reaction_index);
    sub_reaction_from_cluster = reac_from_cluster{reacting_chain}(sub_reaction_index);
    
    
    
else
    sub_reaction_index = min(find(reac_prop{reacting_chain}>p(1)*total_reac_prop(end)));
    sub_reaction_type = reac_type{reacting_chain}(sub_reaction_index);
    sub_reaction_with = reac_with{reacting_chain}(sub_reaction_index);
    sub_reaction_from = reac_from{reacting_chain}(sub_reaction_index);
    sub_reaction_with_cluster = reac_with_cluster{reacting_chain}(sub_reaction_index);
    sub_reaction_from_cluster = reac_from_cluster{reacting_chain}(sub_reaction_index);
    
end

rc = reacting_chain;
rtype = sub_reaction_type;

if rtype == -2
    Val(rc) = Val(rc)+1;
elseif rtype == -1
    Val(rc) = Val(rc)-1;
elseif rtype ==2
    %Identify how to split the cluster
    cluster_reac = Str(rc,rc);
    no_of_clusters =max(diag(Str));
    pos_of_neighbour = randi(bonds(rc));
    bond_choice = Str(rc,:);
    bond_choice(rc) = 0;
    bond_choice = cumsum(bond_choice);
    choice_of_neighbour=find(bond_choice>=pos_of_neighbour,1);
    Str(rc,choice_of_neighbour) = Str(rc,choice_of_neighbour) -1;
    Str(choice_of_neighbour,rc) = Str(choice_of_neighbour,rc) -1;
    %Now, we need to check if this cluster is connected or not!
    chains_in_rc = find(diag(Str)==cluster_reac);
    Adjacency = Str(chains_in_rc,chains_in_rc);
    Adj = sparse(Adjacency);
    [S C] = graphconncomp(Adj);
    if S>1
        %cluster has not broken if S=1, continue with calculations
        new_cluster=chains_in_rc(find(C==2));
        for k=1:length(new_cluster)
            Str(new_cluster(k),new_cluster(k))= no_of_clusters+1;
            
        end
    end
    
elseif rtype ==3
    %Identify how to join the clusters
    joining_cluster = sub_reaction_with_cluster;
    joining_chain = sub_reaction_with;
    if joining_cluster ~= Str(rc,rc)
        Str(joining_chain,rc) = Str(joining_chain,rc)+1;
        Str(rc,joining_chain) = Str(rc,joining_chain)+1;
        chains_in_jc=find(diag(Str)==joining_cluster);
        for p=1:length(chains_in_jc)
            Str(chains_in_jc(p),chains_in_jc(p))=Str(rc,rc);
        end
        chains_to_reduce_cluster = find(diag(Str)>joining_cluster);
        diagonalIdx = (chains_to_reduce_cluster-1) * (length(Str) + 1) + 1;
        Str(diagonalIdx)=Str(diagonalIdx)-1;
    else
        Str(joining_chain,rc) = Str(joining_chain,rc)+1;
        Str(rc,joining_chain) = Str(rc,joining_chain)+1;
    end
end

end