function sequence = alea(n,repetitions)
% S = alea(n,repetitions) returns a row vector consisting of a sequence S of 'n' elements repeted 'repetitions' times each
% ALEA produces a random sequence constructed in a way that reduce any assimilation effect:
% An element E is never followed by itself, nor by any element that has already followed E.
% see Cross, D. V. (1973). “Sequential dependencies and regression in psychophysical judgments,” Percept. Psychophys. 14, 547–552.
%
%
% R.Trapeau Oct.2011

if (nargin < 2)
    error(' There must be at least two args, n_stimuli and repetitions ');
end
nb_it = 0;

fprintf('\nBuilding pseudo-random sequence')

ntrials = n * repetitions;
nblocks = 1;
if repetitions >= n
    seq_max = n * (n-1);
    nblocks = ceil(ntrials / seq_max);
end

while ceil(ntrials/n/nblocks) > 7
    nblocks = nblocks + 1;
end

if nblocks > 1
    blocks_reps(1:nblocks-1) = ceil(ntrials/n/nblocks);
    blocks_reps(nblocks) = (ntrials - n*blocks_reps(1)*(nblocks-1)) / n;
end

sequence = [];
for block = 1:nblocks
    if nblocks == 1
        size_block = n * repetitions;
    else
        size_block = n * blocks_reps(block);
    end
    
    trouve = 0;
    while ~trouve
        tmp_sequence = zeros(1,size_block);
        occ = zeros(n,1);
        test = zeros(n);
        brk = 0;
        
        if isempty(sequence) % first element
            tmp_sequence(1) = randi(n);
        else
            possibilities = find((1:n)~=sequence(end));
            tmp_sequence(1) = possibilities(randi(length(possibilities)));
        end
        occ(tmp_sequence(1)) = occ(tmp_sequence(1)) + 1;
        
        for i = 2:size_block
            if brk == 1
                break;
            end
            flag = 0;
            echec = zeros(n,1);
            echec(tmp_sequence(i-1)) = 1;
            while flag == 0
                if min(echec)
                    brk = 1;
                    break;
                end
                non_echec = find(~echec);
                candidate = non_echec(randi(length(non_echec)));
                if (occ(candidate) == min(occ)) && (test(candidate,tmp_sequence(i-1)) == 0)
                    tmp_sequence(i) = candidate;
                    occ(candidate) = occ(candidate) + 1;
                    test(tmp_sequence(i),tmp_sequence(i-1)) = test(tmp_sequence(i),tmp_sequence(i-1)) + 1;
                    flag = 1;
                    if i == size_block
                        trouve = 1;
                    end
                else
                    echec(candidate) = echec(candidate) + 1;
                end
            end
        end
        nb_it = nb_it + 1;
        fprintf('.')
        if ~mod(nb_it,50)
            fprintf('\n')
        end
        
        if ~trouve
            if nb_it > 1000 % should never happen (but just in case...)
                fprintf('\nMax number of iterations reached > using basic method');
                for i = 1:blocks_reps(block)
                    perm = randperm(n);
                    j = ((i-1) * n) + 1;
                    k = i * n;
                    tmp_sequence(j:k) = perm;
                end
            end
        end
    end
    sequence = [sequence tmp_sequence];
end

fprintf('\nDone after %i iterations\n\n',nb_it)



