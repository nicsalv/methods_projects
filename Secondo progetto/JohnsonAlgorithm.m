function [scheduled] = JohnsonAlgorithm(jobs, M1, M2)
    % Inizializzo il numero di job
    nJob = length(jobs);
    
    % Definisco la lista dei job schedulati
    scheduled = zeros(nJob,1);

    % Definisco gli indici di testa e di coda
    head = 1;
    tail = nJob;

    while(~isempty(jobs))
        A = [M1(jobs) M2(jobs)];
       [m, index] = min(A);

       % if i (position of minimun) is < of 7 -> the job will be scheduled on
       % M1 first
       if(index <= nJob)
           scheduled(head) = jobs(index);

           % Remove job from job's list
           jobs(index) = [];

           % increment head index
           head = head + 1;
       else
           jobIndex = index - nJob;
           scheduled(tail) = jobs(jobIndex);

           % Remove job from job's list
           jobs(jobIndex) = [];

           % decrement head index
           tail = tail - 1;
       end
       nJob = nJob -1;
    end
end