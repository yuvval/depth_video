function [primal_vars, min_val] = solveQuadWithEq(P,q,r,Aeq,beq)
    % solves   min    0.5 xPx + qx +r 
    %          s.t.   Aeq x = b
    %
    % Lior Kirsch 2015
    
    [p,n] = size(Aeq);  % p const , dim n
    assert( n== size(P,1));
    assert( n== size(q,1));
    assert( p== size(beq,1));
    
    KKT_mat = [P, Aeq' ; Aeq, sparse(p,p)];
    KKT_b = [-q;beq];
    
    kkt_opt_var = KKT_mat \ KKT_b;
    
    % if KKT_mat is nonsingular there is a unique solution
    % if kkt is singular but solveable - any solution is an optimal
    % solution.
    % if kkt is non-solvable there is problem is unbounded or infeasible.
    primal_vars = kkt_opt_var(1:n);
    min_val = 0.5*primal_vars'*(P*primal_vars) + q'*primal_vars + r;
end