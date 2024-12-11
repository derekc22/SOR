function [psi, X, Y, iterationCount, err] = getPsi(r, N0, lambda)

    load('data.mat', 'h0', 'U0', 'nx1', 'nx2', 'ny1', 'l1', 'l2', 'ny2', 'lTot', 'hTot', 'epsilon');

    % Calculate the suction flow speed using the suction ratio
    Us = r * (hTot / h0) * U0;

    % Generate the mesh grid
    Nx = (nx1 + nx2) * N0 + 1; % Total points in x-direction
    Ny = (1 + ny1 + ny2) * N0 + 1; % Total points in y-direction
    
    x = linspace(0, lTot, Nx); % x-direction grid
    y = linspace(0, hTot, Ny); % y-direction grid
    [X, Y] = meshgrid(x, y); % mesh

    dx = h0/N0; 
    dy = h0/N0; 
    
    a = 1 / dx^2;
    b = 1 / dy^2;
    c = 2 * (a + b);

    % Initialize psi with Dirichlet boundary conditions
    psi = zeros(Nx, Ny); % Initialize psi array

    % Bottom boundary B1
    psi(:, 1) = 0;

    % Left boundary (x=0)
    for j = 2:Ny-1
        psi(1, j) = U0*y(j);
    end


    % Right boundary B3
    for j = 2:N0
        psi(Nx, j) = Us*y(j);
    end

    % Wall boundary B4
    psi(nx1 * N0 + 1 : (nx1 + nx2) * N0 + 1, N0 + 1 : (1 + ny1) * N0 + 1) = Us * h0;


    % Top boundary B5
    psi(:, Ny) = U0 * hTot;


    iterationCount = 0;
   

    while true
        iterationCount = iterationCount + 1;
        err = 0; % Reset err for convergence check
    
        % Update psi for Domains D1, D2, and D3
        for i = 2:Nx-1
            for j = 2:Ny-1
                % Check if (i, j) belongs to Domain D1, D2, or D3
                if (i <= nx1 * N0) || (i > nx1 * N0 && j <= N0) || (i > nx1 * N0 && j >= (1 + ny1) * N0 + 2)
                    psi_old = psi(i, j);
                    psi_new = (b * (psi(i, j - 1) + psi(i, j + 1)) + a * (psi(i - 1, j) + psi(i + 1, j))) / c;
                    psi(i, j) = psi_old + lambda * (psi_new - psi_old);
                    err = max(abs(psi(i, j) - psi_old), err);
                end
            end
        end
    
        % Update psi for Domain D4
        i = Nx;
        for j = (1 + ny1) * N0 + 2 : Ny - 1
            psi_old = psi(i, j);
            psi_new = (b * (psi(i, j - 1) + psi(i, j + 1)) + 2 * a * psi(i - 1, j)) / c;
            psi(i, j) = psi_old + lambda * (psi_new - psi_old);
            err = max(abs(psi(i, j) - psi_old), err);
        end
    
        % Check for convergence
        if err < epsilon
            break;
        end
    end

    % Display iteration count and convergence error
    fprintf('iterations = %d, err = %g\n', iterationCount, err);
end
