function [u, v, Cp] = getUV(psi, X, Y, N0)
    load('data.mat', 'h0', 'U0', 'nx1', 'nx2', 'ny1', 'l1', 'l2', 'ny2', 'lTot', 'hTot', 'epsilon');
    
    % Initialize grid
    Nx = (nx1 + nx2) * N0 + 1;
    Ny = (ny1 + ny2 + 1) * N0 + 1;
    [u, v] = deal(zeros(Nx, Ny));
    
    % Pad island with NaN
    islandRows = nx1*N0 + 2:(nx1 + nx2)*N0;
    islandCols = N0 + 2:(1 + ny1)*N0;
    u(islandRows, islandCols) = NaN;
    v(islandRows, islandCols) = NaN;
    
    % Grid spacing
    [dx, dy] = deal(h0/N0);
    
    % Interior domain using central differences
    u(2:Nx-1, 2:Ny-1) = (psi(2:Nx-1, 3:Ny) - psi(2:Nx-1, 1:Ny-2)) / (2*dy);
    v(2:Nx-1, 2:Ny-1) = -(psi(3:Nx, 2:Ny-1) - psi(1:Nx-2, 2:Ny-1)) / (2*dx);
    
    % Boundary conditions
    % B1 (left)
    u(1, 1:Ny-1) = (psi(1, 2:Ny) - psi(1, 1:Ny-1)) / dy;
    v(1, 2:Ny-1) = -(psi(2, 2:Ny-1) - psi(1, 2:Ny-1)) / dx;
    
    % B2 (top)
    u(2:Nx-1, Ny) = (psi(2:Nx-1, Ny) - psi(2:Nx-1, Ny-1)) / dy;
    v(2:Nx-1, Ny) = -(psi(3:Nx, Ny) - psi(1:Nx-2, Ny)) / (2*dx);
    
    % B3 (right/suction)
    u(Nx, 2:N0) = (psi(Nx, 3:N0+1) - psi(Nx, 1:N0-1)) / (2*dy);
    v(Nx, 2:N0) = -(psi(Nx, 2:N0) - psi(Nx-1, 2:N0)) / dx;
    
    % B4 (island)
    islandBoundary = nx1*N0 + 1;
    validCols = N0+1:(1+ny1)*N0;
    validCols = validCols(validCols < Ny);
    u(islandBoundary, validCols) = (psi(islandBoundary, validCols+1) - psi(islandBoundary, validCols-1)) / (2*dy);
    v(islandBoundary, validCols) = -(psi(islandBoundary+1, validCols) - psi(islandBoundary-1, validCols)) / (2*dx);
    
    % B5 (bottom)
    u(2:Nx-1, 1) = (psi(2:Nx-1, 2) - psi(2:Nx-1, 1)) / dy;
    v(2:Nx-1, 1) = -(psi(3:Nx, 1) - psi(1:Nx-2, 1)) / (2*dx);
    
    % Calculate Cp if requested
    if nargout > 2
        Cp = 1 - (u.^2 + v.^2) / U0^2;
    end
end