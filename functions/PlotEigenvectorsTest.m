
function PlotEigenvectorsTest(eigenvalues,eigenvectors,non,scalingoff)

scale = eigenvalues(1)*max(abs(eigenvectors(:,1)));

for i_eig=1:55
    first_eig(1:(non+1):(non*non)) = eigenvectors(1:non,i_eig);
    first_eig = reshape(first_eig,[non non]);
    counter = non;
    for k=1:non
        first_eig(k,k+1:end) = eigenvectors(counter+1:counter+non-k,i_eig) ;
        first_eig(k+1:end,k) = eigenvectors(counter+1:counter+non-k,i_eig) ;
        counter = counter + non - k;
    end   
    first_eig = eigenvalues(i_eig)*first_eig./scale;
    scaled_eigenvectors{i_eig}=first_eig;
    scaled_eigenvectors{i_eig}(:,end+1) = zeros(1,non);
    scaled_eigenvectors{i_eig}(end+1,:) = zeros(non+1,1);
end

figure
set(gcf,'Position',[100 200 1400 400])
for i=1:10
    subplot('position',[(.01+(i-1)*.1) .57 .08 .34])
    pcolor(scaled_eigenvectors{i})
    if ~exist('scalingoff')
        caxis([-1 1])
    end
    axis off
    subplot('position',[(.01+(i-1)*.1) .07 .08 .34])
    pcolor(scaled_eigenvectors{i})
    if ~exist('scalingoff')
        caxis([-1 1])
    end
    axis off
end

end