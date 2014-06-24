clc
clear
close all

%% Get info, read images, distort images if needed, find features
folder = 'rad params 032514 cells corrected sc';
nameI1 = file_search('sc_rad_reg_VinTL_Vinculin_\w+_w2TVFRET.TIF',folder);
nameI2 = file_search('sc_rad_reg_VinTL_Vinculin_\w+_w3Teal.TIF',folder);
k_out = ones(1,2);

for i = 1:length(nameI1)
    I1 = imread(nameI1{i});
    I2 = imread(nameI2{i});
%     I2b = lensdistort1(I2,5.7504e-05,1.3782,'ftype',5,'bordertype','fit');
%     imwrite2tif(I2b,[],fullfile(folder,['fit_' nameI2{i}]),'single');
    
    lambda = 1; %length scale of noise to be filtered out; typically 1 pixel
    w = 9; % should be a little greater than the radius of the largest features
    f1 = feature2D(I1,lambda,w);
    f2 = feature2D(I2,lambda,w);
%     f2 = feature2D(I2b,lambda,w);
    f1(:,6) = 1;
    f2(:,6) = 2;
    f1(:,7) = 1;
    f2(:,7) = 2;
    out = vertcat(f1,f2);
    
    %% Particle Tracker
    [lub] = trackmem(out,5,2,2,0);
    x1 = lub(1:2:end,1);
    y1 = lub(1:2:end,2);
    x2 = lub(2:2:end,1);
    y2 = lub(2:2:end,2);
    
    figure
    subplot(2,1,1)
    hist(mod(x1,1));
    subplot(2,1,2)
    hist(mod(x2,1));
    pause()
    close
    
    hold on
    scatter(x1,y1,'.');
    scatter(x2,y2,'.');
    pause();
    close
    hold on
    for i = 1:length(x1)
        plot([x1(i) x2(i)], [y1(i) y2(i)], 'k-')
    end 
    pause();
    close
    
    % Shift x1, y1, x2, y2 to center
    center = round(2039/2);
    x1 = x1-center;
    x2 = x2-center;
    y1 = y1-center;
    y2 = y2-center;
    
    %Convert x-y to r-theta polar coordinates
    [theta1,r1] = cart2pol(x1,y1);
    [theta2,r2] = cart2pol(x2,y2);
    
    %Normalize to R
    R = sqrt(2)*center;
    r1 = r1./R;
    r2 = r2./R;
    w = find(r1>0.8);
    r1(w)=[];
    r2(w)=[];
    
    %Scatter plot where
    % y = s = location of undistorted TVFRET channel (r1)
    % x = r = location of distorted FWCy5 channel (r2)
%     scatter(r2,r1);
%     pause();
%     close
    
    %Scatter plot to fit curves to where
    % y = s-x = TVFRETchannel location - FWCy5 channel location (r1-r2)
    % x = r = location of distorted FWCy5 channel (r2)
    scatter(r2,(r1-r2));
    axis([0 0.8 -.002 0.002]);
    pause();
    close
    
    % Use non linear tool to find parameters for curve fit
    x_in = r2;
    y_in = r1-r2;
    
    % F_type = 1;
    % F = @(k,r)r.*(1./(1+k(1).*r))-r;     % #1 s = r*(1/(1+k*r))-r
    % F_type = 2;
    % F = @(k,r)r.*(1./(1+k(1).*r.^2))-r;  % #2 s = r*(1/(1+k*r^2))-r
    % F_type = 3;
    % F = @(k,r)k(1).*r.^2;               % #3 s = r+k*r^2-r
    % F_type = 4;
    % F = @(k,r)k(1).*r.^3;               % #4 s = r+k*r^3-r
    F_type = 5;
    F = @(k,r)k(1).*r.^k(2);            % Custom (mine) s = r+k*r^b-r
    
    k0 = [1 1];
    % k0 = 1;
    [k,resnorm,~,exitflag,output] = lsqcurvefit(F,k0,x_in,y_in);
    disp(['Ftype: ' num2str(F_type)]);
    disp(['k = ' num2str(k)]);
    disp(['resnorm = ' num2str(resnorm)]);
    hold on
    scatter(x_in,y_in);
    scatter(x_in,F(k,x_in));
    axis([0 0.8 -.002 0.002]);
    pause();
    close
k_out = vertcat(k_out,k);
end
k_out(1,:)=[];
disp(['k = ' num2str(-mean(k_out(:,1)))]);
disp(['b = ' num2str(mean(k_out(:,2)))]);
% save(fullfile(folder,'params_60xTeal_cells.txt'),'-ascii','k_out');
save(fullfile(folder,'params_corrected_sc_60xTeal_cells.txt'),'-ascii','k_out');
