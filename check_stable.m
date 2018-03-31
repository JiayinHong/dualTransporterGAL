function flag = check_stable( y1, y2 )

abs_tol = 10^-4;
rel_tol = 10^-4;

% mode = 'crude';
%
% switch mode
%     case 'precise'
%         abs_tol = 10^-5;
%         rel_tol = 10^-5;
%     case 'crude'
%         abs_tol = 10^-4;
%         rel_tol = 10^-4;
%     case 'more crude'
%         abs_tol = 10^-3;
%         rel_tol = 10^-3;
% end

abs_change = abs(y1-y2);
rel_change = abs_change ./ y1;

flag = all(abs_change < abs_tol | rel_change < rel_tol);

end
