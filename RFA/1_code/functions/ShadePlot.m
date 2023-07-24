function ShadePlot(x1, x2, theColor)
try
	hold on;
	xl = xlim();
	yl = ylim();
	y = [yl(1), yl(1), yl(2), yl(2), yl(1)];
	if theColor ~= 0
		% If the zone is not supposed to be transparent/unshaded/untinted.
		x = [x1, x2, x2, x1, x1];
		fillHandle = fill(x, y, theColor, 'FaceAlpha', 0.1);
		% Put up text labels
		if x1 > 0 && x1 < 255
			text(x1 + 1, 0.95*yl(2), num2str(x1), 'Color', [0, 0.5, 0]);
		end
		if x2 > 0 && x2 < 255
			text(x2 + 1, 0.95*yl(2), num2str(x2), 'Color', [0, 0.5, 0]);
		end
	end
	% Force immediate update of display.
	drawnow;
catch ME
	% Alert the user of the error.
	errorMessage = sprintf('Error in program %s.\nError Message:\n%s', ...
		mfilename, ME.message);
	% Print the error message out to a static text on the GUI.
	fprintf('%s\n', errorMessage);
	uiwait(errordlg(errorMessage));  % Pop up message.
end
return; % from ShadePlot()
end