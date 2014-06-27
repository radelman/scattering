function cool_plot(xr, yr, xlb, ylb, zlb, phis, angles, cr, cm, ax, ay)
	figure();
	line(xr, yr);
	hold('on');
	axis('equal');
	xlim(xr);
	ylim(yr);
	set(gca, 'fontsize', 30);
	xlabel(xlb);
	ylabel(ylb);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	imh = 0;
	for theta = angles
		if (imh == 0)
			imh = imagesc(xr, yr(2 : -1 : 1), real(phis * exp(-1i * theta)));
			if (ax > 0.0 && ay > 0.0)
				rectangle('position', [-ax, -ay, 2.0 * ax, 2.0 * ay], 'curvature', [1.0, 1.0], 'edgecolor', [0.001, 0.001, 0.001], 'facecolor', [0.001, 0.001, 0.001]);
			end
			draw_border();
			colormap(cm);
			cb = colorbar();
			caxis(cr);
			set(cb, 'ticklength', [0.0125, 0.0125]);
			set(cb, 'layer', 'top');
			set(cb, 'fontsize', 20);
			set(get(cb, 'ylabel'), 'fontsize', 30, 'string', zlb);
		else
			set(imh, 'cdata', real(phis * exp(-1i * theta)));
		end
		drawnow();
		pause(0.1);
	end
end
