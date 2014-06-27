function move_legend(legendh, corner, distance)
	gcfp = get(gcf, 'position');
	aspect_ratio = gcfp(3) / gcfp(4);
	gcap = get(gca, 'position');
	legendp = get(legendh, 'position');
	new_legendp = legendp;
	switch (corner)
		case 1
			new_legendp(1) = gcap(1) + distance(1);
			new_legendp(2) = gcap(2) + aspect_ratio * distance(2);
		case 2
			new_legendp(1) = gcap(1) + gcap(3) - distance(1) - legendp(3);
			new_legendp(2) = gcap(2) + aspect_ratio * distance(2);
		case 3
			new_legendp(1) = gcap(1) + distance(1);
			new_legendp(2) = gcap(2) + gcap(4) - aspect_ratio * distance(2) - legendp(4);
		case 4
			new_legendp(1) = gcap(1) + gcap(3) - distance(1) - legendp(3);
			new_legendp(2) = gcap(2) + gcap(4) - aspect_ratio * distance(2) - legendp(4);
	end
	set(legendh, 'position', new_legendp);
end
